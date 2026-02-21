"""
gnomAD integration module for the bioinformatics pipeline.

Provides:
- NCBI Gene ID to Ensembl ID conversion
- gnomAD variant download via GraphQL API
- VCF preparation for bcftools annotation/intersection
"""

import asyncio
import json
import time
from pathlib import Path
from typing import Any

from gql import gql, Client
from gql.transport.aiohttp import AIOHTTPTransport
from gql.transport.exceptions import TransportQueryError, TransportServerError
from loguru import logger

import config

GNOMAD_API_URL = "https://gnomad.broadinstitute.org/api"
GNOMAD_MAX_RETRIES = 3
GNOMAD_RETRY_DELAY_SEC = 20
GNOMAD_HEADERS = {
    # gnomAD currently rejects default aiohttp client signatures with 403.
    "Accept": "application/json",
    "Content-Type": "application/json",
    "Origin": "https://gnomad.broadinstitute.org",
    "Referer": "https://gnomad.broadinstitute.org/",
    "User-Agent": (
        "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
        "AppleWebKit/537.36 (KHTML, like Gecko) "
        "Chrome/122.0.0.0 Safari/537.36"
    ),
}

GNOMAD_QUERY = gql("""
query VariantsInGene($gene_id: String!) {
  gene(gene_id: $gene_id, reference_genome: GRCh38) {
    gene_id
    symbol
    variants(dataset: gnomad_r4) {
      variant_id
      chrom
      pos
      ref
      alt
      consequence
      hgvsc
      hgvsp

      exome {
        af
      }

      genome {
        af
      }

      joint {
        an
        ac
      }
    }
  }
}
""")


def _cache_enabled(key: str, default: bool = True) -> bool:
    cache_cfg = getattr(config, "CACHE_CFG", {}) or {}
    return bool(cache_cfg.get("enabled", True) and cache_cfg.get(key, default))


def _load_ncbi_ensembl_cache() -> dict[str, list[str]]:
    cache_file = getattr(config, "NCBI_TO_ENSEMBL_CACHE_FILE", None)
    if not cache_file:
        return {}
    path = Path(cache_file)
    if not path.exists():
        return {}
    try:
        with open(path) as f:
            raw = json.load(f)
        if not isinstance(raw, dict):
            return {}
        out: dict[str, list[str]] = {}
        for key, value in raw.items():
            if isinstance(value, list):
                out[str(key)] = [str(v) for v in value if v]
        return out
    except Exception:
        return {}


def _save_ncbi_ensembl_cache(cache_data: dict[str, list[str]]) -> None:
    cache_file = getattr(config, "NCBI_TO_ENSEMBL_CACHE_FILE", None)
    if not cache_file:
        return
    path = Path(cache_file)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w") as f:
        json.dump(cache_data, f, indent=2, sort_keys=True)
    tmp.replace(path)


def _is_primary_chr(chr_name: str | None) -> bool:
    """Return True for primary assembly chromosomes (1..22, X, Y, M/MT)."""
    if not chr_name:
        return False
    c = str(chr_name).strip().upper()
    if c.startswith("CHR"):
        c = c[3:]
    return c in {str(i) for i in range(1, 23)} | {"X", "Y", "M", "MT"}


def _fetch_mygene_record(ncbi_gene_id: int) -> dict[str, Any] | None:
    """Fetch minimal MyGene.info record used for Ensembl mapping."""
    import json
    import urllib.request

    url = (
        "https://mygene.info/v3/gene/"
        f"{ncbi_gene_id}?fields=ensembl.gene,genomic_pos_hg38,genomic_pos"
    )
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=30) as response:
        return json.loads(response.read().decode())


def ncbi_to_ensembl_candidates(ncbi_gene_id: int) -> list[str]:
    """
    Convert NCBI Gene ID to ordered Ensembl Gene ID candidates via MyGene.info.

    Returns primary-assembly IDs first, then remaining IDs in source order.
    """
    logger.debug(f"Looking up Ensembl ID for NCBI gene {ncbi_gene_id} via MyGene.info")

    if _cache_enabled("ncbi_to_ensembl"):
        cached = _load_ncbi_ensembl_cache().get(str(ncbi_gene_id))
        if cached:
            logger.debug(
                "Using cached Ensembl candidates for NCBI gene {}: {}",
                ncbi_gene_id,
                ", ".join(cached),
            )
            return cached

    try:
        data = _fetch_mygene_record(ncbi_gene_id)
        if not data:
            logger.warning(f"Empty MyGene.info response for NCBI gene {ncbi_gene_id}")
            return []

        ensembl_data = data.get("ensembl")
        raw_candidates: list[str] = []
        if isinstance(ensembl_data, list):
            raw_candidates = [item.get("gene") for item in ensembl_data if isinstance(item, dict)]
        elif isinstance(ensembl_data, dict):
            raw_candidates = [ensembl_data.get("gene")]

        # Deduplicate while preserving order.
        seen: set[str] = set()
        candidates = []
        for candidate in raw_candidates:
            if not candidate or candidate in seen:
                continue
            seen.add(candidate)
            candidates.append(candidate)

        if not candidates:
            logger.warning(f"No Ensembl ID found for NCBI gene {ncbi_gene_id}")
            return []

        # Build Ensembl -> chromosome map if present in genomic_pos(_hg38).
        ensembl_chr: dict[str, str] = {}
        for field in ("genomic_pos_hg38", "genomic_pos"):
            positions = data.get(field)
            if isinstance(positions, dict):
                positions = [positions]
            if not isinstance(positions, list):
                continue
            for pos in positions:
                if not isinstance(pos, dict):
                    continue
                ens = pos.get("ensemblgene")
                chr_name = pos.get("chr")
                if ens and chr_name and ens not in ensembl_chr:
                    ensembl_chr[ens] = str(chr_name)

        ranked = sorted(
            candidates,
            key=lambda ens: (
                0 if _is_primary_chr(ensembl_chr.get(ens)) else 1,
                candidates.index(ens),
            ),
        )
        logger.debug(
            "Ensembl candidates for NCBI gene {}: {}",
            ncbi_gene_id,
            ", ".join(f"{ens}({ensembl_chr.get(ens, '?')})" for ens in ranked),
        )
        if _cache_enabled("ncbi_to_ensembl"):
            try:
                cache = _load_ncbi_ensembl_cache()
                cache[str(ncbi_gene_id)] = ranked
                _save_ncbi_ensembl_cache(cache)
            except Exception as e:
                logger.warning(
                    "Failed to write NCBI->Ensembl cache for gene {}: {}",
                    ncbi_gene_id,
                    _format_exception(e),
                )
        return ranked

    except Exception as e:
        logger.error(f"Error looking up Ensembl ID via MyGene.info: {e}")
        return []


def ncbi_to_ensembl(ncbi_gene_id: int) -> str | None:
    """
    Convert NCBI Gene ID to a single preferred Ensembl ID.
    """
    candidates = ncbi_to_ensembl_candidates(ncbi_gene_id)
    if not candidates:
        return None
    ensembl_id = candidates[0]
    logger.debug(f"Found Ensembl ID via MyGene.info: {ensembl_id}")
    return ensembl_id


def _select_af(variant: dict) -> float | None:
    """
    Select allele frequency with priority: exome > genome.
    """
    if variant.get("exome") and variant["exome"].get("af") is not None:
        return variant["exome"]["af"]
    if variant.get("genome") and variant["genome"].get("af") is not None:
        return variant["genome"]["af"]
    return None


def _extract_query_error_messages(exc: TransportQueryError) -> list[str]:
    """Extract GraphQL error message strings from gql query exception."""
    errors: Any = getattr(exc, "errors", None)
    if not errors:
        return []

    out: list[str] = []
    if isinstance(errors, list):
        for item in errors:
            if isinstance(item, dict):
                msg = item.get("message")
                if msg:
                    out.append(str(msg))
            elif item is not None:
                out.append(str(item))
    elif isinstance(errors, dict):
        msg = errors.get("message")
        if msg:
            out.append(str(msg))
    else:
        out.append(str(errors))
    return out


def _format_exception(e: Exception) -> str:
    """Return stable, non-empty exception summary for logs."""
    text = str(e).strip()
    if text:
        return f"{type(e).__name__}: {text}"
    return f"{type(e).__name__}: {repr(e)}"


def _write_vcf(variants: list, out_path: Path) -> None:
    """Write variants to VCF format."""
    with open(out_path, "w") as vcf:
        # header
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write("##source=gnomAD_GraphQL_API\n")
        vcf.write("##reference=GRCh38\n")
        vcf.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n')
        vcf.write('##INFO=<ID=MAF,Number=1,Type=Float,Description="Minor Allele Frequency">\n')
        vcf.write('##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence">\n')
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        for v in variants:
            chrom = v["chrom"]
            pos = v["pos"]
            ref = v["ref"]
            alt = v["alt"]
            vid = v["variant_id"]

            af = _select_af(v)
            maf = af  # For common variants, AF ~ MAF

            info_fields = []
            if af is not None:
                info_fields.append(f"AF={af:.6g}")
            if maf is not None:
                info_fields.append(f"MAF={maf:.6g}")
            if v.get("consequence"):
                info_fields.append(f"CSQ={v['consequence']}")

            info = ";".join(info_fields) if info_fields else "."

            vcf.write(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t.\tPASS\t{info}\n")


async def fetch_gnomad_variants(
    ensembl_id: str,
    output_path_full: Path,
 ) -> tuple[int, bool]:
    """
    Download gnomAD variants for a gene via GraphQL API.
    Writes one VCF file with all variants.

    Args:
        ensembl_id: Ensembl Gene ID (e.g., "ENSG00000012048")
        output_path_full: Path for full VCF file (all variants)

    Returns:
        (variant_count, gene_not_found_flag)
    """
    logger.info(f"Fetching gnomAD variants for {ensembl_id}")

    transport = AIOHTTPTransport(url=GNOMAD_API_URL, headers=GNOMAD_HEADERS)
    client = Client(transport=transport, fetch_schema_from_transport=False)

    try:
        data = await client.execute_async(
            GNOMAD_QUERY,
            variable_values={"gene_id": ensembl_id}
        )

        gene = data.get("gene")
        if not gene:
            # gnomAD GraphQL returns gene=None for some Ensembl IDs.
            logger.info(f"gnomAD has no gene entry for {ensembl_id}")
            return 0, True

        # Collect all variants
        all_variants = gene.get("variants", [])

        # Write full VCF only. Filtering is done downstream in analysis.
        _write_vcf(all_variants, output_path_full)

        logger.success(f"Saved gnomAD variants -> {output_path_full} ({len(all_variants)} total)")

        return len(all_variants), False

    except TransportQueryError as e:
        messages = _extract_query_error_messages(e)
        if any("gene not found" in m.lower() for m in messages):
            logger.info(
                "gnomAD has no gene entry for {} (GraphQL: {})",
                ensembl_id,
                "; ".join(messages),
            )
            return 0, True
        details = "; ".join(messages) if messages else _format_exception(e)
        raise RuntimeError(f"gnomAD query error for {ensembl_id}: {details}") from e
    except (TransportServerError, Exception):
        # Let the caller decide retry policy.
        raise

    finally:
        await transport.close()


def download_gnomad_for_gene(ncbi_gene_id: int) -> Path | None:
    """
    Download gnomAD variants for a gene (sync wrapper).
    Uses cache if available.

    Args:
        ncbi_gene_id: NCBI Gene ID

    Returns:
        Full gnomAD VCF path, or None if download failed
    """
    ensembl_candidates = ncbi_to_ensembl_candidates(ncbi_gene_id)
    if not ensembl_candidates:
        logger.warning(f"Cannot download gnomAD: no Ensembl ID for gene {ncbi_gene_id}")
        return None

    use_cache = _cache_enabled("gnomad")
    config.GNOMAD_CACHE_DIR.mkdir(parents=True, exist_ok=True)

    def _cached_vcf_has_variants(path: Path) -> bool:
        try:
            with open(path) as handle:
                return any(not line.startswith("#") for line in handle)
        except Exception:
            return False

    # Check cache for each candidate in order.
    if use_cache:
        for ensembl_id in ensembl_candidates:
            output_vcf_full = config.GNOMAD_CACHE_DIR / f"{ensembl_id}_full.vcf"
            if not output_vcf_full.exists():
                continue
            if _cached_vcf_has_variants(output_vcf_full):
                logger.info(f"Using cached gnomAD data: {output_vcf_full}")
                return output_vcf_full
            logger.warning(f"Ignoring empty cached gnomAD VCF: {output_vcf_full}")

    def _run_fetch_once(ensembl_id: str, output_vcf_full: Path) -> tuple[int, bool]:
        """Run async fetch in both script and already-running-loop contexts."""
        try:
            # If we're already in a running loop (e.g. notebook), run in thread.
            _ = asyncio.get_running_loop()
            import concurrent.futures

            with concurrent.futures.ThreadPoolExecutor() as executor:
                future = executor.submit(
                    asyncio.run,
                    fetch_gnomad_variants(ensembl_id, output_vcf_full),
                )
                return future.result()
        except RuntimeError:
            # No running loop - standard script context.
            return asyncio.run(fetch_gnomad_variants(ensembl_id, output_vcf_full))

    failures: list[str] = []
    for ensembl_id in ensembl_candidates:
        output_vcf_full = config.GNOMAD_CACHE_DIR / f"{ensembl_id}_full.vcf"
        if not use_cache and output_vcf_full.exists():
            output_vcf_full.unlink(missing_ok=True)
        last_error: str | None = None
        failed_by_exception = False

        for attempt in range(1, GNOMAD_MAX_RETRIES + 1):
            try:
                full_count, gene_not_found = _run_fetch_once(ensembl_id, output_vcf_full)
                if gene_not_found:
                    logger.info(f"gnomAD has no data for candidate {ensembl_id}, trying next ID")
                    break
                if full_count == 0:
                    logger.warning(f"gnomAD returned 0 variants for candidate {ensembl_id}, trying next ID")
                    output_vcf_full.unlink(missing_ok=True)
                    break
                return output_vcf_full
            except Exception as e:
                failed_by_exception = True
                last_error = _format_exception(e)
                if attempt < GNOMAD_MAX_RETRIES:
                    logger.warning(
                        "gnomAD fetch attempt {}/{} failed for {}: {}. Retrying in {}s",
                        attempt,
                        GNOMAD_MAX_RETRIES,
                        ensembl_id,
                        last_error,
                        GNOMAD_RETRY_DELAY_SEC,
                    )
                    time.sleep(GNOMAD_RETRY_DELAY_SEC)
                else:
                    failures.append(f"{ensembl_id}: {last_error}")

        if failed_by_exception and last_error:
            logger.warning(
                "gnomAD candidate {} failed after {} attempts: {}",
                ensembl_id,
                GNOMAD_MAX_RETRIES,
                last_error,
            )

    if failures:
        logger.warning(
            "gnomAD fetch failed for all candidate IDs of NCBI gene {}: {}",
            ncbi_gene_id,
            " | ".join(failures[:3]),
        )
    else:
        logger.warning(
            "No gnomAD variants found for NCBI gene {} across candidate IDs: {}",
            ncbi_gene_id,
            ", ".join(ensembl_candidates),
        )
    return None


def prepare_gnomad_vcf(vcf_path: Path, output_dir: Path, suffix: str = "") -> Path | None:
    """
    Prepare gnomAD VCF for bcftools intersection.

    - Copies to output directory
    - Compresses with bgzip
    - Indexes with tabix

    Args:
        vcf_path: Path to source gnomAD VCF
        output_dir: Directory to write prepared VCF
        suffix: Optional suffix for output filename

    Returns:
        Path to compressed/indexed VCF, or None if failed
    """
    if not vcf_path.exists():
        logger.error(f"gnomAD VCF not found: {vcf_path}")
        return None

    output_vcf = output_dir / f"gnomad_gene{suffix}.vcf"
    output_vcf_gz = output_dir / f"gnomad_gene{suffix}.vcf.gz"

    # Copy VCF to output directory
    import shutil
    shutil.copy(vcf_path, output_vcf)

    # Compress and index
    config.run_bio(f"bgzip -c {output_vcf} > {output_vcf_gz}")
    config.run_bio(f"tabix -p vcf {output_vcf_gz}")

    logger.debug(f"Prepared gnomAD VCF: {output_vcf_gz}")
    return output_vcf_gz
