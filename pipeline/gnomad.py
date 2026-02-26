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
GNOMAD_REGION_MIN_WINDOW_BP = 500
GNOMAD_RATE_LIMIT_SLEEP_SEC = 65
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

GNOMAD_GENE_QUERY = gql("""
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

GNOMAD_REGION_QUERY = gql("""
query VariantsInRegion($chrom: String!, $start: Int!, $stop: Int!) {
  region(chrom: $chrom, start: $start, stop: $stop, reference_genome: GRCh38) {
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


def _cached_vcf_has_variants(path: Path) -> bool:
    try:
        with open(path) as handle:
            return any(not line.startswith("#") for line in handle)
    except Exception:
        return False


def _refseq_accession_to_gnomad_chrom(chr_acc: str | None) -> str | None:
    if not chr_acc:
        return None
    c = str(chr_acc).strip()
    if c.startswith("chr"):
        c = c[3:]
    if c in {"X", "Y", "MT", "M"}:
        return "MT" if c == "M" else c
    if c.isdigit():
        value = int(c)
        if value == 23:
            return "X"
        if value == 24:
            return "Y"
        return str(value)

    import re

    match = re.search(r"NC_0+(\d+)\.", c)
    if not match:
        return None
    chrom_num = int(match.group(1))
    if chrom_num == 23:
        return "X"
    if chrom_num == 24:
        return "Y"
    if chrom_num in {12920, 1807}:
        return "MT"
    return str(chrom_num)


def _is_rate_limit_error(messages: list[str]) -> bool:
    joined = " | ".join(m.lower() for m in messages)
    return "rate limit" in joined


def _is_region_split_error(messages: list[str]) -> bool:
    joined = " | ".join(m.lower() for m in messages)
    return (
        "select a smaller region" in joined
        or "too many variants" in joined
        or "temporarily unavailable" in joined
    )


def _variant_key(v: dict[str, Any]) -> tuple[str, int, str, str]:
    return (str(v.get("chrom")), int(v.get("pos")), str(v.get("ref")), str(v.get("alt")))


def _dedupe_variants(variants: list[dict[str, Any]]) -> list[dict[str, Any]]:
    seen: set[tuple[str, int, str, str]] = set()
    out: list[dict[str, Any]] = []
    for v in variants:
        key = _variant_key(v)
        if key in seen:
            continue
        seen.add(key)
        out.append(v)
    return out


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


async def _fetch_region_variants_recursive(
    client: Client,
    chrom: str,
    start: int,
    stop: int,
    *,
    rate_limit_retries: int = 3,
) -> list[dict[str, Any]]:
    try:
        data = await client.execute_async(
            GNOMAD_REGION_QUERY,
            variable_values={"chrom": chrom, "start": int(start), "stop": int(stop)},
        )
    except TransportQueryError as e:
        messages = _extract_query_error_messages(e)
        if _is_rate_limit_error(messages) and rate_limit_retries > 0:
            logger.warning(
                "gnomAD rate limit for region {}:{}-{}. Retrying in {}s ({} retries left)",
                chrom,
                start,
                stop,
                GNOMAD_RATE_LIMIT_SLEEP_SEC,
                rate_limit_retries,
            )
            await asyncio.sleep(GNOMAD_RATE_LIMIT_SLEEP_SEC)
            return await _fetch_region_variants_recursive(
                client,
                chrom,
                start,
                stop,
                rate_limit_retries=rate_limit_retries - 1,
            )

        span = stop - start + 1
        if _is_region_split_error(messages) and span > GNOMAD_REGION_MIN_WINDOW_BP and start < stop:
            mid = (start + stop) // 2
            logger.info(
                "Splitting gnomAD region {}:{}-{} due to API constraint ({}).",
                chrom,
                start,
                stop,
                "; ".join(messages) if messages else "no details",
            )
            left = await _fetch_region_variants_recursive(client, chrom, start, mid)
            right = await _fetch_region_variants_recursive(client, chrom, mid + 1, stop)
            return left + right

        details = "; ".join(messages) if messages else _format_exception(e)
        raise RuntimeError(f"gnomAD region query error for {chrom}:{start}-{stop}: {details}") from e

    region = data.get("region")
    if not region:
        return []
    variants = region.get("variants", [])
    if not isinstance(variants, list):
        return []
    return variants


async def fetch_gnomad_variants_in_region(
    chrom: str,
    start: int,
    stop: int,
    output_path_full: Path,
) -> int:
    """Download gnomAD variants for a genomic region via GraphQL API."""
    logger.info(f"Fetching gnomAD variants for region {chrom}:{start}-{stop}")

    transport = AIOHTTPTransport(url=GNOMAD_API_URL, headers=GNOMAD_HEADERS)
    client = Client(transport=transport, fetch_schema_from_transport=False)

    try:
        variants = await _fetch_region_variants_recursive(client, chrom, int(start), int(stop))
        variants = _dedupe_variants(variants)
        variants.sort(key=lambda v: _variant_key(v))
        _write_vcf(variants, output_path_full)
        logger.success(f"Saved gnomAD variants -> {output_path_full} ({len(variants)} total)")
        return len(variants)
    finally:
        await transport.close()


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
            GNOMAD_GENE_QUERY,
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


def _run_async(coro):
    """Run async coroutine in script and already-running-loop contexts."""
    try:
        _ = asyncio.get_running_loop()
        import concurrent.futures

        with concurrent.futures.ThreadPoolExecutor() as executor:
            future = executor.submit(asyncio.run, coro)
            return future.result()
    except RuntimeError:
        return asyncio.run(coro)


def _download_gnomad_for_gene_via_region(ncbi_gene_id: int, gene_coords: dict[str, Any]) -> Path | None:
    chrom = _refseq_accession_to_gnomad_chrom(gene_coords.get("chr_acc"))
    if not chrom:
        logger.warning(f"Cannot derive gnomAD chromosome from coords for gene {ncbi_gene_id}: {gene_coords}")
        return None

    region_start = int(min(gene_coords["start"], gene_coords["end"])) + 1
    region_stop = int(max(gene_coords["start"], gene_coords["end"])) + 1
    output_vcf_full = (
        config.GNOMAD_CACHE_DIR / f"ncbi_{ncbi_gene_id}_{chrom}_{region_start}_{region_stop}_full.vcf"
    )
    use_cache = _cache_enabled("gnomad")

    if use_cache and output_vcf_full.exists() and _cached_vcf_has_variants(output_vcf_full):
        logger.info(f"Using cached gnomAD region data: {output_vcf_full}")
        return output_vcf_full
    if use_cache and output_vcf_full.exists():
        logger.warning(f"Ignoring empty cached gnomAD region VCF: {output_vcf_full}")
    if not use_cache and output_vcf_full.exists():
        output_vcf_full.unlink(missing_ok=True)

    last_error: str | None = None
    for attempt in range(1, GNOMAD_MAX_RETRIES + 1):
        try:
            full_count = _run_async(
                fetch_gnomad_variants_in_region(chrom, region_start, region_stop, output_vcf_full)
            )
            if full_count == 0:
                logger.warning(
                    "gnomAD region {}:{}-{} returned 0 variants for gene {}",
                    chrom,
                    region_start,
                    region_stop,
                    ncbi_gene_id,
                )
                output_vcf_full.unlink(missing_ok=True)
                return None
            return output_vcf_full
        except Exception as e:
            last_error = _format_exception(e)
            if attempt < GNOMAD_MAX_RETRIES:
                logger.warning(
                    "gnomAD region fetch attempt {}/{} failed for gene {} ({}:{}-{}): {}. Retrying in {}s",
                    attempt,
                    GNOMAD_MAX_RETRIES,
                    ncbi_gene_id,
                    chrom,
                    region_start,
                    region_stop,
                    last_error,
                    GNOMAD_RETRY_DELAY_SEC,
                )
                time.sleep(GNOMAD_RETRY_DELAY_SEC)
            else:
                logger.warning(
                    "gnomAD region fetch failed for gene {} after {} attempts ({}:{}-{}): {}",
                    ncbi_gene_id,
                    GNOMAD_MAX_RETRIES,
                    chrom,
                    region_start,
                    region_stop,
                    last_error,
                )
    return None


def _download_gnomad_for_gene_via_gene_endpoint(ncbi_gene_id: int) -> Path | None:
    ensembl_candidates = ncbi_to_ensembl_candidates(ncbi_gene_id)
    if not ensembl_candidates:
        logger.warning(f"Cannot download gnomAD: no Ensembl ID for gene {ncbi_gene_id}")
        return None

    use_cache = _cache_enabled("gnomad")

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

    failures: list[str] = []
    for ensembl_id in ensembl_candidates:
        output_vcf_full = config.GNOMAD_CACHE_DIR / f"{ensembl_id}_full.vcf"
        if not use_cache and output_vcf_full.exists():
            output_vcf_full.unlink(missing_ok=True)
        last_error: str | None = None
        failed_by_exception = False

        for attempt in range(1, GNOMAD_MAX_RETRIES + 1):
            try:
                full_count, gene_not_found = _run_async(fetch_gnomad_variants(ensembl_id, output_vcf_full))
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


def download_gnomad_for_gene(ncbi_gene_id: int, gene_coords: dict[str, Any] | None = None) -> Path | None:
    """
    Download gnomAD variants for a gene (sync wrapper).

    Strategy:
    1. Preferred: fetch by genomic region from gene_coords (broader and more complete).
    2. Fallback: fetch by Ensembl gene endpoint.
    """
    config.GNOMAD_CACHE_DIR.mkdir(parents=True, exist_ok=True)

    if gene_coords:
        region_vcf = _download_gnomad_for_gene_via_region(ncbi_gene_id, gene_coords)
        if region_vcf:
            return region_vcf
        logger.warning(
            "Falling back to gnomAD gene endpoint for NCBI gene {} after region fetch failure/empty result",
            ncbi_gene_id,
        )

    return _download_gnomad_for_gene_via_gene_endpoint(ncbi_gene_id)


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
