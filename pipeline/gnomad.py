"""
gnomAD integration module for the bioinformatics pipeline.

Provides:
- NCBI Gene ID to Ensembl ID conversion
- gnomAD variant download via GraphQL API
- VCF preparation for bcftools annotation/intersection
"""

import asyncio
from pathlib import Path

from gql import gql, Client
from gql.transport.aiohttp import AIOHTTPTransport
from gql.transport.exceptions import TransportQueryError
from loguru import logger

import config

GNOMAD_API_URL = "https://gnomad.broadinstitute.org/api"

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


def ncbi_to_ensembl(ncbi_gene_id: int) -> str | None:
    """
    Convert NCBI Gene ID to Ensembl Gene ID using MyGene.info API.

    MyGene.info provides reliable cross-references between gene databases.
    See: https://mygene.info/

    Args:
        ncbi_gene_id: NCBI Gene ID (e.g., 672 for BRCA1)

    Returns:
        Ensembl Gene ID (e.g., "ENSG00000012048") or None if not found
    """
    import urllib.request
    import json

    logger.debug(f"Looking up Ensembl ID for NCBI gene {ncbi_gene_id} via MyGene.info")

    try:
        # MyGene.info API endpoint for gene query
        url = f"https://mygene.info/v3/gene/{ncbi_gene_id}?fields=ensembl.gene"

        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=30) as response:
            data = json.loads(response.read().decode())

        # Extract Ensembl gene ID from response
        ensembl_data = data.get("ensembl")

        if ensembl_data:
            # Can be a single dict or a list of dicts
            if isinstance(ensembl_data, list):
                # Take the first one (primary)
                ensembl_id = ensembl_data[0].get("gene")
            else:
                ensembl_id = ensembl_data.get("gene")

            if ensembl_id:
                logger.debug(f"Found Ensembl ID via MyGene.info: {ensembl_id}")
                return ensembl_id

        logger.warning(f"No Ensembl ID found for NCBI gene {ncbi_gene_id}")
        return None

    except Exception as e:
        logger.error(f"Error looking up Ensembl ID via MyGene.info: {e}")
        return None


def _select_af(variant: dict) -> float | None:
    """
    Select allele frequency with priority: exome > genome.
    """
    if variant.get("exome") and variant["exome"].get("af") is not None:
        return variant["exome"]["af"]
    if variant.get("genome") and variant["genome"].get("af") is not None:
        return variant["genome"]["af"]
    return None


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
 ) -> int:
    """
    Download gnomAD variants for a gene via GraphQL API.
    Writes one VCF file with all variants.

    Args:
        ensembl_id: Ensembl Gene ID (e.g., "ENSG00000012048")
        output_path_full: Path for full VCF file (all variants)

    Returns:
        Number of downloaded variants
    """
    logger.info(f"Fetching gnomAD variants for {ensembl_id}")

    transport = AIOHTTPTransport(url=GNOMAD_API_URL)
    client = Client(transport=transport, fetch_schema_from_transport=False)

    try:
        data = await client.execute_async(
            GNOMAD_QUERY,
            variable_values={"gene_id": ensembl_id}
        )

        gene = data.get("gene")
        if not gene:
            logger.warning(f"No gene data returned for {ensembl_id}")
            return 0

        # Collect all variants
        all_variants = gene.get("variants", [])

        # Write full VCF only. Filtering is done downstream in analysis.
        _write_vcf(all_variants, output_path_full)

        logger.success(f"Saved gnomAD variants -> {output_path_full} ({len(all_variants)} total)")

        return len(all_variants)

    except TransportQueryError as e:
        logger.warning(f"gnomAD query failed for {ensembl_id}: {e}")
        return 0

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
    ensembl_id = ncbi_to_ensembl(ncbi_gene_id)
    if not ensembl_id:
        logger.warning(f"Cannot download gnomAD: no Ensembl ID for gene {ncbi_gene_id}")
        return None

    config.GNOMAD_CACHE_DIR.mkdir(parents=True, exist_ok=True)

    output_vcf_full = config.GNOMAD_CACHE_DIR / f"{ensembl_id}_full.vcf"

    # Check cache
    if output_vcf_full.exists():
        logger.info(f"Using cached gnomAD data: {output_vcf_full}")
        return output_vcf_full

    # Download variants - handle both notebook (running loop) and script contexts.
    try:
        # Check if we're in a running event loop (e.g., Jupyter)
        _ = asyncio.get_running_loop()
        # We're in a running loop - use nest_asyncio or run in new thread
        import concurrent.futures
        with concurrent.futures.ThreadPoolExecutor() as executor:
            future = executor.submit(
                asyncio.run,
                fetch_gnomad_variants(ensembl_id, output_vcf_full)
            )
            full_count = future.result()
    except RuntimeError:
        # No running loop - safe to use asyncio.run()
        full_count = asyncio.run(
            fetch_gnomad_variants(ensembl_id, output_vcf_full)
        )

    if full_count == 0:
        logger.warning(f"No gnomAD variants found for {ensembl_id}")
        return None

    return output_vcf_full


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
