"""Process a single gene through the full bioinformatics pipeline.

All functions take an explicit work_dir parameter — no os.chdir().
"""

import gzip
import json
import re
import shutil
from io import StringIO
from pathlib import Path

from Bio import Entrez, SeqIO
from loguru import logger

import config
from bam_filtering import filter_bam_for_gene
from gnomad import download_gnomad_for_gene, prepare_gnomad_vcf
from orthologs import OrthologResult, get_source


def _cache_enabled(cfg: dict, key: str) -> bool:
    cache_cfg = cfg.get("cache", {})
    return bool(cache_cfg.get("enabled", True) and cache_cfg.get(key, True))


def _ortholog_cache_paths(gene_id: int, source_name: str) -> tuple[Path, Path]:
    if source_name == "ncbi_datasets":
        base = config.NCBI_ORTHOLOG_CACHE_DIR
    elif source_name == "blast":
        base = config.BLAST_ORTHOLOG_CACHE_DIR
    else:
        raise ValueError(f"Unknown ortholog cache source: {source_name}")
    return base / f"gene_{gene_id}.fastq.gz", base / f"gene_{gene_id}.meta.json"


def _count_fastq_records(fastq_path: Path) -> int:
    lines = 0
    with open(fastq_path) as fh:
        for _ in fh:
            lines += 1
    return lines // 4


def _load_ortholog_cache(gene_id: int, source_name: str, work_dir: Path) -> tuple[int, int] | None:
    cache_fastq_gz, cache_meta = _ortholog_cache_paths(gene_id, source_name)
    if not cache_fastq_gz.exists() or cache_fastq_gz.stat().st_size == 0:
        return None

    work_fastq = work_dir / "genes_from_coordinates.fastq"
    with gzip.open(cache_fastq_gz, "rb") as src, open(work_fastq, "wb") as dst:
        shutil.copyfileobj(src, dst)

    species_count = 0
    sequence_count = 0
    if cache_meta.exists():
        try:
            with open(cache_meta) as f:
                meta = json.load(f)
            species_count = int(meta.get("species_count", 0))
            sequence_count = int(meta.get("sequence_count", 0))
        except Exception:
            species_count = 0
            sequence_count = 0

    if species_count <= 0 or sequence_count <= 0:
        sequence_count = _count_fastq_records(work_fastq)
        species_count = sequence_count

    if sequence_count <= 0:
        return None

    return species_count, sequence_count


def _save_ortholog_cache(
    gene_id: int,
    source_name: str,
    fastq_path: Path,
    species_count: int,
    sequence_count: int,
) -> None:
    cache_fastq_gz, cache_meta = _ortholog_cache_paths(gene_id, source_name)
    cache_fastq_gz.parent.mkdir(parents=True, exist_ok=True)

    with open(fastq_path, "rb") as src, gzip.open(cache_fastq_gz, "wb", compresslevel=1) as dst:
        shutil.copyfileobj(src, dst)

    payload = {
        "gene_id": int(gene_id),
        "source": source_name,
        "species_count": int(species_count),
        "sequence_count": int(sequence_count),
    }
    with open(cache_meta, "w") as f:
        json.dump(payload, f, indent=2)


def _gene_seq_cache_paths(gene_id: int) -> tuple[Path, Path]:
    return (
        config.GENE_SEQ_CACHE_DIR / f"gene_{gene_id}.fasta",
        config.GENE_SEQ_CACHE_DIR / f"gene_{gene_id}.coords.json",
    )


def _load_gene_seq_cache(gene_id: int, work_dir: Path) -> dict | None:
    cache_fasta, cache_coords = _gene_seq_cache_paths(gene_id)
    if not cache_fasta.exists() or not cache_coords.exists():
        return None

    with open(cache_coords) as f:
        coords = json.load(f)

    shutil.copy(cache_fasta, work_dir / "gene_seq.fasta")
    return {
        "chr_acc": coords["chr_acc"],
        "start": int(coords["start"]),
        "end": int(coords["end"]),
    }


def _save_gene_seq_cache(gene_id: int, work_dir: Path, coords: dict) -> None:
    cache_fasta, cache_coords = _gene_seq_cache_paths(gene_id)
    cache_fasta.parent.mkdir(parents=True, exist_ok=True)

    shutil.copy(work_dir / "gene_seq.fasta", cache_fasta)
    with open(cache_coords, "w") as f:
        json.dump(
            {
                "gene_id": int(gene_id),
                "chr_acc": coords["chr_acc"],
                "start": int(coords["start"]),
                "end": int(coords["end"]),
            },
            f,
            indent=2,
        )


def download_gene_seq(gene_id, work_dir, cfg):
    """Download gene sequence and return genomic coordinates."""
    if _cache_enabled(cfg, "gene_seq"):
        try:
            cached_coords = _load_gene_seq_cache(gene_id, work_dir)
            if cached_coords:
                logger.info(f"Using cached gene sequence for gene {gene_id}")
                return cached_coords
        except Exception as e:
            logger.warning(f"Failed to use cached gene sequence for {gene_id}: {e}")

    logger.info(f"Downloading gene sequence for gene {gene_id}")

    handle = Entrez.esummary(db="gene", id=str(gene_id), retmode="xml")
    gi_list = Entrez.read(handle)["DocumentSummarySet"]["DocumentSummary"][0].get("GenomicInfo", [])
    handle.close()

    if not gi_list:
        logger.warning(f"No genomic info found for gene {gene_id}")
        return None

    coords = gi_list[0]
    chr_acc = coords.get("ChrAccVer")
    start = int(coords["ChrStart"])
    end = int(coords["ChrStop"])

    logger.debug(f"Coordinates: {chr_acc}:{start}-{end}")

    handle = Entrez.efetch(
        db="nuccore",
        id=chr_acc,
        rettype="fasta",
        retmode="text",
        seq_start=min(start, end) + 1,
        seq_stop=max(start, end) + 1,
    )
    fasta_text = handle.read()
    handle.close()

    records = list(SeqIO.parse(StringIO(fasta_text), "fasta"))
    SeqIO.write(records, work_dir / "gene_seq.fasta", "fasta")
    logger.success("Saved gene_seq.fasta")

    coords_out = {
        "chr_acc": chr_acc,
        "start": min(start, end),
        "end": max(start, end),
    }

    if _cache_enabled(cfg, "gene_seq"):
        try:
            _save_gene_seq_cache(gene_id, work_dir, coords_out)
        except Exception as e:
            logger.warning(f"Failed to cache gene sequence for {gene_id}: {e}")

    return coords_out


def generate_pseudoreads(input_fastq, work_dir, phred: int = 30):
    """Create sliding window pseudoreads from ortholog sequences."""
    output_fastq = work_dir / "pseudo_reads.fastq"
    READ_LEN = 75
    STEP = 35

    logger.info(f"Generating pseudoreads (read_len={READ_LEN}, step={STEP}, phred={phred})")

    phred_char = chr(phred + 33)
    total_reads = 0

    with open(output_fastq, "w") as out:
        for record in SeqIO.parse(input_fastq, "fastq"):
            seq = str(record.seq)
            n = len(seq)
            read_index = 1

            for start in range(0, n - READ_LEN + 1, STEP):
                read_seq = seq[start : start + READ_LEN]
                qual = phred_char * len(read_seq)
                header = f"@{record.id}_pseudo_{read_index}_{start+1}-{start+READ_LEN}"

                out.write(f"{header}\n{read_seq}\n+\n{qual}\n")
                read_index += 1
                total_reads += 1

    logger.success(f"Generated {total_reads} pseudoreads -> {output_fastq}")


def align_pseudoreads(work_dir):
    """Run BWA alignment pipeline."""
    logger.info("Aligning pseudoreads with BWA")

    ref = work_dir / "gene_seq.fasta"
    reads = work_dir / "pseudo_reads.fastq"
    sam = work_dir / "aln.sam"
    bam = work_dir / "aln.bam"
    sorted_bam = work_dir / "aln.sorted.bam"

    config.run_bio(f"bwa index {ref}")
    config.run_bio(f"bwa mem {ref} {reads} > {sam}")
    config.run_bio(f"samtools view -bS {sam} > {bam}")
    config.run_bio(f"samtools sort {bam} -o {sorted_bam}")
    config.run_bio(f"samtools index {sorted_bam}")

    logger.success(f"Alignment complete -> {sorted_bam}")


def call_variants(work_dir, bam_path: Path, min_var_freq: float = 0.2):
    """Run VarScan variant calling."""
    logger.info("Calling variants with VarScan")

    ref = work_dir / "gene_seq.fasta"
    sorted_bam = Path(bam_path)
    pileup = work_dir / "gene.mpileup"
    snps = work_dir / "gene_snps.vcf"
    indels = work_dir / "gene_indels.vcf"

    if not sorted_bam.exists():
        raise FileNotFoundError(f"BAM for variant calling not found: {sorted_bam}")

    config.run_bio(f"samtools faidx {ref}")
    config.run_bio(f"samtools mpileup -f {ref} {sorted_bam} > {pileup}")
    config.run_bio(
        f"varscan mpileup2snp {pileup} --min-coverage 8 --min-reads2 2 "
        f"--min-var-freq {min_var_freq} --output-vcf 1 > {snps}"
    )
    config.run_bio(
        f"varscan mpileup2indel {pileup} --min-coverage 8 --min-reads2 2 "
        f"--min-var-freq {min_var_freq} --output-vcf 1 > {indels}"
    )

    logger.success(f"Variant calling complete -> {snps}, {indels}")


def normalize_vcf(gene_coords, work_dir):
    """Convert local VCF coordinates to genomic coordinates."""
    logger.info("Normalizing VCF coordinates")

    chr_acc = gene_coords["chr_acc"]
    region_start = gene_coords["start"] + 1  # Convert to 1-based

    new_chrom = _refseq_accession_to_vcf_chrom(chr_acc)

    logger.debug(f"Converting {chr_acc} -> chr {new_chrom}, start={region_start}")

    input_vcf = work_dir / "gene_snps.vcf"
    output_vcf = work_dir / "gene_snps_normalized.vcf"
    output_vcf_gz = work_dir / "gene_snps_normalized.vcf.gz"

    # Get the region chrom name from the VCF
    region_chrom = None
    with open(input_vcf) as f:
        for line in f:
            if not line.startswith("#"):
                region_chrom = line.split("\t")[0]
                break

    if region_chrom is None:
        logger.warning("No variants found in VCF")
        return

    variant_count = 0
    with open(input_vcf) as fin, open(output_vcf, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                if line.startswith("##contig"):
                    line = f"##contig=<ID={new_chrom}>\n"
                fout.write(line)
                continue

            parts = line.rstrip("\n").split("\t")
            if parts[0] == region_chrom:
                parts[0] = new_chrom
                parts[1] = str(region_start + int(parts[1]) - 1)

            fout.write("\t".join(parts) + "\n")
            variant_count += 1

    logger.info(f"Normalized {variant_count} variants")

    config.run_bio(f"bgzip -c {output_vcf} > {output_vcf_gz}")
    config.run_bio(f"tabix -p vcf {output_vcf_gz}")

    logger.success(f"Created indexed VCF: {output_vcf_gz}")


def _refseq_accession_to_vcf_chrom(chr_acc: str) -> str:
    """Map NCBI genomic accession to VCF chromosome naming used by gnomAD/ClinVar."""
    match = re.search(r"NC_0+(\d+)\.", str(chr_acc))
    if not match:
        return chr_acc

    chrom_num = int(match.group(1))
    if chrom_num == 23:
        return "X"
    if chrom_num == 24:
        return "Y"
    if chrom_num in {12920, 1807}:
        return "MT"
    return str(chrom_num)


def annotate_variants(work_dir, gnomad_vcf_gz=None):
    """Annotate pipeline variants with ClinVar and optionally gnomAD.

    Produces a single annotated VCF via one or two bcftools annotate passes:
    1. Annotate with ClinVar (CLNSIG, CLNDN, CLNREVSTAT)
    2. If gnomad_vcf_gz is provided, annotate the ClinVar-annotated VCF with gnomAD (AF, MAF, CSQ)

    Args:
        work_dir: Gene working directory containing gene_snps_normalized.vcf.gz
        gnomad_vcf_gz: Path to prepared (bgzipped+indexed) gnomAD VCF, or None
    """
    logger.info("Annotating variants with ClinVar and gnomAD")

    sample_vcf = work_dir / "gene_snps_normalized.vcf.gz"
    output_vcf = work_dir / "gene_snps_annotated.vcf"

    if not sample_vcf.exists():
        logger.error(f"Normalized VCF not found: {sample_vcf}")
        return

    if gnomad_vcf_gz:
        # bcftools annotate expects bgzipped+indexed input for a second annotation pass.
        clinvar_vcf_gz = work_dir / "gene_snps_clinvar_annotated.vcf.gz"
        config.run_bio(
            f"bcftools annotate -a {config.CLINVAR_VCF}"
            f" -c INFO/CLNSIG,INFO/CLNDN,INFO/CLNREVSTAT"
            f" -o {clinvar_vcf_gz} -O z"
            f" {sample_vcf}"
        )
        config.run_bio(f"tabix -p vcf {clinvar_vcf_gz}")
        config.run_bio(
            f"bcftools annotate -a {gnomad_vcf_gz}"
            f" -c INFO/AF,INFO/MAF,INFO/CSQ"
            f" -o {output_vcf} -O v"
            f" {clinvar_vcf_gz}"
        )
    else:
        config.run_bio(
            f"bcftools annotate -a {config.CLINVAR_VCF}"
            f" -c INFO/CLNSIG,INFO/CLNDN,INFO/CLNREVSTAT"
            f" -o {output_vcf} -O v"
            f" {sample_vcf}"
        )

    # Count annotations in output
    clinvar_count = 0
    gnomad_count = 0
    total = 0

    with open(output_vcf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            total += 1
            info = line.split("\t")[7] if len(line.split("\t")) > 7 else ""
            if "CLNSIG=" in info:
                clinvar_count += 1
            if "AF=" in info:
                gnomad_count += 1

    logger.success(f"Created annotated VCF: {output_vcf}")
    logger.info(f"Total variants: {total}, with ClinVar: {clinvar_count}, with gnomAD AF: {gnomad_count}")


def _fetch_orthologs(gene_id, work_dir, cfg):
    """Fetch orthologs with cache reuse: NCBI cache -> BLAST cache -> live fetches."""
    output_fastq = Path(work_dir) / "genes_from_coordinates.fastq"

    if _cache_enabled(cfg, "orthologs_ncbi"):
        cached = _load_ortholog_cache(gene_id, "ncbi_datasets", Path(work_dir))
        if cached:
            species_count, sequence_count = cached
            logger.info(f"Using cached NCBI orthologs for gene {gene_id}")

            return OrthologResult(
                fastq_path=output_fastq,
                species_count=species_count,
                sequence_count=sequence_count,
                metadata={"source": "ncbi_datasets_cache", "gene_id": gene_id},
            )

    if _cache_enabled(cfg, "orthologs_blast"):
        cached = _load_ortholog_cache(gene_id, "blast", Path(work_dir))
        if cached:
            species_count, sequence_count = cached
            logger.info(f"Using cached BLAST orthologs for gene {gene_id}")

            return OrthologResult(
                fastq_path=output_fastq,
                species_count=species_count,
                sequence_count=sequence_count,
                metadata={"source": "blast_cache", "gene_id": gene_id},
            )

    ncbi_source = get_source("ncbi_datasets")
    try:
        logger.info("Fetching orthologs with ncbi_datasets source")
        result = ncbi_source.fetch(gene_id=gene_id, output_dir=work_dir)
        if result.sequence_count > 0:
            if _cache_enabled(cfg, "orthologs_ncbi"):
                try:
                    if hasattr(ncbi_source, "save_to_cache"):
                        ncbi_source.save_to_cache(
                            gene_id=gene_id,
                            cache_dir=config.NCBI_ORTHOLOG_CACHE_DIR,
                            source_fastq=result.fastq_path,
                            species_count=result.species_count,
                            sequence_count=result.sequence_count,
                            source_tag="ncbi_datasets",
                        )
                    else:
                        _save_ortholog_cache(
                            gene_id=gene_id,
                            source_name="ncbi_datasets",
                            fastq_path=result.fastq_path,
                            species_count=result.species_count,
                            sequence_count=result.sequence_count,
                        )
                except Exception as e:
                    logger.warning(f"Failed to cache NCBI orthologs for gene {gene_id}: {e}")
            return result
        logger.warning("ncbi_datasets returned 0 sequences, falling back to BLAST")
    except Exception as e:
        logger.warning(f"ncbi_datasets failed: {e}, falling back to BLAST")

    logger.info("Fetching orthologs with blast source")
    blast_source = get_source("blast")
    result = blast_source.fetch(
        gene_id=gene_id,
        output_dir=work_dir,
        hitlist_size=cfg.get("hitlist_size", 5000),
        expect=cfg.get("blast_expect", 10.0),
    )
    if _cache_enabled(cfg, "orthologs_blast") and result.sequence_count > 0:
        try:
            _save_ortholog_cache(
                gene_id=gene_id,
                source_name="blast",
                fastq_path=result.fastq_path,
                species_count=result.species_count,
                sequence_count=result.sequence_count,
            )
        except Exception as e:
            logger.warning(f"Failed to cache BLAST orthologs for gene {gene_id}: {e}")
    return result


def _cleanup_gene_outputs(work_dir: Path, keep_paths: list[Path]) -> None:
    """Delete all files/folders in work_dir except the provided keep paths."""
    keep_resolved = {p.resolve() for p in keep_paths}
    removed_files = 0
    removed_dirs = 0

    for path in work_dir.iterdir():
        if path.resolve() in keep_resolved:
            continue
        if path.is_dir():
            shutil.rmtree(path)
            removed_dirs += 1
        else:
            path.unlink()
            removed_files += 1

    logger.info(
        f"Cleaned intermediates in {work_dir}: removed {removed_files} files and {removed_dirs} dirs"
    )


def run_gene(gene_id, work_dir, cfg):
    """Execute full pipeline for one gene.

    Args:
        gene_id: NCBI Gene ID
        work_dir: Directory for this gene's outputs
        cfg: Config dict with keys like 'hitlist_size', 'blast_expect',
             'pseudo_read_phred', 'min_var_freq', 'bam_filtering'
             and optional 'keep_intermediate_files'
    """
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"{'='*50}")
    logger.info(f"Starting pipeline for gene {gene_id}")
    logger.info(f"Working directory: {work_dir}")
    logger.info(f"{'='*50}")

    # Step 1: Download gene sequence
    gene_coords = download_gene_seq(gene_id, work_dir, cfg)

    # Step 2: Obtain homologous sequences (NCBI Datasets -> BLAST fallback)
    result = _fetch_orthologs(gene_id, work_dir, cfg)

    logger.info(f"Retrieved {result.sequence_count} sequences from {result.species_count} species")

    # Step 3: Generate pseudoreads
    generate_pseudoreads(
        input_fastq=str(result.fastq_path),
        work_dir=work_dir,
        phred=cfg.get("pseudo_read_phred", 30),
    )

    # Step 4: Align pseudoreads
    align_pseudoreads(work_dir)

    # Step 5: Optional LIS BAM filtering (enabled by default via config)
    bam_for_variants = work_dir / "aln.sorted.bam"
    bam_filtering_cfg = cfg.get("bam_filtering", {"enabled": False})
    if bam_filtering_cfg.get("enabled", True):
        logger.info("Running LIS BAM filtering before variant calling")
        filter_result = filter_bam_for_gene(
            work_dir=work_dir,
            filtering_cfg=bam_filtering_cfg,
            verbose=False,
        )
        bam_for_variants = filter_result.output_bam
        logger.info(f"Using filtered BAM for variant calling: {bam_for_variants}")
    else:
        logger.info("BAM filtering disabled by config; using aln.sorted.bam")

    # Step 6: Call variants
    call_variants(
        work_dir,
        bam_for_variants,
        min_var_freq=cfg.get("min_var_freq", 0.2),
    )

    # Step 7: Normalize VCF
    if gene_coords:
        normalize_vcf(gene_coords, work_dir)

    # Step 8: Annotate variants with ClinVar and gnomAD
    gnomad_vcf_full = download_gnomad_for_gene(gene_id, gene_coords=gene_coords)
    gnomad_gz = prepare_gnomad_vcf(gnomad_vcf_full, work_dir) if gnomad_vcf_full else None
    annotate_variants(work_dir, gnomad_gz)

    final_vcf = work_dir / "gene_snps_annotated.vcf"
    if not final_vcf.exists():
        raise FileNotFoundError(f"Expected final output not found: {final_vcf}")

    if not cfg.get("keep_intermediate_files", False):
        _cleanup_gene_outputs(work_dir, [final_vcf])

    logger.success(f"Pipeline completed for gene {gene_id}")
