"""Process a single gene through the full bioinformatics pipeline.

All functions take an explicit work_dir parameter — no os.chdir().
"""

import re
import shutil
from io import StringIO
from pathlib import Path

from Bio import Entrez, SeqIO
from loguru import logger

import config
from gnomad import download_gnomad_for_gene, prepare_gnomad_vcf
from orthologs import get_source


def download_gene_seq(gene_id, work_dir):
    """Download gene sequence and return genomic coordinates."""
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

    return {
        "chr_acc": chr_acc,
        "start": min(start, end),
        "end": max(start, end),
    }


def generate_pseudoreads(input_fastq, work_dir):
    """Create sliding window pseudoreads from ortholog sequences."""
    output_fastq = work_dir / "pseudo_reads.fastq"
    READ_LEN = 75
    STEP = 35
    PHRED = 30

    logger.info(f"Generating pseudoreads (read_len={READ_LEN}, step={STEP})")

    phred_char = chr(PHRED + 33)
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


def call_variants(work_dir):
    """Run VarScan variant calling."""
    logger.info("Calling variants with VarScan")

    ref = work_dir / "gene_seq.fasta"
    sorted_bam = work_dir / "aln.sorted.bam"
    pileup = work_dir / "gene.mpileup"
    snps = work_dir / "gene_snps.vcf"
    indels = work_dir / "gene_indels.vcf"

    config.run_bio(f"samtools faidx {ref}")
    config.run_bio(f"samtools mpileup -f {ref} {sorted_bam} > {pileup}")
    config.run_bio(f"varscan mpileup2snp {pileup} --min-coverage 8 --min-reads2 2 --output-vcf 1 > {snps}")
    config.run_bio(f"varscan mpileup2indel {pileup} --min-coverage 8 --min-reads2 2 --output-vcf 1 > {indels}")

    logger.success(f"Variant calling complete -> {snps}, {indels}")


def normalize_vcf(gene_coords, work_dir):
    """Convert local VCF coordinates to genomic coordinates."""
    logger.info("Normalizing VCF coordinates")

    chr_acc = gene_coords["chr_acc"]
    region_start = gene_coords["start"] + 1  # Convert to 1-based

    match = re.search(r"NC_0+(\d+)\.", chr_acc)
    new_chrom = match.group(1) if match else chr_acc

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
    """Fetch orthologs: try NCBI Datasets first, fall back to BLAST on failure."""
    ncbi_source = get_source("ncbi_datasets")
    try:
        logger.info("Fetching orthologs with ncbi_datasets source")
        result = ncbi_source.fetch(gene_id=gene_id, output_dir=work_dir)
        if result.sequence_count > 0:
            return result
        logger.warning("ncbi_datasets returned 0 sequences, falling back to BLAST")
    except Exception as e:
        logger.warning(f"ncbi_datasets failed: {e}, falling back to BLAST")

    logger.info("Fetching orthologs with blast source")
    blast_source = get_source("blast")
    return blast_source.fetch(
        gene_id=gene_id,
        output_dir=work_dir,
        hitlist_size=cfg.get("hitlist_size", 5000),
    )


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
        cfg: Config dict with keys like 'hitlist_size' and optional 'keep_intermediate_files'
    """
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"{'='*50}")
    logger.info(f"Starting pipeline for gene {gene_id}")
    logger.info(f"Working directory: {work_dir}")
    logger.info(f"{'='*50}")

    # Step 1: Download gene sequence
    gene_coords = download_gene_seq(gene_id, work_dir)

    # Step 2: Obtain homologous sequences (NCBI Datasets -> BLAST fallback)
    result = _fetch_orthologs(gene_id, work_dir, cfg)

    logger.info(f"Retrieved {result.sequence_count} sequences from {result.species_count} species")

    # Step 3: Generate pseudoreads
    generate_pseudoreads(input_fastq=str(result.fastq_path), work_dir=work_dir)

    # Step 4: Align pseudoreads
    align_pseudoreads(work_dir)

    # Step 5: Call variants
    call_variants(work_dir)

    # Step 6: Normalize VCF
    if gene_coords:
        normalize_vcf(gene_coords, work_dir)

    # Step 7: Annotate variants with ClinVar and gnomAD
    gnomad_vcf_full, _ = download_gnomad_for_gene(gene_id)
    gnomad_gz = prepare_gnomad_vcf(gnomad_vcf_full, work_dir) if gnomad_vcf_full else None
    annotate_variants(work_dir, gnomad_gz)

    final_vcf = work_dir / "gene_snps_annotated.vcf"
    if not final_vcf.exists():
        raise FileNotFoundError(f"Expected final output not found: {final_vcf}")

    if not cfg.get("keep_intermediate_files", False):
        _cleanup_gene_outputs(work_dir, [final_vcf])

    logger.success(f"Pipeline completed for gene {gene_id}")
