"""NCBI Datasets CLI-based ortholog source implementation."""

import shutil
import subprocess
import zipfile
from pathlib import Path

from Bio import SeqIO
from loguru import logger

from .base import OrthologResult, OrthologSource
from .registry import register


@register
class NCBIDatasetsOrthologSource(OrthologSource):
    """Retrieve orthologs via NCBI Datasets CLI.

    This source:
    1. Runs `datasets download gene gene-id {id} --ortholog all --include gene`
    2. Reads gene sequences directly from gene.fna (single download vs many Entrez.efetch calls)
    3. Deduplicates sequences (one per organism, preferring NC_ over NW_ accessions)
    4. Outputs sequences in FASTQ format

    Requires the NCBI Datasets CLI to be installed (e.g. via conda/micromamba):
    https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/
    """

    name = "ncbi_datasets"

    @staticmethod
    def _candidate_commands() -> list[str]:
        """Return candidate command paths for the datasets CLI."""
        repo_root = Path(__file__).resolve().parents[2]
        return [
            "datasets",
            str(repo_root / "data" / "datasets"),
        ]

    def _find_datasets_cmd(self) -> str | None:
        """Locate a usable datasets CLI executable."""
        for candidate in self._candidate_commands():
            resolved = shutil.which(candidate)
            if resolved:
                return resolved
            path = Path(candidate).expanduser()
            if path.is_file() and path.stat().st_mode & 0o111:
                return str(path)
        return None

    def is_available(self) -> bool:
        """Check if the datasets CLI is installed and in PATH."""
        return self._find_datasets_cmd() is not None

    def _get_datasets_cmd(self) -> str:
        """Get the datasets CLI command path."""
        datasets_cmd = self._find_datasets_cmd()
        if datasets_cmd:
            return datasets_cmd
        candidates = ", ".join(self._candidate_commands())
        raise RuntimeError(f"datasets CLI not found. Tried: {candidates}")

    def fetch(
        self,
        gene_id: int,
        output_dir: Path,
        phred: int = 30,
    ) -> OrthologResult:
        """Fetch orthologs via NCBI Datasets CLI.

        Args:
            gene_id: NCBI Gene ID
            output_dir: Directory to write output files
            phred: Quality score for FASTQ output

        Returns:
            OrthologResult with path to FASTQ file
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        zip_path = output_dir / "ncbi_dataset.zip"
        data_dir = output_dir / "ncbi_dataset" / "data"
        report_path = data_dir / "data_report.jsonl"
        fastq_path = output_dir / "genes_from_coordinates.fastq"

        # Step 1: Download ortholog data
        self._download_orthologs(gene_id, output_dir, zip_path)

        # Step 2: Extract zip
        self._extract_zip(zip_path, output_dir)

        # Step 3: Read gene sequences from gene.fna and convert to FASTQ
        species_count, sequence_count = self._fetch_sequences(
            data_dir, fastq_path, phred
        )

        return OrthologResult(
            fastq_path=fastq_path,
            species_count=species_count,
            sequence_count=sequence_count,
            metadata={"source": "ncbi_datasets", "gene_id": gene_id},
        )

    def _download_orthologs(
        self, gene_id: int, output_dir: Path, zip_path: Path
    ) -> None:
        """Download ortholog data using datasets CLI."""
        logger.info(f"Downloading orthologs for gene {gene_id} via NCBI Datasets CLI")

        datasets_cmd = self._get_datasets_cmd()
        cmd = [
            datasets_cmd,
            "download",
            "gene",
            "gene-id",
            str(gene_id),
            "--ortholog",
            "all",
            "--include",
            "gene",  # Downloads gene.fna with all sequences
            "--filename",
            str(zip_path),
        ]

        logger.debug(f"Running: {' '.join(cmd)}")

        result = subprocess.run(cmd, capture_output=True, text=True, cwd=output_dir)

        if result.returncode != 0:
            logger.error(f"datasets CLI failed: {result.stderr}")
            raise RuntimeError(f"datasets download failed: {result.stderr}")

        logger.success(f"Downloaded to {zip_path}")

    def _extract_zip(self, zip_path: Path, output_dir: Path) -> None:
        """Extract the downloaded zip file."""
        logger.debug(f"Extracting {zip_path}")

        with zipfile.ZipFile(zip_path, "r") as zf:
            zf.extractall(output_dir)

        logger.debug("Extraction complete")

    def _fetch_sequences(
        self, data_dir: Path, output_path: Path, phred: int
    ) -> tuple[int, int]:
        """Read gene sequences from gene.fna and convert to FASTQ.

        Deduplicates sequences: keeps one per organism, preferring NC_ accessions
        over NW_ accessions, then preferring longer sequences.
        """
        gene_fna = data_dir / "gene.fna"

        if not gene_fna.exists():
            logger.error(f"gene.fna not found at {gene_fna}")
            raise FileNotFoundError(f"gene.fna not found at {gene_fna}")

        logger.info(f"Reading sequences from {gene_fna}")

        # Deduplicate: one sequence per organism, prefer NC_ accessions
        by_organism: dict[str, SeqIO.SeqRecord] = {}
        total_records = 0

        for record in SeqIO.parse(gene_fna, "fasta"):
            total_records += 1
            desc = record.description

            # Extract organism name from description
            if "[organism=" not in desc:
                continue

            org = desc.split("[organism=")[1].split("]")[0]

            # Skip human sequences
            if org.lower() == "homo sapiens":
                continue

            # Deduplication logic: prefer NC_ over NW_, then prefer longer sequences
            if org not in by_organism:
                by_organism[org] = record
            else:
                curr_is_nc = by_organism[org].id.startswith("NC_")
                new_is_nc = record.id.startswith("NC_")

                if new_is_nc and not curr_is_nc:
                    # New record has NC_ accession, current doesn't - prefer new
                    by_organism[org] = record
                elif curr_is_nc == new_is_nc and len(record.seq) > len(by_organism[org].seq):
                    # Same accession type, prefer longer sequence
                    by_organism[org] = record

        logger.info(f"Found {total_records} total records, {len(by_organism)} unique organisms (excluding human)")

        # Write FASTQ
        with open(output_path, "w") as out:
            for org, record in by_organism.items():
                # Add phred quality scores
                record.letter_annotations["phred_quality"] = [phred] * len(record.seq)

                # Extract gene_id from description if available
                if "[GeneID=" in record.description:
                    gene_id = record.description.split("[GeneID=")[1].split("]")[0]
                    record.id = f"gene_{gene_id}"

                record.description = ""
                SeqIO.write(record, out, "fastq")

        logger.success(f"FASTQ saved: {len(by_organism)} sequences from {len(by_organism)} organisms")
        return len(by_organism), len(by_organism)
