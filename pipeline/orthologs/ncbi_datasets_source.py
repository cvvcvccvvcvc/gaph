"""NCBI Datasets CLI-based ortholog source implementation."""

import gzip
import json
import re
import shutil
import subprocess
import tempfile
import zipfile
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from loguru import logger

from .base import OrthologResult, OrthologSource
from .registry import register

SOURCE_FASTQ_PHRED = 30


@register
class NCBIDatasetsOrthologSource(OrthologSource):
    """Retrieve orthologs via NCBI Datasets CLI.

    This source:
    1. Runs `datasets download gene gene-id {id} --ortholog <scope> --include gene`
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

    @staticmethod
    def _normalized_ortholog_scope(scope: str | None) -> str:
        if scope is None:
            return "all"
        normalized = str(scope).strip().lower()
        return normalized if normalized else "all"

    @staticmethod
    def _scope_slug(scope: str | None) -> str:
        normalized = NCBIDatasetsOrthologSource._normalized_ortholog_scope(scope)
        if normalized == "all":
            return "all"
        slug = re.sub(r"[^a-z0-9]+", "_", normalized).strip("_")
        return slug if slug else "all"

    @classmethod
    def _cache_fastq_path(cls, cache_dir: Path, gene_id: int, ortholog_scope: str = "all") -> Path:
        scope_slug = cls._scope_slug(ortholog_scope)
        suffix = "" if scope_slug == "all" else f"__scope_{scope_slug}"
        return Path(cache_dir) / f"gene_{gene_id}{suffix}.fastq.gz"

    @classmethod
    def _cache_meta_path(cls, cache_dir: Path, gene_id: int, ortholog_scope: str = "all") -> Path:
        scope_slug = cls._scope_slug(ortholog_scope)
        suffix = "" if scope_slug == "all" else f"__scope_{scope_slug}"
        return Path(cache_dir) / f"gene_{gene_id}{suffix}.meta.json"

    def has_cached(self, gene_id: int, cache_dir: Path, ortholog_scope: str = "all") -> bool:
        cache_fastq = self._cache_fastq_path(cache_dir, gene_id, ortholog_scope=ortholog_scope)
        return cache_fastq.exists() and cache_fastq.stat().st_size > 0

    def load_cached(
        self,
        gene_id: int,
        cache_dir: Path,
        output_fastq: Path,
        ortholog_scope: str = "all",
    ) -> OrthologResult | None:
        """Restore cached ortholog FASTQ to output_fastq if present."""
        ortholog_scope = self._normalized_ortholog_scope(ortholog_scope)
        cache_fastq = self._cache_fastq_path(cache_dir, gene_id, ortholog_scope=ortholog_scope)
        if not cache_fastq.exists():
            return None

        output_fastq.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(cache_fastq, "rb") as src, open(output_fastq, "wb") as dst:
            shutil.copyfileobj(src, dst)

        meta = self._read_cache_meta(gene_id, cache_dir, ortholog_scope=ortholog_scope)
        if meta:
            species_count = int(meta.get("species_count", 0))
            sequence_count = int(meta.get("sequence_count", 0))
        else:
            sequence_count = self._count_fastq_records(output_fastq)
            species_count = sequence_count

        logger.info(f"Using cached NCBI orthologs for gene {gene_id}: {cache_fastq}")
        return OrthologResult(
            fastq_path=output_fastq,
            species_count=species_count,
            sequence_count=sequence_count,
            metadata={
                "source": "ncbi_datasets_cache",
                "gene_id": gene_id,
                "ortholog_scope": ortholog_scope,
            },
        )

    def save_to_cache(
        self,
        gene_id: int,
        cache_dir: Path,
        source_fastq: Path,
        species_count: int,
        sequence_count: int,
        source_tag: str = "ncbi_datasets",
        ortholog_scope: str = "all",
    ) -> Path:
        """Save an ortholog FASTQ to persistent cache as gzip."""
        ortholog_scope = self._normalized_ortholog_scope(ortholog_scope)
        cache_dir = Path(cache_dir)
        cache_dir.mkdir(parents=True, exist_ok=True)

        cache_fastq = self._cache_fastq_path(cache_dir, gene_id, ortholog_scope=ortholog_scope)
        with open(source_fastq, "rb") as src, gzip.open(cache_fastq, "wb", compresslevel=1) as dst:
            shutil.copyfileobj(src, dst)

        self._write_cache_meta(
            gene_id=gene_id,
            cache_dir=cache_dir,
            source=source_tag,
            species_count=species_count,
            sequence_count=sequence_count,
            ortholog_scope=ortholog_scope,
        )
        return cache_fastq

    def _write_cache_meta(
        self,
        gene_id: int,
        cache_dir: Path,
        source: str,
        species_count: int,
        sequence_count: int,
        ortholog_scope: str = "all",
    ) -> None:
        ortholog_scope = self._normalized_ortholog_scope(ortholog_scope)
        meta_path = self._cache_meta_path(cache_dir, gene_id, ortholog_scope=ortholog_scope)
        payload = {
            "gene_id": int(gene_id),
            "source": source,
            "ortholog_scope": ortholog_scope,
            "species_count": int(species_count),
            "sequence_count": int(sequence_count),
        }
        with open(meta_path, "w") as f:
            json.dump(payload, f, indent=2)

    def _read_cache_meta(self, gene_id: int, cache_dir: Path, ortholog_scope: str = "all") -> dict | None:
        ortholog_scope = self._normalized_ortholog_scope(ortholog_scope)
        meta_path = self._cache_meta_path(cache_dir, gene_id, ortholog_scope=ortholog_scope)
        if not meta_path.exists():
            return None
        try:
            with open(meta_path) as f:
                return json.load(f)
        except Exception:
            return None

    @staticmethod
    def _count_fastq_records(fastq_path: Path) -> int:
        lines = 0
        with open(fastq_path) as fh:
            for _ in fh:
                lines += 1
        return lines // 4

    def fetch(
        self,
        gene_id: int,
        output_dir: Path,
        ortholog_scope: str = "all",
    ) -> OrthologResult:
        """Fetch orthologs via NCBI Datasets CLI.

        Args:
            gene_id: NCBI Gene ID
            output_dir: Directory to write output files

        Returns:
            OrthologResult with path to FASTQ file
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        ortholog_scope = self._normalized_ortholog_scope(ortholog_scope)

        zip_path = output_dir / "ncbi_dataset.zip"
        data_dir = output_dir / "ncbi_dataset" / "data"
        report_path = data_dir / "data_report.jsonl"
        fastq_path = output_dir / "genes_from_coordinates.fastq"

        # Step 1: Download ortholog data
        self._download_orthologs(gene_id, output_dir, zip_path, ortholog_scope=ortholog_scope)

        # Step 2: Extract zip
        self._extract_zip(zip_path, output_dir)

        # Step 3: Read gene sequences from gene.fna and convert to FASTQ
        species_count, sequence_count = self._fetch_sequences(data_dir, fastq_path)

        return OrthologResult(
            fastq_path=fastq_path,
            species_count=species_count,
            sequence_count=sequence_count,
            metadata={
                "source": "ncbi_datasets",
                "gene_id": gene_id,
                "ortholog_scope": ortholog_scope,
            },
        )

    def prefetch_to_cache(
        self,
        gene_ids: list[int],
        cache_dir: Path,
        ortholog_scope: str = "all",
    ) -> set[int]:
        """Batch-download orthologs for many query genes and write per-gene cache files.

        Only writes processed per-gene FASTQ cache files (gzip), no persistent raw batch zip.
        Returns gene IDs for which cached FASTQ was produced.
        """
        ortholog_scope = self._normalized_ortholog_scope(ortholog_scope)
        cache_dir = Path(cache_dir)
        cache_dir.mkdir(parents=True, exist_ok=True)

        dedup_gene_ids: list[int] = []
        seen: set[int] = set()
        for gid in gene_ids:
            gid_int = int(gid)
            if gid_int in seen:
                continue
            seen.add(gid_int)
            dedup_gene_ids.append(gid_int)

        missing = [gid for gid in dedup_gene_ids if not self.has_cached(gid, cache_dir, ortholog_scope=ortholog_scope)]
        if not missing:
            return set()

        requested = {str(gid) for gid in missing}
        logger.info(
            "Batch downloading NCBI orthologs for {} gene(s) (ortholog_scope={})",
            len(missing),
            ortholog_scope,
        )

        with tempfile.TemporaryDirectory(prefix="ncbi_prefetch_", dir=str(cache_dir)) as tmp:
            tmp_dir = Path(tmp)
            input_ids = tmp_dir / "gene_ids.txt"
            batch_zip = tmp_dir / "ncbi_dataset_batch.zip"
            input_ids.write_text("\n".join(str(gid) for gid in missing) + "\n")

            self._download_orthologs_batch(
                input_ids,
                tmp_dir,
                batch_zip,
                ortholog_scope=ortholog_scope,
            )
            self._extract_zip(batch_zip, tmp_dir)

            data_dir = tmp_dir / "ncbi_dataset" / "data"
            grouped = self._collect_records_by_query_gene(data_dir, requested)

            produced: set[int] = set()
            for gid in missing:
                by_organism = grouped.get(str(gid), {})
                if not by_organism:
                    continue

                tmp_fastq = tmp_dir / f"gene_{gid}.fastq"
                species_count, sequence_count = self._write_records_fastq(by_organism, tmp_fastq)
                self.save_to_cache(
                    gene_id=gid,
                    cache_dir=cache_dir,
                    source_fastq=tmp_fastq,
                    species_count=species_count,
                    sequence_count=sequence_count,
                    source_tag=f"ncbi_datasets_batch_{ortholog_scope}",
                    ortholog_scope=ortholog_scope,
                )
                produced.add(gid)

        logger.info(
            "NCBI batch prefetch finished: {} cached, {} missing",
            len(produced),
            len(missing) - len(produced),
        )
        return produced

    def _download_orthologs_batch(
        self,
        gene_ids_file: Path,
        output_dir: Path,
        zip_path: Path,
        ortholog_scope: str = "all",
    ) -> None:
        """Download batch ortholog data using datasets CLI input file."""
        ortholog_scope = self._normalized_ortholog_scope(ortholog_scope)
        datasets_cmd = self._get_datasets_cmd()
        cmd = [
            datasets_cmd,
            "download",
            "gene",
            "gene-id",
            "--inputfile",
            str(gene_ids_file),
            "--ortholog",
            ortholog_scope,
            "--include",
            "gene",
            "--filename",
            str(zip_path),
        ]

        logger.debug(f"Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=output_dir)
        if result.returncode != 0:
            logger.error(f"datasets CLI batch failed: {result.stderr}")
            raise RuntimeError(f"datasets batch download failed: {result.stderr}")

    def _collect_records_by_query_gene(
        self,
        data_dir: Path,
        requested_query_gene_ids: set[str],
    ) -> dict[str, dict[str, SeqRecord]]:
        """Collect and deduplicate FASTA records grouped by query-gene ID."""
        report_path = data_dir / "data_report.jsonl"
        gene_fna = data_dir / "gene.fna"
        if not report_path.exists() or not gene_fna.exists():
            raise FileNotFoundError(f"Expected files not found in {data_dir}")

        ortholog_to_query = self._build_ortholog_to_query_map(report_path, requested_query_gene_ids)
        grouped: dict[str, dict[str, SeqRecord]] = {gid: {} for gid in requested_query_gene_ids}

        for record in SeqIO.parse(gene_fna, "fasta"):
            ortholog_gene_id = self._extract_gene_id(record.description)
            if not ortholog_gene_id:
                continue

            query_gene_id = ortholog_to_query.get(ortholog_gene_id)
            if not query_gene_id:
                continue

            organism = self._extract_organism(record.description)
            if not organism or organism.lower() == "homo sapiens":
                continue

            by_org = grouped.setdefault(query_gene_id, {})
            if organism not in by_org:
                by_org[organism] = record
                continue

            current = by_org[organism]
            curr_is_nc = current.id.startswith("NC_")
            new_is_nc = record.id.startswith("NC_")
            if new_is_nc and not curr_is_nc:
                by_org[organism] = record
            elif curr_is_nc == new_is_nc and len(record.seq) > len(current.seq):
                by_org[organism] = record

        return grouped

    @staticmethod
    def _build_ortholog_to_query_map(report_path: Path, requested_query_gene_ids: set[str]) -> dict[str, str]:
        """Map ortholog GeneID -> requested query GeneID using data_report.jsonl."""
        mapping: dict[str, str] = {}
        with open(report_path) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    row = json.loads(line)
                except json.JSONDecodeError:
                    continue

                ortholog_gene_id = str(row.get("geneId", "")).strip()
                if not ortholog_gene_id:
                    continue

                if ortholog_gene_id in requested_query_gene_ids:
                    mapping[ortholog_gene_id] = ortholog_gene_id
                    continue

                groups = row.get("geneGroups") or []
                query_gene_id = None
                for group in groups:
                    if not isinstance(group, dict):
                        continue
                    group_id = str(group.get("id", "")).strip()
                    if group_id in requested_query_gene_ids:
                        query_gene_id = group_id
                        break

                if query_gene_id:
                    mapping[ortholog_gene_id] = query_gene_id

        return mapping

    @staticmethod
    def _extract_organism(description: str) -> str | None:
        if "[organism=" not in description:
            return None
        return description.split("[organism=")[1].split("]")[0]

    @staticmethod
    def _extract_gene_id(description: str) -> str | None:
        if "[GeneID=" not in description:
            return None
        return description.split("[GeneID=")[1].split("]")[0]

    def _write_records_fastq(
        self,
        by_organism: dict[str, SeqRecord],
        output_fastq: Path,
    ) -> tuple[int, int]:
        """Write deduplicated records to FASTQ (one record per organism)."""
        with open(output_fastq, "w") as out:
            for record in by_organism.values():
                rec = record[:]
                gene_id = self._extract_gene_id(record.description)
                if gene_id:
                    rec.id = f"gene_{gene_id}"
                rec.description = ""
                rec.letter_annotations["phred_quality"] = [SOURCE_FASTQ_PHRED] * len(rec.seq)
                SeqIO.write(rec, out, "fastq")

        count = len(by_organism)
        return count, count

    def _download_orthologs(
        self,
        gene_id: int,
        output_dir: Path,
        zip_path: Path,
        ortholog_scope: str = "all",
    ) -> None:
        """Download ortholog data using datasets CLI."""
        ortholog_scope = self._normalized_ortholog_scope(ortholog_scope)
        logger.info(
            "Downloading orthologs for gene {} via NCBI Datasets CLI (ortholog_scope={})",
            gene_id,
            ortholog_scope,
        )

        datasets_cmd = self._get_datasets_cmd()
        cmd = [
            datasets_cmd,
            "download",
            "gene",
            "gene-id",
            str(gene_id),
            "--ortholog",
            ortholog_scope,
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

    def _fetch_sequences(self, data_dir: Path, output_path: Path) -> tuple[int, int]:
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
                record.letter_annotations["phred_quality"] = [SOURCE_FASTQ_PHRED] * len(record.seq)

                # Extract gene_id from description if available
                if "[GeneID=" in record.description:
                    gene_id = record.description.split("[GeneID=")[1].split("]")[0]
                    record.id = f"gene_{gene_id}"

                record.description = ""
                SeqIO.write(record, out, "fastq")

        logger.success(f"FASTQ saved: {len(by_organism)} sequences from {len(by_organism)} organisms")
        return len(by_organism), len(by_organism)
