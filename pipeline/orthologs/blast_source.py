"""BLAST-based ortholog source implementation."""

import re
import time
from io import StringIO
from pathlib import Path

from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from loguru import logger
from tqdm import tqdm

from .base import OrthologResult, OrthologSource
from .registry import register


@register
class BlastOrthologSource(OrthologSource):
    """Retrieve orthologs via BLASTP search against RefSeq.

    This source:
    1. Downloads the longest protein isoform for the query gene
    2. Runs BLASTP against RefSeq protein database
    3. Maps protein hits to gene IDs
    4. Fetches genomic DNA for each orthologous gene
    5. Outputs sequences in FASTQ format
    """

    name = "blast"

    def is_available(self) -> bool:
        """BLAST source is always available (uses NCBI web service)."""
        return True

    def fetch(
        self,
        gene_id: int,
        output_dir: Path,
        hitlist_size: int = 100,
        phred: int = 30,
    ) -> OrthologResult:
        """Fetch orthologs via BLAST search.

        Args:
            gene_id: NCBI Gene ID
            output_dir: Directory to write output files
            hitlist_size: Maximum number of BLAST hits
            phred: Quality score for FASTQ output

        Returns:
            OrthologResult with path to FASTQ file
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        protein_path = output_dir / "protein_longest_isoform.faa"
        blast_path = output_dir / "blast_result.xml"
        fastq_path = output_dir / "genes_from_coordinates.fastq"

        # Step 1: Download protein
        self._download_longest_protein(gene_id, protein_path)

        # Step 2: BLAST search
        self._blast_search(protein_path, blast_path, hitlist_size)

        # Step 3: Parse BLAST and fetch sequences
        species_count, sequence_count = self._obtain_sequences(
            blast_path, fastq_path, phred
        )

        return OrthologResult(
            fastq_path=fastq_path,
            species_count=species_count,
            sequence_count=sequence_count,
            metadata={"hitlist_size": hitlist_size, "source": "blast"},
        )

    def _download_longest_protein(self, gene_id: int, output_path: Path) -> None:
        """Download the longest protein isoform for a gene with retry logic."""
        logger.info(f"Downloading longest protein isoform for gene {gene_id}")

        links = Entrez.read(
            Entrez.elink(
                dbfrom="gene",
                db="protein",
                id=str(gene_id),
                linkname="gene_protein_refseq",
            )
        )
        protein_ids = [link["Id"] for link in links[0]["LinkSetDb"][0]["Link"]]
        logger.info(f"Found {len(protein_ids)} proteins")

        max_retries = 3
        retry_delay = 5

        for attempt in range(1, max_retries + 1):
            try:
                logger.debug(f"Fetch attempt {attempt}/{max_retries}")

                handle = Entrez.efetch(
                    db="protein",
                    id=",".join(protein_ids),
                    rettype="fasta",
                    retmode="text",
                )
                records = list(SeqIO.parse(handle, "fasta"))
                handle.close()
                logger.info(f"Fetched {len(records)} sequences")

                longest = max(records, key=lambda r: len(r.seq))
                logger.info(f"Longest isoform: {longest.id} ({len(longest.seq)} aa)")

                SeqIO.write(longest, output_path, "fasta")
                logger.success(f"Saved to {output_path}")
                return

            except Exception as e:
                logger.warning(f"Fetch attempt {attempt} failed: {e}")

                if attempt < max_retries:
                    logger.info(f"Retrying in {retry_delay} seconds...")
                    time.sleep(retry_delay)
                    retry_delay *= 2
                else:
                    logger.error(f"Failed to download protein after {max_retries} attempts")
                    raise

    def _blast_search(
        self, protein_path: Path, output_path: Path, hitlist_size: int
    ) -> None:
        """Run BLASTP search against RefSeq with retry logic and progressive hitlist reduction."""
        logger.info(f"Running BLAST search (hitlist_size={hitlist_size})")

        query_protein = SeqIO.read(protein_path, "fasta")
        logger.debug(f"Query protein length: {len(query_protein.seq)} aa")

        result_handle = NCBIWWW.qblast(
            program="blastp",
            database="refseq_protein",
            sequence=str(query_protein.seq),
            expect=10,
            hitlist_size=hitlist_size,
        )

        result = result_handle.read()
        logger.debug(f"BLAST result type: {type(result)}")
        # logger.debug(f'result = {result}')
        # logger.debug(f'Save result to {output_path}')
        if isinstance(result, bytes):
            with open(output_path, "wb") as out:
                out.write(result)
        else:
            with open(output_path, "w") as out:#, encoding="utf-8") as out:
                out.write(result)

        logger.success(f"BLAST results saved to {output_path}")
        return
    
    def _obtain_sequences(
        self, blast_path: Path, output_path: Path, phred: int
    ) -> tuple[int, int]:
        """Parse BLAST results and fetch genomic sequences."""
        logger.info("Parsing BLAST results for homologous sequences")

        with open(blast_path) as fh:
            blast_records = NCBIXML.parse(fh)
            best_by_species = self._get_best_by_species(blast_records)

        protein_ids = [aln.accession for aln, _ in best_by_species.values()]
        protein_ids = list(dict.fromkeys(protein_ids))

        logger.info(
            f"Found {len(best_by_species)} species, {len(protein_ids)} unique protein IDs"
        )

        # Map proteins to genes
        gene_ids = self._batch_protein_to_gene_id(protein_ids)
        logger.info(f"Mapped to {len(gene_ids)} gene IDs")

        # Get gene coordinates
        gene_coords = self._batch_gene_coordinates(gene_ids)

        # Fetch sequences and write FASTQ
        success_count, skip_count = self._fetch_and_write_sequences(
            gene_ids, gene_coords, output_path, phred
        )

        logger.success(f"FASTQ saved: {success_count} sequences, {skip_count} skipped")
        return len(best_by_species), success_count

    def _get_best_by_species(self, blast_records) -> dict:
        """Get best BLAST hit per species (excluding human)."""
        best_by_species = {}

        for record in blast_records:
            for alignment in record.alignments:
                title = alignment.title.lower()

                if "homo sapiens" in title:
                    continue

                m = re.search(r"\[([^\]]+)\]", alignment.title)
                species = m.group(1) if m else "UNKNOWN"

                best_hsp = alignment.hsps[0]
                bits = best_hsp.bits

                if species not in best_by_species:
                    best_by_species[species] = (alignment, bits)
                else:
                    _, prev_bits = best_by_species[species]
                    if bits > prev_bits:
                        best_by_species[species] = (alignment, bits)

        return best_by_species

    def _batch_protein_to_gene_id(
        self, protein_ids: list[str], batch_size: int = 200
    ) -> list[str]:
        """Map protein IDs to gene IDs in batches."""
        gene_ids = set()

        for i in range(0, len(protein_ids), batch_size):
            batch = protein_ids[i : i + batch_size]
            logger.debug(f"Processing batch {i // batch_size + 1}: {len(batch)} proteins")

            with Entrez.elink(
                dbfrom="protein",
                db="gene",
                id=",".join(batch),
                linkname="protein_gene",
            ) as h:
                links = Entrez.read(h)

            for linkset in links:
                if "LinkSetDb" not in linkset:
                    continue
                for db in linkset["LinkSetDb"]:
                    for link in db["Link"]:
                        gene_ids.add(link["Id"])

        return list(gene_ids)

    def _batch_gene_coordinates(self, gene_ids: list[str]) -> dict:
        """Get genomic coordinates for genes."""
        handle = Entrez.esummary(
            db="gene",
            id=",".join(map(str, gene_ids)),
            retmode="xml",
        )
        records = Entrez.read(handle)
        handle.close()

        result = {}

        for doc in records["DocumentSummarySet"]["DocumentSummary"]:
            gid = doc.attributes.get("uid")
            gi_list = doc.get("GenomicInfo", [])

            if not gi_list:
                result[gid] = None
                continue

            coords = gi_list[0]

            try:
                result[gid] = {
                    "chr_acc": coords.get("ChrAccVer"),
                    "start": int(coords["ChrStart"]),
                    "end": int(coords["ChrStop"]),
                }
            except Exception:
                result[gid] = None

        return result

    def _fetch_and_write_sequences(
        self,
        gene_ids: list[str],
        gene_coords: dict,
        output_path: Path,
        phred: int,
    ) -> tuple[int, int]:
        """Fetch genomic sequences and write to FASTQ with retry logic."""
        success_count = 0
        skip_count = 0

        with open(output_path, "w") as out:
            for gid in tqdm(gene_ids, desc="Fetching genes"):
                coords = gene_coords.get(str(gid)) or gene_coords.get(gid)
                if coords is None:
                    skip_count += 1
                    continue

                max_retries = 3
                retry_delay = 2

                for attempt in range(1, max_retries + 1):
                    try:
                        handle = Entrez.efetch(
                            db="nuccore",
                            id=coords["chr_acc"],
                            rettype="fasta",
                            retmode="text",
                            seq_start=min(coords["start"], coords["end"]) + 1,
                            seq_stop=max(coords["start"], coords["end"]) + 1,
                        )
                        fasta_text = handle.read()
                        handle.close()

                        for record in SeqIO.parse(StringIO(fasta_text), "fasta"):
                            seq_len = len(record.seq)
                            if seq_len == 0:
                                continue

                            record.letter_annotations["phred_quality"] = [phred] * seq_len
                            record.id = f"gene_{gid}"
                            record.description = ""

                            SeqIO.write(record, out, "fastq")
                            success_count += 1

                        break  # Success, exit retry loop

                    except Exception as e:
                        if attempt < max_retries:
                            logger.debug(f"Gene {gid} fetch failed (attempt {attempt}), retrying: {e}")
                            time.sleep(retry_delay)
                        else:
                            logger.warning(f"Failed to fetch gene {gid} after {max_retries} attempts: {e}")
                            skip_count += 1

        return success_count, skip_count
