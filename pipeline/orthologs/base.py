"""Base classes for ortholog source implementations."""

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class OrthologResult:
    """Standard output format for all ortholog sources."""

    fastq_path: Path
    species_count: int
    sequence_count: int
    metadata: dict = field(default_factory=dict)


class OrthologSource(ABC):
    """Base class for ortholog retrieval methods.

    Subclasses implement different strategies for obtaining orthologous
    sequences (BLAST, NCBI Datasets CLI, etc.).
    """

    name: str  # Human-readable name, set by subclasses

    @abstractmethod
    def fetch(self, gene_id: int, output_dir: Path, **kwargs) -> OrthologResult:
        """Fetch orthologous sequences for a gene.

        Args:
            gene_id: NCBI Gene ID
            output_dir: Directory to write output files
            **kwargs: Source-specific parameters

        Returns:
            OrthologResult with path to FASTQ and metadata
        """
        pass

    @abstractmethod
    def is_available(self) -> bool:
        """Check if this source's dependencies are available.

        Returns:
            True if the source can be used, False otherwise
        """
        pass
