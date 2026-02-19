"""Modular ortholog source implementations.

This package provides multiple strategies for obtaining orthologous sequences:

- ncbi_datasets: NCBI Datasets CLI (default, requires `datasets` CLI installed)
- blast: BLAST-based search against RefSeq (fallback, always available)

Usage:
    from orthologs import get_default_source, list_sources

    # List available sources
    print(list_sources())  # ['blast', 'ncbi_datasets']

    # Use default source (ncbi_datasets with blast fallback)
    source = get_default_source()
    result = source.fetch(gene_id=672, output_dir=Path("gene_672"))

    # Or use a specific source
    source = get_source("blast")
    result = source.fetch(gene_id=672, output_dir=Path("gene_672"), hitlist_size=5000)

    # Access results
    print(result.fastq_path)      # Path to FASTQ file
    print(result.species_count)   # Number of species
    print(result.sequence_count)  # Number of sequences
"""

from .base import OrthologResult, OrthologSource
from .registry import get_default_source, get_source, list_sources, register

# Import sources to trigger registration
from . import blast_source  # noqa: F401
from . import ncbi_datasets_source  # noqa: F401

__all__ = [
    "OrthologResult",
    "OrthologSource",
    "get_default_source",
    "get_source",
    "list_sources",
    "register",
]
