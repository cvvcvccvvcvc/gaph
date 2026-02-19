"""Registry for ortholog source implementations."""

from typing import Type

from loguru import logger

from .base import OrthologSource

_SOURCES: dict[str, Type[OrthologSource]] = {}


def register(cls: Type[OrthologSource]) -> Type[OrthologSource]:
    """Decorator to register an ortholog source.

    Usage:
        @register
        class MyOrthologSource(OrthologSource):
            name = "my_source"
            ...
    """
    _SOURCES[cls.name] = cls
    return cls


def get_source(name: str) -> OrthologSource:
    """Get an ortholog source instance by name.

    Args:
        name: Source name (e.g., "blast", "ncbi_datasets")

    Returns:
        Instantiated OrthologSource

    Raises:
        KeyError: If source name is not registered
    """
    if name not in _SOURCES:
        available = ", ".join(_SOURCES.keys()) or "(none)"
        raise KeyError(f"Unknown ortholog source: {name!r}. Available: {available}")
    return _SOURCES[name]()


def list_sources() -> list[str]:
    """List all registered ortholog source names."""
    return list(_SOURCES.keys())


def get_default_source(
    primary: str = "ncbi_datasets",
    fallback: str = "blast",
) -> OrthologSource:
    """Get ortholog source with automatic fallback.

    Tries primary source first (ncbi_datasets by default).
    Falls back to secondary if primary is unavailable.

    Args:
        primary: Primary source name (default: "ncbi_datasets")
        fallback: Fallback source name (default: "blast")

    Returns:
        Instantiated OrthologSource
    """
    try:
        source = get_source(primary)
        if source.is_available():
            logger.info(f"Using {primary} ortholog source")
            return source
        logger.warning(f"{primary} unavailable, falling back to {fallback}")
    except KeyError:
        logger.warning(f"{primary} not found, falling back to {fallback}")

    return get_source(fallback)
