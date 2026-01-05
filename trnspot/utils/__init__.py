"""
TRNspot utility functions
=========================

This module contains utility functions for data handling, file I/O,
and common operations used across the package.
"""

from .data import (
    ensure_categorical_obs,
    load_adata,
    save_adata,
    subset_adata_by_cluster,
)
from .io import (
    setup_directories,
    track_files,
    compute_input_hash,
    write_checkpoint,
    check_checkpoint,
)
from .logging import (
    setup_logging,
    log_step,
    log_error,
    get_pipeline_logger,
    get_error_logger,
)

__all__ = [
    # Data utilities
    "ensure_categorical_obs",
    "load_adata",
    "save_adata",
    "subset_adata_by_cluster",
    # I/O utilities
    "setup_directories",
    "track_files",
    "compute_input_hash",
    "write_checkpoint",
    "check_checkpoint",
    # Logging utilities
    "setup_logging",
    "log_step",
    "log_error",
    "get_pipeline_logger",
    "get_error_logger",
]
