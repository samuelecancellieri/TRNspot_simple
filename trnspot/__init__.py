"""
TRNspot - A package for transcriptional regulatory network analysis
===================================================================

TRNspot is a modular Python package for transcriptional regulatory network (TRN)
analysis using single-cell data. It integrates Scanpy, CellOracle, and Hotspot
into a checkpoint-enabled pipeline.

Modules
-------
- utils: Data handling, I/O, and logging utilities
- execution: Pipeline execution functions for preprocessing, CellOracle, Hotspot
- plotting: Visualization functions for QC, embeddings, GRN, and Hotspot
- reporting: HTML and PDF report generation

Legacy Modules (for backward compatibility)
-------------------------------------------
- preprocessing: Scanpy wrappers for QC and normalization
- celloracle_processing: GRN inference with CellOracle
- hotspot_processing: Gene module analysis with Hotspot
- config: Central configuration and parameters
"""

__version__ = "0.1.0"
__author__ = "Samuele Cancellieri"

# Import configuration (always available)
from . import config
from .config import set_random_seed, set_scanpy_settings, get_config, print_config

# Import legacy modules (backward compatibility)
from . import preprocessing
from . import celloracle_processing
from . import hotspot_processing

# Import new modular structure
from . import utils
from . import execution
from . import plotting
from . import reporting

# Convenience re-exports from utils
from .utils import (
    ensure_categorical_obs,
    load_adata,
    save_adata,
    subset_adata_by_cluster,
    stratify_adata,
    setup_directories,
    log_step,
    log_error,
)

# Convenience re-exports from execution
from .execution import (
    run_qc,
    run_normalization,
    run_dimensionality_reduction,
    run_clustering,
    run_full_preprocessing,
    run_celloracle_pipeline,
    run_hotspot_pipeline,
)

# Convenience re-exports from reporting
from .reporting import (
    ReportGenerator,
    generate_report,
    generate_html_report,
    generate_pdf_report,
)

__all__ = [
    # Version info
    "__version__",
    "__author__",
    # Configuration
    "config",
    "set_random_seed",
    "set_scanpy_settings",
    "get_config",
    "print_config",
    # New modular structure
    "utils",
    "execution",
    "plotting",
    "reporting",
    # Legacy modules
    "preprocessing",
    "celloracle_processing",
    "hotspot_processing",
    # Utility functions
    "ensure_categorical_obs",
    "load_adata",
    "save_adata",
    "subset_adata_by_cluster",
    "stratify_adata",
    "setup_directories",
    "log_step",
    "log_error",
    # Execution functions
    "run_qc",
    "run_normalization",
    "run_dimensionality_reduction",
    "run_clustering",
    "run_full_preprocessing",
    "run_celloracle_pipeline",
    "run_hotspot_pipeline",
    # Reporting functions
    "ReportGenerator",
    "generate_report",
    "generate_html_report",
    "generate_pdf_report",
]
