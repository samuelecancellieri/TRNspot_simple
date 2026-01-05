"""
TRNspot execution module
========================

This module contains pipeline execution functions organized by analysis type.
"""

from .preprocessing import (
    run_qc,
    run_normalization,
    run_dimensionality_reduction,
    run_clustering,
    run_full_preprocessing,
)
from .celloracle import (
    run_celloracle_preprocessing,
    run_celloracle_inference,
    run_celloracle_pipeline,
)
from .hotspot import (
    run_hotspot_analysis,
    run_hotspot_pipeline,
)

__all__ = [
    # Preprocessing
    "run_qc",
    "run_normalization",
    "run_dimensionality_reduction",
    "run_clustering",
    "run_full_preprocessing",
    # CellOracle
    "run_celloracle_preprocessing",
    "run_celloracle_inference",
    "run_celloracle_pipeline",
    # Hotspot
    "run_hotspot_analysis",
    "run_hotspot_pipeline",
]
