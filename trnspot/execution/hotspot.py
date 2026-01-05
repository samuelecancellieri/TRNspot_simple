"""
Hotspot execution functions for TRNspot
=======================================

High-level functions for running Hotspot gene module analysis.
"""

from typing import Optional, Any
from anndata import AnnData

from .. import config
from ..utils.logging import log_step, log_error


def run_hotspot_analysis(
    adata: AnnData,
    top_genes: Optional[int] = None,
    layer_key: str = "raw_counts",
    embedding_key: str = "X_pca",
    model: str = "danb",
) -> Any:
    """
    Run Hotspot analysis for gene module identification.

    Parameters
    ----------
    adata : AnnData
        Preprocessed AnnData
    top_genes : int, optional
        Number of top genes to analyze (default from config)
    layer_key : str
        Layer with raw counts
    embedding_key : str
        Embedding for spatial analysis
    model : str
        Statistical model ('danb' = depth-adjusted negative binomial)

    Returns
    -------
    Hotspot
        Hotspot object with analysis results
    """
    from ..hotspot_processing import create_hotspot_object, run_hotspot_analysis as _run

    log_step("Hotspot.Analysis", "STARTED", {"embedding_key": embedding_key})

    try:
        top_genes = top_genes or config.HOTSPOT_TOP_GENES

        # Create Hotspot object
        log_step("Hotspot.CreateObject", "STARTED")
        hs_obj = create_hotspot_object(
            adata,
            top_genes=top_genes,
            layer_key=layer_key,
            model=model,
            embedding_key=embedding_key,
        )
        log_step("Hotspot.CreateObject", "COMPLETED")

        # Run analysis
        log_step("Hotspot.Compute", "STARTED")
        hs_obj = _run(hs_obj)
        log_step("Hotspot.Compute", "COMPLETED")

        log_step("Hotspot.Analysis", "COMPLETED")
        return hs_obj

    except Exception as e:
        log_error("Hotspot.Analysis", e)
        raise


def run_hotspot_pipeline(
    adata: AnnData,
    layer_key: str = "raw_counts",
    embedding_key: str = "X_pca",
    skip: bool = False,
) -> Optional[Any]:
    """
    Run complete Hotspot pipeline.

    Parameters
    ----------
    adata : AnnData
        Input AnnData
    layer_key : str
        Layer with raw counts
    embedding_key : str
        Embedding for analysis
    skip : bool
        Skip Hotspot analysis

    Returns
    -------
    Hotspot or None
        Hotspot object or None if skipped
    """
    if skip:
        log_step("Hotspot", "SKIPPED", {"reason": "skip=True"})
        print("⊘ Skipping Hotspot analysis")
        return None

    log_step("Hotspot.Pipeline", "STARTED")

    try:
        hs_obj = run_hotspot_analysis(
            adata,
            layer_key=layer_key,
            embedding_key=embedding_key,
        )

        log_step("Hotspot.Pipeline", "COMPLETED")
        return hs_obj

    except ImportError:
        log_step("Hotspot", "SKIPPED", {"reason": "Hotspot not installed"})
        print("⊘ Skipping Hotspot (not installed)")
        return None

    except Exception as e:
        log_error("Hotspot.Pipeline", e)
        raise
