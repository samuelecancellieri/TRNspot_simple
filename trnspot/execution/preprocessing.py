"""
Preprocessing execution functions for TRNspot
=============================================

High-level functions for running preprocessing steps.
"""

from typing import Optional, Tuple
import scanpy as sc
from anndata import AnnData

from .. import config
from ..utils.data import ensure_categorical_obs
from ..utils.logging import log_step, log_error


def run_qc(
    adata: AnnData,
    min_genes: Optional[int] = None,
    min_counts: Optional[int] = None,
    max_counts: Optional[int] = None,
    pct_counts_mt_max: Optional[float] = None,
    min_cells: Optional[int] = None,
    save_plots: Optional[str] = None,
) -> AnnData:
    """
    Run quality control on AnnData object.

    This is a wrapper around the preprocessing.perform_qc function
    with logging integration.

    Parameters
    ----------
    adata : AnnData
        Input data
    min_genes, min_counts, max_counts, pct_counts_mt_max, min_cells : optional
        QC thresholds (use config defaults if None)
    save_plots : str, optional
        Base name for saving QC plots

    Returns
    -------
    AnnData
        Filtered AnnData object
    """
    from ..preprocessing import perform_qc

    log_step("QC", "STARTED", {"n_obs": adata.n_obs, "n_vars": adata.n_vars})

    try:
        adata_qc = perform_qc(
            adata,
            min_genes=min_genes,
            min_counts=min_counts,
            max_counts=max_counts,
            pct_counts_mt_max=pct_counts_mt_max,
            min_cells=min_cells,
            save_plots=save_plots,
        )
        log_step(
            "QC", "COMPLETED", {"n_obs": adata_qc.n_obs, "n_vars": adata_qc.n_vars}
        )
        return adata_qc

    except Exception as e:
        log_error("QC", e)
        raise


def run_normalization(
    adata: AnnData,
) -> AnnData:
    """
    Run normalization on AnnData object.

    Parameters
    ----------
    adata : AnnData
        Input data (should be QC'd)

    Returns
    -------
    AnnData
        Normalized AnnData object
    """
    from ..preprocessing import perform_normalization

    log_step("Normalization", "STARTED", {"n_obs": adata.n_obs})

    try:
        adata_norm = perform_normalization(adata)
        log_step("Normalization", "COMPLETED")
        return adata_norm

    except Exception as e:
        log_error("Normalization", e)
        raise


def run_dimensionality_reduction(
    adata: AnnData,
    n_top_genes: Optional[int] = None,
    n_pcs: Optional[int] = None,
    n_neighbors: Optional[int] = None,
) -> AnnData:
    """
    Run dimensionality reduction (HVG selection, PCA, UMAP).

    Parameters
    ----------
    adata : AnnData
        Normalized data
    n_top_genes : int, optional
        Number of highly variable genes
    n_pcs : int, optional
        Number of principal components
    n_neighbors : int, optional
        Number of neighbors for UMAP

    Returns
    -------
    AnnData
        AnnData with embeddings
    """
    log_step("DimensionalityReduction", "STARTED")

    try:
        adata_cc = adata.copy()

        # HVG selection
        n_top_genes = n_top_genes or config.HVGS_N_TOP_GENES
        sc.pp.highly_variable_genes(adata_cc, n_top_genes=n_top_genes, subset=False)
        print(f"Identified top {n_top_genes} highly variable genes")

        # PCA
        n_pcs = n_pcs or config.PCA_N_COMPS
        sc.pp.pca(adata_cc, n_comps=n_pcs, svd_solver=config.PCA_SVDSOLVE)
        print(f"Computed PCA with {n_pcs} components")

        # Neighbors
        n_neighbors = n_neighbors or config.NEIGHBORS_N_NEIGHBORS
        sc.pp.neighbors(
            adata_cc,
            metric=config.NEIGHBORS_METRIC,
            method=config.NEIGHBORS_METHOD,
            n_neighbors=n_neighbors,
            n_pcs=config.NEIGHBORS_N_PCS,
        )
        print(f"Constructed neighborhood graph with {n_neighbors} neighbors")

        # UMAP
        sc.tl.umap(adata_cc)
        print("Computed UMAP embedding")

        log_step("DimensionalityReduction", "COMPLETED")
        return adata_cc

    except Exception as e:
        log_error("DimensionalityReduction", e)
        raise


def run_clustering(
    adata: AnnData,
    resolution: Optional[float] = None,
    algorithm: str = "leiden",
) -> AnnData:
    """
    Run clustering on AnnData object.

    Parameters
    ----------
    adata : AnnData
        AnnData with computed neighbors
    resolution : float, optional
        Clustering resolution
    algorithm : str
        Clustering algorithm ('leiden' or 'louvain')

    Returns
    -------
    AnnData
        AnnData with cluster labels
    """
    log_step("Clustering", "STARTED", {"algorithm": algorithm})

    try:
        adata_cc = adata.copy()

        if algorithm == "leiden":
            resolution = resolution or config.LEIDEN_RESOLUTION
            sc.tl.leiden(adata_cc, resolution=resolution)
        elif algorithm == "louvain":
            resolution = resolution or config.LOUVAIN_RESOLUTION
            sc.tl.louvain(adata_cc, resolution=resolution)
        else:
            raise ValueError(f"Unknown clustering algorithm: {algorithm}")

        print(f"Performed {algorithm} clustering with resolution {resolution}")

        # Ensure categorical
        adata_cc = ensure_categorical_obs(adata_cc, columns=[algorithm])

        n_clusters = len(adata_cc.obs[algorithm].unique())
        log_step("Clustering", "COMPLETED", {"n_clusters": n_clusters})

        return adata_cc

    except Exception as e:
        log_error("Clustering", e)
        raise


def run_full_preprocessing(
    adata: AnnData,
    skip_qc: bool = False,
    save_plots: Optional[str] = None,
) -> AnnData:
    """
    Run full preprocessing pipeline: QC -> Normalization -> DR -> Clustering.

    Parameters
    ----------
    adata : AnnData
        Raw input data
    skip_qc : bool
        Skip QC step (if already done)
    save_plots : str, optional
        Base name for saving plots

    Returns
    -------
    AnnData
        Fully preprocessed AnnData
    """
    log_step("FullPreprocessing", "STARTED")

    try:
        # QC
        if not skip_qc:
            adata = run_qc(adata, save_plots=save_plots)
        else:
            print("Skipping QC (already performed)")

        # Normalization
        adata = run_normalization(adata)

        # Dimensionality reduction
        adata = run_dimensionality_reduction(adata)

        # Clustering
        adata = run_clustering(adata)

        log_step("FullPreprocessing", "COMPLETED")
        return adata

    except Exception as e:
        log_error("FullPreprocessing", e)
        raise
