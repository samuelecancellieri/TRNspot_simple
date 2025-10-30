"""
Data preprocessing module for TRNspot
"""

import os
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import numpy as np
from typing import Optional, Tuple
from anndata import AnnData

from . import config


def perform_qc(
    adata: AnnData,
    min_genes: Optional[int] = None,
    min_counts: Optional[int] = None,
    max_counts: Optional[int] = None,
    pct_counts_mt_max: Optional[float] = None,
    min_cells: Optional[int] = None,
    figsize: Optional[Tuple[int, int]] = None,
    save_plots: Optional[str] = None,
) -> AnnData:
    """
    Perform quality control on AnnData object.

    This function calculates QC metrics, filters cells based on specified thresholds,
    and generates QC visualization plots.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with cells as observations and genes as variables.
    min_genes : int, optional
        Minimum number of genes expressed required for a cell to pass filtering.
        If None, uses config.QC_MIN_GENES (default: 200).
    min_counts : int, optional
        Minimum number of counts required for a cell to pass filtering.
        If None, uses config.QC_MIN_COUNTS (default: 500).
    max_counts : int, optional
        Maximum number of counts allowed for a cell to pass filtering.
        If None, uses config.QC_MAX_COUNTS (default: None, no upper limit).
    pct_counts_mt_max : float, optional
        Maximum percentage of mitochondrial counts allowed for a cell.
        If None, uses config.QC_PCT_MT_MAX (default: 20.0).
    min_cells : int, optional
        Minimum number of cells expressing a gene for the gene to be kept.
        If None, uses config.QC_MIN_CELLS (default: 10).
    plot : bool, default=True
        Whether to generate QC plots.
    figsize : tuple, optional
        Figure size for QC plots (width, height).
        If None, uses config.PLOT_FIGSIZE_LARGE (default: (15, 10)).
    save_plots : str, optional
        Path to save the QC plots. If None, plots are not saved.

    Returns
    -------
    AnnData
        Filtered AnnData object with QC metrics stored in .obs and .var.

    Examples
    --------
    >>> import scanpy as sc
    >>> from trnspot.preprocessing import perform_qc
    >>> from trnspot import config
    >>>
    >>> # Use default config values
    >>> adata = sc.read_h5ad('data.h5ad')
    >>> adata_qc = perform_qc(adata)
    >>>
    >>> # Override specific parameters
    >>> adata_qc = perform_qc(adata, min_genes=300, min_counts=1000)
    >>>
    >>> # Use updated config
    >>> config.update_config(QC_MIN_GENES=500)
    >>> adata_qc = perform_qc(adata)
    """

    # Use config defaults if not specified
    if min_genes is None:
        min_genes = config.QC_MIN_GENES
    if min_counts is None:
        min_counts = config.QC_MIN_COUNTS
    if max_counts is None:
        max_counts = config.QC_MAX_COUNTS
    if pct_counts_mt_max is None:
        pct_counts_mt_max = config.QC_PCT_MT_MAX
    if min_cells is None:
        min_cells = config.QC_MIN_CELLS
    if figsize is None:
        figsize = config.PLOT_FIGSIZE_MEDIUM

    # Make a copy to avoid modifying the original
    adata_cc = adata.copy()

    # Store initial cell and gene counts
    n_cells_initial = adata_cc.n_obs
    n_genes_initial = adata_cc.n_vars

    print(f"Initial data shape: {n_cells_initial} cells × {n_genes_initial} genes")

    # mitochondrial genes, "MT-" for human, "Mt-" for mouse
    adata_cc.var["mt"] = adata_cc.var_names.str.startswith("MT-")
    # ribosomal genes
    adata_cc.var["ribo"] = adata_cc.var_names.str.startswith(("RPS", "RPL"))
    # hemoglobin genes
    adata_cc.var["hb"] = adata_cc.var_names.str.contains("^HB[^(P)]")

    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(
        adata_cc,
        qc_vars=["mt", "ribo", "hb"],
        percent_top=None,
        log1p=True,
        inplace=True,
    )

    # Create violin plot for pre-filtering metrics
    fig, axes = plt.subplots(1, 3, figsize=config.PLOT_FIGSIZE_LARGE)

    sns.violinplot(data=adata_cc.obs, y="n_genes_by_counts", ax=axes[0], inner="box")
    axes[0].set_ylabel("Number of genes")
    axes[0].set_title("Genes per cell")

    sns.violinplot(data=adata_cc.obs, y="total_counts", ax=axes[1], inner="box")
    axes[1].set_ylabel("Total counts")
    axes[1].set_title("Total counts per cell")

    sns.violinplot(data=adata_cc.obs, y="pct_counts_mt", ax=axes[2], inner="box")
    axes[2].set_ylabel("% Mitochondrial counts")
    axes[2].set_title("Mitochondrial percentage")

    if save_plots:
        fig.savefig(
            f"{config.FIGURES_DIR_QC}/violin_pre_filter_{save_plots}.png",
            dpi=300,
            bbox_inches="tight",
        )

    # Create scatter plot for pre-filtering metrics
    fig, ax = plt.subplots(figsize=config.PLOT_FIGSIZE_MEDIUM)
    scatter = ax.scatter(
        adata_cc.obs["total_counts"],
        adata_cc.obs["n_genes_by_counts"],
        c=adata_cc.obs["pct_counts_mt"],
        cmap=config.PLOT_COLOR_PALETTE,
        alpha=0.7,
        s=5,
    )
    ax.set_xlabel("Total counts")
    ax.set_ylabel("Number of genes")
    ax.set_title("Genes vs Total Counts colored by % Mitochondrial")
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label("% Mitochondrial counts")
    if save_plots:
        fig.savefig(
            f"{config.FIGURES_DIR_QC}/scatter_pre_filter_{save_plots}.png",
            dpi=300,
            bbox_inches="tight",
        )

    # Apply filtering
    print("\nApplying filters...")
    print(f"  - Minimum genes per cell: {min_genes}")
    print(f"  - Minimum counts per cell: {min_counts}")
    if max_counts:
        print(f"  - Maximum counts per cell: {max_counts}")
    print(f"  - Maximum mitochondrial percentage: {pct_counts_mt_max}%")

    # Filter cells
    sc.pp.filter_cells(adata_cc, min_genes=min_genes)
    sc.pp.filter_cells(adata_cc, min_counts=min_counts)

    if max_counts:
        sc.pp.filter_cells(adata_cc, max_counts=max_counts)

    adata_cc = adata_cc[adata_cc.obs["pct_counts_mt"] < pct_counts_mt_max, :].copy()

    # Filter genes (keep genes expressed in at least min_cells to preserve rare cell populations)
    sc.pp.filter_genes(adata_cc, min_cells=min_cells)

    # Report filtering results
    n_cells_filtered = adata_cc.n_obs
    n_genes_filtered = adata_cc.n_vars
    cells_removed = n_cells_initial - n_cells_filtered
    genes_removed = n_genes_initial - n_genes_filtered

    print(f"\nFiltering results:")
    print(
        f"  - Cells removed: {cells_removed} ({cells_removed/n_cells_initial*100:.2f}%)"
    )
    print(
        f"  - Genes removed: {genes_removed} ({genes_removed/n_genes_initial*100:.2f}%)"
    )
    print(f"  - Final shape: {n_cells_filtered} cells × {n_genes_filtered} genes")

    # Create seaborn violin plot for post-filtering metrics
    fig, axes = plt.subplots(1, 3, figsize=config.PLOT_FIGSIZE_LARGE)

    sns.violinplot(data=adata_cc.obs, y="n_genes_by_counts", ax=axes[0], inner="box")
    axes[0].set_ylabel("Number of genes")
    axes[0].set_title("Genes per cell")

    sns.violinplot(data=adata_cc.obs, y="total_counts", ax=axes[1], inner="box")
    axes[1].set_ylabel("Total counts")
    axes[1].set_title("Total counts per cell")

    sns.violinplot(data=adata_cc.obs, y="pct_counts_mt", ax=axes[2], inner="box")
    axes[2].set_ylabel("% Mitochondrial counts")
    axes[2].set_title("Mitochondrial percentage")

    if save_plots:
        fig.savefig(
            f"{config.FIGURES_DIR_QC}/violin_post_filter_{save_plots}.png",
            dpi=300,
            bbox_inches="tight",
        )

    # scatter of post-filtering metrics
    fig, ax = plt.subplots(figsize=config.PLOT_FIGSIZE_MEDIUM)
    scatter = ax.scatter(
        adata_cc.obs["total_counts"],
        adata_cc.obs["n_genes_by_counts"],
        c=adata_cc.obs["pct_counts_mt"],
        cmap=config.PLOT_COLOR_PALETTE,
        alpha=0.7,
    )
    ax.set_xlabel("Total counts")
    ax.set_ylabel("Number of genes")
    ax.set_title("Genes vs Total Counts colored by % Mitochondrial")
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label("% Mitochondrial counts")
    if save_plots:
        fig.savefig(
            f"{config.FIGURES_DIR_QC}/scatter_post_filter_{save_plots}.png",
            dpi=300,
            bbox_inches="tight",
        )

    plt.close("all")

    return adata_cc


def perform_normalization(adata: AnnData) -> AnnData:
    """
    Perform normalization on AnnData object.

    This function normalizes the data to a target sum and applies log transformation.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with cells as observations and genes as variables.

    Returns
    -------
    AnnData
        Normalized AnnData object.

    Examples
    --------
    >>> import scanpy as sc
    >>> from trnspot.preprocessing import perform_normalization
    >>> from trnspot import config
    >>>
    >>> adata = sc.read_h5ad('data.h5ad')
    >>> adata_normalized = perform_normalization(adata)
    """

    print("\nPerforming normalization...")
    # Make a copy to avoid modifying the original
    adata_cc = adata.copy()

    # Check if data is already normalized (heuristic: max count > 100 indicates raw counts)
    if adata_cc.X.max() > 100:
        print(
            "Warning: Data appears to be unnormalized (max count > 100). Proceeding with normalization."
        )
    else:
        print(
            "Data appears to be already normalized (max count <= 100). Skipping normalization."
        )
        return adata_cc

    # Saving raw count data
    if "raw_counts" not in adata_cc.layers:
        adata_cc.layers["raw_counts"] = adata_cc.X.copy()
        print("Stored raw counts in layer 'raw_counts'")

    # Normalize to target sum
    sc.pp.normalize_total(adata_cc, target_sum=config.NORMALIZE_TARGET_SUM)
    print(f"Normalized to target sum: {config.NORMALIZE_TARGET_SUM}")

    # Log transform
    sc.pp.log1p(adata_cc)
    print("Applied log1p transformation")

    return adata_cc


# %%


def perform_grn_pre_processing(
    adata: AnnData,
    cluster_key: str,
    cell_downsample: int = 20000,
    top_genes: Optional[int] = None,
    n_neighbors: Optional[int] = None,
    n_pcs: Optional[int] = None,
    svd_solver: Optional[str] = None,
) -> AnnData:
    """
    Perform preprocessing for GRN analysis on AnnData object.

    This function performs QC, normalizes the data, identifies highly variable genes,
    and scales the data in preparation for gene regulatory network analysis.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with cells as observations and genes as variables.
    cluster_key : str, optional
        Key in adata.obs to use for clustering over GRN.
    cell_downsample : int, default=20000
        Number of cells to downsample to (default: 20000).
    top_genes : int, optional
        Number of highly variable genes to select.
        If None, uses config.HVGS_N_TOP_GENES (default: 2000).
    n_neighbors : int, optional
        Number of neighbors for KNN graph.
        If None, uses config.NEIGHBORS_N_NEIGHBORS (default: 15).
    n_pcs : int, optional
        Number of principal components to use.
        If None, uses config.NEIGHBORS_N_PCS (default: 40).
    svd_solver : str, optional
        SVD solver for PCA.
        If None, uses config.PCA_SVDSOLVE (default: 'arpack').

    Returns
    -------
    adata : AnnData
        Preprocessed AnnData object ready for GRN analysis.

    Examples
    --------
    >>> import scanpy as sc
    >>> from trnspot.preprocessing import perform_grn_pre_processing
    >>> from trnspot import config
    >>>
    >>> # Use default config values
    >>> adata = sc.read_h5ad('data.h5ad')
    >>> adata_preprocessed = perform_grn_pre_processing(adata, cluster_key='louvain')
    >>>
    >>> # Override specific parameters
    >>> adata_preprocessed = perform_grn_pre_processing(adata, top_genes=5000, n_pcs=50)
    """

    # Use config defaults if not specified
    if top_genes is None:
        top_genes = config.HVGS_N_TOP_GENES
    if n_neighbors is None:
        n_neighbors = config.NEIGHBORS_N_NEIGHBORS
    if n_pcs is None:
        n_pcs = config.NEIGHBORS_N_PCS
    if svd_solver is None:
        svd_solver = config.PCA_SVDSOLVE

    # Make a copy to avoid modifying the original
    adata_cc = adata.copy()

    # Filter genes based on variance and subset
    sc.pp.highly_variable_genes(adata_cc, n_top_genes=top_genes, subset=True)
    print(f"  Selected {top_genes} highly variable genes")

    # PCA
    sc.tl.pca(adata_cc, svd_solver=svd_solver)
    print(f"  Computed PCA with {svd_solver} solver")

    # Diffusion map
    sc.pp.neighbors(adata_cc, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.diffmap(adata_cc)
    print(f"  Computed diffusion map with {n_neighbors} neighbors and {n_pcs} PCs")

    # Calculate neighbors again based on diffusion map
    sc.pp.neighbors(adata_cc, n_neighbors=n_neighbors, use_rep="X_diffmap")
    print(f"  Recomputed neighbors using diffusion map representation")

    # Clustering
    sc.tl.paga(adata_cc, groups=cluster_key)
    sc.pl.paga(adata_cc, show=None, save=None)
    print(f"  Computed PAGA for cluster key: {cluster_key}")

    sc.tl.draw_graph(adata_cc, init_pos="paga")
    sc.pl.draw_graph(adata_cc, color=cluster_key, show=None, save=None)
    print(f"  Computed draw_graph for cluster key: {cluster_key}")

    if adata_cc.n_obs > cell_downsample:
        print(f"  Downsampling to {cell_downsample} cells for GRN analysis")
        sc.pp.subsample(adata_cc, n_obs=cell_downsample)

    return adata_cc
