"""
Quality Control plotting functions
==================================

Functions for visualizing QC metrics.
"""

from typing import Optional, Tuple
import matplotlib.pyplot as plt
import seaborn as sns
from anndata import AnnData

from .. import config


def plot_qc_violin(
    adata: AnnData,
    metrics: Tuple[str, ...] = ("n_genes_by_counts", "total_counts", "pct_counts_mt"),
    figsize: Optional[Tuple[int, int]] = None,
    save: Optional[str] = None,
) -> plt.Figure:
    """
    Create violin plots for QC metrics.

    Parameters
    ----------
    adata : AnnData
        AnnData with QC metrics computed
    metrics : tuple
        Metrics to plot from adata.obs
    figsize : tuple, optional
        Figure size
    save : str, optional
        Path to save figure

    Returns
    -------
    Figure
        Matplotlib figure object
    """
    figsize = figsize or config.PLOT_FIGSIZE_LARGE
    n_metrics = len(metrics)

    fig, axes = plt.subplots(1, n_metrics, figsize=figsize)
    if n_metrics == 1:
        axes = [axes]

    titles = {
        "n_genes_by_counts": "Genes per cell",
        "total_counts": "Total counts per cell",
        "pct_counts_mt": "Mitochondrial %",
        "pct_counts_ribo": "Ribosomal %",
    }

    for ax, metric in zip(axes, metrics):
        if metric not in adata.obs.columns:
            ax.set_visible(False)
            continue

        sns.violinplot(data=adata.obs, y=metric, ax=ax, inner="box")
        ax.set_ylabel(metric)
        ax.set_title(titles.get(metric, metric))

    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=config.SAVE_DPI, bbox_inches="tight")

    return fig


def plot_qc_scatter(
    adata: AnnData,
    x: str = "total_counts",
    y: str = "n_genes_by_counts",
    color: str = "pct_counts_mt",
    figsize: Optional[Tuple[int, int]] = None,
    save: Optional[str] = None,
) -> plt.Figure:
    """
    Create scatter plot of QC metrics.

    Parameters
    ----------
    adata : AnnData
        AnnData with QC metrics
    x : str
        Column for x-axis
    y : str
        Column for y-axis
    color : str
        Column for coloring points
    figsize : tuple, optional
        Figure size
    save : str, optional
        Path to save figure

    Returns
    -------
    Figure
        Matplotlib figure object
    """
    figsize = figsize or config.PLOT_FIGSIZE_MEDIUM

    fig, ax = plt.subplots(figsize=figsize)

    scatter = ax.scatter(
        adata.obs[x],
        adata.obs[y],
        c=adata.obs[color],
        cmap=config.PLOT_COLOR_PALETTE,
        alpha=0.7,
        s=5,
    )

    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.set_title(f"{y} vs {x} colored by {color}")

    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label(color)

    if save:
        fig.savefig(save, dpi=config.SAVE_DPI, bbox_inches="tight")

    return fig


def plot_qc_summary(
    adata: AnnData,
    save_dir: Optional[str] = None,
    prefix: str = "",
) -> Tuple[plt.Figure, plt.Figure]:
    """
    Create comprehensive QC summary plots.

    Parameters
    ----------
    adata : AnnData
        AnnData with QC metrics
    save_dir : str, optional
        Directory to save figures
    prefix : str
        Prefix for saved files

    Returns
    -------
    tuple
        (violin_figure, scatter_figure)
    """
    violin_save = None
    scatter_save = None

    if save_dir:
        violin_save = (
            f"{save_dir}/qc_violin_{prefix}.png"
            if prefix
            else f"{save_dir}/qc_violin.png"
        )
        scatter_save = (
            f"{save_dir}/qc_scatter_{prefix}.png"
            if prefix
            else f"{save_dir}/qc_scatter.png"
        )

    fig_violin = plot_qc_violin(adata, save=violin_save)
    fig_scatter = plot_qc_scatter(adata, save=scatter_save)

    return fig_violin, fig_scatter
