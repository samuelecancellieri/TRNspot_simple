"""
Embedding visualization functions
=================================

Functions for visualizing PCA, UMAP, and other embeddings.
"""

from typing import Optional, Tuple, List, Union
import matplotlib.pyplot as plt
import scanpy as sc
from anndata import AnnData

from .. import config


def plot_umap(
    adata: AnnData,
    color: Optional[Union[str, List[str]]] = None,
    figsize: Optional[Tuple[int, int]] = None,
    save: Optional[str] = None,
    **kwargs,
) -> plt.Figure:
    """
    Plot UMAP embedding.

    Parameters
    ----------
    adata : AnnData
        AnnData with UMAP coordinates
    color : str or list, optional
        Variables to color by
    figsize : tuple, optional
        Figure size
    save : str, optional
        Path to save figure
    **kwargs
        Additional arguments for sc.pl.umap

    Returns
    -------
    Figure
        Matplotlib figure
    """
    figsize = figsize or config.PLOT_FIGSIZE_MEDIUM

    fig = sc.pl.umap(
        adata,
        color=color,
        show=False,
        return_fig=True,
        **kwargs,
    )

    if save:
        fig.savefig(save, dpi=config.SAVE_DPI, bbox_inches="tight")

    return fig


def plot_pca(
    adata: AnnData,
    color: Optional[Union[str, List[str]]] = None,
    components: str = "1,2",
    figsize: Optional[Tuple[int, int]] = None,
    save: Optional[str] = None,
    **kwargs,
) -> plt.Figure:
    """
    Plot PCA embedding.

    Parameters
    ----------
    adata : AnnData
        AnnData with PCA coordinates
    color : str or list, optional
        Variables to color by
    components : str
        PCA components to plot (e.g., "1,2")
    figsize : tuple, optional
        Figure size
    save : str, optional
        Path to save figure
    **kwargs
        Additional arguments for sc.pl.pca

    Returns
    -------
    Figure
        Matplotlib figure
    """
    figsize = figsize or config.PLOT_FIGSIZE_MEDIUM

    fig = sc.pl.pca(
        adata,
        color=color,
        components=components,
        show=False,
        return_fig=True,
        **kwargs,
    )

    if save:
        fig.savefig(save, dpi=config.SAVE_DPI, bbox_inches="tight")

    return fig


def plot_umap_clusters(
    adata: AnnData,
    cluster_key: str = "leiden",
    figsize: Optional[Tuple[int, int]] = None,
    save: Optional[str] = None,
    **kwargs,
) -> plt.Figure:
    """
    Plot UMAP colored by cluster labels.

    Parameters
    ----------
    adata : AnnData
        AnnData with UMAP and cluster labels
    cluster_key : str
        Key for cluster labels in adata.obs
    figsize : tuple, optional
        Figure size
    save : str, optional
        Path to save figure
    **kwargs
        Additional arguments for sc.pl.umap

    Returns
    -------
    Figure
        Matplotlib figure
    """
    return plot_umap(
        adata,
        color=cluster_key,
        figsize=figsize,
        save=save,
        legend_loc="on data",
        **kwargs,
    )


def plot_pca_variance(
    adata: AnnData,
    n_pcs: int = 50,
    figsize: Optional[Tuple[int, int]] = None,
    save: Optional[str] = None,
) -> plt.Figure:
    """
    Plot PCA variance explained.

    Parameters
    ----------
    adata : AnnData
        AnnData with PCA computed
    n_pcs : int
        Number of PCs to show
    figsize : tuple, optional
        Figure size
    save : str, optional
        Path to save figure

    Returns
    -------
    Figure
        Matplotlib figure
    """
    figsize = figsize or config.PLOT_FIGSIZE_SMALL

    fig, ax = plt.subplots(figsize=figsize)

    variance_ratio = adata.uns["pca"]["variance_ratio"][:n_pcs]
    ax.plot(range(1, len(variance_ratio) + 1), variance_ratio.cumsum(), "o-")
    ax.set_xlabel("Principal Component")
    ax.set_ylabel("Cumulative Variance Explained")
    ax.set_title("PCA Variance")
    ax.grid(True, alpha=0.3)

    if save:
        fig.savefig(save, dpi=config.SAVE_DPI, bbox_inches="tight")

    return fig
