"""
Hotspot visualization functions
===============================

Functions for visualizing Hotspot gene module analysis results.
"""

from typing import Optional, Tuple
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.patches import Patch

from .. import config


def plot_hotspot_modules(
    modules: pd.Series,
    local_correlation_z: pd.DataFrame,
    linkage,
    figsize: Optional[Tuple[int, int]] = None,
    save: Optional[str] = None,
) -> plt.Figure:
    """
    Plot Hotspot gene modules as clustered heatmap.

    Parameters
    ----------
    modules : pd.Series
        Gene module assignments from Hotspot
    local_correlation_z : pd.DataFrame
        Z-scored local correlations
    linkage : array
        Hierarchical clustering linkage
    figsize : tuple, optional
        Figure size
    save : str, optional
        Path to save figure

    Returns
    -------
    Figure
        Matplotlib ClusterGrid
    """
    figsize = figsize or (12, 12)

    # Create module colors
    colors = list(plt.get_cmap("tab10").colors)
    module_colors = {i: colors[(i - 1) % len(colors)] for i in modules.unique()}
    module_colors[-1] = "#ffffff"  # Unassigned genes

    row_colors = pd.Series(
        [module_colors[i] for i in modules],
        index=local_correlation_z.index,
    )

    # Create clustermap
    g = sns.clustermap(
        local_correlation_z,
        row_linkage=linkage,
        col_linkage=linkage,
        row_colors=row_colors,
        cmap="RdBu_r",
        vmin=-8,
        vmax=8,
        xticklabels=False,
        yticklabels=False,
        rasterized=True,
        figsize=figsize,
    )

    # Add legend
    legend_elements = [
        Patch(facecolor=color, edgecolor="k", label=f"Module {module}")
        for module, color in module_colors.items()
        if module != -1
    ]

    g.ax_heatmap.legend(
        handles=legend_elements,
        loc="upper left",
        bbox_to_anchor=(1.05, 1),
        title="Modules",
        frameon=False,
    )

    if save:
        plt.savefig(save, dpi=config.SAVE_DPI, bbox_inches="tight")

    return g


def plot_hotspot_heatmap(
    hs_obj,
    annotate_modules: bool = True,
    enrichment_df: Optional[pd.DataFrame] = None,
    figsize: Optional[Tuple[int, int]] = None,
    save: Optional[str] = None,
):
    """
    Plot Hotspot local correlation heatmap with optional annotations.

    Parameters
    ----------
    hs_obj : Hotspot
        Hotspot object with computed results
    annotate_modules : bool
        Whether to add module annotations
    enrichment_df : pd.DataFrame, optional
        Enrichment results for annotation
    figsize : tuple, optional
        Figure size
    save : str, optional
        Path to save figure

    Returns
    -------
    ClusterGrid
        Seaborn ClusterGrid object
    """
    figsize = figsize or (12, 12)

    # Get module colors
    colors = list(plt.get_cmap("tab10").colors)
    module_colors = {i: colors[(i - 1) % len(colors)] for i in hs_obj.modules.unique()}
    module_colors[-1] = "#ffffff"

    row_colors = pd.Series(
        [module_colors[i] for i in hs_obj.modules],
        index=hs_obj.local_correlation_z.index,
    )

    # Create clustermap
    g = sns.clustermap(
        hs_obj.local_correlation_z,
        row_linkage=hs_obj.linkage,
        col_linkage=hs_obj.linkage,
        row_colors=row_colors,
        cmap="RdBu_r",
        vmin=-8,
        vmax=8,
        xticklabels=False,
        yticklabels=False,
        rasterized=True,
        figsize=figsize,
    )

    # Add module legend
    legend_elements = [
        Patch(facecolor=color, edgecolor="k", label=f"Module {module}")
        for module, color in module_colors.items()
        if module != -1
    ]

    g.ax_heatmap.legend(
        handles=legend_elements,
        loc="upper left",
        bbox_to_anchor=(1.05, 1),
        title="Modules",
        frameon=False,
    )

    if save:
        plt.savefig(save, dpi=config.SAVE_DPI, bbox_inches="tight")

    return g


def plot_module_scores_umap(
    adata,
    module_scores: pd.DataFrame,
    modules_to_plot: Optional[list] = None,
    figsize: Optional[Tuple[int, int]] = None,
    save_dir: Optional[str] = None,
):
    """
    Plot module scores on UMAP embedding.

    Parameters
    ----------
    adata : AnnData
        AnnData with UMAP coordinates
    module_scores : pd.DataFrame
        Module scores from Hotspot (cells x modules)
    modules_to_plot : list, optional
        Specific modules to plot (default: all)
    figsize : tuple, optional
        Figure size
    save_dir : str, optional
        Directory to save figures

    Returns
    -------
    list
        List of figures
    """
    import scanpy as sc

    figsize = figsize or config.PLOT_FIGSIZE_MEDIUM

    if modules_to_plot is None:
        modules_to_plot = module_scores.columns.tolist()

    figs = []
    for module in modules_to_plot:
        if module not in module_scores.columns:
            continue

        # Add module score to adata temporarily
        adata.obs[f"module_{module}"] = module_scores[module].values

        fig = sc.pl.umap(
            adata,
            color=f"module_{module}",
            title=f"Module {module}",
            show=False,
            return_fig=True,
        )

        if save_dir:
            fig.savefig(
                f"{save_dir}/hotspot_module_{module}_umap.png",
                dpi=config.SAVE_DPI,
                bbox_inches="tight",
            )

        figs.append(fig)

        # Clean up
        del adata.obs[f"module_{module}"]

    return figs
