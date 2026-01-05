"""
TRNspot plotting module
=======================

This module contains all visualization functions organized by analysis type.
"""

from .qc import (
    plot_qc_violin,
    plot_qc_scatter,
    plot_qc_summary,
)
from .embeddings import (
    plot_umap,
    plot_pca,
    plot_umap_clusters,
)
from .grn import (
    plot_network_graph,
    plot_heatmap_scores,
    plot_scatter_scores,
    plot_difference_cluster_scores,
    plot_compare_cluster_scores,
)
from .hotspot import (
    plot_hotspot_modules,
    plot_hotspot_heatmap,
)

__all__ = [
    # QC plots
    "plot_qc_violin",
    "plot_qc_scatter",
    "plot_qc_summary",
    # Embedding plots
    "plot_umap",
    "plot_pca",
    "plot_umap_clusters",
    # GRN plots
    "plot_network_graph",
    "plot_heatmap_scores",
    "plot_scatter_scores",
    "plot_difference_cluster_scores",
    "plot_compare_cluster_scores",
    # Hotspot plots
    "plot_hotspot_modules",
    "plot_hotspot_heatmap",
]
