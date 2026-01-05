"""
GRN visualization functions
===========================

Functions for visualizing gene regulatory network analysis results.
"""

from typing import Optional, List, Tuple
from itertools import combinations
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import networkx as nx
from adjustText import adjust_text

from .. import config


def plot_network_graph(
    links_df: pd.DataFrame,
    score_df: pd.DataFrame,
    score: str = "eigenvector_centrality",
    top_percentile: float = 90,
    top_edges: int = 20,
    figsize: Optional[Tuple[int, int]] = None,
    save: Optional[str] = None,
) -> plt.Figure:
    """
    Plot gene regulatory network graph.

    Parameters
    ----------
    links_df : pd.DataFrame
        DataFrame with network edges (source, target, weight)
    score_df : pd.DataFrame
        DataFrame with node scores
    score : str
        Score column for node selection
    top_percentile : float
        Percentile threshold for top genes
    top_edges : int
        Number of top edges per cluster
    figsize : tuple, optional
        Figure size
    save : str, optional
        Path to save figure

    Returns
    -------
    Figure
        Matplotlib figure
    """
    figsize = figsize or config.PLOT_FIGSIZE_LARGE

    # Select top genes by score
    top_genes = pd.DataFrame()
    filtered_links = pd.DataFrame()

    for cluster in score_df["cluster"].unique():
        cluster_scores = score_df.query(f"cluster == '{cluster}'").copy()
        threshold = np.percentile(cluster_scores[score], top_percentile)
        cluster_top = cluster_scores.query(f"{score} > {threshold}")
        cluster_top["gene"] = cluster_top.index
        top_genes = pd.concat([top_genes, cluster_top], axis=0)

        # Filter links
        cluster_links = links_df[links_df["cluster"] == cluster]
        cluster_links = cluster_links[cluster_links["source"].isin(cluster_top.index)]
        cluster_links = cluster_links.nlargest(top_edges, "coef_abs")
        filtered_links = pd.concat([filtered_links, cluster_links], axis=0)

    # Create graph
    graph = nx.from_pandas_edgelist(
        filtered_links,
        source="source",
        target="target",
        create_using=nx.DiGraph(),
    )

    # Plot
    fig, ax = plt.subplots(figsize=figsize)
    pos = nx.spring_layout(graph, seed=42)

    # Node sizes based on degree
    degrees = dict(graph.degree())
    node_sizes = [degrees.get(node, 1) * 100 for node in graph.nodes()]

    nx.draw_networkx_nodes(
        graph,
        pos,
        ax=ax,
        node_size=node_sizes,
        node_color="skyblue",
        alpha=0.7,
        edgecolors="black",
    )
    nx.draw_networkx_edges(
        graph,
        pos,
        ax=ax,
        alpha=0.5,
        arrows=True,
        arrowstyle="-|>",
        arrowsize=10,
    )
    nx.draw_networkx_labels(graph, pos, ax=ax, font_size=8)

    ax.set_title(f"Gene Regulatory Network - Top genes by {score}")
    ax.axis("off")

    if save:
        fig.savefig(save, dpi=config.SAVE_DPI, bbox_inches="tight")

    return fig


def plot_heatmap_scores(
    scores_df: pd.DataFrame,
    score: str = "eigenvector_centrality",
    top_n_genes: int = 10,
    cluster1: Optional[str] = None,
    cluster2: Optional[str] = None,
    figsize: Optional[Tuple[int, int]] = None,
    save: Optional[str] = None,
) -> plt.Figure:
    """
    Plot heatmap of score differences between clusters/stratifications.

    Parameters
    ----------
    scores_df : pd.DataFrame
        DataFrame with scores (must have gene, cluster, stratification columns)
    score : str
        Score column to visualize
    top_n_genes : int
        Number of top genes per stratification
    cluster1, cluster2 : str, optional
        Clusters to compare (uses first two if not specified)
    figsize : tuple, optional
        Figure size
    save : str, optional
        Path to save figure

    Returns
    -------
    Figure
        Matplotlib figure
    """
    # Get clusters
    clusters = scores_df["cluster"].unique()
    if len(clusters) < 2:
        raise ValueError("Need at least 2 clusters to compare")

    cluster1 = cluster1 or clusters[0]
    cluster2 = cluster2 or clusters[1]

    # Melt and compute differences
    heatmap_data = scores_df[scores_df["cluster"].isin([cluster1, cluster2])].copy()
    heatmap_data.loc[heatmap_data["cluster"] == cluster2, score] *= -1

    # Aggregate by gene and stratification
    agg_data = (
        heatmap_data.groupby(["gene", "stratification"])[score].sum().reset_index()
    )

    # Pivot
    pivot_df = agg_data.pivot(
        index="gene", columns="stratification", values=score
    ).fillna(0)

    # Select top genes per stratification
    selected = set()
    for strat in pivot_df.columns:
        selected.update(pivot_df[strat].abs().nlargest(top_n_genes).index)

    pivot_top = pivot_df.loc[list(selected)].reindex(
        pivot_df.loc[list(selected)]
        .abs()
        .sum(axis=1)
        .sort_values(ascending=False)
        .index
    )

    # Plot
    n_genes = len(pivot_top.index)
    n_strats = len(pivot_top.columns)
    fig_width = max(8, n_strats * 0.8)
    fig_height = max(6, n_genes * 0.3)
    figsize = figsize or (fig_width, fig_height)

    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        pivot_top,
        annot=False,
        cmap="RdBu_r",
        center=0,
        cbar_kws={"label": f"{score} (difference)"},
        ax=ax,
    )
    ax.set_ylabel("Gene")
    ax.set_xlabel("Stratification")
    ax.set_title(f"{score}: {cluster1} vs {cluster2}")
    plt.xticks(rotation=90)

    if save:
        fig.savefig(save, dpi=config.SAVE_DPI, bbox_inches="tight")

    return fig


def plot_scatter_scores(
    score_df: pd.DataFrame,
    score: str,
    cluster1: str,
    cluster2: str,
    percentile: float = 99,
    figsize: Optional[Tuple[int, int]] = None,
    save: Optional[str] = None,
) -> Tuple[plt.Figure, List[str]]:
    """
    Plot scatter comparison of scores between two clusters.

    Parameters
    ----------
    score_df : pd.DataFrame
        DataFrame with scores
    score : str
        Score column to compare
    cluster1, cluster2 : str
        Clusters to compare
    percentile : float
        Percentile for highlighting top genes
    figsize : tuple, optional
        Figure size
    save : str, optional
        Path to save figure

    Returns
    -------
    tuple
        (Figure, list of highlighted genes)
    """
    figsize = figsize or config.PLOT_FIGSIZE_SQUARED

    # Prepare data
    res = score_df[score_df["cluster"].isin([cluster1, cluster2])][[score, "cluster"]]
    res = res.reset_index(drop=False)
    piv = pd.pivot_table(res, values=score, columns="cluster", index="gene").fillna(0)

    # Find top genes
    goi1 = piv[piv[cluster1] > np.percentile(piv[cluster1].values, percentile)].index
    goi2 = piv[piv[cluster2] > np.percentile(piv[cluster2].values, percentile)].index
    gois = list(np.union1d(goi1, goi2))

    # Plot
    fig, ax = plt.subplots(figsize=figsize)
    sns.scatterplot(
        x=piv[cluster1],
        y=piv[cluster2],
        s=50,
        edgecolor="black",
        color="#0096FF",
        ax=ax,
    )

    # Diagonal line
    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),
        np.max([ax.get_xlim(), ax.get_ylim()]),
    ]
    ax.plot(lims, lims, "r--", alpha=0.75, zorder=0)
    ax.set_xlim(lims)
    ax.set_ylim(lims)

    # Annotate top genes
    texts = []
    for goi in gois:
        x, y = piv.loc[goi, cluster1], piv.loc[goi, cluster2]
        texts.append(ax.text(x, y, goi))

    adjust_text(texts, arrowprops=dict(arrowstyle="->", color="r", lw=0.5), ax=ax)

    ax.set_title(f"{score}")
    ax.set_xlabel(cluster1)
    ax.set_ylabel(cluster2)
    sns.despine(ax=ax)

    if save:
        fig.savefig(save, dpi=config.SAVE_DPI, bbox_inches="tight")

    return fig, gois


def plot_difference_cluster_scores(
    score_df: pd.DataFrame,
    score: str,
    cluster1: str,
    cluster2: str,
    top_n: int = 5,
    figsize: Optional[Tuple[int, int]] = None,
    save: Optional[str] = None,
) -> plt.Figure:
    """
    Plot rank plot of score differences between clusters.

    Parameters
    ----------
    score_df : pd.DataFrame
        DataFrame with scores
    score : str
        Score column to compare
    cluster1, cluster2 : str
        Clusters to compare
    top_n : int
        Number of top genes to highlight
    figsize : tuple, optional
        Figure size
    save : str, optional
        Path to save figure

    Returns
    -------
    Figure
        Matplotlib figure
    """
    figsize = figsize or config.PLOT_FIGSIZE_SQUARED

    # Get data for each cluster
    data1 = score_df[score_df["cluster"] == cluster1].set_index(
        "gene" if "gene" in score_df.columns else score_df.index.name
    )
    data2 = score_df[score_df["cluster"] == cluster2].set_index(
        "gene" if "gene" in score_df.columns else score_df.index.name
    )

    # Find common genes
    common = list(set(data1.index) & set(data2.index))

    diff_df = pd.DataFrame(
        {
            cluster1: data1.loc[common, score],
            cluster2: data2.loc[common, score],
        }
    )
    diff_df["Difference"] = diff_df[cluster1] - diff_df[cluster2]
    diff_df = diff_df.sort_values("Difference", ascending=False)

    # Plot
    fig, ax = plt.subplots(figsize=figsize)

    colors = ["gray"] * len(diff_df)
    for i in range(top_n):
        colors[i] = "#0096FF"
        colors[-(i + 1)] = "#fb3310"

    ax.scatter(range(len(diff_df)), diff_df["Difference"], c=colors, alpha=0.7, s=20)
    ax.axhline(y=0, color="gray", linestyle="--", alpha=0.5)
    ax.set_xlabel("Rank")
    ax.set_ylabel("Difference")
    ax.set_title(f"{score}: {cluster1} - {cluster2}")
    sns.despine(ax=ax)

    if save:
        fig.savefig(save, dpi=config.SAVE_DPI, bbox_inches="tight")

    return fig


def plot_compare_cluster_scores(
    score_df: pd.DataFrame,
    score: str,
    top_n: int = 10,
    figsize: Optional[Tuple[int, int]] = None,
    save: Optional[str] = None,
) -> plt.Figure:
    """
    Plot comparison of top genes across clusters.

    Parameters
    ----------
    score_df : pd.DataFrame
        DataFrame with scores
    score : str
        Score column to compare
    top_n : int
        Number of top genes per cluster
    figsize : tuple, optional
        Figure size
    save : str, optional
        Path to save figure

    Returns
    -------
    Figure
        Matplotlib figure
    """
    clusters = sorted(score_df["cluster"].unique())
    n_clusters = len(clusters)
    n_cols = int(np.ceil(np.sqrt(n_clusters)))
    n_rows = int(np.ceil(n_clusters / n_cols))

    figsize = figsize or (7 * n_cols, 7 * n_rows)
    fig, axes = plt.subplots(nrows=n_rows, ncols=n_cols, figsize=figsize, squeeze=False)
    axes = axes.flatten()

    for idx, cluster in enumerate(clusters):
        cluster_data = score_df[score_df["cluster"] == cluster]
        top_genes = cluster_data.nlargest(top_n, score)

        ax = axes[idx]
        gene_col = "gene" if "gene" in top_genes.columns else top_genes.index.name
        sns.scatterplot(
            data=top_genes,
            x=score,
            y=gene_col if gene_col else top_genes.index,
            ax=ax,
            s=100,
            edgecolor="black",
        )
        sns.despine(ax=ax)
        ax.set_title(f"{cluster}")
        ax.set_xlabel(score)
        ax.set_ylabel("Gene")

    # Hide unused axes
    for idx in range(len(clusters), len(axes)):
        axes[idx].set_visible(False)

    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=config.SAVE_DPI, bbox_inches="tight")

    return fig
