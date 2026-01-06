"""
GRN Plotting Module for TRNspot
===============================

Functions for generating Gene Regulatory Network visualizations including:
- Network graphs
- Heatmaps of centrality scores
- Scatter plots comparing scores across clusters
- Rank plots showing score differences
"""

import os
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from adjustText import adjust_text
from itertools import combinations
from typing import Optional, List, Dict, Any, Tuple

from .. import config
from .utils import save_plot, plot_exists


# =============================================================================
# Network Graph Plots
# =============================================================================


def plot_network_graph_single(
    graph: nx.DiGraph,
    cluster: str,
    stratification: str,
    score: str,
    node_size_factor: float = 1000.0,
    edge_width_factor: float = 1.0,
    skip_existing: bool = True,
) -> bool:
    """
    Plot a network graph for a single score.

    Parameters
    ----------
    graph : nx.DiGraph
        The directed graph to plot.
    cluster : str
        Cluster name for labeling.
    stratification : str
        Stratification name for labeling.
    score : str
        Score being visualized.
    node_size_factor : float
        Factor to scale node sizes.
    edge_width_factor : float
        Factor to scale edge widths.
    skip_existing : bool
        If True, skip if file already exists.

    Returns
    -------
    bool
        True if plot was generated, False if skipped.
    """
    filepath = (
        f"{config.FIGURES_DIR_GRN}/grn_network_{score}_{stratification}_{cluster}.png"
    )

    if plot_exists(filepath, skip_existing):
        return False

    fig = plt.figure(figsize=config.PLOT_FIGSIZE_WIDE)

    pos = nx.spring_layout(graph, seed=42)

    # Node sizes based on degree
    degrees = dict(graph.degree())
    node_sizes = [degrees[node] * node_size_factor for node in graph.nodes()]

    # Edge widths based on weight
    edge_weights = nx.get_edge_attributes(graph, "weight")
    edge_widths = [
        edge_weights[edge] * edge_width_factor if edge in edge_weights else 1.0
        for edge in graph.edges()
    ]

    nx.draw_networkx_nodes(
        graph,
        pos,
        node_size=node_sizes,
        node_color="skyblue",
        alpha=0.7,
        edgecolors="black",
    )
    nx.draw_networkx_edges(
        graph,
        pos,
        width=edge_widths,
        alpha=0.5,
        arrows=True,
        arrowstyle="-|>",
        arrowsize=10,
    )
    nx.draw_networkx_labels(graph, pos, font_size=8)

    plt.title(f"GRN - Cluster: {cluster}, Stratification: {stratification}")
    plt.axis("off")
    plt.tight_layout()

    return save_plot(
        fig=fig,
        filepath=filepath,
        plot_type="grn",
        metadata={
            "plot_name": "network_graph",
            "score": score,
            "cluster": cluster,
            "stratification": stratification,
        },
        skip_existing=False,
    )


def plot_network_graph(
    score_df: pd.DataFrame,
    links_df: pd.DataFrame,
    scores: Optional[List[str]] = None,
    skip_existing: bool = True,
) -> int:
    """
    Plot network graphs for multiple scores.

    Parameters
    ----------
    score_df : pd.DataFrame
        DataFrame containing scores.
    links_df : pd.DataFrame
        DataFrame containing links.
    scores : list, optional
        List of scores to plot. If None, uses default centrality scores.
    skip_existing : bool
        If True, skip existing plots.

    Returns
    -------
    int
        Number of plots generated.
    """
    if scores is None:
        scores = [
            "degree_all",
            "degree_centrality_all",
            "betweenness_centrality",
            "eigenvector_centrality",
        ]

    plots_generated = 0
    score_df = score_df.copy()
    score_df["cluster"] = score_df["cluster"].astype(str)

    for score in scores:
        top_genes_by_cluster = pd.DataFrame()
        filtered_links_df = pd.DataFrame()

        for cluster in score_df["cluster"].unique():
            cluster_str = str(cluster)
            cluster_scores = score_df.query(f"cluster == '{cluster_str}'").copy()

            if cluster_scores.empty or cluster_scores[score].dropna().empty:
                continue

            percentile = np.percentile(cluster_scores[score].dropna(), 90)
            cluster_scores = cluster_scores.query(f"{score} > {percentile}")

            if cluster_scores.empty:
                continue

            cluster_scores["gene"] = cluster_scores.index
            top_genes_by_cluster = pd.concat(
                [top_genes_by_cluster, cluster_scores], axis=0
            )

            filtered_links_cluster = links_df[links_df["cluster"] == cluster_str]
            filtered_links_cluster = filtered_links_cluster[
                filtered_links_cluster["source"].isin(
                    top_genes_by_cluster.query("cluster==@cluster_str")["gene"]
                )
            ]
            filtered_links_cluster["cluster"] = cluster_str
            filtered_links_cluster = filtered_links_cluster.nlargest(20, "coef_abs")
            filtered_links_df = pd.concat(
                [filtered_links_df, filtered_links_cluster], axis=0
            )

        if top_genes_by_cluster.empty or filtered_links_df.empty:
            continue

        graph = nx.from_pandas_edgelist(
            filtered_links_df,
            source="source",
            target="target",
            create_using=nx.Graph(),
        )

        if graph.number_of_nodes() == 0:
            continue

        components = list(nx.connected_components(graph))
        if not components:
            continue

        # Get stratification name
        stratification = "combined"
        if "stratification" in score_df.columns:
            strats = score_df["stratification"].unique()
            if len(strats) == 1:
                stratification = str(strats[0])

        # Plot for each cluster
        for cluster in score_df["cluster"].unique():
            if plot_network_graph_single(
                graph=graph,
                cluster=str(cluster),
                stratification=stratification,
                score=score,
                skip_existing=skip_existing,
            ):
                plots_generated += 1

    return plots_generated


# =============================================================================
# Heatmap Plots
# =============================================================================


def plot_heatmap_single_score(
    heatmap_data: pd.DataFrame,
    cluster1: str,
    cluster2: str,
    score: str,
    top_n_genes: int = 5,
    skip_existing: bool = True,
) -> bool:
    """
    Plot a heatmap for a single score comparing two clusters.

    Parameters
    ----------
    heatmap_data : pd.DataFrame
        DataFrame with melted score data.
    cluster1 : str
        First cluster name.
    cluster2 : str
        Second cluster name.
    score : str
        Score column to plot.
    top_n_genes : int
        Number of top genes to display per stratification.
    skip_existing : bool
        If True, skip if file already exists.

    Returns
    -------
    bool
        True if plot was generated, False if skipped.
    """
    filepath = (
        f"{config.FIGURES_DIR_GRN}/grn_heatmap_{score}_difference_"
        f"top{top_n_genes}_{cluster1}_vs_{cluster2}.png"
    )

    if plot_exists(filepath, skip_existing):
        return False

    heatmap_data_score = heatmap_data.query("score == @score").copy()

    if heatmap_data_score.empty:
        return False

    # Transform values (negative for cluster2)
    heatmap_transformed = heatmap_data_score.copy()
    heatmap_transformed.loc[heatmap_transformed["cluster"] == cluster2, "value"] *= -1

    # Aggregate by gene and stratification
    heatmap_final = (
        heatmap_transformed.groupby(["gene", "stratification", "score"])["value"]
        .sum()
        .reset_index()
    )

    # Create pivot table
    pivot_df = heatmap_final.pivot(
        index="gene", columns="stratification", values="value"
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

    # Create figure
    n_genes = len(pivot_top.index)
    n_strats = len(pivot_top.columns)
    fig_width = max(8, n_strats * 0.8)
    fig_height = max(6, n_genes * 0.3)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    sns.heatmap(
        pivot_top,
        annot=False,
        fmt=".2f",
        cmap="RdBu_r",
        center=0,
        cbar_kws={"label": f"{score} (difference)"},
        ax=ax,
    )

    ax.set_ylabel("Gene")
    ax.set_xlabel("Stratification")
    plt.xticks(rotation=90)

    return save_plot(
        fig=fig,
        filepath=filepath,
        plot_type="grn",
        metadata={
            "plot_name": "heatmap",
            "score": score,
            "cluster1": cluster1,
            "cluster2": cluster2,
            "top_n_genes": top_n_genes,
        },
        skip_existing=False,
    )


def plot_heatmap_scores(
    scores_df: pd.DataFrame,
    top_n_genes: int = 10,
    scores: Optional[List[str]] = None,
    skip_existing: bool = True,
) -> int:
    """
    Plot heatmaps for multiple scores across cluster combinations.

    Parameters
    ----------
    scores_df : pd.DataFrame
        DataFrame containing scores.
    top_n_genes : int
        Number of top genes to display.
    scores : list, optional
        List of scores to plot. If None, uses default centrality scores.
    skip_existing : bool
        If True, skip existing plots.

    Returns
    -------
    int
        Number of plots generated.
    """
    if scores is None:
        scores = [
            "degree_centrality_all",
            "degree_centrality_in",
            "degree_centrality_out",
            "betweenness_centrality",
            "eigenvector_centrality",
        ]

    clusters = scores_df["cluster"].unique()
    cluster_combinations = list(combinations(clusters, 2))

    heatmap_data_melt = scores_df.melt(
        id_vars=["gene", "cluster", "stratification"],
        value_vars=scores,
        var_name="score",
        value_name="value",
    )

    plots_generated = 0

    for cluster1, cluster2 in cluster_combinations:
        for score in scores:
            filtered_data = heatmap_data_melt[
                (heatmap_data_melt["cluster"].isin([cluster1, cluster2]))
                & (heatmap_data_melt["score"] == score)
            ]

            if plot_heatmap_single_score(
                heatmap_data=filtered_data,
                cluster1=str(cluster1),
                cluster2=str(cluster2),
                score=score,
                top_n_genes=top_n_genes,
                skip_existing=skip_existing,
            ):
                plots_generated += 1

    return plots_generated


# =============================================================================
# Scatter Plots
# =============================================================================


def plot_score_comparison_2d(
    score_df: pd.DataFrame,
    value: str,
    cluster1: str,
    cluster2: str,
    percentile: float = 99,
    dot_color: str = "#0096FF",
    edge_color: Optional[str] = "black",
    fillna_with_zero: bool = True,
) -> Tuple[np.ndarray, plt.Figure, List]:
    """
    Create a 2D scatter plot comparing scores between two clusters.

    Parameters
    ----------
    score_df : pd.DataFrame
        DataFrame with scores.
    value : str
        Score column to compare.
    cluster1 : str
        First cluster (x-axis).
    cluster2 : str
        Second cluster (y-axis).
    percentile : float
        Percentile threshold for annotation.
    dot_color : str
        Color for scatter dots.
    edge_color : str, optional
        Edge color for dots.
    fillna_with_zero : bool
        If True, fill NaN with 0. Otherwise, fill with mean.

    Returns
    -------
    tuple
        (genes_of_interest, figure, text_annotations)
    """
    res = score_df[score_df.cluster.isin([cluster1, cluster2])][[value, "cluster"]]
    res = res.reset_index(drop=False)
    piv = pd.pivot_table(res, values=value, columns="cluster", index="gene")

    if fillna_with_zero:
        piv = piv.fillna(0)
    else:
        piv = piv.fillna(piv.mean(axis=0))

    goi1 = piv[piv[cluster1] > np.percentile(piv[cluster1].values, percentile)].index
    goi2 = piv[piv[cluster2] > np.percentile(piv[cluster2].values, percentile)].index
    gois = np.union1d(goi1, goi2)

    x, y = piv[cluster1], piv[cluster2]

    fig, ax = plt.subplots(figsize=config.PLOT_FIGSIZE_SQUARED)
    sns.scatterplot(
        x=x, y=y, markers="o", s=50, edgecolor=edge_color, color=dot_color, ax=ax
    )
    ax.set_title(f"{value}")

    # Add diagonal line
    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),
        np.max([ax.get_xlim(), ax.get_ylim()]),
    ]
    ax.plot(lims, lims, "r--", alpha=0.75, zorder=0, linewidth=1)
    ax.set_xlim(lims)
    ax.set_ylim(lims)

    # Add annotations
    texts = []
    for goi in gois:
        x_val, y_val = piv.loc[goi, cluster1], piv.loc[goi, cluster2]
        texts.append(ax.text(x_val, y_val, goi))

    adjust_text(
        texts=texts,
        arrowprops=dict(arrowstyle="->", color="r", lw=0.5),
        ax=ax,
        time_lim=2,
    )

    return gois, fig, texts


def plot_scatter_scores(
    score_df: pd.DataFrame,
    scores_list: Optional[List[str]] = None,
    skip_existing: bool = True,
) -> int:
    """
    Plot scatter plots comparing scores between cluster pairs.

    Parameters
    ----------
    score_df : pd.DataFrame
        DataFrame containing scores.
    scores_list : list, optional
        List of scores to plot. If None, uses default centrality scores.
    skip_existing : bool
        If True, skip existing plots.

    Returns
    -------
    int
        Number of plots generated.
    """
    if scores_list is None:
        scores_list = [
            "degree_centrality_all",
            "degree_centrality_in",
            "degree_centrality_out",
            "betweenness_centrality",
            "eigenvector_centrality",
        ]

    clusters = sorted(set(score_df["cluster"].tolist()))
    cluster_combinations = list(combinations(clusters, 2))

    # Ensure output directory exists
    os.makedirs(f"{config.FIGURES_DIR_GRN}/grn_deep_analysis", exist_ok=True)

    plots_generated = 0

    for cluster1, cluster2 in cluster_combinations:
        for score in scores_list:
            cluster1_clean = str(cluster1).replace(" ", "_").strip()
            cluster2_clean = str(cluster2).replace(" ", "_").strip()
            filepath = (
                f"{config.FIGURES_DIR_GRN}/grn_deep_analysis/"
                f"grn_scatter_{score}_{cluster1_clean}_vs_{cluster2_clean}.png"
            )

            if plot_exists(filepath, skip_existing):
                continue

            genes, fig, texts = plot_score_comparison_2d(
                score_df=score_df,
                value=score,
                cluster1=cluster1,
                cluster2=cluster2,
                percentile=99,
                edge_color="black",
                dot_color="#0096FF",
                fillna_with_zero=True,
            )

            fig.axes[0].grid(False)
            sns.despine(ax=fig.axes[0])

            if save_plot(
                fig=fig,
                filepath=filepath,
                plot_type="grn",
                metadata={
                    "plot_name": "scatter",
                    "score": score,
                    "cluster1": cluster1,
                    "cluster2": cluster2,
                },
                skip_existing=False,
            ):
                plots_generated += 1

    return plots_generated


# =============================================================================
# Difference and Comparison Plots
# =============================================================================


def plot_difference_cluster_scores(
    score_df: pd.DataFrame,
    scores: Optional[List[str]] = None,
    skip_existing: bool = True,
) -> int:
    """
    Plot rank plots showing score differences between clusters.

    Parameters
    ----------
    score_df : pd.DataFrame
        DataFrame containing scores.
    scores : list, optional
        List of scores to plot.
    skip_existing : bool
        If True, skip existing plots.

    Returns
    -------
    int
        Number of plots generated.
    """
    if scores is None:
        scores = [
            "degree_centrality_all",
            "degree_centrality_in",
            "degree_centrality_out",
            "betweenness_centrality",
            "eigenvector_centrality",
        ]

    clusters = sorted(set(score_df["cluster"].tolist()))
    stratifications = sorted(set(score_df["stratification"].tolist()))
    stratification_name = (
        stratifications[0] if len(stratifications) == 1 else "combined"
    )

    if len(clusters) < 2:
        print("  Not enough clusters to compare")
        return 0

    cluster_combinations = list(combinations(clusters, 2))
    os.makedirs(f"{config.FIGURES_DIR_GRN}/grn_deep_analysis", exist_ok=True)

    plots_generated = 0

    for cluster_pair in cluster_combinations:
        cluster1, cluster2 = cluster_pair

        cluster_data_1 = score_df[score_df["cluster"] == cluster1]
        cluster_data_2 = score_df[score_df["cluster"] == cluster2]

        tf_overlap = set(cluster_data_1.index.tolist()) & set(
            cluster_data_2.index.tolist()
        )
        common_tfs = list(tf_overlap)

        if not common_tfs:
            continue

        for score in scores:
            filepath = (
                f"{config.FIGURES_DIR_GRN}/grn_deep_analysis/"
                f"top_genes_difference_{cluster1}-{cluster2}_{score}_rankplot.png"
            )

            if plot_exists(filepath, skip_existing):
                continue

            diff_df = pd.DataFrame(
                {
                    cluster1: cluster_data_1.loc[common_tfs, score],
                    cluster2: cluster_data_2.loc[common_tfs, score],
                }
            )
            diff_df["Difference"] = diff_df[cluster1] - diff_df[cluster2]
            diff_df_sorted = diff_df.sort_values("Difference", ascending=False)

            fig, ax = plt.subplots(figsize=config.PLOT_FIGSIZE_SQUARED)

            top_x = 5
            top_positive = diff_df_sorted.head(top_x)
            top_negative = diff_df_sorted.tail(top_x)

            # Color scheme
            colors = ["gray"] * len(diff_df_sorted)
            for i in range(top_x):
                colors[i] = "#0096FF"
            for i in range(len(diff_df_sorted) - top_x, len(diff_df_sorted)):
                colors[i] = "#fb3310"

            ranks = range(len(diff_df_sorted))
            ax.scatter(ranks, diff_df_sorted["Difference"], c=colors, alpha=0.7, s=20)

            # Annotations
            texts = []
            for i, (gene, row) in enumerate(top_positive.iterrows()):
                texts.append(
                    ax.annotate(
                        gene,
                        (i, row["Difference"]),
                        fontsize=8,
                        bbox=dict(
                            boxstyle="round,pad=0.2", facecolor="lightgreen", alpha=0.7
                        ),
                    )
                )

            for i, (gene, row) in enumerate(top_negative.iterrows()):
                rank_pos = len(diff_df_sorted) - top_x + i
                texts.append(
                    ax.annotate(
                        gene,
                        (rank_pos, row["Difference"]),
                        fontsize=8,
                        bbox=dict(
                            boxstyle="round,pad=0.2", facecolor="lightcoral", alpha=0.7
                        ),
                    )
                )

            adjust_text(
                texts, arrowprops=dict(arrowstyle="->", color="black", lw=0.5), ax=ax
            )

            ax.axhline(y=0, color="gray", linestyle="--", alpha=0.5)
            ax.set_xlabel("Rank")
            ax.set_ylabel("Difference")
            ax.set_title(
                f"Rank Plot: {score}\n{cluster1} vs {cluster2} in {stratification_name}"
            )
            sns.despine(ax=ax)
            ax.grid(False)

            if save_plot(
                fig=fig,
                filepath=filepath,
                plot_type="grn",
                metadata={
                    "plot_name": "difference_rankplot",
                    "score": score,
                    "cluster1": cluster1,
                    "cluster2": cluster2,
                },
                skip_existing=False,
            ):
                plots_generated += 1

    return plots_generated


def plot_compare_cluster_scores(
    score_df: pd.DataFrame,
    scores: Optional[List[str]] = None,
    skip_existing: bool = True,
) -> int:
    """
    Plot bar plots comparing top genes across clusters.

    Parameters
    ----------
    score_df : pd.DataFrame
        DataFrame containing scores.
    scores : list, optional
        List of scores to plot.
    skip_existing : bool
        If True, skip existing plots.

    Returns
    -------
    int
        Number of plots generated.
    """
    if scores is None:
        scores = [
            "degree_centrality_all",
            "degree_centrality_in",
            "degree_centrality_out",
            "betweenness_centrality",
            "eigenvector_centrality",
        ]

    clusters = sorted(set(score_df["cluster"].tolist()))
    stratifications = sorted(set(score_df["stratification"].tolist()))
    stratification_name = (
        stratifications[0] if len(stratifications) == 1 else "combined"
    )

    os.makedirs(f"{config.FIGURES_DIR_GRN}/grn_deep_analysis", exist_ok=True)

    plots_generated = 0

    for score in scores:
        filepath = (
            f"{config.FIGURES_DIR_GRN}/grn_deep_analysis/"
            f"grn_barplot_{score}_{stratification_name}_top10.png"
        )

        if plot_exists(filepath, skip_existing):
            continue

        n_clusters = len(clusters)
        n_cols = int(np.ceil(np.sqrt(n_clusters)))
        n_rows = int(np.ceil(n_clusters / n_cols))

        fig, axes = plt.subplots(
            nrows=n_rows,
            ncols=n_cols,
            figsize=(7 * n_cols, 7 * n_rows),
            squeeze=False,
        )
        axes = axes.flatten()

        for idx, cluster in enumerate(clusters):
            cluster_data = score_df[score_df["cluster"] == cluster]
            top10_genes = cluster_data.nlargest(10, score)

            ax = axes[idx]
            sns.scatterplot(
                data=top10_genes,
                x=score,
                y="gene",
                ax=ax,
                s=100,
                edgecolor="black",
                palette="viridis",
            )
            sns.despine(ax=ax)
            ax.grid(False)
            ax.set_title(f"{cluster}")
            ax.set_xlabel(f"{score}")
            ax.set_ylabel("Gene")

        # Hide unused axes
        for idx in range(len(clusters), len(axes)):
            axes[idx].set_visible(False)

        fig.subplots_adjust(wspace=0.2, top=0.9)
        fig.suptitle(
            f"Top 10 Genes by {score} Across Clusters\n{stratification_name}",
            fontsize=16,
        )

        if save_plot(
            fig=fig,
            filepath=filepath,
            plot_type="grn",
            metadata={
                "plot_name": "barplot_comparison",
                "score": score,
                "stratification": stratification_name,
            },
            skip_existing=False,
        ):
            plots_generated += 1

    return plots_generated


# =============================================================================
# Main Entry Point
# =============================================================================


def generate_all_grn_plots(
    score_df: pd.DataFrame,
    links_df: Optional[pd.DataFrame] = None,
    skip_existing: bool = True,
    scores: Optional[List[str]] = None,
) -> Dict[str, int]:
    """
    Generate all GRN plots from score data.

    Parameters
    ----------
    score_df : pd.DataFrame
        DataFrame containing GRN scores.
    links_df : pd.DataFrame, optional
        DataFrame containing GRN links (for network plots).
    skip_existing : bool
        If True, skip existing plots.
    scores : list, optional
        List of scores to plot.

    Returns
    -------
    dict
        Dictionary mapping plot types to counts generated.
    """
    if scores is None:
        scores = [
            "degree_centrality_all",
            "degree_centrality_in",
            "degree_centrality_out",
            "betweenness_centrality",
            "eigenvector_centrality",
        ]

    # Ensure gene column exists
    if "gene" not in score_df.columns:
        score_df = score_df.reset_index()
        if "index" in score_df.columns:
            score_df = score_df.rename(columns={"index": "gene"})

    results = {}

    print("Generating GRN plots...")

    # Scatter plots
    results["scatter"] = plot_scatter_scores(
        score_df, scores_list=scores, skip_existing=skip_existing
    )

    # Difference plots
    results["difference"] = plot_difference_cluster_scores(
        score_df, scores=scores, skip_existing=skip_existing
    )

    # Comparison plots
    results["comparison"] = plot_compare_cluster_scores(
        score_df, scores=scores, skip_existing=skip_existing
    )

    # Heatmap plots
    if "stratification" in score_df.columns:
        results["heatmap"] = plot_heatmap_scores(
            score_df, scores=scores, skip_existing=skip_existing
        )

    # Network plots (if links provided)
    if links_df is not None:
        results["network"] = plot_network_graph(
            score_df, links_df, scores=scores, skip_existing=skip_existing
        )

    total = sum(results.values())
    print(f"  GRN plots: {total} total generated")
    for plot_type, count in results.items():
        if count > 0:
            print(f"    - {plot_type}: {count}")

    return results
