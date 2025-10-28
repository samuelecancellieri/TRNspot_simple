import scanpy as sc
from scipy import cluster
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from adjustText import adjust_text

# from statannotations.Annotator import Annotator
import networkx as nx
import pandas as pd
import marsilea as ma
from . import config
from itertools import combinations
import numpy as np


def plot_heatmap_single_score(
    heatmap_data: pd.DataFrame,
    cluster1: str,
    cluster2: str,
    score: str,
    top_n_genes: int = 5,
):
    """
    Plot a heatmap for a single score and stratification.

    Args:
        scores_df (pd.DataFrame): DataFrame containing scores to plot.
        strat (str): Stratification to filter the DataFrame.
        score (str): Score column to plot.
        top_n_genes (int): Number of top genes to display.

    Returns:
        None
    """

    heatmap_data_melt_score = heatmap_data.query("score == @score").copy()

    heatmap_data_melt_score_transformed = heatmap_data_melt_score.copy()
    heatmap_data_melt_score_transformed.loc[
        heatmap_data_melt_score_transformed["cluster"] == cluster2, "value"
    ] *= -1

    heatmap_data_melt_score_final = (
        heatmap_data_melt_score_transformed.groupby(
            ["gene", "stratification", "score"]
        )["value"]
        .sum()
        .reset_index()
    )

    # Create pivot table: genes x stratifications with the aggregated (signed) eigenvector values
    pivot_df = heatmap_data_melt_score_final.pivot(
        index="gene", columns="stratification", values="value"
    ).fillna(0)

    # Select top genes by absolute aggregated signal across stratifications
    # select top_n_genes per stratification then merge (union)
    selected = set()
    for strat in pivot_df.columns:
        selected.update(pivot_df[strat].abs().nlargest(top_n_genes).index)
    selected = list(selected)
    # keep ordering by descending absolute aggregated signal across stratifications
    pivot_top = pivot_df.loc[selected].reindex(
        pivot_df.loc[selected].abs().sum(axis=1).sort_values(ascending=False).index
    )

    # Plot heatmap
    fig, ax = plt.subplots(figsize=config.PLOT_FIGSIZE_SQUARED_LARGE)
    sns.heatmap(
        pivot_top,
        annot=False,
        fmt=".2f",
        cmap="RdBu_r",
        center=0,
        cbar_kws={"label": f"{score} (difference)"},
        ax=ax,
    )
    # ax.set_title(f"Difference in {score} across stratifications")
    ax.set_ylabel("Gene")
    ax.set_xlabel("Stratification")
    plt.xticks(rotation=90)

    # Save and close
    fig.savefig(
        f"{config.FIGURES_DIR_GRN}/grn_heatmap_{score}_difference_top{top_n_genes}_{cluster1}_vs_{cluster2}.png",
        dpi=config.SAVE_DPI,
        bbox_inches="tight",
    )
    plt.close(fig)

    return 0


def plot_heatmap_scores(
    scores_df: pd.DataFrame,
    top_n_genes: int = 10,
    scores: list[str] = [
        "degree_all",
        "degree_centrality_all",
        "degree_in",
        "degree_centrality_in",
        "degree_out",
        "degree_centrality_out",
        "betweenness_centrality",
        "eigenvector_centrality",
    ],
):
    """
    Plot a heatmap of scores using seaborn.

    Args:
        scores_df (pd.DataFrame): DataFrame containing scores to plot.
        x (str): Column name for x-axis.
        y (str): Column name for y-axis.
        score (str): Column name for score values.
        title (str): Title of the heatmap.
        figsize (tuple[int, int]): Figure size.
        cmap (str): Colormap to use.

    Returns:
        matplotlib.figure.Figure: The generated heatmap figure.
    """

    stratifications = scores_df["stratification"].unique()
    clusters = scores_df["cluster"].unique()
    cluster_combinations = combinations(clusters, 2)

    heatmap_data_melt = scores_df.melt(
        id_vars=["gene", "cluster", "stratification"],
        value_vars=scores,
        var_name="score",
        value_name="value",
    )

    for cluster1, cluster2 in cluster_combinations:
        for score in scores:
            plot_heatmap_single_score(
                heatmap_data=heatmap_data_melt[
                    (heatmap_data_melt["cluster"].isin([cluster1, cluster2]))
                    & (heatmap_data_melt["score"] == score)
                ],
                score=score,
                cluster1=cluster1,
                cluster2=cluster2,
                top_n_genes=top_n_genes,
            )

    return 0


def plot_score_comparison_2D(
    score_df: pd.DataFrame,
    value,
    cluster1,
    cluster2,
    percentile=99,
    dot_color="black",
    edge_color=None,
    annot_shifts=None,
    save=None,
    fillna_with_zero=True,
    plt_show=True,
):
    """
    Make a scatter plot that shows the relationship of a specific network score in two groups.

    Args:
        links (Links object): See network_analisis.Links class for detail.
        value (srt): The network score to be shown.
        cluster1 (str): Cluster nome to analyze. Network scores in the cluste1 are shown as x-axis.
        cluster2 (str): Cluster nome to analyze. Network scores in the cluste2 are shown as y-axis.
        percentile (float): Genes with a network score above the percentile will be shown with annotation. Default is 99.
        annot_shifts ((float, float)): Shift x and y cordinate for annotations.
        save (str): Folder path to save plots. If the folde does not exist in the path, the function create the folder.
            If None plots will not be saved. Default is None.
    """

    # Prepare data for plotting
    res = score_df[score_df.cluster.isin([cluster1, cluster2])][[value, "cluster"]]
    res = res.reset_index(drop=False)
    piv = pd.pivot_table(res, values=value, columns="cluster", index="gene")
    if fillna_with_zero:
        piv = piv.fillna(0)
    else:
        piv = piv.fillna(piv.mean(axis=0))
    piv["sum"] = piv[cluster1] + piv[cluster2]

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
        np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
    ]
    ax.plot(lims, lims, "r--", alpha=0.75, zorder=0, linewidth=1)
    ax.set_xlim(lims)
    ax.set_ylim(lims)

    # if annot_shifts is None:
    #     x_shift, y_shift = (x.max() - x.min()) * 0.03, (y.max() - y.min()) * 0.03
    # else:
    #     x_shift, y_shift = annot_shifts

    texts = list()
    for goi in gois:
        x, y = piv.loc[goi, cluster1], piv.loc[goi, cluster2]
        texts.append(ax.text(x, y, goi))

    adjust_text(
        texts=texts,
        arrowprops=dict(arrowstyle="->", color="r", lw=0.5),
        ax=ax,
        time_lim=2,
    )

    if plt_show:
        plt.show()
        return gois, fig, texts
    else:
        return gois, fig, texts


def plot_scatter_scores(
    score_df: pd.DataFrame,
    scores_list: list[str] = [
        "degree_all",
        "degree_centrality_all",
        "degree_in",
        "degree_centrality_in",
        "degree_out",
        "degree_centrality_out",
        "betweenness_centrality",
        "eigenvector_centrality",
    ],
):
    """
    Plot a scatter plot comparing scores between two stratifications for a given cluster.

    Args:
        scores_df (pd.DataFrame): DataFrame containing scores.
        cluster (str): The cluster to analyze.
        score (str): The score column to compare.
        stratifications (list[str]): List of two stratifications to compare.
    Returns:
        matplotlib.figure.Figure: The generated plot figure.
    """
    clusters = sorted(set(score_df["cluster"].tolist()))
    cluster_combinations = combinations(clusters, 2)

    for cluster1, cluster2 in cluster_combinations:
        for score in scores_list:
            # print(f"Generating scatter plot for {score} - {cluster1} vs {cluster2}")
            genes_t, plot_t, texts_t = plot_score_comparison_2D(
                score_df=score_df,
                value=score,
                cluster1=cluster1,
                cluster2=cluster2,
                percentile=99,
                annot_shifts=None,
                edge_color="black",
                dot_color="#0096FF",  # Use a specific color for the dots
                save=None,
                fillna_with_zero=True,
                plt_show=False,
            )

            cluster1_clean = cluster1.replace(" ", "_").strip()
            cluster2_clean = cluster2.replace(" ", "_").strip()
            plot_t.axes[0].grid(False)
            sns.despine(ax=plot_t.axes[0])

            plot_t.savefig(
                f"{config.FIGURES_DIR_GRN}/grn_deep_analysis/grn_scatter_{score}_{cluster1_clean}_vs_{cluster2_clean}.png",
                dpi=config.SAVE_DPI,
                bbox_inches="tight",
            )
            plt.close(plot_t)

    return 0


def plot_difference_cluster_scores(
    score_df: pd.DataFrame,
    scores: list[str] = [
        "degree_all",
        "degree_centrality_all",
        "degree_in",
        "degree_centrality_in",
        "degree_out",
        "degree_centrality_out",
        "betweenness_centrality",
        "eigenvector_centrality",
    ],
):
    """
    Plot the difference in scores between two stratifications for a given cluster.

    Args:
        scores_df (pd.DataFrame): DataFrame containing scores.
        cluster (str): The cluster to analyze.
        score (str): The score column to compare.
        stratifications (list[str]): List of two stratifications to compare.

    Returns:
        matplotlib.figure.Figure: The generated plot figure.
    """

    # Filter DataFrame for the specific cluster and stratifications
    clusters = sorted(set(score_df["cluster"].tolist()))
    stratifications = sorted(set(score_df["stratification"].tolist()))
    if len(stratifications) == 1:
        stratification_name = stratifications[0]
    else:
        stratification_name = "multiple stratifications"

    if len(clusters) < 2:
        print(f"Not enough clusters to compare for cluster")
        return None
    if len(clusters) > 2:
        print(
            f"More than two clusters found. Will compute differences for all the combinations."
        )

    cluster_combinations = combinations(clusters, 2)

    for cluster_pair in cluster_combinations:
        stratifications = list(cluster_pair)

        cluster_data_1 = score_df[score_df["cluster"] == stratifications[0]]
        cluster_data_2 = score_df[score_df["cluster"] == stratifications[1]]

        tf_treat = set(cluster_data_1.index.tolist())
        tf_untreat = set(cluster_data_2.index.tolist())
        tf_overlap = tf_treat & tf_untreat

        # Prepare DataFrame for plotting differences
        common_tfs = list(tf_overlap)

        for score in scores:
            if common_tfs:
                diff_df = pd.DataFrame(
                    {
                        stratifications[0]: cluster_data_1.loc[common_tfs, score],
                        stratifications[1]: cluster_data_2.loc[common_tfs, score],
                    }
                )
                diff_df["Difference"] = (
                    diff_df[stratifications[0]] - diff_df[stratifications[1]]
                )

            # Create a figure for the plot
            fig, ax = plt.subplots(figsize=config.PLOT_FIGSIZE_SQUARED)

            # Sort by difference for ranking
            diff_df_sorted = diff_df.sort_values("Difference", ascending=False)
            ranks = range(len(diff_df_sorted))

            top_x = 5
            # Get top 5 positive and negative differences
            top_5_positive = diff_df_sorted.head(top_x)
            top_5_negative = diff_df_sorted.tail(top_x)
            # Color the top 5 positive and negative points
            top_5_positive_indices = range(top_x)
            top_5_negative_indices = range(
                len(diff_df_sorted) - top_x, len(diff_df_sorted)
            )

            # Color all points gray first
            colors = ["gray"] * len(diff_df_sorted)

            # Color top 5 positive points
            for i in top_5_positive_indices:
                colors[i] = "#0096FF"

            # Color top 5 negative points
            for i in top_5_negative_indices:
                colors[i] = "#fb3310"

            # Create scatter plot with colored points
            ax.scatter(ranks, diff_df_sorted["Difference"], c=colors, alpha=0.7, s=20)

            # Annotate top 5 positive differences
            texts = []
            for i, (gene, row) in enumerate(top_5_positive.iterrows()):
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

            # Annotate top 5 negative differences
            for i, (gene, row) in enumerate(top_5_negative.iterrows()):
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
                texts,
                arrowprops=dict(arrowstyle="->", color="black", lw=0.5),
                ax=ax,
            )

            ax.axhline(y=0, color="gray", linestyle="--", alpha=0.5)
            ax.set_xlabel("Rank")
            ax.set_ylabel("Difference")
            ax.set_title(
                f"Rank Plot: {score} Difference {stratifications[0]}-{stratifications[1]} in {stratification_name}\nTop genes"
            )
            sns.despine(ax=ax)
            ax.grid(False)
            # plt.tight_layout()
            fig.savefig(
                f"{config.FIGURES_DIR_GRN}/grn_deep_analysis/top_genes_difference_{stratifications[0]}-{stratifications[1]}_{score}_rankplot.png",
                bbox_inches="tight",
                dpi=300,
            )
            # plt.show()

    return 0


def plot_compare_cluster_scores(
    score_df: pd.DataFrame,
    scores: list[str] = [
        "degree_all",
        "degree_centrality_all",
        "degree_in",
        "degree_centrality_in",
        "degree_out",
        "degree_centrality_out",
        "betweenness_centrality",
        "eigenvector_centrality",
    ],
):
    """
    Compare scores across clusters and plot differences.
    """

    clusters = sorted(set(score_df["cluster"].tolist()))
    stratifications = sorted(set(score_df["stratification"].tolist()))
    if len(stratifications) == 1:
        stratification_name = stratifications[0]
    else:
        stratification_name = "multiple stratifications"

    # For each score, create barplots showing top 10 genes per cluster for each stratification
    for score in scores:
        fig, axes = plt.subplots(
            nrows=1,
            ncols=len(clusters),
            figsize=(7 * len(clusters), 7),
            squeeze=False,
        )

        for idx, cluster in enumerate(clusters):
            cluster_data = score_df[score_df["cluster"] == cluster]

            # Get top 10 genes by score for this cluster
            top10_genes = cluster_data.nlargest(10, score)

            # Plot barplot
            ax = axes[0, idx]
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

        fig.subplots_adjust(wspace=0.2, top=0.9)
        fig.suptitle(
            f"Comparison of Top 10 Genes by {score} Across Clusters in {stratification_name}",
            fontsize=16,
        )
        fig.savefig(
            f"{config.FIGURES_DIR_GRN}/grn_deep_analysis/grn_barplot_{score}_{stratification_name}_top10.png",
            dpi=config.SAVE_DPI,
            bbox_inches="tight",
        )
        plt.close(fig)

    return 0


def process_single_score_file(score_file: str) -> pd.DataFrame:
    """
    Process a single score CSV file and return a DataFrame.

    Args:
        score_file (str): Path to the score CSV file.

    Returns:
        pd.DataFrame: Processed DataFrame with scores.
    """

    score_df = pd.read_csv(score_file, index_col=0)
    # score_df.rename(columns={"Unnamed: 0": "gene"}, inplace=True)
    score_df.index.name = "gene"
    stratification_name = score_file.split("/")[-3]
    score_df["stratification"] = stratification_name

    return score_df


def merge_scores(tracked_file: str) -> pd.DataFrame:
    """
    Merge multiple score CSV files into a single DataFrame.
    This function reads multiple CSV files containing scores, adds a stratification
    column based on the parent directory name of each file, and concatenates them
    into a single DataFrame.
    Args:
        tracked_file (str): A file path to a text file containing paths to CSV files with scores.
    Returns:
        pd.DataFrame: A concatenated DataFrame containing all scores from the input
            files, with an additional 'stratification' column indicating the source
            stratification (extracted from the parent directory name of each file).
    Note:
        The stratification name is extracted from the second-to-last component of
        the file path (i.e., file.split("/")[-2]).
        Files are concatenated along axis=0 (rows).
    """

    with open(tracked_file) as f:
        file_list = f.readlines()
    file_list = [file.strip() for file in file_list]

    merged_scores = list()
    for file in file_list:
        if file.count("grn_merged_scores.csv") == 0:
            continue
        print(f"Processing file: {file}")
        score_tmp_df = pd.read_csv(file)
        score_tmp_df.rename(columns={"Unnamed: 0": "gene"}, inplace=True)
        # score_tmp_df.index.name = "gene"
        stratification_name = file.split("/")[-3]
        score_tmp_df["stratification"] = stratification_name
        merged_scores.append(score_tmp_df)

    # concatenate all DataFrames
    merged_scores_df = pd.concat(merged_scores, axis=0)
    merged_scores_df.to_csv(
        f"{config.OUTPUT_DIR}/celloracle/total_merged_scores.csv", index=False
    )
    print(
        f"\n✓ Merged GRN scores saved to: {config.OUTPUT_DIR}/celloracle/total_merged_scores.csv"
    )

    return merged_scores_df
