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


def plot_heatmap_scores(
    scores_df: pd.DataFrame,
    top_n_genes: int = 5,
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

    # Select top N genes based on score
    for score in scores:
        top5_by_cluster = (
            scores_df.sort_values(
                by=["stratification", "cluster", score], ascending=[True, True, False]
            )
            .groupby(by=["stratification", "cluster"])
            .head(top_n_genes)["gene"]
            .unique()
            .tolist()
        )

        # Pivot the DataFrame for heatmap
        heatmap_data = scores_df[scores_df["gene"].isin(top5_by_cluster)]
        heatmap_data = heatmap_data.pivot_table(
            index="gene", columns="stratification", values=score, fill_value=0
        )

        # Plot heatmap
        fig, ax = plt.subplots(figsize=config.PLOT_FIGSIZE_MEDIUM)
        sns.heatmap(
            heatmap_data,
            annot=True,
            fmt=".2f",
            cmap="viridis",
            cbar_kws={"label": "Score"},
            ax=ax,
        )
        ax.set_title(f"Top {top_n_genes} Gene Scores Heatmap - {score}")
        ax.set_ylabel("Gene")
        ax.set_xlabel("Stratification")

        fig.savefig(
            f"{config.FIGURES_DIR_GRN}/grn_heatmap_top{top_n_genes}_{score}.png",
            dpi=config.SAVE_DPI,
            bbox_inches="tight",
        )
        plt.close(fig)


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
        stratification_name = file.split("/")[-3]
        score_tmp_df["stratification"] = stratification_name
        merged_scores.append(score_tmp_df)
    merged_scores_df = pd.concat(merged_scores, axis=0)
    merged_scores_df.to_csv(
        f"{config.OUTPUT_DIR}/celloracle/total_merged_scores.csv", index=False
    )
    print(
        f"\n✓ Merged GRN scores saved to: {config.OUTPUT_DIR}/celloracle/total_merged_scores.csv"
    )

    return merged_scores_df
