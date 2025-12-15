"""
Hotspot processing module for TRNspot
"""

from typing import Optional
import os
import scanpy as sc
import hotspot as hs
from scipy.sparse import csc_matrix
import matplotlib.pyplot as plt
from anndata import AnnData
import pickle
import pandas as pd
import seaborn as sns
import marsilea as ma

from . import config

from . import enrichment_analysis as ea
from matplotlib.patches import Patch


def plot_hotspot_annotation(hs_obj: hs.Hotspot):
    """
    Plot Hotspot gene module annotations on spatial coordinates.

    Parameters:
        hs_obj (hs.Hotspot): An instance of the Hotspot class containing analysis results.

    """

    print("Generating Hotspot local correlation heatmap with annotations...")

    df_enrichment = pd.DataFrame()
    for module in hs_obj.modules.unique():
        genes = hs_obj.modules[hs_obj.modules == module].index.tolist()
        df_module_enrichment = ea.gseapy_ora_enrichment_analysis(genes).results
        df_module_enrichment.columns = [
            x.replace(" ", "_") for x in df_module_enrichment.columns
        ]
        df_module_enrichment["module"] = module
        df_enrichment = pd.concat([df_enrichment, df_module_enrichment])
    df_enrichment.to_csv(
        f"{config.OUTPUT_DIR}/hotspot/hotspot_module_enrichment_results.csv",
        index=False,
    )

    row_colors = None
    colors = list(plt.get_cmap("tab10").colors)
    module_colors = {i: colors[(i - 1) % len(colors)] for i in hs_obj.modules.unique()}
    module_colors[-1] = "#ffffff"

    row_colors1 = pd.Series(
        [module_colors[i] for i in hs_obj.modules],
        index=hs_obj.local_correlation_z.index,
    )

    row_colors = pd.DataFrame(
        {
            "Modules": row_colors1,
        }
    )
    # Get top 3 annotations for each module
    module_annotations = {}
    for module in df_enrichment["module"].unique():
        if module == -1:  # Skip unassigned genes
            continue
        module_df = df_enrichment[df_enrichment["module"] == module]
        if not module_df.empty:
            top3 = module_df.nlargest(3, "Combined_Score")["Term"].tolist()
            module_annotations[module] = "; ".join(top3[:3])
        else:
            module_annotations[module] = "No enrichment"

    # Create annotation labels for genes based on their module
    gene_annotations = pd.Series(
        [module_annotations.get(i, "") if i != -1 else "" for i in hs_obj.modules],
        index=hs_obj.local_correlation_z.index,
        name="Top_Annotations",
    )

    row_colors["Top_Annotations"] = gene_annotations

    # Add heatmap with custom parameters

    # Prepare row colors for clustermap
    # row_colors dataframe needs to be converted to colors if not already
    # In the previous code, row_colors['Modules'] contained hex codes.
    # We need to ensure the index matches the data.

    # Create the clustermap
    g = sns.clustermap(
        hs_obj.local_correlation_z,
        row_linkage=hs_obj.linkage,
        col_linkage=hs_obj.linkage,
        row_colors=row_colors["Modules"],
        cmap="RdBu_r",
        vmin=-8,
        vmax=8,
        xticklabels=False,
        yticklabels=False,
        rasterized=True,
        figsize=(12, 12),
    )

    # Add a legend for the modules
    # Create a list of patches for the legend

    legend_elements = [
        Patch(facecolor=color, edgecolor="k", label=f"Module {module}")
        for module, color in module_colors.items()
        if module != -1
    ]

    # Add the legend to the figure
    g.ax_heatmap.legend(
        handles=legend_elements,
        loc="upper left",
        bbox_to_anchor=(1.05, 1),
        title="Modules",
        frameon=False,
    )

    # Save the figure
    plt.savefig(
        f"{config.FIGURES_DIR_HOTSPOT}/hotspot_local_correlation_heatmap_with_annotations.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()


def save_hotspot_results(
    hotspot_obj: hs.Hotspot,
):
    """
    Save Hotspot analysis results to specified output directory.

    Parameters:
        hotspot_obj (hs.Hotspot): An instance of the Hotspot class containing analysis results.

    """
    # Get results summary
    autocorr_results = hotspot_obj.results

    significant_genes = autocorr_results[
        autocorr_results.FDR < config.HOTSPOT_FDR_THRESHOLD
    ]

    module_scores = hotspot_obj.module_scores
    module_scores.index.name = "cell_id"

    gene_modules = hotspot_obj.modules

    # Save additional results
    results_path = f"{config.OUTPUT_DIR}/hotspot/autocorrelation_results.csv"
    autocorr_results.to_csv(results_path)
    print(f"\n✓ Autocorrelation results saved to: {results_path}")

    significant_path = f"{config.OUTPUT_DIR}/hotspot/significant_genes.csv"
    significant_genes.to_csv(significant_path)
    print(f"✓ Significant genes saved to: {significant_path}")

    module_scores_path = f"{config.OUTPUT_DIR}/hotspot/hotspot_module_scores.csv"
    module_scores.to_csv(module_scores_path)
    print(f"✓ Module scores saved to: {module_scores_path}")

    modules_path = f"{config.OUTPUT_DIR}/hotspot/gene_modules.csv"
    gene_modules.to_csv(modules_path)
    print(f"✓ Gene modules saved to: {modules_path}")

    hotspot_object_path = f"{config.OUTPUT_DIR}/hotspot/hotspot_object.pkl"
    with open(hotspot_object_path, "wb") as f:
        pickle.dump(hotspot_obj, f)
    print(f"✓ Hotspot object saved to: {hotspot_object_path}")


def create_hotspot_object(
    adata: AnnData,
    top_genes: int = config.HOTSPOT_TOP_GENES,
    layer_key: str = "raw_counts",
    model: str = "danb",
    embedding_key: str = "X_pca",
    normalization_key: str = "total_counts",
):
    """
    Creates a Hotspot object for spatial gene module analysis.

    Parameters:
        adata (AnnData): Annotated data object containing gene expression data.
        top_genes (int): Number of top highly variable genes to select. If None, use all genes. Default is 3000.
        layer_key (str): Name of the layer in `adata.layers` that contains the expression data to be used.
        model (str): Statistical model to use for Hotspot analysis. Default is "danb (depth-adjusted negative binomial model)".
        embedding_key (str): Name of the embedding in `adata.obsm` to be used for spatial analysis. Default is "X_pca".
        normalization_key (str): Key in `adata.obs` for normalization. Default is "total_counts".

    Returns:
        hotspot_obj (hs.Hotspot): An instance of the Hotspot class.

    """

    # Create a copy of the AnnData object to avoid modifying the original
    adata_cc = adata.copy()

    if top_genes:
        print(f"Selecting top {top_genes} highly variable genes for Hotspot analysis.")
        sc.pp.highly_variable_genes(adata_cc, n_top_genes=top_genes, subset=True)
    else:
        print("Using all genes for Hotspot analysis.")

    # create a csv matrix from the specified layer or from adata.X
    if layer_key:
        print("Using layer counts for Hotspot analysis.")
        adata_cc.layers[f"csc_{layer_key}"] = csc_matrix(adata_cc.layers[layer_key])
    else:
        print("Using adata.X for Hotspot analysis.")
        adata_cc.layers["csc_X"] = csc_matrix(adata_cc.X)

    # Create Hotspot object
    hotspot_obj = hs.Hotspot(
        adata_cc,
        layer_key=f"csc_{layer_key}" if layer_key else "csc_X",
        model=model,
        latent_obsm_key=embedding_key,
        umi_counts_obs_key=normalization_key,
    )

    return hotspot_obj


def run_hotspot_analysis(hotspot_obj):
    """
    Run Hotspot analysis on the given Hotspot object.

    Parameters:
        hotspot_obj: An instance of the Hotspot class.

    Returns:
        hotspot_obj: The updated Hotspot object with analysis results.
    """

    # Create KNN graph
    hotspot_obj.create_knn_graph(
        weighted_graph=False, n_neighbors=config.HOTSPOT_N_NEIGHBORS
    )
    print("  KNN graph created successfully")

    # Compute autocorrelations
    hs_results = hotspot_obj.compute_autocorrelations(jobs=config.HOTSPOT_N_JOBS)
    print("  Autocorrelations computed successfully")

    # Identify significant genes
    hs_genes = hs_results.loc[
        hs_results.FDR < config.HOTSPOT_FDR_THRESHOLD
    ].index  # Select genes
    local_correlations = hotspot_obj.compute_local_correlations(
        hs_genes, jobs=config.HOTSPOT_N_JOBS
    )
    print(
        f"  Identified {len(hs_genes)} significant genes (FDR < {config.HOTSPOT_FDR_THRESHOLD})"
    )

    modules = hotspot_obj.create_modules(
        min_gene_threshold=config.HOTSPOT_MIN_GENES_PER_MODULE,
        core_only=config.HOTSPOT_CORE_ONLY,
        fdr_threshold=config.HOTSPOT_FDR_THRESHOLD,
    )
    print(f"  Identified {len(modules.unique())} gene modules")

    module_scores = hotspot_obj.calculate_module_scores()
    print("  Module scores calculated")

    plt.close("all")
    hotspot_obj.plot_local_correlations()
    plt.savefig(f"{config.FIGURES_DIR_HOTSPOT}/hotspot_local_correlations.png", dpi=300)
    plt.close()

    # plot_hotspot_annotation(hotspot_obj)

    save_hotspot_results(hotspot_obj)

    return hotspot_obj
