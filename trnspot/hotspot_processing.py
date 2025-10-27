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

from . import config


def save_hotspot_results(
    hotspot_obj: hs.Hotspot,
):
    """
    Save Hotspot analysis results to specified output directory.

    Parameters:
        hotspot_obj (hs.Hotspot): An instance of the Hotspot class containing analysis results.
        output_dir (str): Directory where results will be saved.

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
    print(f"  Identified {len(modules)} gene modules")

    module_scores = hotspot_obj.calculate_module_scores()
    print("  Module scores calculated")

    plt.close("all")
    hotspot_obj.plot_local_correlations()
    plt.savefig(f"{config.FIGURES_DIR_HOTSPOT}/hotspot_local_correlations.png", dpi=300)
    plt.close()

    save_hotspot_results(hotspot_obj)

    return hotspot_obj
