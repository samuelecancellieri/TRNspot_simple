"""
CellOracle processing module for TRNspot
"""

import scanpy as sc
import celloracle as co
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import numpy as np
from typing import Optional, Tuple
from anndata import AnnData
import pickle

from . import config


# Move content from grn_celloracle_processing.py here
def create_oracle_object(
    adata: AnnData,
    cluster_column_name: str,
    embedding_name: str,
    species: str = "human",
    TG_to_TF_dictionary: Optional[str] = None,
    raw_count_layer: Optional[str] = None,
):
    """
    Creates an Oracle object for CellOracle analysis.

    Parameters:
        adata (AnnData): Annotated data object containing gene expression data.
        TG_to_TF_dictionary (str): Path to the file containing the dictionary mapping target genes to transcription factors.
        cluster_column_name (str): Name of the column in `adata.obs` that contains cluster information.
        embedding_name (str): Name of the embedding to be used for analysis.
        raw_count_layer (str, optional): Name of the layer in `adata.layers` that contains raw count data. Defaults to None.

    Returns:
        oracle (Oracle): An instance of the Oracle class.

    """

    # Create a copy of the AnnData object to avoid modifying the original
    adata_cc = adata.copy()

    # Load base GRN based on species
    base_GRN = None
    if species == "human":
        base_GRN = co.data.load_human_promoter_base_GRN()
    elif species == "mouse":
        base_GRN = co.data.load_mouse_promoter_base_GRN()
    else:
        print("if species is not human or mouse no base GRN is loaded")

    # Create Oracle object
    oracle = co.Oracle()

    if raw_count_layer:
        print("Using raw counts layer for Oracle object creation.")
        # use raw counts to build the oracle object
        adata_cc.X = adata_cc.layers[raw_count_layer].copy()
        oracle.import_anndata_as_raw_count(
            adata=adata,
            cluster_column_name=cluster_column_name,
            embedding_name=embedding_name,
        )
    else:
        print("Using normalized counts for Oracle object creation.")
        # use normalized counts to build the oracle object
        oracle.import_anndata_as_normalized_count(
            adata=adata,
            cluster_column_name=cluster_column_name,
            embedding_name=embedding_name,
        )

    # Import base GRN
    if base_GRN is not None:
        oracle.import_TF_data(TF_info_matrix=base_GRN)

    if TG_to_TF_dictionary:
        # Load the TG to TF dictionary
        TG_to_TF_dictionary = pickle.load(
            open(
                TG_to_TF_dictionary,
                "rb",
            )
        )
        # Add the TG to TF dictionary to the oracle object
        oracle.addTFinfo_dictionary(TG_to_TF_dictionary)

    return oracle


def run_PCA(oracle: co.Oracle):
    """
    Perform Principal Component Analysis (PCA) on a CellOracle object and determine optimal number of components.
    This function creates a copy of the input Oracle object, performs PCA analysis, and automatically
    selects the optimal number of principal components based on the explained variance ratio. It visualizes
    the cumulative explained variance and identifies the elbow point to determine where the rate of
    variance explanation significantly decreases.
    Args:
        oracle (co.Oracle): A CellOracle Oracle object containing gene expression data
                           to be analyzed with PCA.
    Returns:
        co.Oracle: A copy of the input Oracle object with PCA analysis performed.
                   The object will have PCA results stored and the optimal number of
                   components determined (capped at maximum 50 components).
    Notes:
        - The function automatically determines the optimal number of components by finding
          the elbow point in the cumulative explained variance curve
        - A plot showing cumulative explained variance is displayed with a vertical line
          indicating the selected number of components
        - The number of components is limited to a maximum of 50
        - The original Oracle object is preserved; a copy is returned with PCA results
    """

    # Copy the oracle object to preserve the original data
    oracle_cc = oracle.copy()

    # Perform PCA
    oracle_cc.perform_PCA()

    # Select important PCs
    fig, ax = plt.subplots(figsize=config.PLOT_FIGSIZE_MEDIUM)
    ax.plot(np.cumsum(oracle_cc.pca.explained_variance_ratio_)[:200])
    n_comps = np.where(
        np.diff(np.diff(np.cumsum(oracle_cc.pca.explained_variance_ratio_)) > 0.002)
    )[0][0]
    ax.axvline(n_comps, c="k")
    fig.savefig(config.FIGURES_DIR + "/pca_explained_variance.png")
    plt.close("all")

    n_comps = min(n_comps, 50)

    return oracle_cc, n_comps


def run_KNN(oracle: co.Oracle, n_comps: int = 50):
    """
    Perform K-Nearest Neighbors (KNN) imputation on a CellOracle object.

    This function creates a copy of the input Oracle object and applies KNN imputation
    to smooth gene expression data based on cellular similarity. The number of neighbors
    is automatically calculated as 2.5% of the total cell count.

    Parameters
    ----------
    oracle : co.Oracle
        A CellOracle Oracle object containing single-cell data and analysis results.
    n_comps : int, optional
        Number of principal components to use for KNN calculation. Default is 50.

    Returns
    -------
    None
        The function modifies the oracle_cc object in-place but does not return it.

    Notes
    -----
    - The function automatically calculates k as 2.5% of the total number of cells
    - The b_sight parameter is set to k * 8
    - The b_maxl parameter is set to k * 4
    - The imputation uses balanced KNN with parallel processing
    - The number of parallel jobs is controlled by config.GRN_N_JOBS

    Examples
    --------
    >>> import celloracle as co
    >>> oracle = co.Oracle()
    >>> # Load and prepare your data
    >>> run_KNN(oracle, n_comps=30)
    cell number is :5000
    Auto-selected k is :125
    """

    oracle_cc = oracle.copy()

    n_cell = oracle_cc.adata.shape[0]
    print(f"cell number is :{n_cell}")
    k = int(0.025 * n_cell)
    print(f"Auto-selected k is :{k}")

    oracle_cc.knn_imputation(
        n_pca_dims=n_comps,
        k=k,
        balanced=True,
        b_sight=k * 8,
        b_maxl=k * 4,
        n_jobs=config.GRN_N_JOBS,
    )

    return oracle_cc


def run_links(
    oracle: co.Oracle,
    cluster_column_name: str,
    p_cutoff: float = 0.001,
):
    """
    Calculate and filter gene regulatory network (GRN) links using CellOracle.
    This function computes regulatory links between transcription factors and target genes
    for each cluster in the dataset, filters them based on statistical significance, and
    generates diagnostic plots and network scores.
    Parameters
    ----------
    oracle : co.Oracle
        A CellOracle Oracle object containing the preprocessed single-cell data and
        base GRN information.
    cluster_column_name : str
        Name of the column in the oracle object that contains cluster assignments
        for cells. Used to compute cluster-specific GRN units.
    p_cutoff : float, optional
        P-value threshold for filtering significant regulatory links. Links with
        p-values above this cutoff will be removed. Default is 0.001.
    Returns
    -------
    links : co.Links
        A CellOracle Links object containing the filtered gene regulatory network
        with methods for further analysis and visualization.
    Notes
    -----
    - The function uses an alpha value of 10 for ridge regression regularization
    - Filtering is based on the absolute value of coefficients ('coef_abs')
    - A degree distribution plot is automatically saved to the figures directory
    - Network scores are computed to assess the quality of the inferred GRN
    """

    # Calculate GRN links
    links = oracle.get_links(
        cluster_name_for_GRN_unit=cluster_column_name,
        alpha=10,
        verbose_level=10,
        n_jobs=config.GRN_N_JOBS,
    )
    # Filter links based on p-value cutoff
    links.filter_links(p=p_cutoff, weight="coef_abs")
    # Calculate network scores like centrality, etc
    links.get_network_score()
    # Plot some stats over the network
    links.plot_degree_distributions(
        plot_model=True, save=config.FIGURES_DIR + "/grn_degree_distribution/"
    )
    links.plot_scores_as_rank(
        cluster=cluster_column_name,
        n_gene=20,
        save=config.FIGURES_DIR + "/grn_ranks/",
    )

    return links
