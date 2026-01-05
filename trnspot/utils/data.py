"""
Data utility functions for TRNspot
==================================

Functions for AnnData manipulation, loading, saving, and subsetting.
"""

import os
from typing import Optional, List, Tuple
import scanpy as sc
import pandas as pd
from anndata import AnnData


def ensure_categorical_obs(
    adata: AnnData,
    columns: Optional[List[str]] = None,
) -> AnnData:
    """
    Convert object/string columns in adata.obs to pandas Categorical type.

    This ensures consistent behavior during stratification, avoiding conflicts
    between string comparisons and numeric values.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    columns : list, optional
        Specific columns to convert. If None, converts all object/string columns
        and common clustering columns (leiden, louvain, cell_type, etc.).

    Returns
    -------
    AnnData
        AnnData with categorical columns in .obs.

    Examples
    --------
    >>> adata = ensure_categorical_obs(adata)
    >>> adata = ensure_categorical_obs(adata, columns=['cell_type', 'batch'])
    """
    # Common stratification/clustering columns to always convert if present
    default_categorical_cols = [
        "leiden",
        "louvain",
        "cell_type",
        "celltype",
        "cluster",
        "clusters",
        "batch",
        "sample",
        "condition",
    ]

    if columns is None:
        # Auto-detect: object dtype columns + known categorical columns
        columns_to_convert = []

        # Add object/string dtype columns
        for col in adata.obs.columns:
            if (
                adata.obs[col].dtype == "object"
                or adata.obs[col].dtype.name == "string"
            ):
                columns_to_convert.append(col)

        # Add default categorical columns if they exist
        for col in default_categorical_cols:
            if col in adata.obs.columns and col not in columns_to_convert:
                columns_to_convert.append(col)
    else:
        columns_to_convert = [c for c in columns if c in adata.obs.columns]

    converted = []
    for col in columns_to_convert:
        if not isinstance(adata.obs[col].dtype, pd.CategoricalDtype):
            adata.obs[col] = adata.obs[col].astype("category")
            converted.append(col)

    if converted:
        print(f"Converted to categorical: {', '.join(converted)}")

    return adata


def load_adata(
    input_path: str,
    verbose: bool = True,
) -> AnnData:
    """
    Load AnnData from various file formats.

    Parameters
    ----------
    input_path : str
        Path to input file (.h5ad, .h5, .loom, etc.)
    verbose : bool
        Whether to print loading information

    Returns
    -------
    AnnData
        Loaded AnnData object

    Raises
    ------
    ValueError
        If file format is not supported
    FileNotFoundError
        If file does not exist
    """
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"File not found: {input_path}")

    if verbose:
        print(f"Loading data from: {input_path}")

    if input_path.endswith(".h5ad"):
        adata = sc.read_h5ad(input_path)
    elif input_path.endswith(".h5"):
        adata = sc.read_10x_h5(input_path)
    elif input_path.endswith(".loom"):
        adata = sc.read_loom(input_path)
    elif input_path.endswith(".csv"):
        adata = sc.read_csv(input_path)
    elif input_path.endswith(".txt") or input_path.endswith(".tsv"):
        adata = sc.read_text(input_path)
    else:
        raise ValueError(
            f"Unsupported file format: {input_path}. "
            "Supported formats: .h5ad, .h5, .loom, .csv, .txt, .tsv"
        )

    if verbose:
        print(f"✓ Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

    return adata


def save_adata(
    adata: AnnData,
    output_path: str,
    verbose: bool = True,
) -> None:
    """
    Save AnnData to file.

    Parameters
    ----------
    adata : AnnData
        AnnData object to save
    output_path : str
        Path to output file (should end with .h5ad)
    verbose : bool
        Whether to print saving information
    """
    # Ensure directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    adata.write(output_path)

    if verbose:
        print(f"✓ Saved: {output_path}")


def subset_adata_by_cluster(
    adata: AnnData,
    cluster_key: str,
    cluster_value: str,
    copy: bool = True,
) -> AnnData:
    """
    Subset AnnData by a specific cluster value.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object
    cluster_key : str
        Column name in adata.obs containing cluster labels
    cluster_value : str
        Value to filter by
    copy : bool
        Whether to return a copy (default True)

    Returns
    -------
    AnnData
        Subset of input AnnData
    """
    if cluster_key not in adata.obs.columns:
        raise ValueError(f"Column '{cluster_key}' not found in adata.obs")

    mask = adata.obs[cluster_key] == cluster_value
    if copy:
        return adata[mask].copy()
    return adata[mask]


def stratify_adata(
    adata: AnnData,
    cluster_key: str,
    clusters: str = "all",
) -> Tuple[List[AnnData], List[str]]:
    """
    Split AnnData into multiple subsets based on cluster annotations.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object
    cluster_key : str
        Column name in adata.obs to stratify by
    clusters : str
        Comma-separated list of clusters, or "all" for all clusters

    Returns
    -------
    tuple
        (list of AnnData subsets, list of cluster names)
    """
    if cluster_key not in adata.obs.columns:
        raise ValueError(f"Column '{cluster_key}' not found in adata.obs")

    # Ensure categorical
    if not isinstance(adata.obs[cluster_key].dtype, pd.CategoricalDtype):
        adata.obs[cluster_key] = adata.obs[cluster_key].astype("category")

    # Get unique clusters
    if clusters == "all":
        unique_clusters = adata.obs[cluster_key].cat.categories
    else:
        requested = set([c.strip() for c in clusters.split(",") if c.strip()])
        unique_clusters = [
            c for c in adata.obs[cluster_key].cat.categories if str(c) in requested
        ]

    adata_list = []
    cluster_names = []

    for cluster in unique_clusters:
        adata_subset = subset_adata_by_cluster(adata, cluster_key, cluster)
        adata_list.append(adata_subset)
        cluster_names.append(str(cluster).replace(" ", "_"))

    return adata_list, cluster_names
