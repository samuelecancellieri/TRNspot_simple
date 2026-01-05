"""
CellOracle execution functions for TRNspot
==========================================

High-level functions for running CellOracle GRN inference.
"""

from typing import Optional, Tuple, Any
from anndata import AnnData

from .. import config
from ..utils.logging import log_step, log_error


def run_celloracle_preprocessing(
    adata: AnnData,
    cluster_key: str = "leiden",
    top_genes: Optional[int] = None,
    cell_downsample: int = 20000,
) -> AnnData:
    """
    Run CellOracle-specific preprocessing.

    This includes HVG selection, diffusion maps, PAGA, and graph drawing.

    Parameters
    ----------
    adata : AnnData
        Preprocessed input data
    cluster_key : str
        Key in adata.obs for cluster labels
    top_genes : int, optional
        Number of top genes for HVG selection
    cell_downsample : int
        Maximum number of cells to use

    Returns
    -------
    AnnData
        Preprocessed AnnData ready for CellOracle
    """
    from ..celloracle_processing import perform_grn_pre_processing

    log_step("CellOracle.Preprocessing", "STARTED", {"cluster_key": cluster_key})

    try:
        adata_prep = perform_grn_pre_processing(
            adata,
            cluster_key=cluster_key,
            top_genes=top_genes,
            cell_downsample=cell_downsample,
        )
        log_step("CellOracle.Preprocessing", "COMPLETED")
        return adata_prep

    except Exception as e:
        log_error("CellOracle.Preprocessing", e)
        raise


def run_celloracle_inference(
    adata: AnnData,
    cluster_key: str = "leiden",
    species: str = "human",
    embedding_name: str = "X_draw_graph_fa",
    raw_count_layer: str = "raw_counts",
    TG_to_TF_dictionary: Optional[str] = None,
) -> Tuple[Any, Any]:
    """
    Run CellOracle GRN inference.

    Parameters
    ----------
    adata : AnnData
        Preprocessed AnnData
    cluster_key : str
        Key in adata.obs for cluster labels
    species : str
        Species for base GRN ('human' or 'mouse')
    embedding_name : str
        Name of embedding in adata.obsm
    raw_count_layer : str
        Name of raw counts layer
    TG_to_TF_dictionary : str, optional
        Path to custom TF dictionary

    Returns
    -------
    tuple
        (Oracle object, Links object)
    """
    from ..celloracle_processing import (
        create_oracle_object,
        run_PCA,
        run_KNN,
        run_links,
        save_celloracle_results,
    )

    log_step("CellOracle.Inference", "STARTED", {"species": species})

    try:
        # Create Oracle object
        log_step("CellOracle.CreateOracle", "STARTED")
        oracle = create_oracle_object(
            adata,
            cluster_column_name=cluster_key,
            embedding_name=embedding_name,
            raw_count_layer=raw_count_layer,
            species=species,
            TG_to_TF_dictionary=TG_to_TF_dictionary,
        )
        log_step("CellOracle.CreateOracle", "COMPLETED")

        # PCA
        log_step("CellOracle.PCA", "STARTED")
        oracle, n_comps = run_PCA(oracle)
        log_step("CellOracle.PCA", "COMPLETED", {"n_comps": n_comps})

        # KNN imputation
        log_step("CellOracle.KNN", "STARTED")
        oracle = run_KNN(oracle, n_comps=n_comps)
        log_step("CellOracle.KNN", "COMPLETED")

        # GRN inference
        log_step("CellOracle.Links", "STARTED")
        links = run_links(oracle)
        log_step("CellOracle.Links", "COMPLETED")

        # Save results
        save_celloracle_results(oracle, links)

        log_step("CellOracle.Inference", "COMPLETED")
        return oracle, links

    except Exception as e:
        log_error("CellOracle.Inference", e)
        raise


def run_celloracle_pipeline(
    adata: AnnData,
    cluster_key: str = "leiden",
    species: str = "human",
    embedding_name: str = "X_draw_graph_fa",
    raw_count_layer: str = "raw_counts",
    TG_to_TF_dictionary: Optional[str] = None,
    skip: bool = False,
) -> Optional[Tuple[Any, Any]]:
    """
    Run complete CellOracle pipeline: preprocessing + inference.

    Parameters
    ----------
    adata : AnnData
        Input AnnData
    cluster_key : str
        Key for cluster labels
    species : str
        Species for base GRN
    embedding_name : str
        Embedding for analysis
    raw_count_layer : str
        Raw counts layer name
    TG_to_TF_dictionary : str, optional
        Path to custom TF dictionary
    skip : bool
        Skip CellOracle analysis

    Returns
    -------
    tuple or None
        (Oracle, Links) objects or None if skipped
    """
    if skip:
        log_step("CellOracle", "SKIPPED", {"reason": "skip=True"})
        print("⊘ Skipping CellOracle analysis")
        return None

    log_step("CellOracle.Pipeline", "STARTED")

    try:
        # Preprocessing
        adata_prep = run_celloracle_preprocessing(adata, cluster_key=cluster_key)

        # Inference
        oracle, links = run_celloracle_inference(
            adata_prep,
            cluster_key=cluster_key,
            species=species,
            embedding_name=embedding_name,
            raw_count_layer=raw_count_layer,
            TG_to_TF_dictionary=TG_to_TF_dictionary,
        )

        log_step("CellOracle.Pipeline", "COMPLETED")
        return oracle, links

    except ImportError:
        log_step("CellOracle", "SKIPPED", {"reason": "CellOracle not installed"})
        print("⊘ Skipping CellOracle (not installed)")
        return None

    except Exception as e:
        log_error("CellOracle.Pipeline", e)
        raise
