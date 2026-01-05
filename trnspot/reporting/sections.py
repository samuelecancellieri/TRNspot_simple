"""
Report section builders for TRNspot
===================================

Functions to create individual report sections.
"""

import os
from typing import Optional, Dict, Any, List
from pathlib import Path
import glob

from .generator import ReportSection


def _find_figures(
    output_dir: str,
    subdirs: List[str],
    patterns: List[str] = ["*.png", "*.jpg", "*.svg"],
) -> List[str]:
    """Find figure files in specified subdirectories."""
    figures = []
    for subdir in subdirs:
        search_dir = os.path.join(output_dir, subdir)
        if not os.path.exists(search_dir):
            continue
        for pattern in patterns:
            figures.extend(glob.glob(os.path.join(search_dir, pattern)))
    return sorted(figures)[:10]  # Limit to 10 figures per section


def create_qc_section(
    adata,
    output_dir: str,
) -> ReportSection:
    """
    Create Quality Control report section.

    Parameters
    ----------
    adata : AnnData
        AnnData object with QC metrics
    output_dir : str
        Output directory with figures

    Returns
    -------
    ReportSection
        QC section for the report
    """
    # Collect metrics
    metrics = {}

    if "n_genes_by_counts" in adata.obs.columns:
        metrics["Median Genes/Cell"] = int(adata.obs["n_genes_by_counts"].median())
    if "total_counts" in adata.obs.columns:
        metrics["Median Counts/Cell"] = int(adata.obs["total_counts"].median())
    if "pct_counts_mt" in adata.obs.columns:
        metrics["Mean MT %"] = round(adata.obs["pct_counts_mt"].mean(), 2)

    metrics["Total Cells"] = adata.n_obs
    metrics["Total Genes"] = adata.n_vars

    # Content
    content = """
    <p>Quality control metrics were calculated for all cells, including the number of 
    genes detected, total counts, and percentage of mitochondrial reads.</p>
    
    <h3>QC Thresholds Applied</h3>
    <ul>
        <li>Cells with too few genes were filtered out</li>
        <li>Cells with too few counts were removed</li>
        <li>Cells with high mitochondrial content were excluded</li>
        <li>Genes expressed in too few cells were removed</li>
    </ul>
    """

    # Find QC figures
    figures = _find_figures(output_dir, ["figures/qc"])

    return ReportSection(
        title="Quality Control",
        content=content,
        figures=figures,
        metrics=metrics,
    )


def create_preprocessing_section(
    adata,
    output_dir: str,
) -> ReportSection:
    """
    Create Preprocessing report section.

    Parameters
    ----------
    adata : AnnData
        Preprocessed AnnData object
    output_dir : str
        Output directory

    Returns
    -------
    ReportSection
        Preprocessing section for the report
    """
    metrics = {
        "Final Cells": adata.n_obs,
        "Final Genes": adata.n_vars,
    }

    # Check for clustering
    if "leiden" in adata.obs.columns:
        n_clusters = len(adata.obs["leiden"].unique())
        metrics["Clusters (Leiden)"] = n_clusters
    elif "louvain" in adata.obs.columns:
        n_clusters = len(adata.obs["louvain"].unique())
        metrics["Clusters (Louvain)"] = n_clusters

    # Check for HVGs
    if "highly_variable" in adata.var.columns:
        n_hvg = adata.var["highly_variable"].sum()
        metrics["Highly Variable Genes"] = int(n_hvg)

    # Check embeddings
    embeddings = []
    if "X_pca" in adata.obsm:
        embeddings.append(f"PCA ({adata.obsm['X_pca'].shape[1]} components)")
    if "X_umap" in adata.obsm:
        embeddings.append("UMAP")
    if "X_draw_graph_fa" in adata.obsm:
        embeddings.append("Force-directed graph")

    content = f"""
    <p>The data was preprocessed using standard single-cell RNA-seq workflows including 
    normalization, highly variable gene selection, dimensionality reduction, and clustering.</p>
    
    <h3>Preprocessing Steps</h3>
    <ol>
        <li><strong>Normalization:</strong> Library size normalization followed by log transformation</li>
        <li><strong>Feature Selection:</strong> Identification of highly variable genes</li>
        <li><strong>Dimensionality Reduction:</strong> {', '.join(embeddings) if embeddings else 'Not computed'}</li>
        <li><strong>Clustering:</strong> Graph-based clustering for cell type identification</li>
    </ol>
    """

    return ReportSection(
        title="Preprocessing",
        content=content,
        metrics=metrics,
    )


def create_celloracle_section(
    celloracle_result,
    output_dir: str,
) -> ReportSection:
    """
    Create CellOracle GRN analysis report section.

    Parameters
    ----------
    celloracle_result : tuple
        (Oracle, Links) objects from CellOracle
    output_dir : str
        Output directory

    Returns
    -------
    ReportSection
        CellOracle section for the report
    """
    metrics = {}
    tables = []

    if celloracle_result is not None:
        oracle, links = celloracle_result

        # Get metrics from links if available
        if hasattr(links, "merged_score") and links.merged_score is not None:
            merged_scores = links.merged_score
            if "cluster" in merged_scores.columns:
                n_clusters = len(merged_scores["cluster"].unique())
                metrics["Clusters Analyzed"] = n_clusters

            # Count TFs
            n_tfs = (
                len(merged_scores.index.unique())
                if hasattr(merged_scores.index, "unique")
                else 0
            )
            metrics["Transcription Factors"] = n_tfs

        if hasattr(links, "filtered_links") and links.filtered_links is not None:
            n_links = sum(
                len(v)
                for v in links.filtered_links.values()
                if isinstance(v, (list, dict))
            )
            metrics["Regulatory Links"] = n_links

    content = """
    <p>Gene regulatory network (GRN) inference was performed using CellOracle, which 
    predicts transcription factor (TF) regulatory relationships from single-cell data.</p>
    
    <h3>Analysis Pipeline</h3>
    <ol>
        <li><strong>Preprocessing:</strong> Diffusion map computation and PAGA graph construction</li>
        <li><strong>Oracle Creation:</strong> Integration of expression data with base GRN</li>
        <li><strong>GRN Inference:</strong> Network inference using machine learning models</li>
        <li><strong>Network Analysis:</strong> Calculation of centrality scores and hub TFs</li>
    </ol>
    
    <h3>Output Files</h3>
    <ul>
        <li><code>celloracle/oracle_object.celloracle.oracle</code> - Oracle object</li>
        <li><code>celloracle/oracle_object.celloracle.links</code> - Links object</li>
        <li><code>celloracle/grn_merged_scores.csv</code> - Network centrality scores</li>
        <li><code>celloracle/grn_filtered_links.pkl</code> - Filtered regulatory links</li>
    </ul>
    """

    # Find GRN figures
    figures = _find_figures(
        output_dir, ["figures/grn", "figures/grn/grn_deep_analysis"]
    )

    return ReportSection(
        title="CellOracle GRN Analysis",
        content=content,
        figures=figures,
        metrics=metrics,
        tables=tables,
    )


def create_hotspot_section(
    hotspot_result,
    output_dir: str,
) -> ReportSection:
    """
    Create Hotspot gene module report section.

    Parameters
    ----------
    hotspot_result : Hotspot
        Hotspot object with results
    output_dir : str
        Output directory

    Returns
    -------
    ReportSection
        Hotspot section for the report
    """
    metrics = {}

    if hotspot_result is not None:
        # Get results
        if hasattr(hotspot_result, "results") and hotspot_result.results is not None:
            n_tested = len(hotspot_result.results)
            n_significant = (hotspot_result.results["FDR"] < 0.05).sum()
            metrics["Genes Tested"] = n_tested
            metrics["Significant Genes"] = int(n_significant)

        if hasattr(hotspot_result, "modules") and hotspot_result.modules is not None:
            n_modules = len([m for m in hotspot_result.modules.unique() if m != -1])
            metrics["Gene Modules"] = n_modules

    content = """
    <p>Hotspot analysis was performed to identify spatially autocorrelated gene modules - 
    groups of genes that show coordinated expression patterns across cells.</p>
    
    <h3>Analysis Steps</h3>
    <ol>
        <li><strong>KNN Graph:</strong> Construction of cell-cell similarity graph</li>
        <li><strong>Autocorrelation:</strong> Testing genes for spatial autocorrelation</li>
        <li><strong>Local Correlations:</strong> Computing pairwise gene correlations</li>
        <li><strong>Module Detection:</strong> Hierarchical clustering of correlated genes</li>
    </ol>
    
    <h3>Output Files</h3>
    <ul>
        <li><code>hotspot/autocorrelation_results.csv</code> - Gene autocorrelation statistics</li>
        <li><code>hotspot/significant_genes.csv</code> - Significantly autocorrelated genes</li>
        <li><code>hotspot/gene_modules.csv</code> - Gene-to-module assignments</li>
        <li><code>hotspot/hotspot_module_scores.csv</code> - Module scores per cell</li>
    </ul>
    """

    # Find Hotspot figures
    figures = _find_figures(output_dir, ["figures/hotspot"])

    return ReportSection(
        title="Hotspot Gene Modules",
        content=content,
        figures=figures,
        metrics=metrics,
    )


def create_summary_section(
    adata=None,
    celloracle_result=None,
    hotspot_result=None,
    output_dir: str = "",
) -> ReportSection:
    """
    Create analysis summary section.

    Parameters
    ----------
    adata : AnnData, optional
        Processed AnnData
    celloracle_result : tuple, optional
        CellOracle results
    hotspot_result : Hotspot, optional
        Hotspot results
    output_dir : str
        Output directory

    Returns
    -------
    ReportSection
        Summary section for the report
    """
    # Build summary content
    analyses_completed = []
    if adata is not None:
        analyses_completed.append("Data preprocessing and clustering")
    if celloracle_result is not None:
        analyses_completed.append("CellOracle GRN inference")
    if hotspot_result is not None:
        analyses_completed.append("Hotspot gene module analysis")

    metrics = {}
    if adata is not None:
        metrics["Final Cells"] = adata.n_obs
        metrics["Final Genes"] = adata.n_vars

    analyses_html = "\n".join(f"<li>{a}</li>" for a in analyses_completed)

    content = f"""
    <h3>Analyses Completed</h3>
    <ul>{analyses_html}</ul>
    
    <h3>Output Directory</h3>
    <p><code>{output_dir}</code></p>
    
    <h3>Key Output Files</h3>
    <ul>
        <li><code>preprocessed_adata.h5ad</code> - Preprocessed AnnData object</li>
        <li><code>clustered_adata.h5ad</code> - Clustered AnnData object</li>
        <li><code>analysis_summary.txt</code> - Text summary of analysis</li>
        <li><code>report.html</code> - This HTML report</li>
    </ul>
    
    <h3>Next Steps</h3>
    <ul>
        <li>Review QC metrics and adjust thresholds if needed</li>
        <li>Explore clustering results and annotate cell types</li>
        <li>Investigate top transcription factors from GRN analysis</li>
        <li>Examine gene modules for biological interpretation</li>
    </ul>
    """

    return ReportSection(
        title="Summary",
        content=content,
        metrics=metrics,
    )


def create_stratification_section(
    stratification_names: List[str],
    output_dir: str,
) -> ReportSection:
    """
    Create stratified analysis summary section.

    Parameters
    ----------
    stratification_names : list
        Names of stratifications analyzed
    output_dir : str
        Output directory

    Returns
    -------
    ReportSection
        Stratification section for the report
    """
    n_strats = len(stratification_names)

    strat_list = "\n".join(f"<li>{name}</li>" for name in stratification_names)

    content = f"""
    <p>The data was stratified into {n_strats} subsets for independent analysis. 
    Each stratification was processed through the full analysis pipeline.</p>
    
    <h3>Stratifications Analyzed</h3>
    <ul>{strat_list}</ul>
    
    <h3>Output Structure</h3>
    <p>Results for each stratification are saved in:</p>
    <code>{output_dir}/stratified_analysis/&lt;stratification_name&gt;/</code>
    """

    return ReportSection(
        title="Stratified Analysis",
        content=content,
        metrics={"Stratifications": n_strats},
    )
