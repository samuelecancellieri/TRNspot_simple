#!/usr/bin/env python
"""
Complete Analysis Pipeline with Configuration
==============================================

This script demonstrates a complete single-cell analysis workflow
using TRNspot with proper configuration management.

Author: TRNspot Team
Date: October 2025
"""

import scanpy as sc
from trnspot import set_random_seed, config
from trnspot.preprocessing import perform_qc, plot_qc_violin, plot_qc_scatter


def main():
    """Run complete analysis pipeline"""

    # ========================================================================
    # 1. Setup and Configuration
    # ========================================================================
    print("=" * 70)
    print("TRNspot Complete Analysis Pipeline")
    print("=" * 70)

    # Set random seed for reproducibility
    print("\n[1/6] Setting up configuration...")
    set_random_seed(42)

    # Optionally customize config
    config.update_config(
        QC_MIN_GENES=200, QC_MIN_COUNTS=500, QC_PCT_MT_MAX=20.0, PLOT_DPI=300
    )

    # Print current configuration
    print("\nCurrent configuration:")
    print(f"  Random seed: {config.RANDOM_SEED}")
    print(f"  Min genes: {config.QC_MIN_GENES}")
    print(f"  Min counts: {config.QC_MIN_COUNTS}")
    print(f"  Max MT%: {config.QC_PCT_MT_MAX}")

    # ========================================================================
    # 2. Load Data
    # ========================================================================
    print("\n[2/6] Loading data...")

    # For demo purposes, use example dataset
    # Replace with: adata = sc.read_h5ad('your_data.h5ad')
    adata = sc.datasets.paul15()
    print(f"  Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

    # ========================================================================
    # 3. Quality Control
    # ========================================================================
    print("\n[3/6] Performing quality control...")

    adata = perform_qc(
        adata,
        min_genes=config.QC_MIN_GENES,
        min_counts=config.QC_MIN_COUNTS,
        pct_counts_mt_max=config.QC_PCT_MT_MAX,
        plot=True,
        save_plots=f"{config.FIGURES_DIR}/qc_histograms.png",
    )

    print(f"  After QC: {adata.n_obs} cells × {adata.n_vars} genes")

    # ========================================================================
    # 4. Additional QC Visualizations
    # ========================================================================
    print("\n[4/6] Generating additional QC plots...")

    # Violin plots
    plot_qc_violin(
        adata,
        figsize=config.PLOT_FIGSIZE_MEDIUM,
        save_path=f"{config.FIGURES_DIR}/qc_violin.png",
    )

    # Scatter plots
    plot_qc_scatter(
        adata,
        color_by="pct_counts_mt",
        figsize=config.PLOT_FIGSIZE_MEDIUM,
        save_path=f"{config.FIGURES_DIR}/qc_scatter.png",
    )

    # ========================================================================
    # 5. Normalization and Preprocessing
    # ========================================================================
    print("\n[5/6] Normalizing and preprocessing...")

    # Normalize to 10,000 counts per cell
    sc.pp.normalize_total(adata, target_sum=config.NORMALIZE_TARGET_SUM)

    # Log transform
    sc.pp.log1p(adata)

    # Identify highly variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes=config.HVGS_N_TOP_GENES)

    print(f"  Identified {adata.var['highly_variable'].sum()} highly variable genes")

    # Save raw data
    adata.raw = adata

    # Subset to HVGs
    adata = adata[:, adata.var["highly_variable"]].copy()

    # Scale data
    sc.pp.scale(adata, max_value=10)

    # PCA
    sc.tl.pca(adata, n_comps=config.PCA_N_COMPS, svd_solver=config.PCA_SVDSOLVE)

    print(f"  Computed {config.PCA_N_COMPS} principal components")

    # ========================================================================
    # 6. Save Results
    # ========================================================================
    print("\n[6/6] Saving results...")

    # Save processed data
    adata.write(f"{config.OUTPUT_DIR}/processed_data.h5ad")
    print(f"  Saved processed data to: {config.OUTPUT_DIR}/processed_data.h5ad")

    # Save configuration for reproducibility
    import json

    with open(f"{config.OUTPUT_DIR}/analysis_config.json", "w") as f:
        json.dump(config.get_config(), f, indent=2, default=str)
    print(f"  Saved configuration to: {config.OUTPUT_DIR}/analysis_config.json")

    # ========================================================================
    # Summary
    # ========================================================================
    print("\n" + "=" * 70)
    print("Analysis Complete!")
    print("=" * 70)
    print(f"Final dataset: {adata.n_obs} cells × {adata.n_vars} genes")
    print(f"Output directory: {config.OUTPUT_DIR}/")
    print(f"Figures directory: {config.FIGURES_DIR}/")
    print("\nFiles generated:")
    print(f"  - {config.OUTPUT_DIR}/processed_data.h5ad")
    print(f"  - {config.OUTPUT_DIR}/analysis_config.json")
    print(f"  - {config.FIGURES_DIR}/qc_histograms.png")
    print(f"  - {config.FIGURES_DIR}/qc_violin.png")
    print(f"  - {config.FIGURES_DIR}/qc_scatter.png")
    print("=" * 70)


if __name__ == "__main__":
    # Create output directories if they don't exist
    import os

    os.makedirs(config.OUTPUT_DIR, exist_ok=True)
    os.makedirs(config.FIGURES_DIR, exist_ok=True)

    # Run pipeline
    main()
