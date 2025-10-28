#!/usr/bin/env python
"""
Complete TRNspot Analysis Pipeline
====================================
This script runs the full analysis workflow:
1. Data loading and quality control
2. Normalization and preprocessing
3. Dimensionality reduction and clustering
4. GRN inference with CellOracle
5. Gene module identification with Hotspot

Usage:
    python complete_pipeline.py --input data/paul15/paul15.h5 --output output --figures figures

Or with default paths:
    python complete_pipeline.py
"""

import os
import shutil
import argparse
import sys
from datetime import datetime
import scanpy as sc
import pickle
import hashlib
import json

# Import TRNspot modules
from trnspot import set_random_seed, set_scanpy_settings, config
from trnspot.preprocessing import (
    perform_qc,
    perform_normalization,
    perform_grn_pre_processing,
)


def compute_input_hash(input_path, **params):
    """Compute hash of input file and parameters for checkpoint verification."""
    hash_obj = hashlib.md5()

    # Hash input file path
    if input_path:
        hash_obj.update(str(input_path).encode())
        # If file exists, hash its modification time
        if os.path.exists(input_path):
            mtime = os.path.getmtime(input_path)
            hash_obj.update(str(mtime).encode())

    # Hash relevant parameters
    for key, value in sorted(params.items()):
        hash_obj.update(f"{key}:{value}".encode())

    return hash_obj.hexdigest()


def write_checkpoint(log_dir, step_name, input_hash, **metadata):
    """Write checkpoint log file for a completed step."""
    os.makedirs(log_dir, exist_ok=True)
    checkpoint_file = os.path.join(log_dir, f"{step_name}.checkpoint")

    checkpoint_data = {
        "step_name": step_name,
        "input_hash": input_hash,
        "timestamp": datetime.now().isoformat(),
        "metadata": metadata,
    }

    with open(checkpoint_file, "w") as f:
        json.dump(checkpoint_data, f, indent=2)


def check_checkpoint(log_dir, step_name, input_hash):
    """Check if step has already been completed with same input."""
    checkpoint_file = os.path.join(log_dir, f"{step_name}.checkpoint")

    if not os.path.exists(checkpoint_file):
        return False

    try:
        with open(checkpoint_file, "r") as f:
            checkpoint_data = json.load(f)

        # Check if input hash matches
        if checkpoint_data.get("input_hash") == input_hash:
            print(f"  ⏭ Checkpoint found - skipping {step_name}")
            print(f"     (completed at {checkpoint_data.get('timestamp')})")
            return True
    except Exception as e:
        print(f"  ⚠ Error reading checkpoint: {e}")

    return False


def setup_directories(output_dir, figures_dir, debug=False):
    """Create necessary output directories."""
    directories = [
        output_dir,
        figures_dir,
        f"{output_dir}/logs",
        f"{output_dir}/celloracle",
        f"{output_dir}/hotspot",
        f"{figures_dir}/qc",
        f"{figures_dir}/grn",
        f"{figures_dir}/hotspot",
    ]

    if debug:
        for directory in directories:
            shutil.rmtree(directory, ignore_errors=True)
            print(f"✓ Removed directory: {directory}")

    for directory in directories:
        os.makedirs(directory, exist_ok=True)

    print(f"✓ Created output directories")


def load_data(input_path):
    """Load data from file or use example dataset."""
    print(f"\n{'='*70}")
    print("STEP 1: Data Loading")
    print(f"{'='*70}")

    if input_path and os.path.exists(input_path):
        print(f"Loading data from: {input_path}")
        if input_path.endswith(".h5ad"):
            adata = sc.read_h5ad(input_path)
        elif input_path.endswith(".h5"):
            adata = sc.read_10x_h5(input_path)
        else:
            raise ValueError(f"Unsupported file format: {input_path}")
    else:
        print("No input file specified or file not found.")
        print("Loading example dataset: Paul et al. 2015 (hematopoietic cells)")
        adata = sc.datasets.paul15()

    print(f"✓ Loaded: {adata.n_obs} cells × {adata.n_vars} genes")
    return adata


def preprocessing_pipeline(adata, name=None, skip_qc=False, log_dir=None):
    """Run preprocessing pipeline."""
    print(f"\n{'='*70}")
    print("STEP 2: Quality Control and Preprocessing")
    print(f"{'='*70}")

    # Check checkpoint for preprocessing
    step_hash = compute_input_hash(
        None,
        n_obs=adata.n_obs,
        n_vars=adata.n_vars,
        min_genes=config.QC_MIN_GENES,
        min_counts=config.QC_MIN_COUNTS,
        pct_mt_max=config.QC_PCT_MT_MAX,
        skip_qc=skip_qc,
    )

    checkpoint_file = os.path.join(config.OUTPUT_DIR, "preprocessed_adata.h5ad")
    if (
        log_dir
        and check_checkpoint(log_dir, "preprocessing", step_hash)
        and os.path.exists(checkpoint_file)
    ):
        print(f"  Loading preprocessed data from: {checkpoint_file}")
        return sc.read_h5ad(checkpoint_file)

    if not skip_qc:
        print("\n[2.1] Performing quality control...")
        adata = perform_qc(
            adata,
            min_genes=config.QC_MIN_GENES,
            min_counts=config.QC_MIN_COUNTS,
            pct_counts_mt_max=config.QC_PCT_MT_MAX,
            save_plots=name,
        )
        print(f"  After QC: {adata.n_obs} cells × {adata.n_vars} genes")
    else:
        print("\n[2.1] Skipping QC (already performed)")

    # Normalization
    print("\n[2.2] Normalizing data...")
    adata = perform_normalization(adata)
    print("  ✓ Normalization complete")

    # Save checkpoint
    if log_dir:
        adata.write(checkpoint_file)
        write_checkpoint(
            log_dir, "preprocessing", step_hash, n_obs=adata.n_obs, n_vars=adata.n_vars
        )

    # Store raw counts for Hotspot
    # if adata.raw is not None:
    #     adata.layers["raw_counts"] = adata.raw.X.copy()
    #     print("  ✓ Stored raw counts in layer 'raw_counts'")

    return adata


def stratification_pipeline(adata, cluster_key_stratification=None):
    """Perform stratification based on specified clustering key."""
    if cluster_key_stratification and cluster_key_stratification in adata.obs.columns:
        print(f"\n{'='*70}")
        print("STEP 2.5: Stratification by Clusters")
        print(f"{'='*70}")

        adata_list = list()
        adata_stratification_list = list()
        unique_clusters = adata.obs[cluster_key_stratification].unique()
        print(
            f"Stratifying data into {len(unique_clusters)} clusters based on '{cluster_key_stratification}'"
        )

        for cluster in unique_clusters:
            adata_cluster = adata[
                adata.obs[cluster_key_stratification] == cluster
            ].copy()
            adata_list.append(adata_cluster)
            adata_stratification_list.append(cluster.replace(" ", "_"))
            print(
                f"  ✓ Cluster '{cluster}': {adata_cluster.n_obs} cells × {adata_cluster.n_vars} genes"
            )

        return adata_list, adata_stratification_list
    else:
        print("\nNo stratification performed (no valid clustering key provided).")
        return [], []  # Return empty lists for compatibility


def dimensionality_reduction_clustering(adata, cluster_key="leiden", log_dir=None):
    """Perform dimensionality reduction and clustering."""
    print(f"\n{'='*70}")
    print("STEP 3: Dimensionality Reduction and Clustering")
    print(f"{'='*70}")

    # Check checkpoint for clustering
    step_hash = compute_input_hash(
        None,
        n_obs=adata.n_obs,
        n_vars=adata.n_vars,
        cluster_key=cluster_key,
        top_genes=config.HVGS_N_TOP_GENES,
        n_neighbors=config.NEIGHBORS_N_NEIGHBORS,
        n_pcs=config.NEIGHBORS_N_PCS,
    )

    checkpoint_file = f"{config.OUTPUT_DIR}/clustered_adata.h5ad"
    if (
        log_dir
        and check_checkpoint(log_dir, "clustering", step_hash)
        and os.path.exists(checkpoint_file)
    ):
        print(f"  Loading clustered data from: {checkpoint_file}")
        return sc.read_h5ad(checkpoint_file)

    # GRN preprocessing (includes HVG, PCA, diffusion map, PAGA)
    print("\n[3.1] Running GRN preprocessing pipeline...")
    adata = perform_grn_pre_processing(
        adata,
        cluster_key=cluster_key,  # Will create leiden clustering
        top_genes=config.HVGS_N_TOP_GENES,
        n_neighbors=config.NEIGHBORS_N_NEIGHBORS,
        n_pcs=config.NEIGHBORS_N_PCS,
    )

    # UMAP for visualization
    print("\n[3.2] Computing UMAP embedding...")
    sc.tl.umap(adata)
    print("  ✓ UMAP complete")

    # Leiden clustering if not already done
    if "leiden" not in adata.obs.columns:
        print("\n[3.3] Performing Leiden clustering...")
        sc.tl.leiden(adata, resolution=1.0)
        n_clusters = len(adata.obs["leiden"].unique())
        print(f"  ✓ Identified {n_clusters} clusters")
    else:
        n_clusters = len(adata.obs["leiden"].unique())
        print(f"\n[3.3] Using existing Leiden clustering ({n_clusters} clusters)")

    # Save checkpoint
    if log_dir:
        adata.write(checkpoint_file)
        write_checkpoint(
            log_dir,
            "clustering",
            step_hash,
            n_obs=adata.n_obs,
            n_vars=adata.n_vars,
            n_clusters=n_clusters,
        )

    return adata


def celloracle_pipeline(
    adata,
    cluster_key="leiden",
    species="human",
    embedding_name="X_draw_graph_fa",
    raw_count_layer="raw_counts",
    TG_to_TF_dictionary=None,
    skip_celloracle=False,
    log_dir=None,
):
    """Run CellOracle GRN inference pipeline."""
    print(f"\n{'='*70}")
    print("STEP 4: CellOracle GRN Inference")
    print(f"{'='*70}")

    if skip_celloracle:
        print("⊘ Skipping CellOracle analysis (--skip-celloracle)")
        return None

    # Check checkpoint for CellOracle
    step_hash = compute_input_hash(
        None,
        n_obs=adata.n_obs,
        cluster_key=cluster_key,
        species=species,
        embedding_name=embedding_name,
    )

    oracle_file = f"{config.OUTPUT_DIR}/celloracle/oracle_object.celloracle.oracle"
    links_file = f"{config.OUTPUT_DIR}/celloracle/oracle_object.celloracle.links"

    if log_dir and check_checkpoint(log_dir, "celloracle", step_hash):
        if os.path.exists(oracle_file) and os.path.exists(links_file):
            try:
                from trnspot.celloracle_processing import load_celloracle_results

                print(f"  Loading CellOracle results from checkpoint...")
                oracle, links = load_celloracle_results(
                    oracle_path=oracle_file, links_path=links_file
                )
                return oracle, links
            except Exception as e:
                print(f"  ⚠ Error loading checkpoint: {e}")
                print(f"  Re-running CellOracle analysis...")

    try:
        from trnspot.celloracle_processing import (
            create_oracle_object,
            run_PCA,
            run_KNN,
            run_links,
            save_celloracle_results,
        )

        print("\n[4.1] Creating Oracle object...")
        oracle = create_oracle_object(
            adata,
            cluster_column_name=cluster_key,
            embedding_name=embedding_name,
            raw_count_layer=raw_count_layer,
            species=species,
            TG_to_TF_dictionary=TG_to_TF_dictionary,
        )
        print("  ✓ Oracle object created")

        print("\n[4.2] Running PCA on Oracle object...")
        oracle, n_comps = run_PCA(oracle)
        print("  ✓ PCA complete")

        print("\n[4.3] Running KNN imputation...")
        oracle = run_KNN(oracle, n_comps=n_comps)
        print("  ✓ KNN imputation complete")

        print("\n[4.4] Inferring GRN links...")
        links = run_links(
            oracle,
            cluster_column_name=cluster_key,
            p_cutoff=0.001,
        )
        print("  ✓ GRN inference complete")

        print("\n[4.5] Saving CellOracle results...")
        save_celloracle_results(oracle, links)
        print("  ✓ CellOracle results saved")

        # Save checkpoint
        if log_dir:
            write_checkpoint(
                log_dir,
                "celloracle",
                step_hash,
                n_obs=adata.n_obs,
                cluster_key=cluster_key,
            )

        return oracle, links

    except ImportError:
        print("\n⚠ CellOracle not installed. Skipping GRN inference.")
        print("  To install: pip install celloracle")
        return None
    except Exception as e:
        print(f"\n⚠ Error during CellOracle analysis: {e}")
        print("  Continuing with remaining analysis...")
        return None


def hotspot_pipeline(
    adata,
    layer_key="raw_counts",
    embedding_key="X_umap",
    skip_hotspot=False,
    log_dir=None,
):
    """Run Hotspot gene module identification pipeline."""
    print(f"\n{'='*70}")
    print("STEP 5: Hotspot Gene Module Identification")
    print(f"{'='*70}")

    if skip_hotspot:
        print("⊘ Skipping Hotspot analysis (--skip-hotspot)")
        return None

    # Check checkpoint for Hotspot
    step_hash = compute_input_hash(
        None,
        n_obs=adata.n_obs,
        top_genes=config.HOTSPOT_TOP_GENES,
        embedding_key=embedding_key,
        fdr_threshold=config.HOTSPOT_FDR_THRESHOLD,
    )

    hotspot_file = f"{config.OUTPUT_DIR}/hotspot/hotspot_object.pkl"
    if (
        log_dir
        and check_checkpoint(log_dir, "hotspot", step_hash)
        and os.path.exists(hotspot_file)
    ):
        print(f"  Loading Hotspot results from checkpoint...")
        with open(hotspot_file, "rb") as f:
            hotspot_obj = pickle.load(f)
        return hotspot_obj

    try:
        from trnspot.hotspot_processing import (
            create_hotspot_object,
            run_hotspot_analysis,
        )

        print("\n[5.1] Creating Hotspot object...")
        hotspot_obj = create_hotspot_object(
            adata,
            top_genes=config.HOTSPOT_TOP_GENES,
            layer_key=layer_key,
            model="danb",
            embedding_key=embedding_key,
            normalization_key="n_counts",
        )
        print("  ✓ Hotspot object created")

        print("\n[5.2] Running Hotspot analysis...")
        print("  (This may take several minutes...)")
        hotspot_obj = run_hotspot_analysis(hotspot_obj)
        print("  ✓ Hotspot analysis complete")

        # Get results summary
        autocorr_results = hotspot_obj.results
        significant_genes = autocorr_results[
            autocorr_results.FDR < config.HOTSPOT_FDR_THRESHOLD
        ]

        print(f"\n  Analysis Summary:")
        print(f"    Total genes analyzed: {len(autocorr_results)}")
        print(f"    Significant genes: {len(significant_genes)}")

        if hasattr(hotspot_obj, "modules"):
            n_modules = len(hotspot_obj.modules.unique())
            print(f"    Gene modules identified: {n_modules}")

        # Save checkpoint
        if log_dir:
            write_checkpoint(
                log_dir,
                "hotspot",
                step_hash,
                n_genes=len(autocorr_results),
                n_significant=len(significant_genes),
            )

        return hotspot_obj

    except ImportError:
        print("\n⚠ Hotspot not installed. Skipping module identification.")
        print("  To install: pip install hotspot-sc")
        return None
    except Exception as e:
        print(f"\n⚠ Error during Hotspot analysis: {e}")
        print("  Continuing with remaining analysis...")
        return None


def generate_summary(adata, celloracle_result, hotspot_result, start_time, output_dir):
    """Generate analysis summary report."""
    print(f"\n{'='*70}")
    print("STEP 6: Generating Analysis Summary")
    print(f"{'='*70}")

    end_time = datetime.now()
    duration = end_time - start_time

    summary = []
    summary.append("=" * 70)
    summary.append("TRNspot Complete Pipeline - Analysis Summary")
    summary.append("=" * 70)
    summary.append(f"\nAnalysis completed at: {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    summary.append(f"Total runtime: {duration}")
    summary.append(f"\n{'='*70}")
    summary.append("Dataset Information")
    summary.append("=" * 70)
    summary.append(f"Final dataset: {adata.n_obs} cells × {adata.n_vars} genes")

    if "leiden" in adata.obs.columns:
        n_clusters = len(adata.obs["leiden"].unique())
        summary.append(f"Clusters identified: {n_clusters}")

    # CellOracle results
    summary.append(f"\n{'='*70}")
    summary.append("CellOracle GRN Inference")
    summary.append("=" * 70)
    if celloracle_result:
        oracle, links = celloracle_result
        summary.append(f"Status: ✓ Completed")
        summary.append(
            f"Oracle object: {output_dir}/celloracle/oracle_object.celloracle.oracle"
        )
        summary.append(f"GRN links: {output_dir}/celloracle/grn_links.celloracle.links")
    else:
        summary.append(f"Status: ⊘ Skipped or failed")

    # Hotspot results
    summary.append(f"\n{'='*70}")
    summary.append("Hotspot Module Identification")
    summary.append("=" * 70)
    if hotspot_result:
        autocorr_results = hotspot_result.results
        significant = autocorr_results[
            autocorr_results.FDR < config.HOTSPOT_FDR_THRESHOLD
        ]
        summary.append(f"Status: ✓ Completed")
        summary.append(f"Genes analyzed: {len(autocorr_results)}")
        summary.append(f"Significant genes: {len(significant)}")
        if hasattr(hotspot_result, "modules"):
            n_modules = len(hotspot_result.modules.unique())
            summary.append(f"Modules identified: {n_modules}")
        summary.append(f"Results: {output_dir}/hotspot/")
    else:
        summary.append(f"Status: ⊘ Skipped or failed")

    # Output files
    summary.append(f"\n{'='*70}")
    summary.append("Output Files")
    summary.append("=" * 70)
    summary.append(f"Preprocessed data: {output_dir}/preprocessed_adata.h5ad")
    summary.append(f"Figures directory: {config.FIGURES_DIR}/")
    summary.append(f"Results directory: {output_dir}/")

    # summary.append(f"\n{'='*70}")
    # summary.append("Next Steps")
    # summary.append("=" * 70)
    # summary.append("1. Explore cell clusters and marker genes")
    # summary.append("2. Analyze gene regulatory networks (CellOracle results)")
    # summary.append("3. Examine gene modules and their functions (Hotspot results)")
    # summary.append("4. Perform downstream analysis (GO enrichment, pathway analysis)")
    # summary.append("5. Generate publication-quality visualizations")

    summary_text = "\n".join(summary)
    print("\n" + summary_text)

    # Save summary to file
    summary_path = f"{output_dir}/analysis_summary.txt"
    with open(summary_path, "w") as f:
        f.write(summary_text)
    print(f"\n✓ Summary saved to: {summary_path}")

    return summary_text


def track_files(output_dir):
    """Track output files generated during the pipeline."""
    tracked_files = []
    for root, dirs, files in os.walk(output_dir):
        for file in files:
            tracked_files.append(os.path.join(root, file))

    tracked_files_path = os.path.join(output_dir, "tracked_files.txt")
    if os.path.exists(tracked_files_path):
        os.remove(tracked_files_path)
    with open(tracked_files_path, "w") as f:
        for file in tracked_files:
            f.write(f"{file}\n")
    print(f"\n✓ Tracked output files saved to: {tracked_files_path}")

    return tracked_files


def grn_deep_analysis_pipeline(grn_score_path):
    """Run GRN deep analysis using tracked output files."""
    from trnspot.grn_deep_analysis import (
        process_single_score_file,
        plot_scatter_scores,
        plot_compare_cluster_scores,
        plot_difference_cluster_scores,
    )

    print(f"\n{'='*70}")
    print("STEP 7: GRN Deep Analysis")
    print(f"{'='*70}")

    os.makedirs(f"{config.OUTPUT_DIR}/grn_deep_analysis", exist_ok=True)
    os.makedirs(f"{config.FIGURES_DIR_GRN}/grn_deep_analysis", exist_ok=True)

    score_df = process_single_score_file(grn_score_path)
    print("  ✓ Processed GRN score file")
    plot_compare_cluster_scores(score_df)
    print("  ✓ Generated cluster comparison plots")
    plot_difference_cluster_scores(score_df)
    print("  ✓ Generated cluster difference plots")
    plot_scatter_scores(score_df)
    print("  ✓ Generated scatter plots")


def main():
    """Main pipeline execution."""
    parser = argparse.ArgumentParser(
        description="Run complete TRNspot analysis pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run with example dataset
  python complete_pipeline.py
  
  # Run with custom input
  python complete_pipeline.py --input data/my_data.h5ad
  
  # Skip specific analyses
  python complete_pipeline.py --skip-celloracle --skip-hotspot
  
  # Custom configuration
  python complete_pipeline.py --seed 123 --n-jobs 16
        """,
    )

    # Input/Output arguments
    parser.add_argument(
        "--input",
        "-i",
        type=str,
        default=None,
        help="Input data file (.h5ad or .h5). If not provided, uses example dataset.",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default="output",
        help="Output directory (default: output)",
    )
    parser.add_argument(
        "--name",
        "-n",
        type=str,
        default="test_run",
        help="Name of the analysis run (optional)",
    )

    # Analysis parameters
    parser.add_argument(
        "--species",
        "-s",
        type=str,
        default="human",
        help="Species for GRN inference (default: human)",
    )
    parser.add_argument(
        "--cluster-key",
        type=str,
        default="leiden",
        help="Clustering column name (default: leiden)",
    )
    parser.add_argument(
        "--cluster-key-stratification",
        type=str,
        default=None,
        help="Clustering column name to perform stratification on (default: None)",
    )
    parser.add_argument(
        "--embedding-grn",
        type=str,
        default="X_draw_graph_fa",
        help="Embedding name for GRN inference (default: X_draw_graph_fa)",
    )
    parser.add_argument(
        "--embedding-hotspot",
        type=str,
        default="X_umap",
        help="Embedding name for Hotspot analysis (default: X_umap)",
    )
    parser.add_argument(
        "--raw-count-layer",
        type=str,
        default="raw_counts",
        help="Layer name for raw counts (default: raw_counts)",
    )
    parser.add_argument(
        "--tf-dictionary",
        type=str,
        default=None,
        help="Path to TF to target gene dictionary pickle file (default: None)",
    )

    # Quality control parameters
    parser.add_argument(
        "--min-genes",
        type=int,
        default=200,
        help="Minimum genes per cell for QC (default: 200)",
    )
    parser.add_argument(
        "--min-counts",
        type=int,
        default=500,
        help="Minimum counts per cell for QC (default: 500)",
    )

    # Computational parameters
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed (default: 42)",
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=8,
        help="Number of parallel jobs (default: 8)",
    )

    # Pipeline control flags
    parser.add_argument(
        "--skip-qc",
        action="store_true",
        default=False,
        help="Skip quality control (use if already performed)",
    )
    parser.add_argument(
        "--skip-celloracle",
        action="store_true",
        default=False,
        help="Skip CellOracle GRN inference",
    )
    parser.add_argument(
        "--skip-hotspot",
        action="store_true",
        default=False,
        help="Skip Hotspot module identification",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        default=False,
        help="Enable debug mode",
    )

    args = parser.parse_args()

    # Print header
    print("\n" + "=" * 70)
    print("TRNspot Complete Analysis Pipeline")
    print("=" * 70)
    print(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    start_time = datetime.now()

    # Setup configuration
    print(f"\n{'='*70}")
    print("STEP 0: Configuration Setup")
    print(f"{'='*70}")

    # Create directories
    setup_directories(
        output_dir=args.output,
        figures_dir=os.path.join(args.output, "figures"),
        debug=args.debug,
    )

    # Set random seed and scanpy settings
    set_random_seed(args.seed)
    set_scanpy_settings()

    sc.settings.logfile = os.path.join(args.output, "logs", "scanpy_log.txt")

    # Update TRNspot configuration
    config.update_config(
        OUTPUT_DIR=args.output,
        FIGURES_DIR=os.path.join(args.output, "figures"),
        FIGURES_DIR_QC=os.path.join(args.output, "figures", "qc"),
        FIGURES_DIR_GRN=os.path.join(args.output, "figures", "grn"),
        FIGURES_DIR_HOTSPOT=os.path.join(args.output, "figures", "hotspot"),
        GRN_N_JOBS=args.n_jobs,
        HOTSPOT_N_JOBS=args.n_jobs,
        QC_MIN_GENES=args.min_genes,
        QC_MIN_COUNTS=args.min_counts,
    )

    print(f"Random seed: {args.seed}")
    print(f"Output directory: {args.output}")
    print(f"Figures directory: {os.path.join(args.output, 'figures')}")
    print(f"Parallel jobs: {args.n_jobs}")
    print(f"QC thresholds: min_genes={args.min_genes}, min_counts={args.min_counts}")

    try:
        # Run pipeline
        log_dir = os.path.join(args.output, "logs")

        adata = load_data(args.input)
        adata_preprocessed = preprocessing_pipeline(
            adata, args.name, skip_qc=args.skip_qc, log_dir=log_dir
        )
        adata_list, adata_stratification_list = stratification_pipeline(
            adata_preprocessed, args.cluster_key_stratification
        )

        if args.cluster_key_stratification:
            # Process each stratified adata
            for idx, (adata_l, adata_stratification) in enumerate(
                zip(adata_list, adata_stratification_list)
            ):
                # Create stratified folder INSIDE the main output directory
                stratified_output_dir = os.path.join(
                    args.output, str(adata_stratification)
                )
                stratified_figures_dir = os.path.join(stratified_output_dir, "figures")
                stratified_log_dir = os.path.join(stratified_output_dir, "logs")

                print(f"\n{'='*70}")
                print(f"Processing stratified dataset: {adata_stratification}")
                print(f"{'='*70}")
                print(f"  Output directory: {stratified_output_dir}")

                # Create directories for this stratification
                setup_directories(stratified_output_dir, stratified_figures_dir)

                # Update config for this stratification
                config.update_config(
                    OUTPUT_DIR=stratified_output_dir,
                    FIGURES_DIR=stratified_figures_dir,
                    FIGURES_DIR_QC=os.path.join(stratified_figures_dir, "qc"),
                    FIGURES_DIR_GRN=os.path.join(stratified_figures_dir, "grn"),
                    FIGURES_DIR_HOTSPOT=os.path.join(stratified_figures_dir, "hotspot"),
                )

                # Process each stratified adata
                adata_l = dimensionality_reduction_clustering(
                    adata_l, cluster_key=args.cluster_key, log_dir=stratified_log_dir
                )
                celloracle_result = celloracle_pipeline(
                    adata_l,
                    cluster_key=args.cluster_key,
                    species=args.species,
                    raw_count_layer=args.raw_count_layer,
                    embedding_name=args.embedding_grn,
                    TG_to_TF_dictionary=args.tf_dictionary,
                    skip_celloracle=args.skip_celloracle,
                    log_dir=stratified_log_dir,
                )
                hotspot_result = hotspot_pipeline(
                    adata_l,
                    layer_key=args.raw_count_layer,
                    embedding_key=args.embedding_hotspot,
                    skip_hotspot=args.skip_hotspot,
                    log_dir=stratified_log_dir,
                )

                # Generate summary
                generate_summary(
                    adata_l,
                    celloracle_result,
                    hotspot_result,
                    start_time,
                    stratified_output_dir,
                )

                # Run GRN deep analysis on the tracked files
                grn_deep_analysis_pipeline(
                    os.path.join(
                        stratified_output_dir, "celloracle", "grn_merged_scores.csv"
                    )
                )
        else:
            # No stratification; process the original adata
            adata_l = dimensionality_reduction_clustering(
                adata_preprocessed, cluster_key=args.cluster_key, log_dir=log_dir
            )
            celloracle_result = celloracle_pipeline(
                adata_l,
                cluster_key=args.cluster_key,
                species=args.species,
                raw_count_layer=args.raw_count_layer,
                embedding_name=args.embedding_grn,
                TG_to_TF_dictionary=args.tf_dictionary,
                skip_celloracle=args.skip_celloracle,
                log_dir=log_dir,
            )
            hotspot_result = hotspot_pipeline(
                adata_l,
                layer_key=args.raw_count_layer,
                embedding_key=args.embedding_hotspot,
                skip_hotspot=args.skip_hotspot,
                log_dir=log_dir,
            )

            # Generate summary
            generate_summary(
                adata_l,
                celloracle_result,
                hotspot_result,
                start_time,
                args.output,
            )

            # Run GRN deep analysis on the tracked files
            grn_deep_analysis_pipeline(
                os.path.join(
                    stratified_output_dir, "celloracle", "grn_merged_scores.csv"
                )
            )

        # Print final summary for all stratifications
        print(f"\n{'='*70}")
        print("Pipeline completed successfully! ✓")
        print(f"{'='*70}")

        if adata_stratification_list:
            print(f"\nProcessed {len(adata_stratification_list)} stratified datasets:")
            for stratification in adata_stratification_list:
                print(f"  - {os.path.join(args.output, str(stratification))}/")
        else:
            print(f"\nResults saved to: {args.output}/")

        # Track output files to keep record of all the files generated by the pipeline
        track_files(args.output)

        # Reset config back to main output directory
        config.update_config(
            OUTPUT_DIR=args.output,
            FIGURES_DIR=os.path.join(args.output, "figures"),
            FIGURES_DIR_QC=os.path.join(args.output, "figures", "qc"),
            FIGURES_DIR_GRN=os.path.join(args.output, "figures", "grn"),
            FIGURES_DIR_HOTSPOT=os.path.join(args.output, "figures", "hotspot"),
        )

        from trnspot.grn_deep_analysis import merge_scores, plot_heatmap_scores

        total_merged_scores = merge_scores(
            os.path.join(args.output, "tracked_files.txt")
        )
        plot_heatmap_scores(total_merged_scores)
        print("  ✓ Generated overall GRN deep analysis heatmap")

        return 0

    except Exception as e:
        print(f"\n{'='*70}")
        print(f"Pipeline failed with error: {e}")
        print(f"{'='*70}\n")
        import traceback

        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
