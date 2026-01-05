# TRNspot AI Coding Instructions

**TRNspot** is a Python package for transcriptional regulatory network (TRN) analysis using single-cell data. It integrates Scanpy, CellOracle, and Hotspot into a modular, checkpoint-enabled pipeline.

## Architecture

### Data Flow

`AnnData (.h5ad)` → **Preprocessing** (QC/normalization) → **Clustering** → **CellOracle** (GRN inference) + **Hotspot** (gene modules) → **Deep Analysis** (visualization) → **Reporting** (HTML/PDF)

### Key Components

| File                               | Purpose                                           |
| ---------------------------------- | ------------------------------------------------- |
| `trnspot/config.py`                | **SINGLE SOURCE OF TRUTH** - all parameters       |
| `examples/complete_pipeline.py`    | `PipelineController` class - central orchestrator |
| `run_complete_analysis.py`         | Entry point wrapper                               |
| `trnspot/preprocessing.py`         | Scanpy wrappers (QC, normalize, cluster)          |
| `trnspot/celloracle_processing.py` | GRN inference with CellOracle                     |
| `trnspot/hotspot_processing.py`    | Spatial autocorrelation modules                   |
| `trnspot/grn_deep_analysis.py`     | Network visualization (NetworkX, Marsilea)        |
| `trnspot/reporting/`               | HTML/PDF report generation module                 |

## Critical Patterns

### 1. Never Hardcode Parameters

```python
# ✅ CORRECT - use config
from trnspot import config
sc.pp.filter_cells(adata, min_genes=config.QC_MIN_GENES)
plt.figure(figsize=config.PLOT_FIGSIZE_MEDIUM)
plt.savefig(f"{config.FIGURES_DIR_GRN}/plot.png", dpi=config.SAVE_DPI)

# ❌ WRONG - hardcoded values
sc.pp.filter_cells(adata, min_genes=200)
```

### 2. Extend PipelineController for New Steps

```python
# In examples/complete_pipeline.py - add new method:
def run_step_new_analysis(self, adata, log_dir=None):
    """Execute new analysis step."""
    log_step("Controller.NewAnalysis", "STARTED")
    try:
        result = my_analysis_function(adata)
        log_step("Controller.NewAnalysis", "COMPLETED")
        return result
    except Exception as e:
        log_error("Controller.NewAnalysis", e)
        raise
```

**NEVER** create standalone scripts - always integrate into the controller.

### 3. Logging Pattern

```python
from examples.complete_pipeline import log_step, log_error

log_step("MyStep", "STARTED", {"n_cells": adata.n_obs})
log_step("MyStep", "COMPLETED", {"result_count": len(results)})
log_error("MyStep", exception)  # Logs to error.log with traceback
```

### 4. Function Signatures with Config Defaults

```python
def my_function(
    adata: AnnData,
    param: Optional[int] = None,  # Allow override
) -> AnnData:
    if param is None:
        param = config.MY_PARAM  # Fall back to config
```

## Developer Workflows

### Run Pipeline

```bash
# Full run with example data
python run_complete_analysis.py

# With custom data
python run_complete_analysis.py --input data.h5ad --output results

# Specific steps only (checkpoints auto-resume)
python examples/complete_pipeline.py --steps load preprocessing clustering

# Stratified parallel analysis (by cell type)
python examples/complete_pipeline.py --cluster-key-stratification celltype --parallel --n-jobs 4
```

### Testing

```bash
pytest tests/                    # Run all tests
pytest tests/test_config.py     # Config tests specifically
```

When adding new config parameters, add corresponding tests to `tests/test_config.py`.

### Output Structure

```
output/
├── preprocessed_adata.h5ad
├── report.html                  # Interactive HTML report
├── report.pdf                   # PDF report (requires weasyprint)
├── logs/{pipeline.log, error.log, *.checkpoint}
├── figures/{qc/, grn/, hotspot/}
├── celloracle/{grn_merged_scores.csv, grn_filtered_links.pkl}
├── hotspot/{autocorrelation_results.csv, gene_modules.csv}
└── stratified_analysis/<ClusterName>/  # Per-cluster results
```

## Reporting Module

Generate comprehensive HTML/PDF reports with the `trnspot.reporting` module:

```python
from trnspot.reporting import generate_report

# Generate both HTML and PDF reports
outputs = generate_report(
    output_dir="results/",
    title="My Analysis Report",
    subtitle="Sample Dataset Analysis",
    adata=adata,
    celloracle_result=(oracle, links),
    hotspot_result=hs,
    log_file="results/logs/pipeline.log",
)
print(f"HTML: {outputs['html']}")
print(f"PDF: {outputs.get('pdf', 'Not generated')}")

# For smaller HTML file size (images as relative paths)
outputs = generate_report(
    output_dir="results/",
    title="My Analysis Report",
    embed_images=False,  # Use relative paths instead of base64
)
```

### Report Contents

- **Input Data Summary**: Dataset dimensions, annotations, embeddings
- **Configuration Settings**: All parameters used in the analysis
- **QC Section**: Quality control metrics and filtering results
- **Preprocessing**: Normalization and HVG selection details
- **Clustering**: Cluster sizes and visualization
- **CellOracle GRN**: Network analysis results and TF rankings
- **Hotspot Modules**: Gene module analysis
- **Operations Log**: Complete pipeline execution log
- **Plot Gallery**: Interactive gallery with lightbox navigation

### Report Options

| Parameter      | Default           | Description                                                                                                   |
| -------------- | ----------------- | ------------------------------------------------------------------------------------------------------------- |
| `embed_images` | `True`            | If `True`, embed images as base64 (self-contained but larger). If `False`, use relative paths (smaller file). |
| `formats`      | `["html", "pdf"]` | List of output formats to generate                                                                            |

## AnnData Conventions

- **Metrics**: `.obs` (per-cell), `.var` (per-gene)
- **Embeddings**: `.obsm['X_pca']`, `.obsm['X_umap']`
- **Raw counts**: Store in `.layers['raw_count']` before normalization
- **Import**: `import scanpy as sc`, access as `adata`

## Adding New Config Parameters

1. Add to `trnspot/config.py` with docstring
2. Add to `get_config()` return dict in same file
3. Add test in `tests/test_config.py`
