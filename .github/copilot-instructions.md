# TRNspot AI Coding Instructions

**TRNspot** is a Python package for transcriptional regulatory network (TRN) analysis using single-cell data. It integrates Scanpy, CellOracle, and Hotspot into a modular, checkpoint-enabled pipeline.

## Architecture

### Data Flow

`AnnData (.h5ad)` → **Preprocessing** (QC/normalization) → **Clustering** → **CellOracle** (GRN inference) + **Hotspot** (gene modules) → **Deep Analysis** (visualization) → **Reporting** (HTML/PDF)

### Package Structure (Modular)

```
trnspot/
├── __init__.py              # Main exports
├── config.py                # SINGLE SOURCE OF TRUTH - all parameters
├── utils/                   # Utility functions
│   ├── data.py              # AnnData loading, saving, subsetting
│   ├── io.py                # File I/O, checkpoints, directories
│   └── logging.py           # Logging setup, log_step, log_error
├── execution/               # Pipeline execution
│   ├── preprocessing.py     # QC, normalization, clustering
│   ├── celloracle.py        # GRN inference pipeline
│   └── hotspot.py           # Gene module analysis
├── plotting/                # Visualization
│   ├── qc.py                # QC violin/scatter plots
│   ├── embeddings.py        # UMAP, PCA plots
│   ├── grn.py               # Network graphs, heatmaps
│   └── hotspot.py           # Module heatmaps, UMAP overlays
├── reporting/               # Report generation
│   ├── generator.py         # ReportGenerator class
│   └── sections.py          # Section builders
└── (legacy modules)         # Backward compatibility
    ├── preprocessing.py
    ├── celloracle_processing.py
    ├── hotspot_processing.py
    └── grn_deep_analysis.py
```

### Key Components

| Module                          | Purpose                                           |
| ------------------------------- | ------------------------------------------------- |
| `trnspot.config`                | **SINGLE SOURCE OF TRUTH** - all parameters       |
| `trnspot.utils`                 | Data handling, I/O, logging utilities             |
| `trnspot.execution`             | Pipeline execution functions                      |
| `trnspot.plotting`              | Visualization functions                           |
| `trnspot.reporting`             | HTML/PDF report generation                        |
| `examples/complete_pipeline.py` | `PipelineController` class - central orchestrator |

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

### 2. Use Modular Imports

```python
# ✅ CORRECT - use new modular structure
from trnspot.utils import load_adata, save_adata, ensure_categorical_obs
from trnspot.execution import run_full_preprocessing, run_celloracle_pipeline
from trnspot.plotting import plot_umap, plot_network_graph
from trnspot.reporting import generate_report

# ✅ Also valid - convenience imports from top level
from trnspot import (
    load_adata,
    run_full_preprocessing,
    generate_report,
)
```

### 3. Extend PipelineController for New Steps

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

### 4. Logging Pattern

```python
from trnspot.utils import log_step, log_error

log_step("MyStep", "STARTED", {"n_cells": adata.n_obs})
log_step("MyStep", "COMPLETED", {"result_count": len(results)})
log_error("MyStep", exception)  # Logs to error.log with traceback
```

### 5. Categorical Data Handling

```python
# Always ensure categorical columns before stratification/grouping
from trnspot.utils import ensure_categorical_obs

adata = ensure_categorical_obs(adata, columns=["leiden", "celltype"])
```

### 6. Function Signatures with Config Defaults

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

### Generate Reports

```python
from trnspot.reporting import generate_report

# Generate both HTML and PDF reports
outputs = generate_report(
    output_dir="results/",
    title="My Analysis Report",
    adata=adata,
    celloracle_result=(oracle, links),
    hotspot_result=hs,
)
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
├── report.html                    # HTML report
├── report.pdf                     # PDF report (if weasyprint installed)
├── logs/{pipeline.log, error.log, *.checkpoint}
├── celloracle/{grn_merged_scores.csv, grn_filtered_links.pkl}
├── hotspot/{autocorrelation_results.csv, gene_modules.csv}
├── figures/{qc/, grn/, hotspot/}  # All generated plots
└── stratified_analysis/<ClusterName>/  # Per-cluster results
```

## AnnData Conventions

- **Metrics**: `.obs` (per-cell), `.var` (per-gene)
- **Embeddings**: `.obsm['X_pca']`, `.obsm['X_umap']`
- **Raw counts**: Store in `.layers['raw_count']` before normalization
- **Categorical columns**: Use `pd.CategoricalDtype` for cluster/stratification keys
- **Import**: `import scanpy as sc`, access as `adata`

## Adding New Config Parameters

1. Add to `trnspot/config.py` with docstring
2. Add to `get_config()` return dict in same file
3. Add test in `tests/test_config.py`
