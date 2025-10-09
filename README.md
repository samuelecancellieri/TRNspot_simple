# TRNspot

A Python package for transcriptional regulatory network analysis.

## Installation

### From source

```bash
# Clone the repository
git clone https://github.com/yourusername/trnspot.git
cd trnspot

# Create and activate virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install the package in development mode
pip install -e .

# Or install with development dependencies
pip install -e ".[dev]"
```

## Usage

### Configuration and Reproducibility

Set random seed and configure default parameters:

```python
from trnspot import set_random_seed, config

# Set random seed for reproducibility
set_random_seed(42)

# View all configuration parameters
config.print_config()

# Update specific parameters
config.update_config(
    QC_MIN_GENES=300,
    QC_MIN_COUNTS=1000,
    PLOT_DPI=600
)
```

### Quality Control

Perform comprehensive quality control on single-cell RNA-seq data:

```python
import scanpy as sc
from trnspot import config
from trnspot.preprocessing import perform_qc, plot_qc_violin, plot_qc_scatter

# Load your data
adata = sc.read_h5ad('your_data.h5ad')

# Perform QC - uses config defaults automatically
adata_qc = perform_qc(adata)
# Equivalent to: min_genes=200, min_counts=500, pct_counts_mt_max=20.0, min_cells=10

# Or override specific parameters
adata_qc = perform_qc(
    adata,
    min_genes=300,      # Override
    min_counts=1000     # Override
    # Other params use config defaults
)
```

### Complete Workflow Example

```python
import scanpy as sc
from trnspot import set_random_seed, config
from trnspot.preprocessing import perform_qc, perform_grn_pre_processing

# 1. Set up reproducibility
set_random_seed(42)

# 2. Optionally customize config
config.update_config(QC_MIN_GENES=300, QC_MIN_COUNTS=1000)

# 3. Load and process data
adata = sc.read_h5ad('your_data.h5ad')
adata = perform_qc(adata)  # Uses config defaults

# 4. Normalize
sc.pp.normalize_total(adata, target_sum=config.NORMALIZE_TARGET_SUM)
sc.pp.log1p(adata)

# 5. GRN preprocessing
adata = perform_grn_pre_processing(adata)  # Uses config defaults
```

### CellOracle Integration

Perform gene regulatory network inference using CellOracle:

```python
from trnspot.celloracle_processing import (
    create_oracle_object,
    run_PCA,
    run_KNN,
    run_links
)

# Note: Requires CellOracle installation
# pip install celloracle

# Create Oracle object
oracle = create_oracle_object(
    adata=adata,
    cluster_column_name='leiden',
    embedding_name='X_umap',
    raw_count_layer='raw_counts'
)

# Perform PCA and KNN imputation
oracle = run_PCA(oracle)
run_KNN(oracle, n_comps=50)

# Infer regulatory links
links = run_links(
    oracle,
    cluster_column_name='leiden',
    p_cutoff=0.001
)

# Save results
oracle.to_hdf5('oracle_object.celloracle.oracle')
links.to_hdf5('grn_links.celloracle.links')
```

### Running Examples

```bash
# Activate environment
source venv/bin/activate

# Run QC example
python examples/example_qc.py

# Run configuration example
python examples/config_example.py

# Run CellOracle workflow (requires celloracle)
python examples/celloracle_workflow.py

# Run quick demo
python examples/quick_demo.py

# Test config integration
python examples/test_config_integration.py
```

## Features

- **Configuration Management**: Centralized configuration for reproducibility
  - Global random seed setting
  - Default parameters for all analysis steps
  - Easy parameter updates
  - Configuration profiles for different analysis types
- **Quality Control**: Comprehensive QC with multiple visualization options
  - Cell filtering based on gene count, total counts, and mitochondrial percentage
  - Automated QC metrics calculation
  - Before/after filtering comparison plots
  - Violin and scatter plots for detailed inspection
- **Data Preprocessing**: Complete preprocessing pipeline
  - Normalization and scaling for single-cell RNA-seq data
  - Highly variable genes selection
  - PCA and dimensionality reduction
  - Neighborhood graph construction
- **CellOracle Integration**: Gene regulatory network inference
  - Oracle object creation with raw or normalized counts
  - Automated PCA component selection
  - KNN imputation for noise reduction
  - Regulatory link inference with statistical filtering
  - Network visualization and quality metrics
- **Gene Regulatory Network Analysis**: Network construction and analysis tools
  - TF-target gene relationship inference
  - Network topology analysis
  - Cluster-specific GRN construction

## Documentation

- **[Configuration Guide](docs/CONFIG.md)** - Complete configuration documentation
- **[QC Functions](docs/QC_FUNCTIONS.md)** - Quality control functions guide
- **[CellOracle Processing](docs/CELLORACLE_PROCESSING.md)** - CellOracle integration guide
- **[Package Structure](docs/PACKAGE_STRUCTURE.md)** - Package organization overview
- **[Preprocessing Updates](docs/PREPROCESSING_CONFIG_UPDATE.md)** - Config integration details

## Development

### Running Tests

```bash
# Run all tests
pytest tests/

# Run specific test file
pytest tests/test_config.py -v
pytest tests/test_preprocessing.py -v
pytest tests/test_celloracle.py -v

# Run with coverage
pytest tests/ --cov=trnspot --cov-report=html
```

### Code Quality

```bash
# Format code
black trnspot/

# Check linting
flake8 trnspot/

# Type checking
mypy trnspot/
```

### Testing Notes

- CellOracle tests use mocking when CellOracle is not installed
- Some tests are skipped if CellOracle is not available
- Use `pytest -v` for verbose output
- Tests cover all major functionality with both unit and integration tests

## License

MIT License

## Authors

Samuele Cancellieri
