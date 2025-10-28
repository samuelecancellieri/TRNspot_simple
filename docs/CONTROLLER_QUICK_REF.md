# Pipeline Controller - Quick Reference

## New Command-Line Flags

### `--parallel`

Enable parallel execution for stratified analyses.

```bash
python examples/complete_pipeline.py \
    --cluster-key-stratification celltype \
    --parallel \
    --n-jobs 4
```

### `--steps [step1 step2 ...]`

Run only specific pipeline steps.

**Available steps:**

- `load` - Load input data
- `preprocessing` - QC and normalization
- `stratification` - Split by clustering key
- `clustering` - Dimensionality reduction and clustering
- `celloracle` - GRN inference with CellOracle
- `hotspot` - Gene module identification
- `grn_analysis` - Deep GRN analysis
- `summary` - Generate summary report

```bash
# Run only preprocessing and clustering
python examples/complete_pipeline.py \
    --steps load preprocessing clustering

# Run only analysis steps (assumes preprocessing done)
python examples/complete_pipeline.py \
    --steps celloracle hotspot
```

## PipelineController API

### Core Methods

```python
from examples.complete_pipeline import PipelineController

controller = PipelineController(args, start_time)

# Individual steps
controller.run_step_load()
controller.run_step_preprocessing()
controller.run_step_stratification()
controller.run_step_clustering(adata)
controller.run_step_celloracle(adata)
controller.run_step_hotspot(adata)
controller.run_step_grn_analysis(csv_path)

# Execution modes
controller.run_complete_pipeline(steps=None, parallel=False)
controller.run_stratified_pipeline_sequential()
controller.run_stratified_pipeline_parallel(n_jobs=4)
controller.process_single_stratification(adata, name)
```

## Common Use Cases

### 1. Run Full Pipeline

```bash
python examples/complete_pipeline.py --input data.h5ad --output results
```

### 2. Stratified Analysis (Parallel)

```bash
python examples/complete_pipeline.py \
    --input data.h5ad \
    --output results \
    --cluster-key-stratification celltype \
    --parallel \
    --n-jobs 4
```

### 3. Run Specific Steps

```bash
# Only preprocessing
python examples/complete_pipeline.py \
    --input data.h5ad \
    --output results \
    --steps load preprocessing

# Only analysis (skip preprocessing)
python examples/complete_pipeline.py \
    --input preprocessed.h5ad \
    --output results \
    --steps clustering celloracle hotspot
```

### 4. Resume from Checkpoint

```bash
# First run - completes preprocessing
python examples/complete_pipeline.py \
    --input data.h5ad \
    --output results \
    --steps load preprocessing clustering

# Second run - skips completed steps via checkpoints
python examples/complete_pipeline.py \
    --input data.h5ad \
    --output results \
    --steps load preprocessing clustering celloracle
```

### 5. Programmatic Control

```python
from examples.complete_pipeline import PipelineController
from argparse import Namespace
from datetime import datetime

args = Namespace(
    input="data.h5ad",
    output="results",
    cluster_key_stratification="celltype",
    parallel=True,
    n_jobs=4,
    # ... other args
)

controller = PipelineController(args, datetime.now())

# Run preprocessing only
controller.run_complete_pipeline(
    steps=['load', 'preprocessing', 'stratification']
)

# Inspect results
print(f"Found {len(controller.adata_list)} stratifications")

# Run remaining steps in parallel
controller.run_stratified_pipeline_parallel(n_jobs=8)
```

## Performance Guide

### Memory Optimization

```bash
# Limit parallel workers for large datasets
python examples/complete_pipeline.py \
    --cluster-key-stratification celltype \
    --parallel \
    --n-jobs 2  # Reduce memory usage

# Or run sequentially
python examples/complete_pipeline.py \
    --cluster-key-stratification celltype
    # No --parallel flag
```

### CPU Optimization

```bash
# Maximum parallelization
python examples/complete_pipeline.py \
    --cluster-key-stratification celltype \
    --parallel \
    --n-jobs 16  # Use all cores
```

### Checkpoint Optimization

```bash
# Save time by using checkpoints
python examples/complete_pipeline.py \
    --input data.h5ad \
    --output results
    # Automatically skips completed steps
```

## Workflow Examples

### Development Workflow

```bash
# 1. Test preprocessing
python examples/complete_pipeline.py \
    --steps load preprocessing

# 2. Test clustering
python examples/complete_pipeline.py \
    --steps load preprocessing clustering

# 3. Run full pipeline
python examples/complete_pipeline.py
```

### Production Workflow

```bash
# 1. Preprocess all data
python examples/complete_pipeline.py \
    --input data.h5ad \
    --output production \
    --steps load preprocessing stratification

# 2. Run stratified analysis in parallel
python examples/complete_pipeline.py \
    --input data.h5ad \
    --output production \
    --cluster-key-stratification celltype \
    --parallel \
    --n-jobs 8
```

## See Also

- [Full Controller Guide](PIPELINE_CONTROLLER_GUIDE.md)
- [Complete Pipeline Guide](COMPLETE_PIPELINE_GUIDE.md)
- [Usage Examples](../examples/controller_usage_examples.py)
