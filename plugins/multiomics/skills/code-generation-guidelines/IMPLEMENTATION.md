# Code Generation — Implementation

Concrete patterns and copy-paste snippets for common analysis tasks. See [SKILL.md](./SKILL.md) for import conventions and pipeline order, and [REFERENCE.md](./REFERENCE.md) for the full library listing.

---

## Table of Contents

1. [Performance Optimization](#performance-optimization)
2. [Code Structure Patterns](#code-structure-patterns)
3. [Common Patterns](#common-patterns)
4. [Best Practices](#best-practices)

## Performance Optimization

### UMAP on Large Datasets (>5,000 samples)

Use optimized parameters to prevent timeouts:

```python
import umap

reducer = umap.UMAP(
    n_neighbors=30,      # Faster than default 15
    min_dist=0.3,        # Faster than default 0.1
    low_memory=True,     # Memory-efficient mode
    verbose=True,        # Shows progress
    random_state=42
)
embedding = reducer.fit_transform(data)

# Or with Scanpy
sc.tl.umap(adata, n_neighbors=30, min_dist=0.3)
```

### Very Large Datasets (>20,000 samples)

Consider random subsampling for exploratory analysis:

```python
# Subsample for initial exploration
sample_idx = np.random.choice(len(df), size=10000, replace=False)
df_sample = df.iloc[sample_idx]

# For AnnData
sc.pp.subsample(adata, n_obs=10000)
```

### Sparse Matrix Operations

Many AVITI24 count matrices are sparse. Handle appropriately:

```python
import scipy.sparse as sp

# Check if sparse
if sp.issparse(adata.X):
    # Use sparse-compatible operations
    total_counts = np.asarray(adata.X.sum(axis=1)).ravel()

    # For dense operations, convert only when necessary
    X_dense = adata.X.toarray()

    # Better: use sparse methods
    from scipy.sparse import csr_matrix
    X_normalized = csr_matrix(X_dense / total_counts[:, np.newaxis])
```

### Parallel Processing

Use multiprocessing for independent operations:

```python
# For squidpy spatial analysis
sq.gr.spatial_autocorr(adata, mode='moran', genes=genes, n_jobs=4)

# For PyDESeq2
dds = DeseqDataSet(adata=pdata, design_factors="condition")
dds.deseq2(n_cpus=4)

# For custom operations
from joblib import Parallel, delayed

results = Parallel(n_jobs=4)(
    delayed(process_well)(well_data)
    for well_data in wells
)
```

### Memory Management

For large datasets, manage memory explicitly:

```python
import gc

# Process one well at a time
for well in wells:
    well_adata = adata[adata.obs['Well'] == well].copy()
    # Process...
    del well_adata
    gc.collect()

# Use generators for large iterations
def process_cells_generator(adata, batch_size=1000):
    for i in range(0, adata.n_obs, batch_size):
        yield adata[i:i+batch_size]

# Clear intermediate results
del intermediate_result
gc.collect()
```

---

## Code Structure Patterns

### End-to-End Cell Sequence

Run these as ordered notebook cells (not a standalone `.py` script), writing outputs to a local directory you control:

```python
# === Imports ===
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import elembio_spatialdata_analysis as ebsd
from elembio_spatialdata_analysis.load import Loader

# === Configuration ===
DATA_PATH = "/path/to/run.zarr"
OUTPUT_DIR = "/path/to/output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# === Discover (always first) ===
loader = Loader(DATA_PATH)
tables = set(loader.available_tables)
print("tables:", sorted(tables))
print("types:", sorted(loader.available_target_types))
# Run QC next (run-quality-overview) — well/segmentation + present modalities only
# Then cell QC below for what exists

# === Example: 3′ path when present ===
if "ThreePrimeUntargeted" not in tables:
    raise ValueError("This template expects ThreePrimeUntargeted; branch on available_tables instead")

adata = loader.load_tables(table="ThreePrimeUntargeted", lazy=False)
ebsd.pp.compute_qc_metrics(adata, "ThreePrimeUntargeted")
ebsd.pp.filter_cells(adata, target_type="ThreePrimeUntargeted")  # MT default 0.10
ebsd.pp.filter_genes(adata, "ThreePrimeUntargeted", min_cells_per_gene=10)
ebsd.pp.normalize(adata, "ThreePrimeUntargeted")  # never normalize imaging/OPS

sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=False)
sc.pp.pca(adata, n_comps=50)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5, flavor="igraph", n_iterations=2, directed=False)

sc.pl.umap(adata, color=["leiden", "Well"], save="_clusters.png")
adata.write(os.path.join(OUTPUT_DIR, "processed.h5ad"))
print(f"Results saved to {OUTPUT_DIR}")
```

Imaging / OPS branches (only if present): load → `compute_qc_metrics` → `filter_cells` with `Imaging` / `OPS`; **skip normalize**. Resolve SpecializedTargeted DISS via `run.table_suffixes_for_target_type(("OPS", "SpecializedTargeted"))`, then keep OPS assay only (drop `QC_Panel` / gene panels).

### Reproducible Output Pattern

```python
# Always set random state for reproducibility
np.random.seed(42)
sc.settings.seed = 42

# Include version information
import sys
print(f"Python: {sys.version}")
print(f"Scanpy: {sc.__version__}")
print(f"AnnData: {ad.__version__}")

# Log parameters used
params = {
    "leiden_resolution": 0.5,
    "n_neighbors": 15,
    "n_pcs": 30,
    "min_cells_per_gene": 10,  # RNA gene filter; not assignment_rate
}
print(f"Parameters: {params}")
```

---

## Common Patterns

### Saving Figures

```python
# Always specify output directory, dpi, and bbox_inches
plt.savefig(os.path.join(OUTPUT_DIR, 'figure_name.png'),
            dpi=150, bbox_inches='tight')
plt.close()  # Close to free memory
```

### Summary Statistics

```python
# Export summary as JSON for programmatic access
import json

summary = {
    'n_cells': int(adata.n_obs),
    'n_genes': int(adata.n_vars),
    'n_clusters': int(adata.obs['leiden'].nunique()),
    'wells': adata.obs['Well'].unique().tolist()
}

with open(os.path.join(OUTPUT_DIR, 'summary.json'), 'w') as f:
    json.dump(summary, f, indent=2)
```

### Error Handling

```python
from elembio_spatialdata_analysis.load import Loader

loader = Loader(DATA_PATH)  # same store as in Load Data section

# Branch on what exists — never invent modalities
present = set(loader.available_tables)
if not present:
    raise ValueError("No tables in store")

if "ThreePrimeUntargeted" in present:
    rna = loader.load_tables(table="ThreePrimeUntargeted", lazy=False)
else:
    rna = None

img_key = "imaging" if "imaging" in present else ("CellPaint" if "CellPaint" in present else None)
if img_key:
    imaging = loader.load_tables(table=img_key, lazy=False)
else:
    imaging = None

run = next(iter(loader._runs.values()))
specialized = run.table_suffixes_for_target_type(("OPS", "SpecializedTargeted"))
ops_suffixes = [t for t in specialized if "qc_panel" not in t.lower()]  # OPS assay; skip QC/gene panels
```

---

## Best Practices

1. **Use absolute paths** for output files to ensure they're found
2. **Close plots** after saving to free memory
3. **Copy AnnData objects** before modification: `adata.copy()`
4. **Check sparse matrices** before using numpy operations
5. **Set random seeds** for reproducibility
6. **Document parameters** used in analysis
7. **Save intermediate results** for long-running analyses
8. **Use appropriate dtypes** (e.g., category for cluster labels)
