# SpatialData Loading and Access — Accessing AnnData Slots

Open once a table is loaded and you need `.X`, `.obs`, `.var`, `.obsm` controls or `.layers`, or need to subset by well or feature.

## Contents
- [Access the expression matrix (.X)](#access-the-expression-matrix-x)
- [Access cell metadata (.obs)](#access-cell-metadata-obs)
- [Access feature metadata (.var)](#access-feature-metadata-var)
- [Access control and background counts (.obsm)](#access-control-and-background-counts-obsm)
- [Access alternative matrices (.layers)](#access-alternative-matrices-layers)
- [Subset by well or feature](#subset-by-well-or-feature)

---

## Access the expression matrix (.X)

```python
import numpy as np
import scipy.sparse as sp

X = adata.X
if sp.issparse(X):
    print(f"Sparse, sparsity={1 - X.nnz / (X.shape[0] * X.shape[1]):.2%}")
    X_dense = X.toarray()   # caution: memory
else:
    X_dense = X

# Counts for a specific gene
if "GAPDH" in adata.var_names:
    gene_idx = adata.var_names.get_loc("GAPDH")
    print(X_dense[:, gene_idx].max())

# Total counts per cell (use obsm['total_counts'] for per-batch totals)
total = X.sum(axis=1).A1 if sp.issparse(X) else X.sum(axis=1)
print(f"Median counts/cell: {np.median(total):.0f}")
```

---

## Access cell metadata (.obs)

```python
obs = adata.obs

# Tile column varies across store versions
tile_col = "tile_id" if "tile_id" in obs.columns else ("Tile" if "Tile" in obs.columns else None)

# Area: prefer obsm['cell_geometry'] (what filter_cells reads)
if "cell_geometry" in adata.obsm and "Area" in adata.obsm["cell_geometry"].columns:
    area = adata.obsm["cell_geometry"]["Area"]
    print(f"Area (obsm): median={area.median():.1f}")
elif "Area" in obs.columns:
    print(f"Area (obs fallback): median={obs['Area'].median():.1f}")

# To use area in obs explicitly:
# adata.obs["Area"] = adata.obsm["cell_geometry"]["Area"]

# Per-batch assigned totals live in obsm['total_counts'], not obs columns
if "total_counts" in adata.obsm:
    print("Batch columns:", adata.obsm["total_counts"].columns.tolist())

well_a1 = adata[adata.obs["Well"] == "A1"]
print(f"Cells in A1: {well_a1.n_obs}")
```

---

## Access feature metadata (.var)

```python
var = adata.var

# Analysis targets vs controls
analysis_targets = var[var["control_type"].fillna("") == ""]
pos_ctrl = var[var["control_type"] == "Positive Control 1"]
print(f"Analysis targets: {len(analysis_targets)}, Positive controls: {len(pos_ctrl)}")
# Positive controls (e.g. GAPDH) are rows in .X/.var — not a separate obsm column

# Features per batch
for batch in var["batch_name"].dropna().unique():
    print(f"  {batch}: {(var['batch_name'] == batch).sum()} features")
```

---

## Access control and background counts (.obsm)

Which keys a counts table carries depends on its type; imaging tables have none of them. Check each key before reading it, and report only those present — do not narrate NSB/GAPDH/assignment-rate as general QC.

- `total_counts` — every counts table (3′, OPS, barcoding).
- `control_variables` — negative-control targets kept out of `.X`: NSB probes (barcoding) or negative-control guides (OPS), when the panel defines them. Not on 3′.
- `unassigned_counts` — barcoding only; the basis for assignment rate. New format uses unified DataFrames (not `control_variables_*` per-batch keys).

```python
import numpy as np

# Only when present (barcoding or OPS)
if "control_variables" in adata.obsm:
    ctrl = adata.obsm["control_variables"]
    print("Control columns:", ctrl.columns.tolist())

if "unassigned_counts" in adata.obsm:
    print("Unassigned columns:", adata.obsm["unassigned_counts"].columns.tolist())

if "total_counts" in adata.obsm and "unassigned_counts" in adata.obsm:
    assigned = adata.obsm["total_counts"].sum(axis=1).values
    unassigned = adata.obsm["unassigned_counts"].sum(axis=1).values
    rate = assigned / (assigned + unassigned + 1e-9)
    print(f"Assignment rate (approx): {rate.mean():.1%}")

if "control_type" in adata.var.columns:
    mask = (adata.var["control_type"] == "Positive Control 1").values
    if mask.any():
        g = adata.X[:, mask]
        gapdh = np.asarray(g.sum(axis=1)).ravel()
        print(f"GAPDH median/cell: {np.median(gapdh):.1f}")
```

---

## Access alternative matrices (.layers)

```python
import scipy.sparse as sp

if "nuclear" in adata.layers:
    nuclear = adata.layers["nuclear"]
    print(f"Nuclear layer: {nuclear.shape}")

# Layers added during preprocessing:
# 'raw'         — snapshot before pp.normalize
# 'normalized'  — post-normalize counts
# 'raw_counts'  — legacy name in some older stores
```

---

## Subset by well or feature

```python
# Specific wells (cells axis), table already in memory
wells_mask = adata.obs["Well"].isin(["A1", "A2", "B1", "B2"])
subset = adata[wells_mask].copy()

# Lazy table: slice first, then .to_memory() on the slice (.copy() keeps it lazy)
subset = adata[wells_mask].to_memory()

# Imaging: pick targets or a kind of feature (columns absent on some older stores)
if "feature_class" in adata.var:
    base = (adata.var["feature_class"] == "base").to_numpy()    # per-target intensity stats, not CellProfiler texture/shape
    whole_cell = base & (adata.var["suffix"] == "").to_numpy()  # whole-cell intensity only
    proteins = (adata.var["target_type"] == "ProteinIF").to_numpy()  # IF markers only
    picked = adata[:, whole_cell & proteins].to_memory()
named = adata[:, adata.var_names.isin(["<Target1>", "<Target2>"])]  # whole-cell feature of named targets
```
