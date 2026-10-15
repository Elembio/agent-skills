---
name: raw-data-ingestion
description: Loading non-zarr data into AnnData — including .h5ad files, .zarr AnnData directories, AVITI24 RawCellStats .parquet files, and generic CSV or parquet count matrices with accompanying cell/gene metadata — for use with all downstream analysis skills. Use when the user points to a .h5ad, .parquet, .csv, or saved AnnData file, or asks how to reload previously saved data or ingest a generic tabular count matrix.
metadata:
  version: 0.2.3
  author: elembio
---

# Raw Data Ingestion

## When to Use This Skill

Use this skill (not `spatialdata-loading-and-access`) when the input is:

| Input format | Typical source |
|---|---|
| `*.h5ad` | AnnData saved by a previous session (`adata.write_h5ad`) |
| `*.zarr` (AnnData, not SpatialData) | AnnData saved with `adata.write_zarr` |
| `*RawCellStats.parquet` | AVITI24 raw cell stats upstream of the zarr builder |
| `*.csv` / `*.parquet` count matrix | Generic tabular counts (cells x genes or genes x cells), optionally with a cell metadata file and/or a gene metadata file |

If the input is a SpatialData `.zarr` or `.zarr.zip` store, use [spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md) instead.

---

## Quick Start

### Reload a saved AnnData (h5ad or zarr)

```python
import anndata as ad

# From h5ad (most common — single file)
adata = ad.read_h5ad("/path/to/transcript.h5ad")

# From AnnData zarr directory
adata = ad.read_zarr("/path/to/transcript.zarr")

print(adata)
# Verify expected structure
print(adata.obs.columns.tolist())    # Well, tile_id, X_tile, Y_tile, …
print(adata.var.columns.tolist())    # target_name, batch_name, …
print(list(adata.obsm.keys()))       # cell_geometry, qc_metrics, …
print(list(adata.layers.keys()))     # raw (if previously normalized)
```

### Load from RawCellStats parquet

```python
import pandas as pd
import anndata as ad
import numpy as np

df = pd.read_parquet("/path/to/RawCellStats.parquet")

# Parquet column names (RawCellStats). Rename to Loader slots below.
PARQUET_OBS = ["Cell", "Well", "WellLabel", "Tile", "X", "Y", "Area"]
MORPH_PREFIXES = ("AreaShape_", "Granularity_", "Texture_", "Intensity_", "Location_", "RadialDistribution_")

def is_barcoding_rna(c):  return any(c.endswith(f".B0{i}") for i in range(1, 7))
def is_barcoding_protein(c): return any(c.endswith(f".P0{i}") for i in range(1, 7)) and "_Nuclear" not in c
def is_imaging(c):  # RawCellStats still uses .CPxx suffixes for imaging channels + CP morphology
    return (any(c.endswith(f".CP0{i}") for i in range(1, 9)) and "_Nuclear" not in c) or c.startswith(MORPH_PREFIXES)

# Prefer SpatialData zarr + Loader when available — this path is for pre-zarr parquet only
transcript_cols = [c for c in df.columns if is_barcoding_rna(c)]      # legacy barcoding
protein_cols    = [c for c in df.columns if is_barcoding_protein(c)] # legacy barcoding
imaging_cols    = [c for c in df.columns if is_imaging(c)]           # imaging (CP-group ± IF in zarr era)

obs = df[[c for c in PARQUET_OBS if c in df.columns]].copy()
obs = obs.rename(columns={"Tile": "tile_id", "X": "X_tile", "Y": "Y_tile"})
obs.index = obs["Cell"].astype(str)

def make_adata(feature_cols, dtype=np.float32):
    a = ad.AnnData(
        X=df[feature_cols].values.astype(dtype),
        obs=obs.drop(columns=["Area"], errors="ignore").copy(),
        var=pd.DataFrame(index=feature_cols),
    )
    if "Area" in obs.columns:
        a.obsm["cell_geometry"] = obs[["Area"]].copy()
    return a

transcript_adata = make_adata(transcript_cols)           # legacy
protein_adata    = make_adata(protein_cols)              # legacy
imaging_adata    = make_adata(imaging_cols, dtype=np.float64)  # → imaging-pipeline MAD-z (not pp.normalize)
```

### Load from a generic CSV or parquet count matrix

For tabular count data from any source (public atlases, other platforms, Kallisto/Salmon outputs, etc.).  The minimal input is a **counts file**; a **cell metadata file** (obs) and **gene metadata file** (var) are optional.

```python
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix

# 1. Counts: rows = cells, columns = genes.  First column is the cell id.
#    Use read_parquet for .parquet instead.
counts = pd.read_csv("/path/to/counts.csv", index_col=0)

# If the file is genes x cells instead, transpose:
# counts = counts.T

# 2. Optional: cell metadata (obs) indexed by the same cell ids as counts rows
obs = pd.read_csv("/path/to/cell_metadata.csv", index_col=0)
obs = obs.loc[counts.index]    # keep order aligned with counts

# 3. Optional: gene metadata (var) indexed by gene ids matching counts columns
var = pd.read_csv("/path/to/gene_metadata.csv", index_col=0)
var = var.loc[counts.columns]

# 4. Build AnnData.  Store integer counts in layers["raw"] so pseudobulk DE works.
X = csr_matrix(counts.values)
adata = ad.AnnData(X=X, obs=obs, var=var)
adata.layers["raw"] = adata.X.copy()

# 5. If the CSV carries spatial coordinates, populate obsm["spatial"] for spatial-statistics
_xy = None
if {"X", "Y"}.issubset(adata.obs.columns):
    _xy = adata.obs[["X", "Y"]].to_numpy(dtype=float)
elif {"x", "y"}.issubset(adata.obs.columns):
    _xy = adata.obs[["x", "y"]].to_numpy(dtype=float)
elif {"X_tile", "Y_tile"}.issubset(adata.obs.columns):
    _xy = adata.obs[["X_tile", "Y_tile"]].to_numpy(dtype=float)
if _xy is not None:
    adata.obsm["spatial"] = _xy

print(adata)
```

Downstream skills then work unchanged: QC / normalization / HVG / PCA / UMAP / clustering / pseudobulk DE / spatial statistics all operate on this AnnData the same way they do on AVITI24-derived AnnData.

---

## After Loading

All formats produce AnnData objects compatible with downstream skills. Verify structure before proceeding:

```python
# Check whether normalization has already been applied
if "raw" in adata.layers:
    print("Data is normalized — adata.X is normalized, adata.layers['raw'] has counts")
else:
    print("Data is raw counts — safe to normalize")

# Check whether QC metrics exist
if "qc_metrics" in adata.obsm:
    print("QC already computed")
```

---

## Detail Files

- [IMPLEMENTATION.md](./IMPLEMENTATION.md) — full parquet conversion with `obsm["cell_geometry"]`, `var` metadata, protein nuclear compartment, validation checks, and save patterns.
- [REFERENCE.md](./REFERENCE.md) — parquet column naming convention, batch ID extraction, expected shapes, `obs`/`var` cross-reference.

## Related Skills

- [spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md) — preferred when a SpatialData zarr is available
- [imaging-pipeline](../imaging-pipeline/SKILL.md) — imaging / IF after ingest (MAD-z; not `pp.normalize`)
- [cell-quality-control](../cell-quality-control/SKILL.md)
- [multimodal-data-preparation](../multimodal-data-preparation/SKILL.md)

## Disclaimer

This content is provided "as-is" without warranty of any kind. AI-generated output may contain
inaccuracies. Review and verify all outputs before relying on them.
