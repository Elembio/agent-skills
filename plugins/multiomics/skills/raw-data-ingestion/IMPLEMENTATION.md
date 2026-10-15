# Raw Data Ingestion — Implementation

## Table of Contents

1. [Overview](#overview)
2. [Load and Inspect](#1-load-and-inspect)
3. [Column Selectors](#2-column-selectors)
4. [Build obs](#3-build-obs)
5. [Extract batch ID from column names](#4-extract-batch-id-from-column-names)
6. [Build var](#5-build-var)
7. [Build per-modality AnnData](#6-build-per-modality-anndata)
8. [Protein nuclear compartment (optional)](#7-protein-nuclear-compartment-optional)
9. [Validation checks](#8-validation-checks)
10. [Save to h5ad (optional)](#9-save-to-h5ad-optional)
11. [Next steps](#10-next-steps)
12. [Pitfalls](#pitfalls)

## Overview

This document covers the full parquet → AnnData conversion path. For loading previously saved `.h5ad` or `.zarr` AnnData files, see the quick start in `SKILL.md`.

A `*RawCellStats.parquet` file is the upstream wide table from which the SpatialData zarr builder creates per-modality global tables. Converting it directly lets you run all downstream skills (QC, normalization, etc.) before a zarr store is built.

---

## 1. Load and Inspect

```python
import pandas as pd
import anndata as ad
import numpy as np
import re

path = "/path/to/RawCellStats.parquet"
df = pd.read_parquet(path)
print(df.shape)          # (n_cells, n_features)
print(df.columns[:10])   # sanity-check column names
```

---

## 2. Column Selectors

```python
# Parquet names in RawCellStats
PARQUET_OBS = ["Cell", "Well", "WellLabel", "Tile", "X", "Y", "Area"]
MORPH_PREFIXES = (
    "AreaShape_", "Granularity_", "Texture_",
    "Intensity_", "Location_", "RadialDistribution_",
)

def is_transcript(c):
    return any(c.endswith(f".B0{i}") for i in range(1, 7))

def is_protein_whole(c):
    return any(c.endswith(f".P0{i}") for i in range(1, 7)) and "_Nuclear" not in c

def is_protein_nuclear(c):
    return "_Nuclear" in c and any(f".P0{i}" in c.split("_Nuclear")[1] for i in range(1, 7))

def is_imaging(c):
    # raw channel (no morphology prefix, no _Nuclear)
    is_channel = any(c.endswith(f".CP0{i}") for i in range(1, 9)) and "_Nuclear" not in c
    # per-channel morphology (e.g. Intensity_MeanIntensity_WGA.CP01)
    is_morph   = c.startswith(MORPH_PREFIXES) and "_Nuclear" not in c
    return is_channel or is_morph

transcript_cols = [c for c in df.columns if is_transcript(c)]       # legacy barcoding
protein_cols = [c for c in df.columns if is_protein_whole(c)]       # legacy barcoding
protein_nuc_cols = [c for c in df.columns if is_protein_nuclear(c)]
imaging_cols = [c for c in df.columns if is_imaging(c)]
```

---

## 3. Build obs

```python
obs = df[[c for c in PARQUET_OBS if c in df.columns]].copy()
obs = obs.rename(columns={"Tile": "tile_id", "X": "X_tile", "Y": "Y_tile"})
obs.index = obs["Cell"].astype(str)
obs.index.name = None
# Keep WellLabel — matches Loader modality tables (sample-sheet name). Well is the plate position.
```

---

## 4. Extract batch ID from column names

The batch suffix (`.B01`, `.P02`, etc.) encodes the AVITI cycle/batch:

```python
def batch_from_col(col: str) -> str:
    m = re.search(r"\.(B\d+|P\d+|CP\d+)$", col)
    return m.group(1) if m else ""

def target_from_col(col: str) -> str:
    # strip trailing batch suffix
    return re.sub(r"\.(B\d+|P\d+|CP\d+)$", "", col)
```

---

## 5. Build var

```python
def make_var(cols, target_type: str) -> pd.DataFrame:
    var = pd.DataFrame(index=cols)
    var["target_name"] = [target_from_col(c) for c in cols]
    var["batch_name"]  = [batch_from_col(c)  for c in cols]
    var["target_type"] = target_type
    # Mark unassigned controls
    var["control_type"] = var["target_name"].apply(
        lambda n: "Unassigned" if n.lower() == "unassigned" else ""
    )
    return var
```

---

## 6. Build per-modality AnnData

```python
def make_adata(feature_cols, target_type, dtype=np.float32, obs_df=None):
    obs_df = obs_df if obs_df is not None else obs
    X = df[feature_cols].values.astype(dtype)
    adata = ad.AnnData(
        X=X,
        obs=obs_df.drop(columns=["Area"], errors="ignore").copy(),
        var=make_var(feature_cols, target_type),
    )
    # Match Loader: Area in cell_geometry; tile pixels already on obs as X_tile/Y_tile
    if "Area" in obs_df.columns:
        adata.obsm["cell_geometry"] = obs_df[["Area"]].copy()
    return adata

transcript_adata = make_adata(transcript_cols, "Transcript")  # legacy barcoding
protein_adata    = make_adata(protein_cols,    "Protein")     # legacy barcoding
imaging_adata    = make_adata(imaging_cols,    "Imaging")     # imaging table stand-in
```

---

## 7. Protein nuclear compartment (optional)

If you need whole-cell vs. nuclear protein side by side:

```python
protein_nuc_adata = make_adata(protein_nuc_cols, "Protein_Nuclear")

# Or store nuclear matrix in obsm of the whole-cell AnnData:
nuc_X = df[protein_nuc_cols].values.astype(np.float32)
protein_adata.obsm["protein_nuclear"] = pd.DataFrame(
    nuc_X, index=protein_adata.obs_names, columns=protein_nuc_cols
)
```

---

## 8. Validation checks

```python
assert transcript_adata.n_obs == len(df), "Row count mismatch"
assert not protein_adata.obs_names.duplicated().any(), "Duplicate cell IDs"
assert "cell_geometry" in transcript_adata.obsm
print("Imaging var sample:")
print(imaging_adata.var.head())
```

---

## 9. Save to h5ad (optional)

```python
out_dir = "/path/to/output/"
import os; os.makedirs(out_dir, exist_ok=True)

transcript_adata.write_h5ad(f"{out_dir}/transcript.h5ad")
protein_adata.write_h5ad(f"{out_dir}/protein.h5ad")
imaging_adata.write_h5ad(f"{out_dir}/imaging.h5ad")
```

---

## 10. Next steps

Prefer a SpatialData zarr + `Loader` when available. After parquet conversion:

1. **QC** — `ebsd.pp.compute_qc_metrics` / `filter_cells` with the right `target_type` (`ThreePrimeUntargeted` / `Imaging` / legacy `Transcript` / `Protein`) → [cell-quality-control](../cell-quality-control/SKILL.md)
2. **Counts** — `ebsd.pp.normalize` → [normalization-strategies](../normalization-strategies/SKILL.md)
3. **Imaging** — MAD-z feature matrix → [imaging-pipeline](../imaging-pipeline/SKILL.md) (**not** `pp.normalize`)
4. **Multimodal** — [multimodal-data-preparation](../multimodal-data-preparation/SKILL.md)

---

## Pitfalls

- **Imaging columns in RawCellStats still use `.CPxx` suffixes** — treat as the imaging modality, not a separate “CellPaint product.” Prefer zarr `table="imaging"` when the store exists.
- **Do not `pp.normalize` imaging floats** — use imaging-pipeline MAD-z.
- **`_Nuclear` columns are a separate compartment.** Do not mix them into whole-cell protein matrices.
- **Cell ID type.** Use `obs["Cell"].astype(str)` as the index.
- **Counts are integers** for barcoding Transcript/Protein; imaging channels / morphology are floats.
- **Batch suffix is the AVITI cycle**, not a biological batch — experimental batch is in `Well` / `WellLabel`.
