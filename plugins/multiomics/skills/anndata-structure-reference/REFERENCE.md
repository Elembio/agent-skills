# AnnData Structure Reference — Full Slot Inventory

Detailed reference for **current** AVITI24 **`global-table-*`** AnnData objects. See [SKILL.md](./SKILL.md) for the short overview and quick-start patterns.

---

## Table of Contents

1. [AnnData Tree Diagram](#1-anndata-tree-diagram)
2. [.X — Expression / Feature Matrix](#2-x--expression--feature-matrix)
3. [.obs — Cell Metadata](#3-obs--cell-metadata)
   - [Standard columns](#standard-columns)
   - [Spatial coordinate columns](#spatial-coordinate-columns)
   - [Morphology / area (current build)](#morphology--area-current-build)
   - [Per-batch counts (`obsm["total_counts"]`)](#per-batch-counts-obsmtotal_counts)
4. [.var — Feature Metadata](#4-var--feature-metadata)
5. [.obsm — Multi-dimensional Cell Annotations](#5-obsm--multi-dimensional-cell-annotations)
   - [Core unified keys](#core-unified-keys)
   - [QC metrics](#qc-metrics)
   - [Filter flags](#filter-flags)
   - [OPS-specific](#ops-specific)
   - [Embeddings (added during analysis)](#embeddings-added-during-analysis)
6. [.layers — Alternative Matrices](#6-layers--alternative-matrices)
7. [.uns — Unstructured Annotations](#7-uns--unstructured-annotations)
8. [.obsp — Pairwise Cell Annotations](#8-obsp--pairwise-cell-annotations)
9. [Modality-Specific Structure](#9-modality-specific-structure)
10. [Common Access Patterns](#10-common-access-patterns)
11. [Legacy Key Migration Notes](#11-legacy-key-migration-notes)
12. [Validation Function](#12-validation-function)

---

## 1. AnnData Tree Diagram

```
AnnData object
    ├── .X                    # Counts (3′/Transcript/Protein/OPS) or dense imaging features
    ├── .obs                  # Cell metadata (Well, tile_id, …; OPS Confidence when merged)
    ├── .var                  # Feature metadata (target_type, batch_name, …; imaging: feature_class)
    ├── .obsm                 # Multi-dimensional cell annotations — discover keys at runtime
    │   ├── cell_geometry         # DataFrame — Area, NuclearArea, Xum_tile, Yum_tile, … (common)
    │   ├── total_counts          # DataFrame — per-batch assigned (Transcript/Protein; may exist on 3′)
    │   ├── unassigned_counts     # DataFrame — per-batch unassigned (barcoding; not for imaging QC)
    │   ├── control_variables     # DataFrame — NSB* (Transcript / some Protein — not 3′/imaging/OPS)
    │   ├── filter_flags          # DataFrame — written by filter_cells
    │   ├── qc_metrics            # DataFrame — written by compute_qc_metrics
    │   ├── target_assignment     # OPS only (also mergeable into obs)
    │   ├── spatial               # DataFrame — centroids X, Y in the region pixel frame (build)
    │   ├── spatial_artifacts     # boolean artifact membership (build; optional)
    │   └── X_pca / X_umap        # added during analysis
    ├── .layers               # nuclear, raw (post-normalize snapshot), …
    └── .uns                  # spatialdata_attrs, metadata_attrs, filter_params, spatial_artifact_params, …
```

Keys under `obsm` are **modality-dependent** — always `print(list(adata.obsm.keys()))`.

---

## 2. .X — Expression / Feature Matrix

Shape `(n_cells, n_features)`. **Not always integer counts.**

### By modality

| Modality | Table suffix | Content | Data Type |
|----------|--------------|---------|-----------|
| **ThreePrimeUntargeted** | `ThreePrimeUntargeted` | Gene counts / depth (DISS; batches **summed** into one matrix) | Sparse CSR int |
| **Imaging** | `imaging` (alias `CellPaint`) | Morphology (± IF) features | Dense float |
| **OPS** | custom `{batch}-{panel}` | OPS library-member counts (DISS) | Sparse CSR int |
| **Transcript** (legacy) | `Transcript` | Targeted transcript counts | Sparse CSR int |
| **Protein** (legacy) | `Protein` | Targeted protein counts | Sparse CSR int |

After `pp.normalize` on **count** modalities (3′, Transcript, Protein), pre-norm → `layers["raw"]` and `X` becomes normalized. **Do not normalize imaging or OPS.**

### Access patterns

```python
import numpy as np, scipy.sparse as sp

is_sparse = sp.issparse(adata.X)
# For count modalities: row sums = depth. For imaging: sum of feature intensities (not counts/depth).
row_sum = np.asarray(adata.X.sum(axis=1)).ravel() if is_sparse else np.asarray(adata.X).sum(axis=1)

if "GAPDH" in adata.var_names:  # often Transcript; may be absent on 3′/imaging
    gene_idx = adata.var_names.get_loc("GAPDH")
    gene_expr = adata.X[:, gene_idx].toarray().ravel() if is_sparse else adata.X[:, gene_idx]
```

---

## 3. .obs — Cell Metadata

Per-cell metadata with shape `(n_cells, n_columns)`. Index contains cell identifiers.

### Standard columns

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `Cell` | int/str | Cell identifier. On well-keyed stores, matches pixel values in `{well}-labels-Cell` | — |
| `cell_id_local` | int | ID unique within a tile / reused across wells — **not** the well-mask pixel value | `123` |
| `cell_id_global` | int | Same ids as `Cell` on current well-keyed stores | — |
| **`tile_id`** | str | **Tile identifier — preferred in current builds** | `"L1R01C01S1"` |
| `Tile` | str | Legacy / alternate tile name (some tables) | same role |
| `Well` | str | Plate position | `"A1"` |
| `WellLabel` | str | Sample-sheet name on that well | `"Sample1_DMSO_1"` |
| `region_key` | str | SpatialData region reference | — |
| `well_id` | str | Well ID string (optional) | — |

**`Well` vs `WellLabel` (do not mix):**

| Column | Means | Where |
|---|---|---|
| `Well` | Plate position (`A1`) | `adata.obs` |
| `WellLabel` | Sample-sheet string | `adata.obs` |
| `well_label` | Same strings as `WellLabel` | `run.wells` / `run.tiles` only |

On cells, color and groupby use `WellLabel` or `condition`. There is no `obs["well_label"]`. Parse design from the cell table when it is loaded ([experiment-design](../experiment-design/SKILL.md)).

Safe tile-column lookup:

```python
tile_col = "tile_id" if "tile_id" in adata.obs.columns else "Tile"
tiles = adata.obs[tile_col].unique()
```

### Spatial coordinate columns

Discover present names (`print(adata.obs.columns)`, `print(adata.obsm["spatial"].columns)` if present). Typical layout:

| Location | Columns | Units | Use |
|----------|---------|-------|-----|
| `obsm["spatial"]` | `X`, `Y` | region pixels | squidpy / centroids (convert to ndarray) |
| `obs` | `X_tile`, `Y_tile` | tile pixels | tile-local position |
| `obsm["cell_geometry"]` | `Xum_tile`, `Yum_tile` (and `Area`, …) | tile-local microns | micron neighborhoods within a well |

Do not invent `obs["X"]` / `obs["Y"]` or `obs["Xum"]` / `obs["Yum"]`.

### Morphology / area (current build)

**Cell area and related morphology live in `obsm["cell_geometry"]`** as a DataFrame with columns such as `Area`, `NuclearArea`, `AreaUm`, `Xum_tile`, `Yum_tile`, `AreaShape_*` (print the columns). This DataFrame is **not** copied into `obs` automatically.

- `pp.filter_cells` resolves area from `obsm["cell_geometry"]` first, then falls back to `obs["Area"]`. Default loads without copying are valid for QC.
- To use area in `obs` explicitly: `adata.obs["Area"] = adata.obsm["cell_geometry"]["Area"]`
- `load_tables(..., include_cell_geometry=True)` merges geometry columns into `.X` and `.var` (not into `obs`). Rarely needed; emits a warning for non-imaging tables.

```python
def cell_area_series(adata):
    if "cell_geometry" in adata.obsm and "Area" in adata.obsm["cell_geometry"].columns:
        return adata.obsm["cell_geometry"]["Area"]
    if "Area" in adata.obs.columns:
        return adata.obs["Area"]
    raise KeyError("No Area in obsm['cell_geometry'] or obs['Area']")
```

### Per-batch counts (`obsm["total_counts"]`)

When present (typical **Transcript / Protein**), per-batch assigned totals are a DataFrame in **`obsm["total_counts"]`** (columns = batch names) — not `obs['total_counts_B*']`.

For **ThreePrimeUntargeted**, expression in `.X` is already **summed across sequencing batches**. Cell QC count gating uses summed `.X` row sums (`passes_counts`), not per-batch filters. Do not invent a per-batch filter story for 3′ just because `obsm["total_counts"]` exists.

```python
if "total_counts" in adata.obsm:
    b01 = adata.obsm["total_counts"]["B01"]   # only if that column exists
```

### Access patterns

```python
# All columns
print(adata.obs.columns.tolist())

# Filter by well
well_a1 = adata[adata.obs["Well"] == "A1"]
```

---

## 4. .var — Feature Metadata

Per-feature metadata with shape `(n_features, n_columns)`. Index is feature names (genes / proteins).

### Standard columns

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `target_name` | str | Feature name (same as index) | `"GAPDH"` |
| `target_id` | int | Numeric ID within batch | `42` |
| `target_type` | str | Feature class | `"ThreePrimeUntargeted"`, `"Imaging"` / `"CellPaint"`, `"ProteinIF"`, `"OPS"`, legacy `"Transcript"` / `"Protein"` |
| `batch_name` | str | Molecular / imaging batch | `"B01"`, `"CP01"`, `"P01"`, … |
| `control_type` | str | `""` = analysis target; `"Positive Control 1"` = spike-in in `.X` (legacy Transcript GAPDH, etc.) | — |

### Imaging-oriented columns (when present)

| Column | Type | Description |
|--------|------|-------------|
| `feature_class` | str | e.g. `"base"`, `"cellprofiler"` |
| `suffix` | str | Compartment suffix (`""`, `"_Nuclear"`, …) |
| `is_background` | bool | Background channel flag |

### Additional columns (often barcoding)

| Column | Type | Description |
|--------|------|-------------|
| `target_sequence` | str | Barcode sequence |
| `probe_concentration` | float | Probe concentration used |
| `gene_id` | str | Ensembl gene ID (if available) |

**NSB / unassigned (legacy Transcript, some Protein):** NSB → `obsm["control_variables"]` (`NSB*`); unassigned → `obsm["unassigned_counts"]`. **Not** part of 3′ / imaging / OPS QC narratives.

### Access patterns

```python
print(adata.var.columns.tolist())

# Genes in a specific batch
batch_b01_genes = adata.var[adata.var["batch_name"] == "B01"].index

# Subset table to one batch
adata_b01 = adata[:, adata.var["batch_name"] == "B01"]

# Analysis targets only (exclude positive controls)
analysis_genes = adata.var[adata.var["control_type"].fillna("") == ""].index

# Positive controls (e.g. GAPDH)
pos_ctrl = adata.var[adata.var["control_type"] == "Positive Control 1"]
```

---

## 5. .obsm — Multi-dimensional Cell Annotations

Dictionary of arrays / DataFrames with shape `(n_cells, *)`. **Current build** uses unified keys.

### Core unified keys

| Key | Type | Content | Source |
|-----|------|---------|--------|
| **`cell_geometry`** | DataFrame | `Area`, `NuclearArea`, `AreaUm`, `Xum_tile`, `Yum_tile`, `AreaShape_*` — primary area source | build (common) |
| **`total_counts`** | DataFrame | Assigned counts per **batch** per cell | barcoding / some RNA tables — guard |
| **`unassigned_counts`** | DataFrame | Unassigned / background per batch | barcoding — guard; not 3′ QC |
| **`control_variables`** | DataFrame | NSB* columns | **Transcript** (legacy) — omit for 3′/imaging/OPS narratives |
| `data_missing` | DataFrame | Optional missingness flags | optional |
| **`spatial`** | DataFrame | Cell centroids `X`, `Y` in the region pixel frame | build (common) |
| **`spatial_artifacts`** | DataFrame | Boolean artifact membership. Rollups: `any-imaging`, `any-ThreePrimeUntargeted`, `any-OPS`, plus `any`. Per-unit columns exist (`imaging-Actin`, …). Absent on tiles / skip-artifact / older stores. When no spatial artifacts are present, all columns are `False`. | build (optional) |

```python
# Transcript / barcoding only — guard keys
if "total_counts" in adata.obsm and "unassigned_counts" in adata.obsm:
    import numpy as np
    tc, ua = adata.obsm["total_counts"], adata.obsm["unassigned_counts"]
    assigned, unassigned = tc.sum(axis=1).values, ua.sum(axis=1).values
    denom = assigned + unassigned
    assignment_rate = np.divide(assigned, denom, where=denom > 0, out=np.zeros(adata.n_obs))

if "control_variables" in adata.obsm:
    nsb_cols = [c for c in adata.obsm["control_variables"].columns if "NSB" in c.upper()]
```

### QC metrics

| Key | Type | Content | Source |
|-----|------|---------|--------|
| `qc_metrics` | DataFrame | Computed QC metrics — absent until `compute_qc_metrics` runs | `pp.compute_qc_metrics` |

QC metric columns are **modality-specific and data-dependent**. Imaging often has a sparse
metric set (intensity filtering uses `filter_flags`, not MT/NSB). Discover with
`adata.obsm["qc_metrics"].columns` and guard every access.

Common columns (when present):

- `total_counts` — depth / intensity total metric when written
- `n_genes_by_counts` — genes with counts > 0 (3′, Transcript)
- `pct_counts_mt` — mitochondrial % in **`qc_metrics`** after `compute_qc_metrics` (3′) — not a guaranteed `obs` column

See [cell-quality-control](../cell-quality-control/SKILL.md) for the full per-modality breakdown.

### Filter flags

| Key | Type | Content | Source |
|-----|------|---------|--------|
| `filter_flags` | DataFrame | Boolean per-cell filter flags | `pp.filter_cells` |

Filter flag columns (discover at runtime — set is modality-specific):

- `passes_nan_filter`
- `passes_empty_filter`
- `passes_counts` — count-depth gate (3′: summed `.X` per well; Transcript/Protein: per batch). Legacy reports may still show `passes_counts_per_batch`
- `passes_area`
- `passes_imaging_intensity` — Imaging
- `passes_mt_outlier` — 3′
- `passes_nsb_outlier` / `passes_gapdh_outlier` / `passes_assignment_rate` — Transcript only when applied
- `passes_spatial_artifacts` — written always; all-True unless `filter_spatial_artifacts=True` (gates on the modality rollup of `obsm["spatial_artifacts"]`)
- `passes_all_filters`

`filter_cells` also writes **`uns["filter_params"]`** — the resolved threshold dict used in that run. Use it to audit or reproduce filtering:

```python
if "filter_params" in adata.uns:
    print(adata.uns["filter_params"])
    # e.g. {'min_quant_count': 0.025, 'max_quant_nsb': 0.99, 'max_quant_gapdh': None, ...}
```

### OPS-specific

Present only for OPS tables (table name is custom, e.g. `"B01-<PanelName>"`). Loaded by `load_tables(..., include_target_assignment=True)`.

| Key | Type | Content |
|-----|------|---------|
| `target_assignment` | DataFrame (cells × 5) | Per-cell OPS identity call — columns: `MostAbundantTarget` (str, cell-line name or `"Undetermined"`), `MostAbundantTargetCount` (int, winning library-member polony count), `TotalCount` (int, total polony reads), `Confidence` (float [0–1], `MostAbundantTargetCount / (TotalCount + 1)`), `target_gene` (str, gene-level annotation) |

When `include_target_assignment=True` is passed to `load_tables()`, these five columns are also merged into `obs`:

```python
run = next(iter(loader._runs.values()))
ops_tables = run.table_suffixes_for_target_type(("OPS", "SpecializedTargeted"))
ops = loader.load_tables(table=ops_tables[0], include_target_assignment=True, lazy=False)
# Never table="OPS" — suffixes are custom, e.g. "B01-<PanelName>"

# Available directly in obs after include_target_assignment=True
ops.obs["MostAbundantTarget"]      # cell-line call (str), or "Undetermined"
ops.obs["MostAbundantTargetCount"] # winning library-member count (int)
ops.obs["TotalCount"]              # total reads across all barcodes (int)
ops.obs["Confidence"]              # MostAbundantTargetCount / (TotalCount + 1) (float 0–1)
ops.obs["target_gene"]             # gene-level annotation (str)

# Standard filter
assigned = ops.obs["Confidence"] >= 0.5
print(f"Assigned: {assigned.sum():,} / {ops.n_obs:,}")
```

**`uns["metadata_attrs"]` for OPS tables:**

```python
ops.uns["metadata_attrs"]
# {
#     "batch_class":      "DISSBatch",
#     "target_type":      "OPS",           # identifies this as an OPS table
#     "batch_name":       ["B01"],
#     "target_site_name": "<custom>",      # custom panel name from run manifest — varies per run
# }
```

### Embeddings (added during analysis)

| Key | Type | Shape | Content |
|-----|------|-------|---------|
| `X_pca` | ndarray | (n_cells, n_pcs) | PCA embedding |
| `X_umap` | ndarray | (n_cells, 2) | UMAP embedding |
| `X_tsne` | ndarray | (n_cells, 2) | t-SNE embedding |

```python
print(list(adata.obsm.keys()))
```

---

## 6. .layers — Alternative Matrices

Dictionary of matrices with the same shape as `.X`.

| Key | Content | When Present |
|-----|---------|--------------|
| `nuclear` | Nuclear compartment counts | Nuclear segmentation available |
| `raw` | Pre-normalization counts snapshot | After `pp.normalize()` |
| `log_normalized` | Log-transformed normalized values | After log1p transform |

```python
print(list(adata.layers.keys()))

# Nuclear fraction
if "nuclear" in adata.layers:
    import scipy.sparse as sp
    total  = adata.X.toarray() if sp.issparse(adata.X) else adata.X
    nuclear_frac = adata.layers["nuclear"] / (total + 1e-10)

# Save snapshot before normalization (pp.normalize uses "raw" by default)
adata.layers["raw"] = adata.X.copy()

# Access raw snapshot after normalization
raw = adata.layers["raw"]
```

---

## 7. .uns — Unstructured Annotations

Dictionary for arbitrary metadata.

### SpatialData linkage

| Key | Type | Content |
|-----|------|---------|
| `spatialdata_attrs` | dict | Linkage to SpatialData elements |

Attributes inside `spatialdata_attrs`:

- `region` — list of label element names this table annotates
- `region_key` — column in `.obs` specifying region (usually `"region_key"`)
- `instance_key` — column linking to label instances. Well-keyed masks match `obs["Cell"]`.

```python
sd_attrs = adata.uns["spatialdata_attrs"]
# region list looks like ['A1-labels-Cell', ...]; paint well masks with obs['Cell']
```

### Run-level metadata

| Key | Type | Content |
|-----|------|---------|
| `metadata_attrs` | dict | Metadata attribute names |
| `spatial_artifact_params` | dict | Present with `obsm["spatial_artifacts"]`: `version`, `spatial_method`, `detection_enabled`, `detection_disabled_reason`, `rollup_columns`, per-modality params |
| `control_variables_metadata` | DataFrame | Metadata about control targets — columns: `target_name`, `target_id`, `control_type`, `batch_name` |

```python
ctrl_meta = adata.uns["control_variables_metadata"]
print(ctrl_meta)
```

### Analysis results (added during analysis)

| Key | Type | Content |
|-----|------|---------|
| `filter_params` | dict | Resolved filter thresholds written by `filter_cells` / `compute_filter_diagnostics` |
| `pca` | dict | PCA results (`variance`, `variance_ratio`, `components`) |
| `neighbors` | dict | Neighbor graph parameters |
| `umap` | dict | UMAP parameters |
| `leiden` | list | Leiden clustering info |
| `rank_genes_groups` | dict | Marker gene results |
| `moranI` | DataFrame | Moran's I spatial autocorrelation |
| `{cluster}_nhood_enrichment` | dict | Neighborhood enrichment |

```python
print(list(adata.uns.keys()))

if "rank_genes_groups" in adata.uns:
    import scanpy as sc
    markers = sc.get.rank_genes_groups_df(adata, group="0")

if "pca" in adata.uns:
    variance_ratio = adata.uns["pca"]["variance_ratio"]
```

---

## 8. .obsp — Pairwise Cell Annotations

Sparse matrices of shape `(n_cells, n_cells)` for cell-cell relationships.

| Key | Type | Content |
|-----|------|---------|
| `connectivities` | sparse | Neighbor graph edge weights |
| `distances` | sparse | Neighbor graph distances |
| `spatial_connectivities` | sparse | Spatial neighbor graph |
| `spatial_distances` | sparse | Spatial neighbor distances |

```python
print(list(adata.obsp.keys()))

if "connectivities" in adata.obsp:
    graph = adata.obsp["connectivities"]
    print(f"Neighbor graph: {graph.nnz} edges")
```

---

## 9. Modality-Specific Structure

### ThreePrimeUntargeted (`global-table-ThreePrimeUntargeted`)

```python
adata.X                                   # Sparse int gene counts — batches already summed
adata.n_vars                              # Large (thousands+)
adata.obsm["cell_geometry"]["Area"]
# After compute_qc_metrics:
adata.obsm["qc_metrics"]["pct_counts_mt"]  # prefer this over obs
adata.obsm["qc_metrics"]["n_genes_by_counts"]
# After filter_cells: passes_counts (summed depth), passes_mt_outlier, passes_area
# No NSB / GAPDH gates — do not narrate them
# filter_genes before normalize; pass min_cells_per_gene=10 explicitly (package default may be None)
```

### Imaging (`global-table-imaging`, legacy list name `CellPaint`)

```python
adata.X                                   # Dense float features (not counts)
adata.var["target_type"]                  # "CellPaint" stains and/or "ProteinIF"
adata.var["batch_name"]                   # CP* morphology vs P* IF
# Often: feature_class, suffix, is_background
adata.obsm["cell_geometry"]
# filter_cells(..., "Imaging") → passes_imaging_intensity (CP* only); no filter_genes; no normalize
# Load with table="imaging" even if available_tables lists "CellPaint"
```

### OPS (custom `global-table-{batch}-{panel}`)

```python
# Discover: run.table_suffixes_for_target_type(("OPS", "SpecializedTargeted"))
ops = loader.load_tables(table=name, include_target_assignment=True, lazy=False)
ops.obs["Confidence"]                     # primary call quality (not assignment_rate)
ops.obs["MostAbundantTarget"]
# filter_cells(..., "OPS") — area-focused; count quantiles off
# No normalize / no filter_genes for expression-style workflows
```

### Legacy Transcript (`global-table-Transcript`)

```python
adata.X                                   # Sparse int counts (positive-control rows included)
adata.obsm["total_counts"]               # per-batch assigned — used by passes_counts
adata.obsm["unassigned_counts"]
adata.obsm["control_variables"]          # NSB*
adata.var["control_type"]               # "" vs "Positive Control 1" (GAPDH in X)
```

### Legacy Protein (`global-table-Protein`)

```python
# Similar barcoding obsm to Transcript; var target_type "Protein"
# NSB when control_variables present; no MT
```

---

## 10. Common Access Patterns

### Describe AnnData contents

```python
import scipy.sparse as sp

def describe_adata(adata):
    """Print summary of AnnData contents."""
    print(f"Shape: {adata.n_obs} cells × {adata.n_vars} features")
    print(f"\n.obs columns:  {list(adata.obs.columns)}")
    print(f"\n.var columns:  {list(adata.var.columns)}")
    print(f"\n.obsm keys:    {list(adata.obsm.keys())}")
    print(f"\n.layers keys:  {list(adata.layers.keys())}")
    print(f"\n.uns keys:     {list(adata.uns.keys())}")
    print(f"\n.obsp keys:    {list(adata.obsp.keys())}")
    if sp.issparse(adata.X):
        print(f"\n.X is sparse with {adata.X.nnz} non-zero entries")
    else:
        print(f"\n.X is dense")

describe_adata(adata)
```

### Per-batch assignment rate (Transcript / barcoding only)

Skip unless both keys exist. Not a standard 3′ filter; not OPS confidence.

```python
import numpy as np

if "total_counts" in adata.obsm and "unassigned_counts" in adata.obsm:
    tc = adata.obsm["total_counts"]
    ua = adata.obsm["unassigned_counts"]
    assigned   = tc.sum(axis=1).values
    unassigned = ua.sum(axis=1).values
    denom = assigned + unassigned
    assignment_rate = np.divide(assigned, denom, where=denom > 0, out=np.zeros(adata.n_obs))
```

### Spatial coordinates for squidpy

```python
import numpy as np

sp = adata.obsm["spatial"]
adata.obsm["spatial"] = (
    sp[["X", "Y"]].to_numpy(dtype=float) if hasattr(sp, "columns") else np.asarray(sp, dtype=float)
)
# Then subset obs["Well"] before sq.gr.spatial_neighbors
```

### Save / restore raw counts snapshot

```python
# Before normalization
adata.layers["raw"] = adata.X.copy()

# Restore (e.g. to re-normalize with different method)
adata.X = adata.layers["raw"].copy()
```

---

## 11. Legacy Key Migration Notes

These patterns appear in older notebooks and pipeline versions. Do **not** use them for new work.

| Legacy pattern | Current equivalent |
|----------------|--------------------|
| `adata.obs["total_counts_B01"]` | `adata.obsm["total_counts"]["B01"]` |
| `adata.obsm["control_variables_B01"]` (per-batch key) | `adata.obsm["control_variables"]` (unified DataFrame, columns = probe names) |
| `adata.obsm["control_variables"]["Unassigned"]` | `adata.obsm["unassigned_counts"]` |
| `adata.obs["Area"]` as primary area source | `adata.obsm["cell_geometry"]["Area"]` |
| `adata.obs["Tile"]` | `adata.obs["tile_id"]` (with fallback to `"Tile"`) |

---

## 12. Validation Function

```python
import scipy.sparse as sp

def validate_aviti24_adata(adata, modality: str = "ThreePrimeUntargeted") -> bool:
    """Lightweight checks. modality: ThreePrimeUntargeted | Imaging | OPS | Transcript | Protein."""
    errors = []

    if "Well" not in adata.obs.columns:
        errors.append("Missing .obs column: Well")
    if "tile_id" not in adata.obs.columns and "Tile" not in adata.obs.columns:
        errors.append("Missing tile column: need tile_id or Tile")
    if "cell_geometry" not in adata.obsm and "Area" not in adata.obs.columns:
        errors.append("No cell area: need obsm['cell_geometry'] or obs['Area']")

    if modality in ("Transcript", "Protein"):
        for col in ("batch_name", "target_type"):
            if col not in adata.var.columns:
                errors.append(f"Missing .var column: {col}")
        if "total_counts" not in adata.obsm:
            errors.append("Missing obsm['total_counts'] (expected for barcoding tables)")

    if modality == "Imaging" and sp.issparse(adata.X):
        errors.append("Imaging X unexpectedly sparse — expected dense features")

    if errors:
        print("Validation FAILED:")
        for e in errors:
            print(f"  - {e}")
        return False
    print("Validation passed")
    return True

validate_aviti24_adata(adata, "ThreePrimeUntargeted")
```
