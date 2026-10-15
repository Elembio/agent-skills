# SpatialData Loading and Access — Reference

## Contents

1. [AnnData Structure Overview](#1-anndata-structure-overview)
2. [obsm Key Roles](#2-obsm-key-roles)
3. [Multi-Run Concatenation](#3-multi-run-concatenation)
4. [export_tables provenance (optional)](#4-export_tables-provenance-optional)
5. [Old patterns: legacy barcoding stores](#5-old-patterns-legacy-barcoding-stores)

---

## 1. AnnData Structure Overview

| Attribute | Content | Shape |
|---|---|---|
| `.X` | Count matrix for 3′ / OPS / legacy barcoding; continuous floats for imaging | `(n_cells, n_features)` |
| `.obs` | Cell metadata: `Well`, tile (`tile_id`), tile-local pixels (`X_tile`, `Y_tile`), cell IDs | `(n_cells, n_columns)` |
| `.var` | Feature metadata: `target_name` / imaging schema (`target_type`, `feature_class`, `suffix`, `is_background`, `batch_name`), … | `(n_features, n_columns)` |
| `.obsm` | Per-cell auxiliary DataFrames: QC inputs, geometry, flags | `str → DataFrame` |
| `.layers` | Alternative matrices (e.g. `"nuclear"` on some count tables; imaging has none) | `str → matrix` |
| `.uns` | Unstructured metadata (`metadata_attrs`, control tables, `export_provenance` after export, …) | `dict` |

**Table names for `load_tables(table=...)`:** use suffixes from `loader.available_tables` (or the literals `"ThreePrimeUntargeted"`, `"imaging"`, custom OPS `{batch}-{panel}`). Do not hand-roll modality→table maps — and never pass `"OPS"` as a table name.

**Imaging discovery:** `available_tables` may show `"imaging"` or legacy `"CellPaint"`. Both mean the imaging modality. Prefer loading with `table="imaging"` (Loader aliases either direction); `table="CellPaint"` also works when that suffix is listed.

---

## 2. obsm Key Roles

These keys follow the unified contract required by `pp.filter_cells` and `pp.normalize`. Do **not** use per-batch `control_variables_*` keys from older stores.

| Key | Written by | Role |
|---|---|---|
| `total_counts` | Store build | DataFrame, columns = batch names. Assigned counts per batch per cell. Used by count-modality filters / QC. |
| `unassigned_counts` | Store build | DataFrame, columns = batch names. Unassigned counts per batch; used with `total_counts` for assignment rate. |
| `control_variables` | Store build | DataFrame, columns like `NSB*`. Non-specific binding probes (count modalities). |
| `cell_geometry` | Store build | DataFrame: `Area`, `NuclearArea`, `AreaUm`, `Xum_tile`, `Yum_tile`, `AreaShape_*`, etc. **Not copied to `obs` by default.** Primary geometry source for imaging and for `filter_cells` area gates. |
| `spatial` | Store build | DataFrame of cell centroids `X`, `Y` in the region pixel frame. Convert to an `(n, 2)` ndarray before squidpy. |
| `data_missing` | Store build (optional) | Boolean DataFrame. Flags cells or batches where data is absent. |
| `qc_metrics` | `compute_qc_metrics` | Written before filtering and diagnostic plots. |
| `filter_flags` | `filter_cells` / diagnostics | Boolean DataFrame recording which filter each cell failed. |

### Notes on specific keys

**`total_counts`** — columns = batch names. Summing across columns gives per-cell assigned totals; this differs from `X.sum(axis=1)`, which includes control features.

**`control_variables`** — NSB-style probes. Positive-control signal (when present) lives as rows in `.X` where `var["control_type"]` matches the control label.

**`cell_geometry`** — not copied into `obs` automatically. `load_tables(..., include_cell_geometry=True)` merges columns into `.X`/`.var` and is intended for **imaging** tables (warns for other modalities).

**`qc_metrics` / `filter_flags`** — absent on freshly loaded tables; written after `pp.compute_qc_metrics` / `pp.filter_cells`.

---

## 3. Multi-Run Concatenation

### Use Loader, not `anndata.concat`

```python
loader = Loader(["/path/to/run1.zarr", "/path/to/run2.zarr"])
adata = loader.load_tables(table="ThreePrimeUntargeted", lazy=False)
# concatenated=True is the default — one stacked AnnData across both runs
# Pass concatenated=False for dict[str, AnnData] keyed by run
```

Do **not** hand-concatenate with `anndata.concat` — it bypasses `var`/`varm` merge semantics and provenance columns.

### Provenance columns added by Loader

When `concatenated=True` across multiple stores:

| Column | Content |
|---|---|
| `run_name` | Human-readable run identifier derived from the store path |
| `source_key` | Unique key identifying the source store within the loader |
| `source_store` | Full path (or URI) to the originating Zarr store |

These columns are absent when loading a single store.

### var / varm merge semantics

`Loader` applies merge rules for `.var` and `.varm`: features present in all runs are retained; features present in only a subset follow the `merge` policy (default: outer join with `NaN` fill). Differing `.uns` entries may be dropped on concat.

---

## 4. export_tables provenance (optional)

**Not part of the default workflow.** Only use when the user asks to write tables to disk. Normal analysis stays on `load_tables` without saving.

`Loader.export_tables(out_dir, tables=None, fmt="h5ad"|"zarr", ...)` writes selected modality tables as `{run}__{table}.{ext}` and stamps `uns["export_provenance"]`. Imagery is not exported.

```python
paths = loader.export_tables(
    "/path/to/out",
    tables=["ThreePrimeUntargeted", "imaging"],
    fmt="h5ad",
)
```

---

## 5. Old patterns: legacy barcoding stores

Older SDOs from barcoding runs carry these tables. Current runs do not emit them; load them only when they appear in `available_tables`. Analysis for these stores is in [legacy-targeted-pipeline](../legacy-targeted-pipeline/SKILL.md).

| Modality | `load_tables(table=...)` | Notes |
|---|---|---|
| Transcript | `"Transcript"` | Deprecated barcoding RNA — current runs do not emit |
| Protein | `"Protein"` | Deprecated barcoding protein — not the same as imaging IF (`ProteinIF`) |
