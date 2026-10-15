# SpatialData Loading and Access — Backends and Loading

Open for edge cases and variant details the SKILL.md workflow leaves out: store backends beyond a local `.zarr` or `.zarr.zip` (S3, several stores at once, an in-memory `SpatialData`), example discovery outputs, `RunInfo` inspection, per-modality load variants, checking a loaded table, and common pitfalls.

## Contents
- [Store backends and path formats](#store-backends-and-path-formats)
- [Inspect run metadata (RunInfo)](#inspect-run-metadata-runinfo)
- [Load modality tables (common entry point)](#load-modality-tables-common-entry-point)
- [Inspect a loaded AnnData](#inspect-a-loaded-anndata)
- [Load each modality](#load-each-modality)
- [Validation](#validation)
- [Common pitfalls](#common-pitfalls)

---

## Store backends and path formats

`Loader` auto-selects a backend based on the input type — no configuration needed beyond the path string.

> The `s3://` backends below need AWS credentials in the environment (`storage_options`). Without them (the ElemBio Cloud sandbox has none), mount or download the store first, as in [SKILL.md step 1](./SKILL.md#1-entry-a-confirmed-store-as-a-local-path).

| Input | Backend | Notes |
|---|---|---|
| `"/data/run.zarr"` | Local dir | Standard case |
| `"/data/run.zarr.zip"` | Local zip | Standard delivered format; reads parquet entries directly from the zip |
| `"s3://bucket/run.zarr"` | S3 dir | `s3fs` + `fsspec`; requires `storage_options` for auth |
| `"s3://bucket/run.zarr.zip"` | S3 zip | Selective remote extraction to local temp dir; only tables/metadata extracted (~1% of zip size) |
| `SpatialData` object | In-memory | Wrap an already-loaded `SpatialData`; rasters stay on `sdata` |

```python
from elembio_spatialdata_analysis.load import Loader

# Local dir
loader = Loader("/path/to/run.zarr")

# Local zip
loader = Loader("/path/to/run.zarr.zip")

# S3 dir — uses the default credential chain (env vars / instance role / ~/.aws/config)
loader = Loader("s3://bucket/path/to/run.zarr")

# S3 zip — tables extracted to local temp dir; points read via remote range requests
loader = Loader("s3://bucket/path/to/run.zarr.zip")

# Named profile override (only when the default chain doesn't have the right credentials)
loader = Loader("s3://bucket/path/to/run.zarr", storage_options={"profile": "your-profile"})

# SpatialData already in memory (e.g. after spatialdata.read_zarr)
import spatialdata as sd
sdata = sd.read_zarr("/path/to/run.zarr")
loader = Loader(sdata)
```

### Common `storage_options` keys for S3

| Key | Purpose | Example |
|---|---|---|
| `profile` | AWS CLI named profile | `"your-profile"` |
| `region` | AWS region override | `"us-west-2"` |

`storage_options` is ignored for local and in-memory backends. Prefer the environment’s default credential chain when a named profile is not required.

### Multi-run and mixed backends

<!-- canonical: spatialdata-loading-and-access memory guard. Copy verbatim; edit only in that skill's SKILL.md -->
**Memory guard — before any `lazy=False`, `.to_memory()` or `.compute()`:**

1. Load lazily first and run `check_before` from `scripts/check_load.py` in `spatialdata-loading-and-access`. Run the script unchanged in the Python session that holds the data (`exec` its contents); don't rewrite it. It compares each table's estimated in-memory size with the memory this session can use, lists large session objects you could `del`, and raises on a FAIL.
2. Materialize only the tables and runs the plan needs. With several runs, load `concatenated=False` so it checks each run and their total. A WARN says what to `del` or to load one run at a time. A FAIL means subset before materializing: propose a subset to the user (runs, wells, or a set of features), then slice the lazy table and call `.to_memory()` on the slice, for example `adata[adata.obs["Well"].isin(wells)].to_memory()`. `.copy()` keeps a lazy table lazy. Feature subsets by modality: [ACCESS.md — Subset](../spatialdata-loading-and-access/ACCESS.md#subset-by-well-or-feature).
3. Keep `.X` in the form it was stored (`scipy.sparse.issparse(adata.X)`). Sparse tables (3′, OPS) are mostly zeros: never densify a whole one (`.toarray()`, `.todense()`, `np.asarray(adata.X)`); densify only a subset of cells or features. Dense tables (imaging, legacy barcoding) are float32 with `NaN` where a batch has no data: do not convert them to sparse, which turns missing into zero.
<!-- /canonical -->

Pass a list of paths (any mix of local, zip, S3) — `Loader` handles each with the appropriate backend:

```python
loader = Loader([
    "s3://bucket/path/to/run1.zarr",
    "/path/to/run2.zarr",
])
# Pass storage_options only if S3 needs a named profile:
# loader = Loader([...], storage_options={"profile": "your-profile"})

runs = loader.load_tables(table="ThreePrimeUntargeted", concatenated=False)  # lazy, one per run
check_before(runs)  # each run and their total; raises on FAIL
adata = loader.load_tables(table="ThreePrimeUntargeted", lazy=False)
# adata.obs["run_name"] identifies the source run for each cell
```

---

## Inspect run metadata (RunInfo)

`Loader._runs` is a `dict[str, RunInfo]` — one entry per store. `RunInfo` carries all run-level metadata as lazily cached properties. Inspect it before loading to confirm what's available and catch mismatches early.

```python
loader = Loader("/path/to/run.zarr")

# Human-readable summary (run name, tables, batches, wells, tiles, points mode, build version)
print(loader.summary())

# Access the single RunInfo for a one-run loader
run = next(iter(loader._runs.values()))

# Or by run name (key = run name with __N dedup suffix if needed)
# run = loader._runs["<run-name>"]

# What's available — current product examples (not all are always present)
print(run.available_tables)        # e.g. {'ThreePrimeUntargeted', 'imaging', 'B01-<TargetSite>'}
print(run.available_target_types)  # e.g. {'ThreePrimeUntargeted', 'imaging', 'OPS'}
# Older SDOs may still list 'CellPaint' (load with table="imaging" — aliased)
# and/or legacy {'Transcript', 'Protein', 'imaging'}. More example sets below.

# Well / tile / batch structure
print(run.well_names)   # ['A1', 'A2', 'B1', 'B2']
print(run.batch_names)  # e.g. DISS batches + 'CP01' (+ optional 'P01' IF)
print(run.batch_types)  # batch_name → batch_type / class

# Full DataFrames from the metadata table
print(run.wells)    # indexed by well name (A1): well_index, well_label, cell_type, ...
# Cell tables use obs["Well"] (position) and obs["WellLabel"] (sample-sheet name)
print(run.tiles)    # DataFrame indexed by tile_name: tile_index, x_mm, y_mm, well_name, ...
print(run.batches)  # DataFrame: batch_name, batch_class, batch_type, num_targets, ...

# Target panel info — DISS vs legacy barcoding are not both present
print(run.diss_targets_metadata)  # current DISS / OPS runs
print(run.targets_imaging)        # imaging targets (CellPaint stains ± optional ProteinIF)
# Legacy barcoding SDOs only:
# print(run.targets_barcoding)

# Build version
print(run.build_provenance)    # {'pipeline_version': '<version>', ...}
```

**Example `available_tables` sets** (table suffixes; OPS tables use custom names):

- `{'ThreePrimeUntargeted', 'imaging'}`
- `{'B01-<TargetSite>', 'imaging'}` — OPS custom suffix + imaging
- `{'ThreePrimeUntargeted', 'B01-...', 'CellPaint'}` — older stores list `CellPaint` instead of `imaging`
- Legacy barcoding SDOs: `{'Transcript', 'Protein', 'CellPaint'}` or `{'Transcript', 'Protein', 'imaging'}` — see [REFERENCE.md — Old patterns](./REFERENCE.md#5-old-patterns-legacy-barcoding-stores)

`available_target_types` may include `'ThreePrimeUntargeted'`, `'imaging'`, and `'OPS'`. These are target types, not table names: with the default build, an OPS table is named `{batch_name}-{target_site}`.

**Note:** `RunInfo` is a private attribute (`_runs`). The public surface is `loader.available_tables`, `loader.available_target_types`, and `loader.summary()`. Use `_runs` when you need per-run metadata like well lists, tile layouts, or the target panel.

---

## Load modality tables (common entry point)

For any workflow that only needs AnnData (pp, QC, normalise, Scanpy): use `Loader`. No need for `spatialdata.read_zarr`.

```python
from elembio_spatialdata_analysis.load import Loader

loader = Loader("/path/to/run.zarr")
print(loader.summary())
print("Available tables:", sorted(loader.available_tables))
print("Target types:", sorted(loader.available_target_types))
# e.g. current: ['B01-<TargetSite>', 'ThreePrimeUntargeted', 'imaging']
# e.g. older imaging table name listed as 'CellPaint' — still load with table="imaging"

# Defaults: lazy=True, concatenated=True — single stacked AnnData, lazy-backed .X
adata = loader.load_tables(table="ThreePrimeUntargeted")

# Pass lazy=False before pp / Scanpy (sc.pp.pca, neighbors, UMAP, HVG, pp.filter_cells, …)
adata = loader.load_tables(table="ThreePrimeUntargeted", lazy=False)
```

**`concatenated`:** default `True` returns one stacked `AnnData`. Pass `concatenated=False` only when you need a `dict[str, AnnData]` keyed by run.

**`lazy`:** default `True` is preferred for large tables. Pass `lazy=False` when downstream code requires a fully materialised `.X` — lazy `AnnData` silently returns delayed arrays that break `sc.pp.*`, neighbors, UMAP, and most `pp` functions.

### Export tables without imagery (optional)

**Not required for analysis.** Do not call `export_tables` on every load. Keep working from `load_tables` in memory unless the user explicitly wants table files on disk (e.g. pull only AnnData from a remote/zip store without copying imagery).

When needed, `Loader.export_tables` writes selected modality tables to local `h5ad` or `zarr` and stamps `uns["export_provenance"]`:

```python
paths = loader.export_tables(
    "/path/to/out",
    tables=["ThreePrimeUntargeted", "imaging"],  # default: all available_tables
    fmt="h5ad",
)
```

---

## Inspect a loaded AnnData

```python
adata = loader.load_tables(table="ThreePrimeUntargeted")  # lazy=True default; add lazy=False before pp/Scanpy

print(f"Shape: {adata.shape}")         # (n_cells, n_features)
print(f"Features: {adata.var_names[:10].tolist()}")

print("\n.obs columns (Well, tile_id, X_tile, Y_tile, cell IDs):")
print(adata.obs.columns.tolist())

print("\n.var columns:")
print(adata.var.columns.tolist())

print("\n.obsm keys:")
print(list(adata.obsm.keys()))
# Common: cell_geometry, spatial (region-pixel X/Y), spatial_artifacts; count tables may add total_counts
# After QC: also qc_metrics, filter_flags

print("\n.layers keys:")
print(list(adata.layers.keys()))   # e.g. 'nuclear' on some count tables; imaging has none

print("\n.uns keys:")
print(list(adata.uns.keys()))      # 'metadata_attrs', control metadata, etc.
```

---

## Load each modality

### 3′ transcriptome (standard RNA)

3′ untargeted; select HVGs before dimensionality reduction. See [3prime-untargeted-pipeline](../3prime-untargeted-pipeline/SKILL.md).

```python
if "ThreePrimeUntargeted" in loader.available_tables:
    rna = loader.load_tables(table="ThreePrimeUntargeted", lazy=False)
    print(f"3′ shape: {rna.shape}")
```

### OPS

OPS tables have **custom names**: one table per batch and target site, `{batch_name}-{target_site}` by default. Find them from table metadata rather than assuming a name. Always pass `include_target_assignment=True`. Prefer pairing with 3′ RNA when both are present.

Discover suffixes from table metadata (accepts both `"OPS"` and `"SpecializedTargeted"` aliases — do not use `startswith` on batch names):

```python
run = next(iter(loader._runs.values()))
ops_tables = run.table_suffixes_for_target_type(("OPS", "SpecializedTargeted"))
# e.g. ["B01-<TargetSite>", "B01-<ControlSite>"]
```

**Control sites.** When there is more than one OPS table, read each one's control-site flag. A lazy load reads metadata, not counts:

```python
control = {
    name: loader.load_tables(table=name, lazy=True).uns.get("metadata_attrs", {}).get("is_control_site", False)
    for name in ops_tables
}
assay_tables = [name for name, is_control in control.items() if not is_control]
```

`True` marks a control site (for example a housekeeping panel), not the assay. Stores built from older panels read `False` for every table, so control sites cannot be told apart there. When more than one table remains in `assay_tables`, the confirmation gate asks which one to analyze.

```python
if assay_tables:
    ops_table_name = assay_tables[0]  # the only one, or the one the user chose at the gate
    print(f"OPS table: {ops_table_name!r}")

    ops_adata = loader.load_tables(
        table=ops_table_name,
        lazy=False,
        include_target_assignment=True,
    )
    print(ops_adata.obs[["MostAbundantTarget", "TotalCount", "Confidence"]].head())

    # Merge OPS identity into paired 3′ RNA (before RNA filtering) when present
    if "ThreePrimeUntargeted" in loader.available_tables:
        rna = loader.load_tables(table="ThreePrimeUntargeted", lazy=False)
        assert ops_adata.obs_names.equals(rna.obs_names), "obs_names mismatch!"
        rna.obs["OPS_line"] = ops_adata.obs["MostAbundantTarget"].values
        rna.obs["OPS_Confidence"] = ops_adata.obs["Confidence"].values
```

See [ops-pipeline](../ops-pipeline/SKILL.md) for the full workflow.

### Imaging

Imaging `.X` is continuous intensities / CellProfiler morphology features — **not** counts. Do not call count-style normalize on it.

Prefer `table="imaging"`. Loader aliases `imaging` ↔ `CellPaint` either direction, so one call works for current (`global-table-imaging`) and legacy (`global-table-CellPaint`) stores.

**Discovery gotcha:** `available_tables` reports the store’s real suffix. Legacy SDOs list `"CellPaint"`, not `"imaging"`. That entry **is** imaging — do not treat missing `"imaging"` as “no imaging modality.” Downstream code should map either listed name to the imaging path.

```python
has_imaging = bool({"imaging", "CellPaint"} & loader.available_tables)
if has_imaging:
    # Prefer the canonical name; Loader resolves CellPaint stores via alias
    img = loader.load_tables(table="imaging", lazy=False)
    # Equally valid on a legacy store:
    # img = loader.load_tables(table="CellPaint", lazy=False)
    print(f"Imaging shape: {img.shape}")
    print(img.var[["target_type", "feature_class", "batch_name", "suffix"]].head())

    # Morphology stains (CP*) vs optional IF markers (P* / ProteinIF) — IF may be absent
    print(img.var["target_type"].value_counts(dropna=False))
    print(img.var["batch_name"].value_counts())
```

Analysis: [imaging-pipeline](../imaging-pipeline/SKILL.md). Overlays: [spatial-multimodal-visualization](../spatial-multimodal-visualization/SKILL.md). Cell QC: `target_type="Imaging"`.

### Legacy barcoding (Transcript / Protein)

Only on older SDOs. Current runs do not emit these tables. Barcoding `Protein` ≠ imaging `ProteinIF`.

```python
if "Transcript" in loader.available_tables:
    transcript = loader.load_tables(table="Transcript", lazy=False)
if "Protein" in loader.available_tables:
    protein = loader.load_tables(table="Protein", lazy=False)
```

---

## Validation

```python
assert len(sdata.images) > 0, "No images"
assert len(sdata.tables) > 0, "No tables"
assert any(k.startswith("global-table-") for k in sdata.tables), "No global tables"

def validate_adata(adata, modality="ThreePrimeUntargeted"):
    assert adata.n_obs > 0 and adata.n_vars > 0
    assert "Well" in adata.obs.columns
    assert "batch_name" in adata.var.columns
    if modality in ("ThreePrimeUntargeted", "Transcript", "Protein"):
        assert "total_counts" in adata.obsm
    elif modality in ("imaging", "CellPaint"):
        assert "cell_geometry" in adata.obsm
        # Optional IF: check var, do not require ProteinIF
        print(adata.var["target_type"].value_counts(dropna=False))
    print(f"OK: {modality} shape {adata.shape}")
```

---

## Common pitfalls

1. **Assuming modalities are always present** — Inspect `available_tables` / `available_target_types` first. A run may have any subset of 3′, OPS, and imaging.

2. **`lazy=True` in pipelines** — Fine for exploration. Pass `lazy=False` before `pp`, `qc`, or Scanpy (`sc.pp.pca`, neighbors, UMAP).

3. **`Loader(path)` has no rasters** — Path-backed `Loader` only implements tables + points + metadata. For images or labels, call `sd.read_zarr(path)` first, then `Loader(sdata)`.

4. **Imaging table name** — Prefer `table="imaging"`. Loader aliases `imaging` ↔ `CellPaint`. If discovery lists only `"CellPaint"`, that is imaging (legacy suffix) — still load it; do not conclude imaging is missing.

5. **IF Protein is optional** — Imaging may be morphology-only (`CP*` / CellPaint stains). Check `var["target_type"]` / `batch_name` before assuming `ProteinIF` / `P*` batches exist.

6. **Do not normalize imaging as counts** — `pp.normalize_imaging` is unimplemented; fluorescence is not count data.

7. **OPS is never `table="OPS"`** — Use `run.table_suffixes_for_target_type(("OPS", "SpecializedTargeted"))`; pass `include_target_assignment=True`.

8. **`cell_geometry` not in `obs`** — Lives in `obsm["cell_geometry"]`. `filter_cells` reads there first.

9. **Legacy vs DISS metadata** — Current DISS runs use `diss_targets_metadata` (+ `targets_imaging`); legacy barcoding uses `targets_barcoding`. They are not both present.

10. **Multi-run concatenation** — Use `Loader([path1, path2, ...])`, not `anndata.concat`.

11. **Points are lazy** — `load_points` returns a Polars `LazyFrame`; call `.collect()` before pandas/numpy ops.

12. **Missing elements after `filter_by_coordinate_system`** — Only elements registered to that coordinate system are returned. Global tables may not be included.

13. **`read_zarr` can't open a `.zarr.zip`.** It needs a *directory* store; passing a `.zarr.zip` (or a `zarr.ZipStore` group) fails on zarr 3.x. FUSE-mount the zip with `ratarmount` first, see [Reading a `.zarr.zip` with `read_zarr` (FUSE mount)](./SPATIAL.md#reading-a-zarrzip-with-read_zarr-fuse-mount). `Loader` reads a zip natively for tables/points; this only bites when you need a `SpatialData` object from a zip.
