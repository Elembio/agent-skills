---
name: spatialdata-loading-and-access
description: Opens a confirmed AVITI24 SpatialData store (.zarr or .zarr.zip) for analysis — the per-modality AnnData tables, the polony points, or the full SpatialData object with its images, segmentation labels and shapes. Chooses between Loader and spatialdata.read_zarr, finds which modality tables the run actually contains, and loads only what the next step needs. Use whenever code opens a store or calls Loader or load_tables, after locate-spatialdata-store has validated the store path, including a path the user gave. Not for finding a run's store (locate-spatialdata-store) or the full AnnData slot inventory (anndata-structure-reference).
metadata:
  version: 0.5.8
  author: elembio
---

# SpatialData Loading and Access

Turn a confirmed AVITI24 store into the per-modality AnnData tables (and, when needed, the rasters) that the next skill works on.

## Workflow

Copy this checklist and tick it off as you go:

```
- [ ] 1. Entry    — confirmed store(s) from locate-spatialdata-store, as local paths
- [ ] 2. Discover — summary, available_tables, map each table to a modality
- [ ] 3. Plan     — method, tables per run, rasters; confirm once with the user
- [ ] 4. Load     — size check first, then only the planned tables
- [ ] 5. Verify   — check what loaded, then hand off to the next skill
```

### 1. Entry: a confirmed store, as a local path

Start from the store(s) that [locate-spatialdata-store](../locate-spatialdata-store/SKILL.md) confirmed with the user. Its handoff can take more than one form. Make each store a local `.zarr` or `.zarr.zip` path before step 2:

| You have | Do |
|---|---|
| A local or mounted `.zarr` / `.zarr.zip` path | Use it as is. |
| An `s3://` URI, or a run or execution plus the store's path in it | Mount it, or download only the store, with [elembio-cloud-data-access](../elembio-cloud-data-access/SKILL.md). For a raw `s3://` URI, first [resolve its storage connection](../elembio-cloud-data-access/SKILL.md#resolving-a-storage-uri). Do not pass `s3://` to `Loader` unless AWS credentials are set; the ElemBio Cloud sandbox has none. |
| No store, a run name, a directory, or a path that is not a `.zarr` / `.zarr.zip` | Go to locate-spatialdata-store. Do not search for or pick a store here. |

`Loader` auto-selects from the path — local dir (`.zarr`), local zip (`.zarr.zip`), or an existing `SpatialData` object. For several stores at once, or `s3://` with credentials, open [IMPLEMENTATION.md — Store backends](./IMPLEMENTATION.md#store-backends-and-path-formats). Requires Python 3.12+ (matches `elembio-spatialdata-analysis`).

### 2. Discover what the run has

A run carries **whatever combination it was configured for** — none of 3′, OPS, or imaging is universal. Always inspect before loading:

```python
from elembio_spatialdata_analysis.load import Loader

loader = Loader("/path/to/run.zarr")
print(loader.summary())
print("Available:", sorted(loader.available_tables))
print("Types:", sorted(loader.available_target_types))
```

Map each table to a modality with [Modality table names](#modality-table-names). Two names need care:

- **`"CellPaint"`** — if `available_tables` shows `"CellPaint"` and not `"imaging"`, that **is** the imaging table (legacy name) — do not conclude imaging is absent, and do not wait for a separate `"imaging"` entry. Load with `table="imaging"` (or `table="CellPaint"`; both work). Treat `"CellPaint"` in the table list as imaging everywhere downstream (QC, viz, pipelines).
- **OPS** — discover suffixes from metadata (handles both `"OPS"` and `"SpecializedTargeted"` aliases):

  ```python
  run = next(iter(loader._runs.values()))
  ops_tables = run.table_suffixes_for_target_type(("OPS", "SpecializedTargeted"))
  # e.g. ["B01-<TargetSite>", "B01-<ControlSite>"]
  ```

  `loader._runs` is private, but it is the access path the library itself documents for per-run `RunInfo`, and no public alternative exists yet. Use it as shown.

  With more than one OPS table, find the control sites before planning: [IMPLEMENTATION.md — OPS](./IMPLEMENTATION.md#ops).

**Imaging content varies by run:**

- **Morphology-only** — typically `CP*` batches / `var["target_type"] == "CellPaint"` (e.g. four structural stains). IF is not assumed.
- **Morphology + IF** — `CP*` stains plus `P*` antibody markers (`var["target_type"] == "ProteinIF"`). Same table; inspect `var` / `batch_name`.

### 3. Plan the load, then confirm once

Choose the method from what the request needs:

| You need | Use |
|---|---|
| AnnData tables — QC, normalize, Scanpy DR, run analysis | **`Loader(path).load_tables(table=...)`** (default) |
| Polony points without rasters | `Loader(path).load_points(batch=..., well=..., target_type=...)` |
| Images, labels, shapes, napari, `filter_by_coordinate_system` | `spatialdata.read_zarr(path)` on a **directory** store |
| Both rasters AND table/point API | `sdata = sd.read_zarr(dir_path)` → `loader = Loader(sdata)` |
| Export tables only (no imagery) — **optional** | `Loader(path).export_tables(out_dir, fmt="h5ad")` |

**Default for “analyze this run”:** `Loader` → discover tables → load modalities. Do **not** start with `read_zarr` unless the user needs imagery.

**`Loader(path)` does not give you `sdata.images`** — path backends only implement tables + points + metadata. For rasters, call `sd.read_zarr` on a directory store (or a FUSE-mounted zip).

**`.zarr.zip`:** `Loader` opens zip natively for tables/points. `read_zarr` **cannot** open a zip directly — when rasters from a zip are required, FUSE-mount it read-only with `ratarmount` (no extraction) and `read_zarr` the nested `.zarr` inside the mount. `ratarmount` needs FUSE (Linux `fuse3` or macOS macFUSE) and may not be installed — check `shutil.which("ratarmount")` first. Without it, read only the elements you need, and extract only with the user's agreement. Recipe and fallbacks: [SPATIAL.md — FUSE mount](./SPATIAL.md#reading-a-zarrzip-with-read_zarr-fuse-mount).

**Confirm the plan with the user in one message before you load anything for analysis.** The gate applies to every load, whatever its size, and also when the request already names the modality. State the plan: the method, which tables from which runs, whether rasters are needed (and for a `.zarr.zip`, that this means a FUSE mount and an eager `read_zarr`), and how each table will load. Say why, and name any assumption you made about the request. Then wait.

- **Ask only about intent.** Do not ask about anything discovery already settled (`CellPaint` is imaging; OPS tables have custom names), and never re-ask which store — locate-spatialdata-store confirmed it.
- **Name only what the plan uses.** Do not list modalities the run lacks unless the user asked for one.
- **Several OPS tables:** plan the ones not flagged as control sites (see [Discover](#2-discover-what-the-run-has)). Ask which to analyze only when more than one remains.
- **A modality the user asked for is absent:** say what the run does have. Do not substitute another modality.
- **No user to ask** (subagent, headless run): write the plan and its assumptions into your output, then proceed.
- Plan only the load. Reports, figures and decks belong to the skills after this one.

### 4. Load

<!-- canonical: spatialdata-loading-and-access memory guard. Copy verbatim; edit only in that skill's SKILL.md -->
**Memory guard — before any `lazy=False`, `.to_memory()` or `.compute()`:**

1. Load lazily first and run `check_before` from `scripts/check_load.py` in `spatialdata-loading-and-access`. Run the script unchanged in the Python session that holds the data (`exec` its contents); don't rewrite it. It compares each table's estimated in-memory size with the memory this session can use, lists large session objects you could `del`, and raises on a FAIL.
2. Materialize only the tables and runs the plan needs. With several runs, load `concatenated=False` so it checks each run and their total. A WARN says what to `del` or to load one run at a time. A FAIL means subset before materializing: propose a subset to the user (runs, wells, or a set of features), then slice the lazy table and call `.to_memory()` on the slice, for example `adata[adata.obs["Well"].isin(wells)].to_memory()`. `.copy()` keeps a lazy table lazy. Feature subsets by modality: [ACCESS.md — Subset](../spatialdata-loading-and-access/ACCESS.md#subset-by-well-or-feature).
3. Keep `.X` in the form it was stored (`scipy.sparse.issparse(adata.X)`). Sparse tables (3′, OPS) are mostly zeros: never densify a whole one (`.toarray()`, `.todense()`, `np.asarray(adata.X)`); densify only a subset of cells or features. Dense tables (imaging, legacy barcoding) are float32 with `NaN` where a batch has no data: do not convert them to sparse, which turns missing into zero.
<!-- /canonical -->

- `lazy=True` — default; fine for exploration. Pass `lazy=False` before `pp`, Scanpy (`sc.pp.pca`, neighbors, UMAP, HVG), or filtering.
- `concatenated=True` — default; single stacked `AnnData`. Pass `False` for `dict[str, AnnData]` keyed by run.

```python
planned = ["ThreePrimeUntargeted", "imaging"]  # exactly the tables confirmed in step 3
tables = {}
for name in planned:
    runs = loader.load_tables(table=name, concatenated=False)  # lazy, one AnnData per run
    check_before(runs)  # raises on FAIL -> subset first; WARN -> do what it says
    tables[name] = loader.load_tables(table=name, lazy=False)

# When you also need images / labels:
import spatialdata as sd
sdata = sd.read_zarr("/path/to/run.zarr")
loader = Loader(sdata)
```

Rarely needed: `include_cell_geometry` (imaging only) is described under `cell_geometry` in [REFERENCE.md — obsm key roles](./REFERENCE.md#2-obsm-key-roles); `validate_points` under [SPATIAL.md — points](./SPATIAL.md#access-points-polonies-via-loader).

### 5. Verify, then hand off

Run `check_after` on every table you loaded, before passing it on:

```python
check_after(tables)  # expects materialized; use expect="lazy" for exploration
```

- **`.X` still lazy** where you meant to materialize: reload with `lazy=False`.
- **`.X` in memory where you asked for `lazy=True`:** the lazy read failed and the library read the whole table eagerly, without a warning. Count that memory before you load anything else.
- **Opened the full SpatialData object:** `print(sdata)` and confirm the images, labels, points or shapes you planned are present.
- **No cells:** stop and report it; the store or the filters are wrong.
- Tell the next skill what you loaded: tables, runs, lazy or materialized, and any assumption the user confirmed.

Next: [run-quality-overview](../run-quality-overview/SKILL.md) → [cell-quality-control](../cell-quality-control/SKILL.md) → normalize count modalities only → pipeline(s) for what is present. For image or spatial requests, [spatial-multimodal-visualization](../spatial-multimodal-visualization/SKILL.md). Code order: [code-generation-guidelines](../code-generation-guidelines/SKILL.md).

## Modality table names

| Modality | `load_tables(table=...)` | Notes |
|---|---|---|
| 3′ transcriptome | `"ThreePrimeUntargeted"` | Standard RNA; 3′ untargeted counts; select HVGs before DR |
| OPS / specialized | **custom** — e.g. `"B01-<TargetSite>"` | One table per batch and target site, named `{batch_name}-{target_site}` by default. Find them with `table_suffixes_for_target_type` (see [Discover](#2-discover-what-the-run-has)) rather than assuming a name. Pass `include_target_assignment=True`. |
| Imaging | `"imaging"` | Morphology (CellPaint stains) ± optional IF (`ProteinIF`). Continuous floats, not counts |
| Legacy barcoding (older stores only) | `"Transcript"`, `"Protein"` | Current runs do not emit these. See [REFERENCE.md — Old patterns](./REFERENCE.md#5-old-patterns-legacy-barcoding-stores) |

## Detail Files

- [IMPLEMENTATION.md](./IMPLEMENTATION.md) — open for edge cases and variant details this workflow leaves out: backends beyond a local `.zarr` or `.zarr.zip` (S3, several stores at once, an in-memory `SpatialData`), example discovery outputs, `RunInfo` inspection, per-modality load variants, and common pitfalls
- [ACCESS.md](./ACCESS.md) — open once a table is loaded and you need `.X`, `.obs`, `.var`, `.obsm` controls or `.layers`, or a well or feature subset
- [SPATIAL.md](./SPATIAL.md) — open when the request needs images, labels, points, coordinate systems or spatial subsetting, including rasters from a `.zarr.zip` (FUSE)
- [REFERENCE.md](./REFERENCE.md) — open when you need what an `obsm` key holds, how multi-run concatenation merges, `export_tables` provenance, or a legacy barcoding store (`Transcript` / `Protein`)

## Related Skills

- [locate-spatialdata-store](../locate-spatialdata-store/SKILL.md) — upstream: resolves and confirms the store this skill starts from
- [elembio-cloud-data-access](../elembio-cloud-data-access/SKILL.md) — mounts or downloads a confirmed cloud store so `Loader` can open it
- [metadata-and-target-mapping](../metadata-and-target-mapping/SKILL.md) — wells, tiles, batches, targets from RunInfo / DISS metadata
- [run-quality-overview](../run-quality-overview/SKILL.md) — run / well / tile QC after loading
- [cell-quality-control](../cell-quality-control/SKILL.md) — per-cell filtering after run QC
- [3prime-untargeted-pipeline](../3prime-untargeted-pipeline/SKILL.md) — standard RNA workflow when 3′ is present
- [ops-pipeline](../ops-pipeline/SKILL.md) — OPS pooled screens
- [imaging-pipeline](../imaging-pipeline/SKILL.md) — morphology / IF analysis after load
- [spatial-multimodal-visualization](../spatial-multimodal-visualization/SKILL.md) — image / spatial overlays
- [anndata-structure-reference](../anndata-structure-reference/SKILL.md) — full slot inventory

## Disclaimer

This content is provided "as-is" without warranty of any kind. AI-generated output may contain
inaccuracies. Review and verify all outputs before relying on them.
