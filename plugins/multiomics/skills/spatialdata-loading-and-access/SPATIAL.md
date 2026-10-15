# SpatialData Loading and Access — Rasters, Points and Spatial

Open when the request needs images, labels, polony points, coordinate systems or spatial subsetting — including rasters from a `.zarr.zip` (FUSE mount).

## Contents
- [Load with full SpatialData (images / labels / spatial)](#load-with-full-spatialdata-images--labels--spatial)
- [Reading a .zarr.zip with read_zarr (FUSE mount)](#reading-a-zarrzip-with-read_zarr-fuse-mount)
- [Access images](#access-images)
- [Access labels (segmentation masks)](#access-labels-segmentation-masks)
- [Access points (polonies) via Loader](#access-points-polonies-via-loader)
- [Coordinate systems](#coordinate-systems)
- [Bounding box spatial subsetting](#bounding-box-spatial-subsetting)
- [Linking elements across modalities](#linking-elements-across-modalities)
- [Memory management for large stores](#memory-management-for-large-stores)

---

## Load with full SpatialData (images / labels / spatial)

Only when you need rasters, spatial filtering, or napari.

```python
import spatialdata as sd

# Instant — reads metadata only; data is dask-backed and lazy
sdata = sd.read_zarr("/path/to/run.zarr")

print(f"Images: {list(sdata.images.keys())}")
print(f"Labels: {list(sdata.labels.keys())}")
print(f"Points: {list(sdata.points.keys())}")
print(f"Tables: {list(sdata.tables.keys())}")

# Then wrap with Loader if you also need load_tables / load_points
from elembio_spatialdata_analysis.load import Loader
loader = Loader(sdata)
adata = loader.load_tables(table="ThreePrimeUntargeted", lazy=False)
# Rasters stay on sdata; Loader does not add a separate raster API.
# Guard on available_tables — load imaging / OPS / legacy only if present.
```

Selective loading (skip expensive element types):

```python
sdata = sd.read_zarr(path, selection=("images", "tables"))
sdata = sd.read_zarr(path, on_bad_files="warn")   # "error" | "warn" | "skip"
```

---

## Reading a `.zarr.zip` with `read_zarr` (FUSE mount)

`spatialdata.read_zarr` reads a **directory** store. It cannot open a `.zarr.zip` directly, and the no-FUSE workarounds fail on the current stack (zarr 3.x / spatialdata 0.7.x):

- `zarr.storage.ZipStore` → `open_group(path=...)` → `sd.read_zarr(group)` raises `AttributeError: module 'zarr.storage' has no attribute 'ConsolidatedMetadataStore'`.
- Bare `fsspec ZipFileSystem` + `zarr FsspecStore` fails the zarr-v3 async-filesystem check.
- `Loader` opens a `.zarr.zip` for **tables/points** but returns AnnData/points, not a `SpatialData` (no images/labels/shapes).

So when you need a `SpatialData` object from a zip, FUSE-mount the zip read-only with `ratarmount` (it indexes the zip's central directory via range reads — **no extraction**) and point `read_zarr` at the mountpoint. Conversely, once you hold the full `SpatialData`, `Loader(sdata)` exposes the same `load_tables`/`load_points` API — so the mounted full object is not a separate table/point source, and Loader-tailored downstream tools run on it unchanged.

### Mount the zip (FUSE + `ratarmount`)

`ratarmount` and a working FUSE layer (Linux `fuse3`; macOS macFUSE) are **not part of the analysis dependency stack** — they are not guaranteed to be present. Confirm the capability first; if absent, install it or fall back to [Targeted reads](#targeted-reads--skip-read_zarr-when-you-only-need-a-few-elements).

```python
import shutil
assert shutil.which("ratarmount"), "ratarmount missing — `pip install ratarmount`, or use targeted reads"
```

```bash
pip install ratarmount
# macOS arm64 only — point the binding at the universal libfuse, or the dlopen fails:
export FUSE_LIBRARY_PATH=/usr/local/lib/libfuse.2.dylib
ratarmount /path/to/run.zarr.zip /mnt/store         # /mnt/store/run.zarr now browsable
```

```python
import spatialdata as sd
sdata = sd.read_zarr("/mnt/store/run.zarr")          # NESTED store path, not the mount root
```

Unmount: `umount /mnt/store` (macOS) or `fusermount -u /mnt/store` (Linux). If the store is remote, localize or mount it first ([`elembio-cloud-data-access`](../elembio-cloud-data-access/SKILL.md)); a managed runtime may pre-provide `ratarmount` — follow its guidance when present.

### The `read_zarr` open is eager and expensive

`read_zarr` is **eager on open** — it walks metadata across every element (measured ~130 s cold / ~40 s warm on a 35 GB, ~47k-entry store; ~4 GB peak RAM) even though array data stays lazy (dask) afterward. Both the open time and the peak RAM **scale with the number of elements** — for images that is wells × channels × pyramid levels, plus tiles, labels, points, and tables — so a larger store costs proportionally more than the benchmark above.

- **If you need only part of the store, do not open the full `SpatialData`** — use targeted reads (next section). The full open is the costly step; most questions don't require it.
- **In a managed runtime with a tool-call timeout or memory ceiling** (e.g. a code-interpreter sandbox), the eager open can exceed the timeout or approach the memory limit, and crossing a hard ceiling may not fail gracefully. Follow that environment's guidance for backgrounding the open and budgeting RAM.

### Targeted reads — skip `read_zarr` when you only need a few elements

Opening the full `SpatialData` is the costly part. If you only need one image or a per-gene count, read directly off the mount (each is fast and pulls only the bytes touched):

```python
import zarr
arr = zarr.open("/mnt/store/run.zarr/images/<well>-images-cell_paint/5")[:]   # coarsest pyramid level
```

For a per-gene point count, a full dask compute can OOM — use a synchronous scheduler, project the single column, and loop per part with gc:

```python
import dask, gc, glob, pyarrow.parquet as pq
dask.config.set(scheduler="synchronous")
n = 0
for part in sorted(glob.glob("/mnt/store/run.zarr/points/A1-points-B0*/points.parquet/part.*.parquet")):
    t = pq.read_table(part, columns=["target_name"])          # gene column is `target_name`
    n += t.column("target_name").to_pylist().count("GAPDH")
    gc.collect()
```

### Troubleshooting

| Symptom | Cause / fix |
|---|---|
| `GroupNotFoundError` / empty dir | Pointed at the mount root; the store is at `<mount>/<name>.zarr` (nested). |
| macOS `dlopen … need arm64` / `No module named fuse` | Binding loaded x86 `libosxfuse`; `export FUSE_LIBRARY_PATH=/usr/local/lib/libfuse.2.dylib` before mounting. |
| `read_zarr` hangs / slow on first open | Slow-cold (~minutes) is expected — see "The `read_zarr` open is eager and expensive" above. |
| Remount: "already running (pid …)" after unmount | Stale `ratarmount` lock; mount at a fresh path. |
| `points.parquet` "is a directory" | It's a directory of `part.N.parquet` — glob the parts. |

### Mount first; extract only with the user's agreement

- **Mount by default.** `ratarmount` range-reads the zip in place: nothing is copied, and only the bytes you touch are read.
- **When FUSE is unavailable** and cannot be installed, try [targeted reads](#targeted-reads--skip-read_zarr-when-you-only-need-a-few-elements) first. Extract only if the request needs a full `SpatialData` object.
- **Ask before extracting.** A store is often close to 100 GB, so extraction takes a long time and needs a second full-size copy on disk. State the size and the free disk, and wait for the user to agree. After extraction, `read_zarr` the nested `<name>.zarr` directory inside the extracted folder, the same nesting as the mount.

---

## Access images

```python
image = sdata.images["L1R01C01S1-images-cell_paint"]
print(f"Shape: {image.shape}")           # (channels, y, x), dask-backed
print(f"Channels: {image.coords['c'].values}")

channel = image.sel(c="Mitochondria")    # select by channel name
subset  = image[:, 100:200, 100:200].compute()   # load partial region
```

---

## Access labels (segmentation masks)

```python
import numpy as np

cell_labels = sdata.labels["L1R01C01S1-labels-Cell"]
labels_array = cell_labels.to_numpy()
print(f"Cells: {len(np.unique(labels_array)) - 1}")  # exclude background (0)

if "L1R01C01S1-labels-Nuclear" in sdata.labels:
    nuclear_labels = sdata.labels["L1R01C01S1-labels-Nuclear"]
```

---

## Access points (polonies) via Loader

`Loader.load_points` returns a Polars `LazyFrame` (or `dict` when `concatenated=False`). This is **separate from table loading** and does not affect `obsm` on any AnnData.

**`validate_points=False` (default)** — schema validation is skipped at construction for speed (slow on many-batch or remote/zip stores). Pass `validate_points=True` to surface missing required columns eagerly rather than at `load_points` call time.

```python
# Eager schema check at construction (useful when debugging a new store)
loader = Loader("/path/to/run.zarr", validate_points=True)
```

### Filters

All filters compose as AND. Pass any combination:

```python
# Filter by a target type present in this run (check available_target_types first)
points = loader.load_points(target_type="ThreePrimeUntargeted")
# Legacy barcoding SDOs may still use target_type="Transcript"

# Specific batch
points = loader.load_points(batch="B01")

# Batch + well
points = loader.load_points(batch="B01", well="A1")

# Batch + tile (tile-mode runs only; returns [] in well-mode runs)
points = loader.load_points(batch="B01", tile="L1R01C01S1")

# Multiple wells or batches — pass a list
points = loader.load_points(batch=["B01", "B02"], well=["A1", "A2"])

# Explicit point element name (bypasses all other filters)
points = loader.load_points(point_name="A1-points-B01")

# No filters — all points in the run
points = loader.load_points()

points_df = points.collect()   # materialise from LazyFrame
```

### Points schema (Polars columns)

Key columns after `.collect()`:

| Column | Dtype | Description |
|---|---|---|
| `x`, `y` | `Int64` | Pixel coordinates (global) |
| `Xum`, `Yum` | `Float32` | Coordinates in micrometers |
| `target_name` | `String` | Gene / protein name |
| `tile_id` | `String` | Tile identifier |
| `well_id` | `String` | Well identifier |
| `cell` | `UInt16` | Local cell assignment within tile |
| `is_nuclear` | `Boolean` | Whether point is in nucleus |
| `label_Cell` | `Int64` | Cell label (local, per tile) |
| `label_Cell_global` | `Int32` | Global cell label |
| `batch_name` | added by Loader | Batch name |
| `run_name` | added by Loader | Run name (multi-run only) |

`label_*` column names vary by segmentation mask (e.g. `label_Nuclei`, `label_Nuclei_global`).

### Tile-mode vs well-mode runs

Point element names differ between run types:

| Mode | Element name pattern | Example |
|---|---|---|
| Well mode | `{well_id}-points-{batch}` | `A1-points-B01` |
| Tile mode | `{tile_name}-points-{batch}` | `L1R01C01S1-points-B01` |

`Loader` auto-detects mode per run. In tile mode, `load_points(well="A1")` still works — it resolves through the tile→well mapping automatically.

```python
# Check mode for a run
run = next(iter(loader._runs.values()))
print(run._point_mode)          # "well" or "tile"
print(run.point_names[:5])      # list all point element names
```

### Via SpatialData (dask-backed)

```python
points_dask = sdata.points["L1R01C01S1-points-B01"]
print(points_dask.columns.tolist())
points_pd = points_dask.compute()
```

---

## Coordinate systems

Elements are spatially linked via coordinate systems:

- `tile_{tile_id}` — local pixel coordinates
- `well_{well_id}` — well-level coordinates (translated from tiles)
- `global` — run-level coordinates (aggregates all wells)

```python
print(sdata.coordinate_systems)

well_a1_sdata = sdata.filter_by_coordinate_system("well_A1")
tile_sdata    = sdata.filter_by_coordinate_system("tile_L1R01C01S1")
```

**Note:** Use exact names like `"well_A1"` or `"tile_L1R01C01S1"` — not bare `"A1"` or `"L1R01C01S1"`.

---

## Bounding box spatial subsetting

```python
from spatialdata import bounding_box_query

subset = bounding_box_query(
    sdata,
    axes=("y", "x"),
    min_coordinate=[1000, 1000],
    max_coordinate=[2000, 2000],
    target_coordinate_system="tile_L1R01C01S1"
)
```

---

## Linking elements across modalities

Tables annotate labels via `region_key` and `cell_id_local`.

```python
table = sdata.tables["L1R01C01S1-table-B01"]
print(table.obs["region_key"].unique())   # e.g. ['L1R01C01S1-labels-Cell']

cell_id = 123
labels = sdata.labels["L1R01C01S1-labels-Cell"]
cell_mask = labels.to_numpy() == cell_id

cell_row = table[table.obs["cell_id_local"] == cell_id]
# Tile-local pixels on obs; region-pixel centroids in obsm["spatial"]
cell_x = float(cell_row.obs["X_tile"].iloc[0])
cell_y = float(cell_row.obs["Y_tile"].iloc[0])
# Or: xy = table.obsm["spatial"].loc[cell_row.obs_names[0], ["X", "Y"]]
```

**Global vs local cell IDs:** `cell_id_global = cell_id_local + (2^16 * tile_index)`. Use `cell_id_local` when linking to labels within a tile.

---

## Memory management for large stores

```python
# Process one well at a time
for well in ["A1", "A2", "B1", "B2"]:
    well_sdata = sdata.filter_by_coordinate_system(f"well_{well}")
    # ... process ...
    del well_sdata

# Chunked point processing
for chunk in sdata.points["L1R01C01S1-points-B01"].to_delayed():
    chunk_df = chunk.compute()
    # ... process chunk ...
```
