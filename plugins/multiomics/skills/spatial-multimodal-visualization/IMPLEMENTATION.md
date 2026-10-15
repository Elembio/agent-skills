# Spatial Multimodal Visualization — Implementation

## Table of Contents

1. [Open rasters](#1-open-rasters)
2. [Plot the well image](#2-plot-the-well-image)
3. [Channels and colors](#3-channels-and-colors)
4. [Segmentation outlines](#4-segmentation-outlines)
5. [Fill cells from a table](#5-fill-cells-from-a-table)
6. [Plot one cell](#6-plot-one-cell)
7. [Common pitfalls](#7-common-pitfalls)

---

## 1. Open rasters

Directory store:

```python
import spatialdata as sd
import elembio_spatialdata_analysis as ebsd

sdata = sd.read_zarr("/path/to/run.zarr")
```

`.zarr.zip` — mount first, then `read_zarr` the **nested** `{run}.zarr` (not the mount root). Full FUSE recipe: [spatialdata-loading-and-access](../spatialdata-loading-and-access/SPATIAL.md#reading-a-zarrzip-with-read_zarr-fuse-mount).

```python
# after: ratarmount /path/to/run.zarr.zip /mnt/store
sdata = sd.read_zarr("/mnt/store/run.zarr")
```

`read_zarr` is eager on open (can take a minute on a large store). Array data stays lazy.

Tables-only work can still use `Loader("/path/to/run.zarr.zip")` without mounting.

Cell tables for ids and `color=` are on the same store:

```python
print(list(sdata.tables))
adata = sdata.tables["global-table-imaging"]
```

`obsm["spatial"]` is not an input to `plot_image` / `plot_cell`. Use it in
[spatial-statistics](../spatial-statistics/SKILL.md).

---

## 2. Plot the well image

```python
print(sorted(sdata.images))
print(ebsd.vis.list_image_channels(sdata, well="A1"))
ebsd.vis.plot_image(sdata, well="A1")
ebsd.vis.plot_image(sdata)  # plate: every well that has an image
# ebsd.vis.plot_image(sdata, wells=["A1", "B1"])
```

Default channels are Nucleus (red) / Actin (green) / Cell-Membrane (blue),
overlaid as one image. The plate uses the same tints on every well.

Layout is letter columns × number rows (12-well is two rows: A1–F1 then
A2–F2). Wells without an image stay empty so the grid stays rectangular.
`separate=True` and a list of `color=` columns are single-well only.

---

## 3. Channels and colors

Morphology stains and protein IF targets live on the **same**
`{well}-images-cell_paint` element as extra `c` channels. The default overlay
is Nucleus (red) / Actin (green) / Cell-Membrane (blue). Anything else that
`list_image_channels` returns — including IF protein targets — is valid
`channels=` input for `plot_image` and `plot_cell`.

```python
names = ebsd.vis.list_image_channels(sdata, well="A1")
print(names)  # morphology + IF targets present on this run; do not invent names

ebsd.vis.plot_image(sdata, well="A1", channels="Nucleus")
ebsd.vis.plot_image(sdata, well="A1", channels="Nucleus", channel_colors="cyan")
ebsd.vis.plot_image(
    sdata,
    well="A1",
    channels=["Nucleus", "Actin"],
    channel_colors=["blue", "green"],
)
# Protein IF: pick names from `names` (not Bkg_* / IgG unless asked)
# ebsd.vis.plot_image(sdata, well="A1", channels=["Nucleus", names[i]])
# ebsd.vis.plot_image(sdata, channels=names[i])  # plate, same IF on every well
# ebsd.vis.plot_cell(sdata, cell_id, channels=names[i], channel_colors="magenta")
ebsd.vis.plot_image(
    sdata,
    well="A1",
    channels=["Nucleus", "Actin"],
    separate=True,  # one panel per channel
)
```

Do not invent channels. Tint IF with `channel_colors` (one color, a list, or a
`{channel: color}` dict). Extra channels without a named default cycle a
fluorescence palette.

`scale=0` is full resolution (~15k × 16k × all channels). Do not use it for a
whole-well preview. `plot_image` rejects `scale=0` when more than one well
is drawn.

Background / IgG channels (`Bkg_*`, `Mouse_IgG`, `Rabbit_IgG`) are diagnostics —
do not use them as the default overlay.

Raster IF (`channels=`) is the pixel stain. A table fill (`color=`) is a
per-cell number (for example intensity in `obs` / `qc_metrics`). Do not
substitute one for the other.

---

## 4. Segmentation outlines

```python
ebsd.vis.plot_image(sdata, well="A1", outlines=True)
ebsd.vis.plot_image(sdata, well="A1", channels="Nucleus", outlines=True)
# ebsd.vis.plot_image(sdata, well="A1", outlines=True, outline_width=0.5)

# Labels only: do not load or draw fluorescence
ebsd.vis.plot_image(
    sdata,
    well="A1",
    show_image=False,
    outlines=True,
)
```

This reads `{well}-labels-Cell` (mask borders, not the stored polygons).
No AnnData table is required. `plot_cell` is the helper that outlines one
cell from `{well}-shapes-Cell`.
`show_image=False` uses a black background so white outlines remain visible.
It requires `outlines=True` and/or `color=`.

---

## 5. Fill cells from a table

Optional. Paints each mask by a cell-level value (`obs`, `var_names`, or
`qc_metrics`) — not by raw channel pixels.

```python
adata = sdata.tables["global-table-imaging"]
# Area lives in obsm["cell_geometry"] — copy if you want it as color=
adata.obs["Area"] = adata.obsm["cell_geometry"]["Area"]

ebsd.vis.plot_image(sdata, well="A1", adata=adata, color="Area", outlines=True)

# Fill and outlines only, with no fluorescence pixels loaded
ebsd.vis.plot_image(
    sdata,
    adata=adata,
    color="Area",
    show_image=False,
    outlines=True,
    alpha=1,
)
```

Discover fill names from `adata.obs`, `adata.var_names`, and `qc_metrics`.
The helper subsets to `Well == well` and maps through `obs["Cell"]`. On the
plate, each panel uses that well's rows and a continuous `color=` shares one
scale across wells. A column that is constant within a well (for example
`WellLabel`) is a poor plate fill — every cell in a panel gets the same color.

With `show_image=True` (default), `alpha=0.55` lets fluorescence remain
visible below a fill. With `show_image=False`, use `alpha=1` for solid fill.
Image-only options (`channels`, `channel_colors`, `pct`, `separate`) are
invalid when `show_image=False`; the helper raises instead of ignoring them.

---

## 6. Plot one cell

`plot_cell` crops the well image around one cell id and draws that cell's
`{well}-shapes-Cell` polygon as a vector contour. Neighbors stay visible as
context but are not outlined.

```python
ebsd.vis.plot_cell(sdata, cell_id)
ebsd.vis.plot_cell(sdata, cell_id, channels="Nucleus", channel_colors="cyan")
# ebsd.vis.plot_cell(sdata, cell_id, channels=if_name)  # if_name from list_image_channels
ebsd.vis.plot_cell(sdata, cell_id, scale=1)  # same FOV, coarser pixels
ebsd.vis.plot_cell(sdata, cell_id, crop_size=128)  # tighter field, still full res
ebsd.vis.plot_cell(sdata, cell_id, smooth=0)  # exact traced boundary, no rounding
# ebsd.vis.plot_cell(sdata, cell_id, well="A1")  # skip the table lookup
```

The well comes from a table on the store (`obs["Cell"]` and `obs["Well"]`),
so no separate `Loader` table is needed. Default `scale=0` is full
resolution; `crop_size` is the field of view in full-resolution pixels
(default 256), so changing `scale` only changes sharpness. Crops at the
well edge are padded so the cell stays centered. The contour is smoothed
(`smooth=6`) by default.

**"A representative high-area cell" is a selection step, not an argument.**
Pick the id first, then plot it:

```python
import numpy as np

img = sdata.tables["global-table-imaging"]
area = img.obsm["cell_geometry"]["Area"].to_numpy()
pick = int(np.argmin(np.abs(area - np.percentile(area, 90))))  # nearest the 90th pct
ebsd.vis.plot_cell(sdata, int(img.obs["Cell"].iloc[pick]))
```

Use a quantile-nearest cell rather than the maximum — the maximum is usually a
segmentation artifact, not a representative cell.

---

## 7. Common pitfalls

- **`Loader(path)` has no images.** Mount + `read_zarr`, or pass a directory store.
- **`read_zarr` on the zip file fails.** Mount with `ratarmount`, then open the nested `.zarr`.
- **Painting with `cell_id_local`.** Well-keyed masks match `obs["Cell"]`.
- **Scale 0 OOM / hang.** Leave `scale=None` unless the user asked for full resolution.
  Never use `scale=0` for the plate. `plot_cell` defaults to `scale=0` but
  reads only its crop window, so that is cheap.
- **Asking `plot_cell` to choose a cell.** It takes an id. Select from the table first.
