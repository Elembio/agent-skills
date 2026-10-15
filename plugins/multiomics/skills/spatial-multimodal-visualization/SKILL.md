---
name: spatial-multimodal-visualization
description: >-
  Plot AVITI24 well fluorescence: morphology stains and protein IF targets
  on the same well image, the whole flow-cell plate, or a single-cell
  close-up. Use when the user asks to show the well image, immunofluorescence,
  IF / protein targets as channels, the plate of well images, overlay
  segmentation, fill cells by a table column, or show a representative
  cell crop. Not for UMAP (embedding-visualization) or squidpy neighborhood
  stats (spatial-statistics).
metadata:
  version: 0.3.3
  author: elembio
---

# Spatial Multimodal Visualization

`plot_image` for wells and the plate. `plot_cell` for one cell close up.

```text
sd.read_zarr (directory store)
  → ebsd.vis.list_image_channels(sdata, well=…)
  → ebsd.vis.plot_image(sdata, well=…)             # one well
  → ebsd.vis.plot_image(sdata)                     # plate
  → ebsd.vis.plot_image(..., show_image=False, outlines=True)
  → ebsd.vis.plot_cell(sdata, cell_id)             # one cell, full-res crop
```

## Hard rules

- **Rasters need `SpatialData`.** `Loader` loads tables/points only. Open images
  with `spatialdata.read_zarr` on a **directory** `.zarr`. A `.zarr.zip` must be
  FUSE-mounted first — [spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md).
- **Well-keyed names.** Images `{well}-images-cell_paint`. Labels
  `{well}-labels-Cell`. Shapes `{well}-shapes-Cell`. Print `list(sdata.images)`
  before plotting.
- **Use `ebsd.vis.plot_image` and `ebsd.vis.plot_cell`.** Do not hand-roll
  `imshow` recipes. If they are missing, the analysis package is too old.
- **IF protein targets are extra channels on the same image.** Default overlay
  is Nucleus / Actin / Cell-Membrane. Protein IF (and other stains present on
  the run) are selected with `channels=` after
  `list_image_channels(sdata, well)`. Do not invent target names. Skip
  `Bkg_*` / IgG unless the user asked for those diagnostics.
  `color=` is a table fill (per-cell numbers), not the IF raster.
- **Scale 0 is a whole-well image, not a cell crop.** `plot_image` defaults to
  coarsest (`scale=None`). Never pass `scale=0` for the plate. `plot_cell`
  defaults to `scale=0` and only reads a 256 px full-res window — that is cheap.
  `crop_size` is the field of view in full-resolution pixels; changing `scale`
  does not zoom.
- **Use `show_image=False` for segmentation alone.** Fills and/or outlines on
  black; fluorescence pixels are not loaded. Requires `color=` and/or
  `outlines=True`. Do not pass `channels` / `channel_colors` / `separate`.
- **Label pixels and shape ids match `obs["Cell"]`**, not `cell_id_local`.
- **`plot_cell` takes a cell id, not a rule.** Pick the id from a table on
  `sdata` first. Do not pass a separate Loader AnnData.
- **Tables for ids and fills live on `sdata`.** Prefer
  `sdata.tables["global-table-imaging"]` (discover keys with
  `list(sdata.tables)`). `obsm["spatial"]` is the cell-centroid table for
  [spatial-statistics](../spatial-statistics/SKILL.md), not for these image
  plotters. Use `Loader` only when you need tables without opening rasters.

## Prerequisites

- [spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md)
- A directory `.zarr` (or FUSE-mounted zip) with images / labels
- For `plot_image(..., color=)`: a cell table — same
  `sdata.tables["global-table-imaging"]` (or `Loader` if rasters are not open)

## Quick start

```python
import numpy as np
import spatialdata as sd
import elembio_spatialdata_analysis as ebsd

sdata = sd.read_zarr("/path/to/run.zarr")  # directory; zip → ratarmount first
print(list(sdata.tables))
print(sorted(sdata.images))
print(ebsd.vis.list_image_channels(sdata, well="A1"))
img = sdata.tables["global-table-imaging"]

ebsd.vis.plot_image(sdata, well="A1")
ebsd.vis.plot_image(sdata, well="A1", channels=["Nucleus", "Actin"], outlines=True)
# IF: names come from list_image_channels — do not invent them
# ebsd.vis.plot_image(sdata, well="A1", channels=["Nucleus", "YOUR_IF_TARGET"])
# ebsd.vis.plot_cell(sdata, cell_id, channels="YOUR_IF_TARGET", channel_colors="magenta")
ebsd.vis.plot_image(sdata)  # plate: letter columns × number rows
ebsd.vis.plot_image(sdata, show_image=False, outlines=True)

area = img.obsm["cell_geometry"]["Area"].to_numpy()
pick = int(np.argmin(np.abs(area - np.percentile(area, 90))))
ebsd.vis.plot_cell(sdata, int(img.obs["Cell"].iloc[pick]))
```

## Detail files

- [IMPLEMENTATION.md](./IMPLEMENTATION.md) — plate, channels, outlines, fills, `plot_cell`
- [REFERENCE.md](./REFERENCE.md) — element names, ID mapping, argument tables

## Related skills

- [spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md) — `read_zarr` vs Loader, zip / FUSE
- [cell-quality-control](../cell-quality-control/SKILL.md) — filter before coloring by QC
- [imaging-pipeline](../imaging-pipeline/SKILL.md) — morphology / IF analysis (not image display)
- [embedding-visualization](../embedding-visualization/SKILL.md) — UMAP / PCA, not rasters
- [spatial-statistics](../spatial-statistics/SKILL.md) — squidpy on coordinates, not images

## Disclaimer

This content is provided "as-is" without warranty of any kind. AI-generated output may contain
inaccuracies. Review and verify all outputs before relying on them.
