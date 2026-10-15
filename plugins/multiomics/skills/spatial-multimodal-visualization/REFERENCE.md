# Spatial Multimodal Visualization — Reference

## Element names (current / well-keyed)

| Element | Pattern | Example |
|---|---|---|
| Image | `{well}-images-cell_paint` | `A1-images-cell_paint` |
| Cell labels | `{well}-labels-Cell` | `A1-labels-Cell` |
| Nuclear labels | `{well}-labels-Nuclear` | `A1-labels-Nuclear` |
| Shapes | `{well}-shapes-Cell` | `A1-shapes-Cell` |

## Label ID mapping

On well-keyed stores, `{well}-labels-Cell` pixel values match `adata.obs["Cell"]` (same ids as `cell_id_global`). They do **not** match `cell_id_local`. `{well}-shapes-Cell` uses the same ids as its index.

`plot_image(..., adata=..., color=...)` maps through `obs["Cell"]`. `plot_cell(sdata, cell_id)` looks the same id up in the store's own table.

## Pyramid

Levels `0`–`5` on current well images. Level `0` is full resolution (`c × ~14688 × ~16488`). Level `5` is the `plot_image` default (`~459 × 515`). `plot_cell` defaults to level `0` but reads only `crop_size` (256 full-res pixels).

Image and label arrays at the same level are the same `y × x` (or off by one pixel — the helper crops to the overlap).

## Default channel tints

| Channel | Color |
|---|---|
| `Nucleus` | red |
| `Actin` | green |
| `Cell-Membrane` | blue |

Further channels — including protein IF targets — cycle a fluorescence
palette. Override with `channel_colors`. List names with
`list_image_channels`. Do not invent targets. Do not use IgG / `Bkg_*` as
the default overlay.

## `plot_image`

| Argument | Default | Notes |
|---|---|---|
| `well` | `None` (plate) | One well id, e.g. `"A1"`. Omit for every well that has an image |
| `wells` | all with images | Optional subset for the plate. Do not combine with `well=` |
| — | `list_image_channels(sdata, well)` | Channel names on that well image |
| `channels` | Nucleus / Actin / Cell-Membrane | Morphology **or protein IF** names from `list_image_channels`; tinted overlay. Shared on the plate |
| `channel_colors` | see table above | One color, a list, or `{channel: color}` |
| `outlines` | `False` | White cell borders from `{well}-labels-Cell` (mask strokes, not shapes) |
| `outline_width` | `1` | Line weight (points). ``0.5`` is thinner; ``0`` draws none |
| `adata`, `color` | off | Fill masks by `obs` / `var_names` / `qc_metrics`. Plate: one shared color scale |
| `separate` | `False` | One panel per channel (single well only) |
| `show_image` | `True` | Set `False` for fill/outlines on black without loading fluorescence |
| `scale` | `None` (coarsest) | `0` = full resolution; rejected when more than one well is drawn |
| `alpha` | `0.55` | Opacity of a table fill |
| `cmap` | `"viridis"` | Continuous table fill |

With `show_image=False`, provide `color=` and/or `outlines=True`. Do not pass
image-only options: `channels`, `channel_colors`, `pct`, or `separate`.

## `plot_cell`

| Argument | Default | Notes |
|---|---|---|
| `cell_id` | required | Value from `obs["Cell"]`; unique across the store |
| `well` | resolved from a store table | Pass it to skip the lookup |
| `scale` | `0` | Full resolution. Only the crop window is read. Coarser levels stay the same FOV |
| `crop_size` | `256` | Square field in **full-resolution** pixels; edge crops are padded |
| `outline` | `True` | Vector contour of `{well}-shapes-Cell` for this cell only |
| `smooth` | `6` | Simplification in full-res pixels, then corner-cut. `0` is the exact traced boundary |
| `outline_color`, `outline_width` | `"white"`, `1.5` | Contour style |
| `channels`, `channel_colors` | as `plot_image` | Morphology or protein IF; same tint rules |

Stored shapes are pixel-boundary polygons at full resolution. The default
`smooth=6` rounds that staircase for figures.

## Common pitfalls

**Zip without mount** — `spatialdata.read_zarr` cannot open `.zarr.zip`. Use `ratarmount`, then `read_zarr` on the nested store.

**`read_zarr` open time** — metadata walk is eager; a large well store can take a minute. Pixel data stays lazy until a scale is read.

**Y-axis** — `imshow` already uses image coordinates (row 0 at the top). Do not invert Y on these helpers.

**Memory** — never `to_numpy()` scale 0 of a whole well for a preview.
