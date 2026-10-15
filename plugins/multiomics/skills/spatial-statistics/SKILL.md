---
name: spatial-statistics
description: >-
  Spatial neighborhood analysis with squidpy on AVITI24 cell tables:
  neighbors, Moran's I / Geary's C, co-occurrence, enrichment, Ripley's L.
  Use when the user asks about spatial patterns, proximity, or cell-cell
  neighborhoods. Coordinates come from obsm["spatial"] on the cell table.
  Not for well fluorescence images (spatial-multimodal-visualization).
metadata:
  version: 0.3.1
  author: elembio
---

# Spatial Statistics

Squidpy on **cell centroids**, not rasters. Build a neighborhood graph
**per well**, then run autocorrelation / enrichment / Ripley on that graph.

```text
Loader or sdata.tables
  → xy from obsm["spatial"] as an (n, 2) array
  → subset Well
  → sq.gr.spatial_neighbors
  → Moran's I / nhood_enrichment / Ripley
```

## Hard rules

- **Use the store's spatial table.** Current tables carry `obsm["spatial"]`
  (`X`, `Y` in the region pixel frame). Convert to a numpy array before
  squidpy. Do not invent `obs["X"]` / `obs["Y"]` or `obs["Xum"]` / `obs["Yum"]`.
  Tile-local pixels are `obs["X_tile"]` / `obs["Y_tile"]`.
- **Do not mix wells in one graph.** Wells are separate physical spaces.
  Subset `obs["Well"]` first; rebuild neighbors per well.
- **Tables only — no `read_zarr`.** `Loader.load_tables` or
  `sdata.tables["global-table-…"]` is enough. Rasters are
  [spatial-multimodal-visualization](../spatial-multimodal-visualization/SKILL.md).
- **Discover the table.** `print(loader.available_tables)` or
  `print(list(sdata.tables))`. 3′ uses `ThreePrimeUntargeted`; imaging uses
  `imaging`.
- **Units match the coordinates.** `obsm["spatial"]` is pixels. A squidpy
  `radius` is then pixels, not microns. Tile-local microns live in
  `obsm["cell_geometry"]` as `Xum_tile` / `Yum_tile` when present.

## Prerequisites

- [spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md)
- [clustering-and-annotation](../clustering-and-annotation/SKILL.md) before neighborhood enrichment
- An AnnData with `obs["Well"]` and `obsm["spatial"]` (or `cell_geometry`)

## Quick start

```python
import numpy as np
import squidpy as sq
from elembio_spatialdata_analysis.load import Loader

loader = Loader("/path/to/run.zarr")  # zip is fine — tables only
print(loader.available_tables)
adata = loader.load_tables(table="ThreePrimeUntargeted", lazy=False)

sp = adata.obsm["spatial"]
adata.obsm["spatial"] = (
    sp[["X", "Y"]].to_numpy(dtype=float) if hasattr(sp, "columns") else np.asarray(sp, dtype=float)
)

well = adata[adata.obs["Well"].astype(str) == "A1"].copy()
sq.gr.spatial_neighbors(well, coord_type="generic", n_neighs=10, spatial_key="spatial")
sq.gr.spatial_autocorr(well, mode="moran", genes=list(well.var_names[:100]), n_perms=100)
print(well.uns["moranI"].sort_values("I", ascending=False).head(10))
```

If `read_zarr` is already open, the same table is
`sdata.tables["global-table-ThreePrimeUntargeted"]` (or `global-table-imaging`).

## Detail files

- [IMPLEMENTATION.md](./IMPLEMENTATION.md) — neighbors, Moran's / Geary's, enrichment, Ripley, per-well loop
- [REFERENCE.md](./REFERENCE.md) — parameters, output slots, pitfalls

## Related skills

- [clustering-and-annotation](../clustering-and-annotation/SKILL.md)
- [spatial-multimodal-visualization](../spatial-multimodal-visualization/SKILL.md) — well images, not squidpy
- [spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md)

## Disclaimer

This content is provided "as-is" without warranty of any kind. AI-generated output may contain
inaccuracies. Review and verify all outputs before relying on them.
