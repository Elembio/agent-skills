# Raw Data Ingestion — Reference

## Column Naming Convention

```
<TargetName>.<BatchID>          → whole-cell measurement
<TargetName>_Nuclear.<BatchID>  → nuclear compartment measurement
```

### Batch ID suffixes by modality

| Modality | Suffix pattern | Example |
|---|---|---|
| Transcript (legacy barcoding) | `.B01` – `.B06` | `GAPDH.B01` |
| Protein (legacy barcoding) | `.P01` – `.P06` | `AKT1.P02` |
| Imaging channel (RawCellStats `.CPxx`) | `.CP01`, `.CP02`, … | `Phalloidin.CP01` |
| Imaging CellProfiler | `AreaShape_*`, `Granularity_*_<stain>.CP01`, … | `Intensity_MeanIntensity_WGA.CP01` |

### Unassigned control marker

Any target whose name is exactly `"Unassigned"` (case-insensitive) is treated as an unassigned / background control. Set `var["control_type"] = "Unassigned"` for those features.

---

## Expected Shapes (legacy Transcript + Protein example)

| AnnData | `.n_obs` | `.n_vars` | Notes |
|---|---|---|---|
| `transcript_adata` | 325 815 | 748 | 748 transcripts × 1 batch cycle |
| `protein_adata` | 325 815 | 56 | 56 proteins × 1 batch cycle, whole-cell |
| `protein_nuc_adata` | 325 815 | 56 | same proteins, nuclear compartment |
| `imaging_adata` | 325 815 | 168 | 6 raw channels + 162 per-channel morphology features |

---

## obs Column Reference

RawCellStats parquet still uses `Tile` / `X` / `Y` / `Area`. Map them to Loader slots when building AnnData:

| Parquet column | Loader equivalent | Description |
|---|---|---|
| `Cell` | `obs_names` / `obs["Cell"]` | Unique integer cell ID |
| `Well` | `obs["Well"]` | Plate position (`A1`) |
| `WellLabel` | `obs["WellLabel"]` | Sample-sheet name. `run.wells` stores the same strings as `well_label`. |
| `Tile` | `obs["tile_id"]` | Tile / FOV identifier |
| `X` | `obs["X_tile"]` | Tile-local centroid X (pixels) |
| `Y` | `obs["Y_tile"]` | Tile-local centroid Y (pixels) |
| `Area` | `obsm["cell_geometry"]["Area"]` | Cell area (pixels²) |

Region-pixel centroids for squidpy live in `obsm["spatial"]` (`X`, `Y`) on zarr builds. RawCellStats usually only has tile-local `X`/`Y` — rename to `X_tile`/`Y_tile`; do not put them in `cell_geometry`.

---

## Imaging columns in RawCellStats (`.CPxx` naming)

RawCellStats parquet still uses `.CPxx` suffixes for imaging channels. In the SpatialData product these map to the **imaging** table (`imaging` / `CellPaint` alias).

| Type | Column pattern | Notes |
|---|---|---|
| Base intensity | `<stain>.CPxx` | Whole-cell (and related) intensities |
| Per-channel CellProfiler | `<Prefix>_<stain>.CPxx` | Texture / Granularity / … |
| Whole-cell geometry | `AreaShape_*` | Often cell-level; zarr keeps geometry in `obsm["cell_geometry"]` |

> Do **not** count-normalize imaging. Use [imaging-pipeline](../imaging-pipeline/SKILL.md) MAD-z, with the reference (control wells / global / per-well) chosen from the design.

---

## var Column Reference

| Column | Values | Description |
|---|---|---|
| `target_name` | e.g. `"GAPDH"`, stain name | Target / gene / stain (batch suffix stripped) |
| `batch_name` | e.g. `"B01"`, `"P02"`, `"CP01"` | Cycle or imaging batch |
| `target_type` | `"Imaging"`, `"Transcript"`, `"Protein"`, … | Prefer `"Imaging"` for imaging columns; legacy parquet demos may still say `"CellPaint"` |
| `control_type` | `"Unassigned"` or `""` | Panel controls when present |
