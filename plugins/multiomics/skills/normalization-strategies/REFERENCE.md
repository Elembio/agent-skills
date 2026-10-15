# Normalization Strategies — Reference

## Contents

1. [normalize() parameter table](#1-normalize-parameter-table)
2. [Expected outputs by modality](#2-expected-outputs-by-modality)
3. [Log1p primitive](#3-log1p-primitive)

---

## 1. normalize() parameter table

### `normalize(adata, target_type, **kwargs)`

| Parameter | Type | Default | Description |
|---|---|---|---|
| `adata` | `AnnData` | — | AnnData to normalize |
| `target_type` | `str` | — | `"Transcript"`, `"ThreePrimeUntargeted"`, `"Protein"` — not `"Imaging"` / `"CellPaint"` (`NotImplementedError`; use imaging-pipeline) |
| `batch_colname` | `str` | `"batch_name"` | `.var` column identifying batch |
| `well_colname` | `str` | `"Well"` | `.obs` column identifying well |
| `run_colname` | `str \| None` | `None` | `.obs` column for run; auto-detects `run_name`/`run_id` when `None` |
| `layer` | `str \| None` | `None` | Input layer; `None` uses `.X` |
| `log_base` | `float \| None` | `None` | Log base for ebsd.pp.log1p step; `None` = natural log |
| `raw_layer_name` | `str` | `"raw"` | Name for the pre-normalization raw snapshot |
| `agg_norm_layer_name` | `str \| None` | `None` | If set, writes intermediate agg-norm result to this layer |
| `mor_norm_layer_name` | `str \| None` | `None` | If set, writes intermediate MoR result to this layer (Transcript only) |
| `inplace` | `bool` | `True` | Modify `adata` in place; returns `None`. If `False`, returns a normalized copy. |
| `target_counts` | `float \| None` | `None` | Fixed scaling target for agg-norm step; `None` uses per-well per-batch median |
| `apply_agg_norm` | `bool \| None` | `None` | **Protein only**: run aggregated-count scaling before ebsd.pp.log1p; `None` → `True` |

---

## 2. Expected outputs by modality

### Transcript / ThreePrimeUntargeted

| Output | Description |
|---|---|
| `adata.X` | Log1p-normalized, bias-corrected values |
| `adata.layers["raw"]` | Integer counts snapshot before any transformation |
| `adata.layers["norm_agg_counts"]` | After agg-count scaling, before MoR (only if `agg_norm_layer_name` set) |
| `adata.layers["norm_med_of_ratios"]` | After MoR, before ebsd.pp.log1p (only if `mor_norm_layer_name` set) |

Typical `.X` range after Transcript normalization: **0 – 10** (natural log scale). Values above 15 may indicate normalization did not run or data contains outlier genes.

### Protein

| Output | Description |
|---|---|
| `adata.X` | ebsd.pp.log1p-normalized values |
| `adata.layers["raw"]` | Integer counts snapshot |

### Imaging (alias CellPaint)

`ebsd.pp.normalize(..., "Imaging")` raises `NotImplementedError`. See [imaging-pipeline](../imaging-pipeline/SKILL.md): base + CellProfiler + `AreaShape_*` → MAD-z against a **control-well, global, or per-well** reference depending on whether the biology is constant within a well; **exclude Bkg / IgG / Cleave** from the embedding matrix; keep raw `X` for S/B and localization.

---

## 3. Log1p primitive

For ad-hoc ebsd.pp.log1p on a specific layer without running the full pipeline:

```python
import elembio_spatialdata_analysis as ebsd

ebsd.pp.log1p(adata, layer=None, base=None)   # base=None → natural log
```
