---
name: anndata-structure-reference
description: Canonical AnnData slot layout for AVITI24 tables — where counts, imaging features, QC metrics, geometry, OPS assignments, and filter flags live. Use when asking where a field is stored, why a column is missing, or which slots 3′, OPS, imaging, or legacy Transcript/Protein use. Discover present modalities first; NSB/GAPDH apply to Transcript only.
metadata:
  version: 0.2.6
  author: elembio
---

# AnnData Structure Reference

## Purpose

Canonical AnnData slot layout for AVITI24 `global-table-*` objects as consumed by `elembio_spatialdata_analysis`. Look up `obs` / `var` / `obsm` / `layers` / `uns` by **modality actually present**. Legacy layouts (`control_variables_{batch}`, `obs['total_counts_*']`) are obsolete.

```text
Current tables: ThreePrimeUntargeted → OPS (custom names) → imaging (CellPaint alias)
Legacy: Transcript / Protein
```

## Prerequisites

- [spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md)
- [cell-quality-control](../cell-quality-control/SKILL.md) for QC that writes `qc_metrics` / `filter_flags`

## Key Concepts

1. **`X` depends on modality** — 3′ / Transcript / Protein / OPS: raw integer counts (sparse). **Imaging:** dense continuous features (not counts). After `pp.normalize` on count modalities, pre-norm snapshot → `layers["raw"]`, `X` becomes normalized.
2. **Cell area** lives in `obsm["cell_geometry"]` (`Area`, …) — not copied into `obs` by default. `filter_cells` reads geometry first, then `obs["Area"]`.
3. **Discover `obsm` keys** — do not assume `total_counts`, `unassigned_counts`, `control_variables`, or `spatial_artifacts` exist on every table.
4. **NSB / GAPDH are Transcript (legacy barcoding) concepts** — NSB in `obsm["control_variables"]` (`NSB*`); GAPDH as rows in `X` with `var["control_type"] == "Positive Control 1"`. Do **not** narrate for 3′, imaging, or OPS.
5. **3′ count depth for QC** uses summed `.X` (batches already summed into the matrix). Per-batch `obsm["total_counts"]` is a Transcript/Protein pattern when present.
6. **`tile_id`** preferred over `Tile`; safe fallback: `tile_col = "tile_id" if "tile_id" in adata.obs.columns else "Tile"`.
7. **`Well` vs `WellLabel`:** `obs["Well"]` is the plate position (`A1`). `obs["WellLabel"]` is the sample-sheet string. Color/groupby on cells uses `WellLabel` (or `condition`). `run.wells` stores the same strings as `well_label` — there is no `obs["well_label"]`.
8. **`obsm["qc_metrics"]`** exists only after `compute_qc_metrics`. **`obsm["filter_flags"]`** after `filter_cells` (e.g. `passes_counts`, `passes_area`, `passes_all_filters`; imaging → `passes_imaging_intensity`; 3′ → `passes_mt_outlier`).
9. **OPS** tables are DISS (`SpecializedTargeted`) with custom names — discover via `run.table_suffixes_for_target_type(("OPS", "SpecializedTargeted"))`, then drop QC/gene panels before Confidence workflows; load assay tables with `include_target_assignment=True`.
10. **Centroids for squidpy** live in `obsm["spatial"]` (`X`, `Y` pixels, often a DataFrame). Convert to an ndarray. Tile-local pixels are `obs["X_tile"]` / `obs["Y_tile"]`. Image plots use rasters + `obs["Cell"]`, not the spatial slot.

## Quick-Start Access Patterns

```python
import numpy as np
import scipy.sparse as sp

# Discover what this AnnData actually has
print("obsm:", list(adata.obsm.keys()))
print("var cols:", adata.var.columns.tolist())

# .X — sparse counts or dense imaging features
if sp.issparse(adata.X):
    total = np.asarray(adata.X.sum(axis=1)).ravel()
else:
    total = np.asarray(adata.X).sum(axis=1)

# Cell area
if "cell_geometry" in adata.obsm and "Area" in adata.obsm["cell_geometry"].columns:
    area = adata.obsm["cell_geometry"]["Area"]
elif "Area" in adata.obs.columns:
    area = adata.obs["Area"]
else:
    area = None

tile_col = "tile_id" if "tile_id" in adata.obs.columns else "Tile"

# Transcript-only when present — skip on 3′ / imaging / OPS
if "control_variables" in adata.obsm:
    nsb_cols = [c for c in adata.obsm["control_variables"].columns if "NSB" in c.upper()]

# After QC
if "qc_metrics" in adata.obsm:
    print(adata.obsm["qc_metrics"].columns.tolist())
if "filter_flags" in adata.obsm:
    passes = adata.obsm["filter_flags"]["passes_all_filters"]
```

## Detail Files

- [REFERENCE.md](./REFERENCE.md) — full slot inventory, per-modality sections, validation

## Related Skills

- [spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md)
- [run-quality-overview](../run-quality-overview/SKILL.md)
- [cell-quality-control](../cell-quality-control/SKILL.md)
- [imaging-pipeline](../imaging-pipeline/SKILL.md) — imaging / IF analysis orchestration (geometry in `obsm["cell_geometry"]`)

## Disclaimer

This content is provided "as-is" without warranty of any kind. AI-generated output may contain
inaccuracies. Review and verify all outputs before relying on them.
