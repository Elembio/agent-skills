---
name: legacy-targeted-pipeline
description: End-to-end orchestration for legacy AVITI24 Transcript / Protein barcoding panels — load → QC → normalize → PCA/UMAP → cluster → pseudobulk DE. Use only when available_tables include Transcript and/or Protein barcoding. Not for current 3′ / OPS / imaging-first runs (use 3prime-untargeted-pipeline, ops-pipeline, or imaging-pipeline). Not for CellPaint/imaging (use imaging-pipeline).
metadata:
  version: 0.2.5
  author: elembio
---

# Legacy Targeted Pipeline

> **Confirm the SDO first.** Before running, make sure you are using the SDO the user selected in [locate-spatialdata-store](../locate-spatialdata-store/SKILL.md) — a run may have several (instrument-direct, `cells2stats`, `elembio-multiomics-nf`).

## Purpose

Orchestrator for **legacy** Transcript / Protein barcoding panels. Links to specialist skills for each step. Current-product runs (3′, OPS, imaging) should **not** use this skill.

## Prerequisites

- AVITI24 spatialdata Zarr store with `Transcript` and/or `Protein` tables
- Python **3.12+** with `elembio_spatialdata_analysis`

### Check data type before starting

```python
from elembio_spatialdata_analysis.load import Loader
loader = Loader("/path/to/run.zarr")
print(sorted(loader.available_tables))
```

| Present | Go to |
|---|---|
| `ThreePrimeUntargeted` | [3prime-untargeted-pipeline](../3prime-untargeted-pipeline/SKILL.md) |
| OPS / `SpecializedTargeted` custom table | [ops-pipeline](../ops-pipeline/SKILL.md) |
| `imaging` / `CellPaint` | [imaging-pipeline](../imaging-pipeline/SKILL.md) |
| `Transcript` and/or `Protein` (barcoding) | **continue here** |

If unsure, ask which modality to analyze. Never invent absent modalities.

## Key Concepts

### Loader vs `read_zarr`

| Need | Use |
|------|-----|
| Tables, metadata, analysis | `Loader(STORE_PATH)` |
| Images, labels, spatial overlays | `spatialdata.read_zarr`, then optionally `Loader(sdata)` |

**Materialization:** Always `lazy=False` for this pipeline. Lazy `X` conflicts with PCA / neighbors / most `pp` steps.

### Pipeline Step Decision Table

| Step | Modality / Decision | Guidance |
|------|---------------------|----------|
| **Load** | Transcript and/or Protein | Check `loader.available_tables` |
| **QC filter** | Default or custom thresholds | Start with defaults; adjust if >30 % cells rejected |
| **Normalization** | Transcript / Protein | `pp.normalize` — see [normalization-strategies](../normalization-strategies/SKILL.md) |
| **HVG** | Transcript (targeted) | Skip — panels are small |
| **PCA / UMAP** | After normalize | [dimensionality-reduction](../dimensionality-reduction/SKILL.md) |
| **Batch correction** | Wells dominate PCA | [batch-effect-assessment](../batch-effect-assessment/SKILL.md) |
| **Clustering** | Leiden resolution | Start 0.5 |
| **DE** | Pseudobulk DESeq2 | ≥ 2 biological replicates per condition |

### Pipeline Order

```
Step 0  Setup
Step 1  Loader(STORE_PATH)          ← spatialdata-loading-and-access
Step 2  RunInfo metadata            ← metadata-and-target-mapping
Step 3  load_tables(table=..., lazy=False)  ← spatialdata-loading-and-access
Step 4  Experimental metadata join  ← metadata-and-target-mapping
Step 5  compute_qc_metrics +
        filter_cells                ← cell-quality-control
Step 6  normalize                   ← normalization-strategies
Step 7  PCA + batch assessment      ← batch-effect-assessment
Step 8  neighbors + UMAP            ← dimensionality-reduction
Step 9  UMAP visualization          ← embedding-visualization
Step 10 Leiden + rank_genes_groups  ← clustering-and-annotation
Step 11 Pseudobulk DE (DESeq2)      ← pseudobulk-differential-expression
Step 12 Save results

Optional: [cluster-validation](../cluster-validation/SKILL.md) — invoke deliberately, not by default.
```

## Step Index

| Step | Name | Specialist skill |
|------|------|-----------------|
| 1 | Open store with `Loader` | [spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md) |
| 2 | Parse `RunInfo` metadata | [metadata-and-target-mapping](../metadata-and-target-mapping/SKILL.md) |
| 3 | Load Transcript / Protein tables | [spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md) |
| 4 | Add experimental metadata | [metadata-and-target-mapping](../metadata-and-target-mapping/SKILL.md) |
| 5 | QC filtering | [cell-quality-control](../cell-quality-control/SKILL.md) |
| 6 | Normalize | [normalization-strategies](../normalization-strategies/SKILL.md) |
| 7 | PCA + batch effect assessment | [batch-effect-assessment](../batch-effect-assessment/SKILL.md) |
| 8 | Neighbors + UMAP | [dimensionality-reduction](../dimensionality-reduction/SKILL.md) |
| 9 | Embedding visualization | [embedding-visualization](../embedding-visualization/SKILL.md) |
| 10 | Clustering + marker genes | [clustering-and-annotation](../clustering-and-annotation/SKILL.md) |
| 11 | Pseudobulk DE | [pseudobulk-differential-expression](../pseudobulk-differential-expression/SKILL.md) |
| *(opt)* | Cluster validation | [cluster-validation](../cluster-validation/SKILL.md) |

## Detail Files

- **[IMPLEMENTATION.md](./IMPLEMENTATION.md)** — notebook-style code for the legacy barcoding path

## Related Skills

- [3prime-untargeted-pipeline](../3prime-untargeted-pipeline/SKILL.md) — current-product RNA
- [ops-pipeline](../ops-pipeline/SKILL.md) — OPS pooled screening
- [imaging-pipeline](../imaging-pipeline/SKILL.md) — imaging / CellPaint / IF
- [spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md)
- [normalization-strategies](../normalization-strategies/SKILL.md)
- [dimensionality-reduction](../dimensionality-reduction/SKILL.md)
- [embedding-visualization](../embedding-visualization/SKILL.md)
- [clustering-and-annotation](../clustering-and-annotation/SKILL.md)
- [pseudobulk-differential-expression](../pseudobulk-differential-expression/SKILL.md)

## Disclaimer

This content is provided "as-is" without warranty of any kind. AI-generated output may contain
inaccuracies. Review and verify all outputs before relying on them.
