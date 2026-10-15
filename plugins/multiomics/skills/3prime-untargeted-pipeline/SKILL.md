---
name: 3prime-untargeted-pipeline
description: End-to-end workflow for AVITI24 3′ untargeted DISS transcriptome data (ThreePrimeUntargeted table). Covers HVG selection, mitochondrial QC, filter_genes, and de novo cell type discovery. Use when analyzing 3′ untargeted data or when the table is ThreePrimeUntargeted — not for legacy Transcript/Protein barcoding (use legacy-targeted-pipeline) or imaging (use imaging-pipeline). Prefer depth / total_counts language, not UMIs.
metadata:
  version: 0.1.9
  author: elembio
---

# 3′ Untargeted Pipeline

> **Confirm the SDO first.** Before running, make sure you are using the SDO the user selected in [locate-spatialdata-store](../locate-spatialdata-store/SKILL.md) — a run may have several (instrument-direct, `cells2stats`, `elembio-multiomics-nf`). If 3′ is declared in the panel but `ThreePrimeUntargeted` isn't in the store's tables, that SDO lacks the aligned 3′ counts — stop and resolve it via locate-spatialdata-store ([If aligned 3′ is missing](../locate-spatialdata-store/SKILL.md#if-aligned-3-is-missing--call-out-and-confirm)) before running this pipeline.

## Purpose

This skill provides an end-to-end workflow for AVITI24 3′ untargeted **DISS** transcriptome data. Unlike targeted panels, it requires `ebsd.pp.filter_genes`, highly variable gene (HVG) selection, mitochondrial percentage filtering, and well-aware normalization via `ebsd.pp.normalize` before dimensionality reduction and clustering. Cell types are discovered de novo rather than assumed from panel design. Counts are DISS depth / `total_counts` — not UMIs.

## Prerequisites

- AVITI24 SpatialData Zarr store with a `ThreePrimeUntargeted` table
- Python **3.12+** with `elembio-spatialdata-analysis`

**Loading:** Use `Loader` + `load_tables(table="ThreePrimeUntargeted")`; see [spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md). Use `read_zarr` only when you need images or labels.

**Filtering:** Use `elembio_spatialdata_analysis.pp` — `ebsd.pp.compute_qc_metrics` → `ebsd.pp.filter_cells` → `ebsd.pp.filter_genes` with `"ThreePrimeUntargeted"` (same contract as [cell-quality-control](../cell-quality-control/SKILL.md)). `ebsd.pp.compute_qc_metrics` may call Scanpy internally for MT / ribo / HB metrics when those columns are missing; `ebsd.pp.filter_cells` reads `obsm["qc_metrics"]` for `max_pct_mt`, not ad-hoc `obs` slices.

## Key Concepts

### Targeted vs 3′ Untargeted

For HVG method comparison and reference database options, see [REFERENCE.md](./REFERENCE.md).

| Aspect | Targeted | 3′ Untargeted |
|--------|----------|---------------|
| **Genes** | Panel of 100–500 | Thousands detected |
| **HVG selection** | Not needed | Required (~2000 genes) |
| **MT% QC** | Usually not needed | Important — use `max_pct_mt` |
| **Normalization** | Optional / minimal | `ebsd.pp.normalize` (well-aware agg-count scaling → median-of-ratios → log1p) |
| **Cell typing** | Often known from panel | De novo discovery |
| **Reference DBs** | Panel design | CellMarker, PanglaoDB, HCA |

### Pipeline Order

For the full step-by-step pipeline code, expected outputs, and common pitfalls, see [IMPLEMENTATION.md](./IMPLEMENTATION.md).

1. Load `ThreePrimeUntargeted` table
2. `ebsd.pp.compute_qc_metrics` → visualize QC
3. `ebsd.pp.filter_cells` with `max_pct_mt` → `ebsd.pp.filter_genes`
4. `ebsd.pp.normalize(adata, "ThreePrimeUntargeted")` → well-aware agg-count scaling → median-of-ratios → log1p
5. HVG selection (`highly_variable_genes`)
6. Scale → PCA (on HVG subset)
7. Assess batch effects
8. `neighbors` → UMAP
9. Leiden clustering at multiple resolutions
10. Marker gene identification (`rank_genes_groups`)
11. Cell type annotation (manual or reference-based)
12. *(Optional)* Cell cycle scoring
13. *(Optional)* Pseudobulk DE
14. Save results

## Quick Start

<!-- canonical: spatialdata-loading-and-access memory guard. Copy verbatim; edit only in that skill's SKILL.md -->
**Memory guard — before any `lazy=False`, `.to_memory()` or `.compute()`:**

1. Load lazily first and run `check_before` from `scripts/check_load.py` in `spatialdata-loading-and-access`. Run the script unchanged in the Python session that holds the data (`exec` its contents); don't rewrite it. It compares each table's estimated in-memory size with the memory this session can use, lists large session objects you could `del`, and raises on a FAIL.
2. Materialize only the tables and runs the plan needs. With several runs, load `concatenated=False` so it checks each run and their total. A WARN says what to `del` or to load one run at a time. A FAIL means subset before materializing: propose a subset to the user (runs, wells, or a set of features), then slice the lazy table and call `.to_memory()` on the slice, for example `adata[adata.obs["Well"].isin(wells)].to_memory()`. `.copy()` keeps a lazy table lazy. Feature subsets by modality: [ACCESS.md — Subset](../spatialdata-loading-and-access/ACCESS.md#subset-by-well-or-feature).
3. Keep `.X` in the form it was stored (`scipy.sparse.issparse(adata.X)`). Sparse tables (3′, OPS) are mostly zeros: never densify a whole one (`.toarray()`, `.todense()`, `np.asarray(adata.X)`); densify only a subset of cells or features. Dense tables (imaging, legacy barcoding) are float32 with `NaN` where a batch has no data: do not convert them to sparse, which turns missing into zero.
<!-- /canonical -->

```python
from elembio_spatialdata_analysis.load import Loader
import elembio_spatialdata_analysis as ebsd
import scanpy as sc

loader = Loader("/path/to/run.zarr")
check_before(loader.load_tables(table="ThreePrimeUntargeted", concatenated=False))  # memory guard: lazy, raises on FAIL
adata = loader.load_tables(table="ThreePrimeUntargeted", lazy=False)

ebsd.pp.compute_qc_metrics(adata, "ThreePrimeUntargeted")

ebsd.pp.filter_cells(adata, "ThreePrimeUntargeted", inplace=True)  # max_pct_mt=0.10 default
ebsd.pp.filter_genes(adata, "ThreePrimeUntargeted", min_cells_per_gene=10)  # always inplace

print(f"After QC: {adata.shape}")
```

## Detail Files

- [IMPLEMENTATION.md](./IMPLEMENTATION.md) — Full step-by-step pipeline code (steps 1–15), Expected Outputs, Common Pitfalls
- [REFERENCE.md](./REFERENCE.md) — HVG method comparison and reference database options

## Related Skills

- [cell-quality-control](../cell-quality-control/SKILL.md) — QC metrics and filtering contract
- [dimensionality-reduction](../dimensionality-reduction/SKILL.md) — HVG selection and PCA details
- [clustering-and-annotation](../clustering-and-annotation/SKILL.md) — Cell type annotation strategies
- [cell-cycle-analysis](../cell-cycle-analysis/SKILL.md) — Cell cycle scoring
- [pseudobulk-differential-expression](../pseudobulk-differential-expression/SKILL.md) — Pseudobulk DE
- [batch-effect-assessment](../batch-effect-assessment/SKILL.md) — Batch correction guidance

## Disclaimer

This content is provided "as-is" without warranty of any kind. AI-generated output may contain
inaccuracies. Review and verify all outputs before relying on them.
