---
name: index
description: Index of all spatial-data analysis skills with one-line "Use when…" triggers. Read first to pick the right specialist before generating analysis code, then open and follow the matching SKILL.md. For run analysis use Loader (not read_zarr on zip), ebsd cell QC (not Scanpy calculate_qc_metrics), and DISS depth vocabulary (not UMIs).
metadata:
  version: 0.3.7
  author: elembio
---

# Spatial Data Skills — Index

**Read this first.** Pick the skill(s) whose "Use when…" line matches the user's request, then **open and follow** the matching `SKILL.md` before writing analysis code.

**Analysis package:** `elembio_spatialdata_analysis` **≥ 0.16.0** (`import elembio_spatialdata_analysis as ebsd`). Skills in this stack are validated against this version — older versions may be missing APIs referenced here.

### Platform defaults (AVITI24 / DISS)

These override generic single-cell habits:

| Task | Use | Do not |
|---|---|---|
| Tables, QC, normalize, embeddings | `Loader` + `load_tables` | `spatialdata.read_zarr` on a `.zarr.zip` for table work |
| Full images / labels from a zip | FUSE-mount then `read_zarr` (only if rasters are needed) | Assume `read_zarr` opens zip natively |
| Per-cell QC / filtering | `ebsd.pp.compute_qc_metrics` → `filter_cells` ([cell-quality-control](../cell-quality-control/SKILL.md)) | `sc.pp.calculate_qc_metrics` as the AVITI24 QC path |
| Count depth language | depth / `total_counts` / polony counts (DISS) | “UMI” / droplet-seq vocabulary |
| Pipeline order | load → run-quality-overview → cell-quality-control → normalize (counts only) | Skip specialist skills and improvise with Scanpy defaults |

---

## Quick start — orienting yourself on a new analysis

If the user asks to analyze a run or is new to the stack, read in this order **before** generating code:

0. [locate-spatialdata-store](../locate-spatialdata-store/SKILL.md) — resolve a run name, directory, or path to a validated `.zarr` store; discover available tables. **When a run has several SDOs (instrument-direct, `cells2stats`, `elembio-multiomics-nf`), confirm which one the user wants before proceeding.**
1. [spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md) — open the `.zarr` store with `Loader`; discover `available_tables` (never invent modalities)
2. [run-quality-overview](../run-quality-overview/SKILL.md) — run / well / tile health for **present** modalities (no cell drops)
3. [cell-quality-control](../cell-quality-control/SKILL.md) — per-cell `filter_cells` / `filter_modalities` on what exists
4. [normalization-strategies](../normalization-strategies/SKILL.md) — normalize **count** modalities only (not imaging / not OPS)
5. Modality pipeline(s) for what is present — [3prime-untargeted-pipeline](../3prime-untargeted-pipeline/SKILL.md), [ops-pipeline](../ops-pipeline/SKILL.md), [imaging-pipeline](../imaging-pipeline/SKILL.md)
6. [code-generation-guidelines](../code-generation-guidelines/SKILL.md) — imports, order, and platform defaults when writing code

Optional orientation: [anndata-structure-reference](../anndata-structure-reference/SKILL.md) (slots by present modality).

**Shortcuts — skip step 0 when:**
- You already have a `.zarr` / `.zarr.zip` path → go straight to step 1.
- Input is not a `.zarr` store (`.h5ad`, parquet, CSV, …) → swap step 1 for [raw-data-ingestion](../raw-data-ingestion/SKILL.md).

---

## Pipelines (end-to-end orchestration)

- [3prime-untargeted-pipeline](../3prime-untargeted-pipeline/SKILL.md) — Use when: data has a `ThreePrimeUntargeted` table (3′ untargeted DISS transcriptome).
- [ops-pipeline](../ops-pipeline/SKILL.md) — Use when: data has a true OPS (Optical Pooled Screening) identity library — DISS `SpecializedTargeted`, not QC/housekeeping gene panels that also label as OPS.
- [imaging-pipeline](../imaging-pipeline/SKILL.md) — Use when: imaging / morphology / CellPaint / IF (ProteinIF) — inventory, S/B, localization, MAD-z feature matrix and its design-dependent scaling reference (control wells / global / per-well), well-level treatment contrasts. Then DR + embedding-visualization for UMAP. Not legacy barcoding Protein.
- [legacy-targeted-pipeline](../legacy-targeted-pipeline/SKILL.md) — Use when: **legacy** Transcript / Protein barcoding panels only — not for current 3′ / OPS / imaging-first runs.

## Data access

- [locate-spatialdata-store](../locate-spatialdata-store/SKILL.md) — Use when: user provides a run name, a directory, or asks to analyze a run; resolves to a validated `.zarr` path — and confirms which SDO to use when a run has several — before loading.
- [spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md) — Use when: opening AVITI24 `.zarr` stores; discovering and loading present modalities (3′, OPS, imaging/CellPaint alias, legacy Transcript/Protein).
- [raw-data-ingestion](../raw-data-ingestion/SKILL.md) — Use when: input is `.h5ad`, AnnData `.zarr`, AVITI24 `RawCellStats.parquet`, or generic CSV / parquet count matrices.
- [elembio-cloud-data-access](../elembio-cloud-data-access/SKILL.md) — Use when: pulling user data from ElemBio Cloud — listing or downloading runs, executions, storage; FUSE-mounting via the `elembio` CLI.

## Orientation & metadata

- [anndata-structure-reference](../anndata-structure-reference/SKILL.md) — Use when: you need to know where counts, QC metrics, geometry, controls, or filter flags live (canonical `obs` / `var` / `obsm` / `layers` / `uns` slot layout).
- [metadata-and-target-mapping](../metadata-and-target-mapping/SKILL.md) — Use when: listing wells, tiles, batches, panels, or targets from `RunInfo` / `global-table-metadata` (`uns`); print `Well` + `WellLabel` as recorded. Prefer `Loader` over raw `uns` parsing. Not for inferring condition groupings.
- [experiment-design](../experiment-design/SKILL.md) — Use when: user asks what conditions, replicates, cell lines, treatments, or doses are in the run (inferred from `WellLabel`); also reports whether condition is confounded with `Well`, which decides the imaging scaling reference.

## QC & preprocessing

- [run-quality-overview](../run-quality-overview/SKILL.md) — Use when: cells per well / tile, segmentation quality, well-to-well consistency, or overall run quality — first QC step after loading; not per-cell filtering.
- [cell-quality-control](../cell-quality-control/SKILL.md) — Use when: filtering low-quality cells after run QC — 3′ (MT + `filter_genes`), Imaging (`passes_imaging_intensity`), OPS (area-only), or legacy panels; only for modalities present in the store.
- [normalization-strategies](../normalization-strategies/SKILL.md) — Use when: normalizing **counts** (agg / median-of-ratios / log1p). Imaging / IF scaling → [imaging-pipeline](../imaging-pipeline/SKILL.md) (MAD-z; not `pp.normalize`).

## Dimensionality reduction & batch

- [dimensionality-reduction](../dimensionality-reduction/SKILL.md) — Use when: computing PCA/UMAP/neighbors (canonical). Counts after normalize; imaging after the imaging-pipeline MAD-z matrix (check which scaling reference it used).
- [batch-effect-assessment](../batch-effect-assessment/SKILL.md) — Use when: detecting or quantifying batch effects across wells, runs, or molecular batches; deciding whether to apply Harmony / scVI.
- [cell-cycle-analysis](../cell-cycle-analysis/SKILL.md) — Use when: scoring cell-cycle phase (S / G2M) or regressing out proliferation effects on normalized AnnData before downstream steps.

## Clustering, annotation & validation

- [clustering-and-annotation](../clustering-and-annotation/SKILL.md) — Use when: Leiden clustering, marker-gene identification (Wilcoxon), resolution selection, cell-type annotation.
- [cluster-validation](../cluster-validation/SKILL.md) — Use when: validating cluster robustness — silhouette scores, classifier-based reconstruction, batch-contamination checks. Invoke deliberately, not by default in the standard pipeline.

## Differential expression

- [pseudobulk-differential-expression](../pseudobulk-differential-expression/SKILL.md) — Use when: pseudobulk DE with PyDESeq2, aggregating by sample / well / replicate (requires ≥ 2 biological replicates per condition).
- [de-visualization-and-interpretation](../de-visualization-and-interpretation/SKILL.md) — Use when: volcano / MA plots, heatmaps of top DEGs, GSEA, or pathway analysis on a DE results table.

## Visualization

- [embedding-visualization](../embedding-visualization/SKILL.md) — Use when: plotting/coloring UMAP or PCA (canonical) — genes, clusters, `Well` (plate position), `WellLabel` / condition (sample name), or imaging `sb_*` / `nucfrac_*`.
- [spatial-multimodal-visualization](../spatial-multimodal-visualization/SKILL.md) — Use when: well fluorescence (morphology stains **or protein IF channels**), the flow-cell plate, overlaying segmentation, filling cells by a table column, or a single / representative cell close-up (`plot_image` / `plot_cell`).

## Spatial statistics

- [spatial-statistics](../spatial-statistics/SKILL.md) — Use when: neighborhood analysis, spatial autocorrelation (Moran's I, Geary's C), neighborhood enrichment, Ripley's statistics (squidpy on `obsm["spatial"]`, per well).

## Multimodal

- [multimodal-data-preparation](../multimodal-data-preparation/SKILL.md) — Use when: aligning present modalities (3′ / OPS / imaging) before joint analysis — per-modality transforms (count normalize vs imaging MAD-z).
- [integrated-analysis](../integrated-analysis/SKILL.md) — Use when: joint analysis of aligned multimodal data — concatenated PCA / UMAP, modality-weighted integration, CCA, MOFA+.

## Code conventions

- [code-generation-guidelines](../code-generation-guidelines/SKILL.md) — Use when: writing elembio analysis code for local or interactive (notebook) use — discover → run QC → cell QC → normalize-only-counts order, imports, performance.

---

## Disclaimer

This content is provided "as-is" without warranty of any kind. AI-generated output may contain
inaccuracies. Review and verify all outputs before relying on them.
