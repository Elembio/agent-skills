# Run Quality Overview — Reference

## Modality patterns

| Pattern | Typical modalities |
|---|---|
| Full multimodal | 3′ + OPS + imaging (often listed as `CellPaint`) |
| Transcript + imaging | 3′ + imaging (`CellPaint`) |
| Legacy targeted | Transcript + Protein + CellPaint |
| Morphology-only | OPS + `imaging` |

Skills must discover per run — do not hardcode one flagship.

## uns key map (`run._metadata`)

Try newer prefix first.

| Rollup | Newer | Older |
|---|---|---|
| Per-well cytometry | `c2s_runstats_cyto_wells` | `runstats_cyto_wells` |
| Per-well raw reads | `c2s_runstats_rawread_wells` | `runstats_rawread_wells` |
| Per-batch cytometry | `c2s_runstats_cyto_batches` | `runstats_cyto_batches` |
| 3′ sample metrics (optional) | `f2p_3p_runstats_*` | — |
| Target-site stats (optional) | `c2s_runstats_targetstats_*` | — |

## Cytometry columns (`*runstats_cyto_wells`)

| Column | Meaning |
|---|---|
| `WellLocation` | Physical well (`A1`) — same role as `obs["Well"]` |
| `CellCount` | Segmented cells (pre per-cell QC) |
| `PercentNucleatedCells` | Nucleation % (0–100) |
| `PercentConfluency` | Confluency % (0–100) |
| `MedianCellDiameter` | Median diameter (px) |
| `MeanAssignedBarcodingCountPerCell` | Barcoding-assignment metric, not an OPS health metric; use AnnData `Confidence` / `TotalCount` |
| `AssignedCountsPerMM2` | Barcoding-assignment metric, not an OPS health metric |

No sample-sheet column on cyto. Backfill `WellLabel` from `run.wells` (`well_label`). Preserve a native `WellLabel` if a later store adds one.

## OPS assignment columns (`include_target_assignment=True`)

| Column | Meaning |
|---|---|
| `MostAbundantTarget` | Winning barcode / `"Undetermined"` |
| `MostAbundantTargetCount` | Reads for winner |
| `TotalCount` | Total barcode reads for cell |
| `Confidence` | `MostAbundantTargetCount / (TotalCount + 1)` (+1 pseudocount; 0 if no reads) |

**OPS assay vs targeted DISS QC/gene panel:** both are DISS (`DISSBatch` / `SpecializedTargeted`). Metadata `target_type=OPS` covers OPS identity libraries **and** targeted gene/QC panels (e.g. `B01-QC_Panel` = GAPDH+ACTB). **Confidence≥0.5 is only for OPS assay tables.** For QC/gene panels, report depth/detection only — never headline Confidence. This is not legacy `BarcodingBatch`.

Discover suffixes: `run.table_suffixes_for_target_type(("OPS", "SpecializedTargeted"))`, then split OPS assay vs `qc_panel` / housekeeping before any OPS metrics.

## Soft gate bands (guidance only)

From `elembio_spatialdata_analysis/qc/config/qc_gates.yml`. Guidance only, **not** hard filters in this skill:

| Metric | Warn guidance |
|---|---|
| Nucleation (fraction) | below ~0.85 (convert cyto `%` ÷ 100) |
| Confluency % | below ~20 or above ~70 |
| Cell count / well | below ~2000 (context-dependent) |
| OPS assay assigned-cells % | `ops_cells_assigned_rate` warn below ~50 (assay sites) |
| OPS run assignment % | `run_assignment_rate_ops` warn below ~5 (different scale — see YAML comments) |
| QC control assigned % | **inverted** — high single-assign can be bad (`ops_control_assigned_cells`) |

Prefer relative well-to-well comparison over absolute cutoffs.

## Cyto run-QC checklist (core checks)

Cytoprofiling run health is commonly framed as a fixed set of checks. This table
maps each to its metric and gate in
`elembio_spatialdata_analysis/qc/config/qc_gates.yml`.

| Check | Metric | Current gate (warn / fail) | Notes |
|---|---|---|---|
| Assigned rate | `assignment_rate` (per-cell) or `run_assignment_rate_protein` (run-level) | 50% / 40% (per-cell) — 80% / 60% (run-level, Protein) | Two candidate metrics at different scales; state which one you're citing |
| Mismatch rate | `mismatch_assignment_rate` | >30% / >50% | Percent of assigned polonies with ≥1 base mismatch |
| Extracellular ratio | *(no gate)* | — | Metric exists (`qc.aggregate` well output, `ExtraCellularRatioBeforeTrimming`) but has no threshold in `qc_gates.yml` |
| Cell count / well | `cell_count_per_well` | <2000 / <1000 | Context-dependent — compare wells, don't treat as a hard floor |
| Confluency % | `confluency_pct` | <20 or >70 / <15 or >80 | |
| Percent nucleated | `nucleation_rate` | <85% / <70% | Cyto reports 0-100; gate is 0-1 fraction — convert |
| GAPDH (positive control) | `passes_gapdh_outlier` (`pp.filter`) | opt-in, quantile-based (`max_quant_gapdh`, default `None`) | Per-cell relative outlier flag, not an absolute density gate — a different design than a fixed counts/mm² threshold, not a gap |
| NSB (non-specific binding) | `passes_nsb_outlier` (`pp.filter`) | quantile-based (`max_quant_nsb`) | Same design as GAPDH — relative, not absolute |
| Assigned counts / mm² | `counts_per_mm2_per_batch` | <100,000 / <50,000 | Shared across Protein/Transcript/3′ — confirm this is the metric you mean before citing a specific number |

**Do not cite a Q30 threshold as part of this checklist.** `overall_q30` exists in
`qc_gates.yml` but is a **sequencing chemistry** metric (bases2fastq), not a
cytoprofiling run-health check — there is no Q30 gate for cyto/cell metrics.

## `qc.aggregate` helpers

```python
from elembio_spatialdata_analysis.qc import per_tile_cell_count, per_tile_total_counts, per_well_overview

per_tile_cell_count(adata) -> dict[str, int]
```

## Scope boundary

| This skill | Specialist |
|---|---|
| Run/well/modality **summaries** | — |
| Soft interpret / flag | — |
| Per-cell `filter_cells` | [cell-quality-control](../cell-quality-control/SKILL.md) |
| OPS merge / DE / confidence filtering | [ops-pipeline](../ops-pipeline/SKILL.md) |
| Imaging intensity gate | [cell-quality-control](../cell-quality-control/SKILL.md) (`passes_imaging_intensity`) |
