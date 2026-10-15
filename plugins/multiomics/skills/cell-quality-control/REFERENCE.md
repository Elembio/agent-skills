# Cell Quality Control — Reference

Defaults match `DEFAULT_BOUNDS_PER_MODALITY` in `elembio_spatialdata_analysis.pp.filter`.

## Table of Contents

1. [Defaults by Modality](#defaults-by-modality)
2. [`qc_metrics` Columns](#qc_metrics-columns)
3. [Filter Flags Columns](#filter-flags-columns)
4. [Modality-Specific Behavior](#modality-specific-behavior)

## Defaults by Modality

### ThreePrimeUntargeted (standard RNA)

```python
filter_cells(adata, "ThreePrimeUntargeted")
# min_quant_count=0.025, max_quant_count=0.99
# min_quant_area=0.01, max_quant_area=0.99
# max_pct_mt=0.10          # ON by default for 3′
# No NSB / GAPDH / assignment-rate gates — do not narrate them for 3′
# Count gate (passes_counts): summed adata.X totals per well — not per-batch
# filter_spatial_artifacts=False  # opt in to gate on any-ThreePrimeUntargeted
```

### Imaging (alias: CellPaint)

```python
filter_cells(adata, "Imaging")   # or "CellPaint" — normalized to imaging
# min_quant_count=0.025, max_quant_count=0.99
# min_quant_area=0.01, max_quant_area=0.99
# passes_nan_filter, passes_empty_filter, passes_imaging_intensity all scoped
# to CP* primary morphology channels (IF / Bkg / CellProfiler / compartment excluded).
# NaN confined to those excluded columns does NOT fail a cell.
# qc_metrics['total_counts'] is also morphology-channel-only for imaging.
# filter_spatial_artifacts=False  # opt in to gate on any-imaging
```

### OPS

```python
filter_cells(adata, "OPS")
# min_quant_count=None, max_quant_count=None   # count quantiles OFF
# min_quant_area=0.01, max_quant_area=0.99     # area ON
# No NSB / GAPDH / MT
# Confidence filtering is NOT this gate — see ops-pipeline / run-quality-overview
# filter_spatial_artifacts=False  # opt in to gate on any-OPS
```

### Legacy Transcript (only place NSB / GAPDH / assignment-rate apply)

```python
filter_cells(adata, "Transcript")
# min_quant_count=0.025, max_quant_count=0.99
# area 0.01 / 0.99; max_quant_nsb=0.99 when control_variables present
# min_assignment_rate=None; max_quant_gapdh=None; max_pct_mt=None
# filter_spatial_artifacts=False  # no Transcript rollup — opt-in is a no-op
# Narrate NSB/GAPDH/assignment only when this table is present
```

### Legacy Protein

```python
filter_cells(adata, "Protein")
# count + area quantiles; NSB when present; no MT
```

### Opt-in examples

```python
# Transcript only:
filter_cells(adata, "Transcript", min_assignment_rate=0.5)
filter_cells(adata, "Transcript", max_quant_gapdh=0.99)
# 3′ only:
filter_cells(adata, "ThreePrimeUntargeted", max_pct_mt=0.2)  # looser than default 0.10 (keeps cells with up to 20% mitochondrial reads)
# Any artifact-bearing modality (Imaging / 3′ / OPS):
filter_cells(adata, "Imaging", filter_spatial_artifacts=True)
```

### `filter_genes` (RNA only — 3′ / Transcript)

```python
filter_genes(adata, "ThreePrimeUntargeted")                          # default is 10 for 3′
filter_genes(adata, "ThreePrimeUntargeted", min_cells_per_gene=100)  # common tighten
# Skip for Imaging / OPS / Protein — their default is None, so a bare call no-ops
```

| `min_cells_per_gene` | Guidance |
|---|---|
| **10** | Default for 3′ and a sensible floor |
| 50–100 | Reasonable escalate when sparsity / noisy embeddings persist |
| 100–500 | Large `n_obs` only; prefer ~1–5% of cells; ask before applying |

Not a substitute for HVG selection (after normalize).

### Quantile overrides

`QuantileSpec` is a type alias (`float | Mapping | None`), not a constructor class. Override with floats or well/batch mappings:

```python
filter_cells(
    adata,
    "ThreePrimeUntargeted",
    min_quant_area=0.01,
    max_quant_area=0.99,
)
```

---

## `qc_metrics` Columns

Written by `compute_qc_metrics`. **Discover at runtime** — do not hardcode.

| Column | 3′ | Imaging | OPS | Transcript | Protein |
|---|:---:|:---:|:---:|:---:|:---:|
| `total_counts` | ✓ | ✓ | ✓ | ✓ | ✓ |
| `n_genes_by_counts` | ✓ | — | — | ✓ | — |
| `total_assignment_rate` | △ | — | — | △ | △ |
| `pct_counts_mt` | △ | — | — | — | — |
| `total_nsb_counts` / `total_gapdh_counts` | — | — | — | ✓* | — |
| `nuclear_frac` | △ | — | — | △ | △ |

✓ always when modality applies; △ when source data present; \* may be zeros if controls absent.

Imaging `total_counts` is the sum of **primary morphology channels** (`CP*` whole-cell base) only — aligned with the NaN / empty / intensity gates so Bkg / IF NaNs do not poison the metric. Other modalities sum the full matrix or `obsm['total_counts']`.

---

## Filter Flags Columns

| Flag | Typical modality | Meaning |
|---|---|---|
| `passes_all_filters` | all | Conjunction of applied flags |
| `passes_nan_filter` | all | NaN check. **Imaging:** scoped to `CP*` primary morphology only — NaN in IF / Bkg / CellProfiler / compartment columns is ignored. Other modalities: full matrix. |
| `passes_empty_filter` | all | Empty-cell check. **Imaging:** row sum of `CP*` primary morphology only. Other modalities: full matrix. |
| `passes_counts` | 3′, Transcript, Protein | Count-depth quantile gate. **3′:** summed `.X` per well. **Transcript/Protein:** per batch via `obsm['total_counts']`. (Legacy alias `passes_counts_per_batch` may appear in old reports.) |
| `passes_area` | all with area gate | Area quantile |
| `passes_imaging_intensity` | Imaging | Morphology (`CP*`) intensity gate |
| `passes_imaging_counts` | Imaging (legacy name) | Alias of intensity gate — prefer `passes_imaging_intensity` |
| `passes_mt_outlier` | 3′ | MT % gate |
| `passes_nsb_outlier` | **Transcript only** | NSB gate — omit from non-Transcript narratives |
| `passes_gapdh_outlier` | **Transcript only** (opt-in) | GAPDH gate |
| `passes_assignment_rate` | **Transcript only** (opt-in) | Assignment-rate gate (not OPS confidence) |
| `passes_spatial_artifacts` | Imaging / 3′ / OPS (opt-in) | Not in `obsm["spatial_artifacts"]` modality rollup. Default off (all-True; report hides). Transcript / Protein: no-op. |

---

## Modality-Specific Behavior

| Aspect | 3′ | Imaging | OPS | Transcript | Protein |
|---|---|---|---|---|---|
| Count quantiles | Yes (`passes_counts`, summed `.X`) | Intensity (`passes_imaging_intensity`) | **No** | Yes (`passes_counts`, per batch) | Yes (`passes_counts`, per batch) |
| Area | Yes | Yes | Yes | Yes | Yes |
| MT | Default 0.10 | — | — | — | — |
| Primary specialty flag | `passes_mt_outlier` | `passes_imaging_intensity` | — | NSB | — |
| Normalize after? | Yes | **No** | **No** | Yes | Yes |

`filter_modalities`: area owner is Transcript or ThreePrimeUntargeted when present; only pass modalities that exist in the run.
