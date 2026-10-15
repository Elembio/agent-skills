---
name: cell-quality-control
description: Cell- and gene-level QC filtering for AVITI24 data after run-quality-overview — compute_qc_metrics, filter_cells / filter_modalities for modalities actually present, and filter_genes for RNA (3′ / Transcript) before normalize. Use when filtering low-quality cells or genes, setting QC thresholds, or generating filter reports. Not for run/well rollups. Never invent modalities absent from the store.
metadata:
  version: 0.2.13
  author: elembio
---

# Cell Quality Control

## Purpose

QC with `elembio_spatialdata_analysis.pp`: compute metrics, filter low-quality **cells**, and for **RNA** also drop near-zero **genes** (`filter_genes`) before normalize. Comes **after** [run-quality-overview](../run-quality-overview/SKILL.md).

```text
Load → run-quality-overview → cell-quality-control
  RNA: compute_qc_metrics → filter_cells → filter_genes → normalize
  Imaging / OPS: compute_qc_metrics → filter_cells  (no filter_genes, no normalize)
```

**Hard rules:**
- Discover `available_tables` / types first; filter **only** modalities that exist.
- Use **`ebsd.pp.compute_qc_metrics` → `ebsd.pp.filter_cells`** (this skill). Do **not** use `sc.pp.calculate_qc_metrics` as the AVITI24 cell-QC path — Scanpy metrics omit Element geometry, MT defaults, imaging intensity, and OPS/area rules.
- Narrate **only** filters that apply to the modality in hand. Do **not** mention NSB, GAPDH, assignment-rate, MT, or imaging-intensity unless that modality uses them and the inputs exist.

## Prerequisites

- [spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md)
- Prefer [run-quality-overview](../run-quality-overview/SKILL.md) first (well/run health)
- AnnData from `load_tables(..., lazy=False)` for each present modality

## Key Concepts

### What `filter_cells` does

Modality-aware per-cell tests; cell must pass all applicable flags (`passes_all_filters`). Use the **specialty** row for the present modality only — ignore the rest.

**Common (most modalities):**

| Test | Catches | Default |
|---|---|---|
| Min/max counts (quantile) → `passes_counts` | Empty / extreme depth cells | On for 3′ (summed `.X`) and legacy Transcript/Protein (per batch); **off for OPS** |
| Min/max cell area | Debris / clumps | On (shared segmentation) |

**Specialty — only if that modality is present:**

| Modality | Specialty gate | Default |
|---|---|---|
| **ThreePrimeUntargeted** | Mitochondrial % (`passes_mt_outlier`); count gate `passes_counts` on **summed** `.X` per well (not per-batch) | MT **on** at `max_pct_mt=0.10` |
| **Imaging** | NaN, empty, and intensity gates all scoped to `CP*` primary morphology (`passes_nan_filter`, `passes_empty_filter`, `passes_imaging_intensity` via `select_morphology_channels`). `total_counts` in `qc_metrics` is also morphology-only. IF (ProteinIF, e.g. `P01`) columns are not evaluated by these gates, so their NaNs remain; handle them in [imaging-pipeline](../imaging-pipeline/SKILL.md) before building the analysis matrix. | On; NaN/empty in IF / background / CellProfiler / compartment columns do **not** fail cells. Alias `target_type="CellPaint"` OK |
| **OPS** | *(none beyond area)* | Count quantiles **off**; confidence is not this gate — see run-QC / [ops-pipeline](../ops-pipeline/SKILL.md) |
| **Legacy Transcript** | NSB / optional GAPDH / optional assignment rate | See subsection below — **do not discuss on 3′, imaging, or OPS** |
| **Legacy Protein** | Optional NSB when `control_variables` present | Counts + area; no MT |

#### Legacy Transcript only (skip unless `Transcript` is present)

NSB, GAPDH, and assignment-rate are **not** part of 3′ / imaging / OPS QC. Mention them only for Transcript:

| Test | Default |
|---|---|
| NSB outlier (`passes_nsb_outlier`) | On when `control_variables` present |
| GAPDH outlier (`passes_gapdh_outlier`) | **Off** — opt in via `max_quant_gapdh` |
| Assignment rate (`passes_assignment_rate`) | **Off** — opt in via `min_assignment_rate` (not OPS confidence) |

### Always `compute_qc_metrics` before `filter_cells`

Writes `obsm["qc_metrics"]` for plots/reports. Column set is modality-specific — guard with `if col in ...columns`.

### RNA: `filter_genes` is required (not optional)

For **ThreePrimeUntargeted** and legacy **Transcript**, after `filter_cells` and **before** normalize:

```python
# 3′ defaults to min_cells_per_gene=10. Legacy Transcript has no default — pass it there.
ebsd.pp.filter_genes(adata, "ThreePrimeUntargeted")
```

Skip for Imaging, OPS, and Protein. Skipping gene filter on 3′ leaves ultra-sparse genes that often make PCA/UMAP noisy even when cells look fine. Gene filter ≠ HVG (HVG is later, after normalize).

**Start at 10, then escalate with evidence.** Print `n_obs` / `n_vars` before and after; only raise thresholds if distributions or a bad embedding warrant it — pause for approval before large jumps.

| `min_cells_per_gene` | When |
|---|---|
| **10** | Recommended starting floor (pass explicitly) |
| ~50–100 | Common tighten when many near-zero genes remain or UMAP looks noisy/splotchy |
| ~100–500 | Only on large `n_obs`; think as a **fraction of cells** (e.g. ~1–5%), not a magic absolute. Confirm with the user before applying |

### Critical defaults

For the complete defaults table per modality (source of truth), see [REFERENCE.md](./REFERENCE.md).

| Parameter | Default | Notes |
|---|---|---|
| `inplace` | `True` | Use `False` to keep the unfiltered copy |
| `diagnostics` | `False` | With `inplace=False`, `True` → `(adata, FilterDiagnostics)` |
| `max_pct_mt` (3′ only) | **`0.10`** | Do not discuss for non-3′ modalities |
| `min_cells_per_gene` (RNA) | pass **`10`** | Recommended start; pass explicitly (package bound may be `None`) |
| `max_quant_gapdh` / `min_assignment_rate` | `None` (off) | **Transcript only** — do not pass or narrate elsewhere |
| `filter_spatial_artifacts` | `False` (off) | Opt in to drop cells in `obsm["spatial_artifacts"]` modality rollup (`any-imaging` / `any-ThreePrimeUntargeted` / `any-OPS`). Report hides the row when off. |

### Filter flags (report only what ran)

- Most modalities: `passes_nan_filter`, `passes_empty_filter` (full matrix), `passes_counts`, `passes_area`, `passes_all_filters`
- **Imaging**: `passes_nan_filter`, `passes_empty_filter`, `passes_imaging_intensity` — all three scoped to `CP*` primary morphology channels only (IF / background / CellProfiler / compartment NaNs do not fail cells). IF-column NaNs remain; handle them in [imaging-pipeline](../imaging-pipeline/SKILL.md) before matrix construction.
- 3′: `passes_mt_outlier`
- Transcript only: `passes_nsb_outlier`; optional `passes_gapdh_outlier` / `passes_assignment_rate`
- Optional (any artifact-bearing modality): `passes_spatial_artifacts` — only narrate when `filter_spatial_artifacts=True`

### Load / lazy

Use `Loader.load_tables` with `lazy=False`. Do not `read_zarr` for table-only QC. Area lives in `obsm["cell_geometry"]["Area"]` — do not need `include_cell_geometry=True`.

## Quick Start

For complete per-modality workflows, filter diagnostics, and common pitfalls, see [IMPLEMENTATION.md](./IMPLEMENTATION.md).

<!-- canonical: spatialdata-loading-and-access memory guard. Copy verbatim; edit only in that skill's SKILL.md -->
**Memory guard — before any `lazy=False`, `.to_memory()` or `.compute()`:**

1. Load lazily first and run `check_before` from `scripts/check_load.py` in `spatialdata-loading-and-access`. Run the script unchanged in the Python session that holds the data (`exec` its contents); don't rewrite it. It compares each table's estimated in-memory size with the memory this session can use, lists large session objects you could `del`, and raises on a FAIL.
2. Materialize only the tables and runs the plan needs. With several runs, load `concatenated=False` so it checks each run and their total. A WARN says what to `del` or to load one run at a time. A FAIL means subset before materializing: propose a subset to the user (runs, wells, or a set of features), then slice the lazy table and call `.to_memory()` on the slice, for example `adata[adata.obs["Well"].isin(wells)].to_memory()`. `.copy()` keeps a lazy table lazy. Feature subsets by modality: [ACCESS.md — Subset](../spatialdata-loading-and-access/ACCESS.md#subset-by-well-or-feature).
3. Keep `.X` in the form it was stored (`scipy.sparse.issparse(adata.X)`). Sparse tables (3′, OPS) are mostly zeros: never densify a whole one (`.toarray()`, `.todense()`, `np.asarray(adata.X)`); densify only a subset of cells or features. Dense tables (imaging, legacy barcoding) are float32 with `NaN` where a batch has no data: do not convert them to sparse, which turns missing into zero.
<!-- /canonical -->

```python
import elembio_spatialdata_analysis as ebsd
from elembio_spatialdata_analysis.load import Loader

loader = Loader("/path/to/run.zarr")
tables = set(loader.available_tables)
print("tables:", sorted(tables))

# --- Example: 3′ when present (cells → genes → then normalize elsewhere) ---
if "ThreePrimeUntargeted" in tables:
    check_before(loader.load_tables(table="ThreePrimeUntargeted", concatenated=False))  # memory guard: lazy, raises on FAIL
    rna = loader.load_tables(table="ThreePrimeUntargeted", lazy=False)
    ebsd.pp.compute_qc_metrics(rna, "ThreePrimeUntargeted")
    if "filter_flags" not in rna.obsm:          # guard: skip if already filtered
        ebsd.pp.filter_cells(rna, target_type="ThreePrimeUntargeted")  # max_pct_mt=0.10 default
    ebsd.pp.generate_filter_report(rna)
    print("before genes:", rna.shape)
    ebsd.pp.filter_genes(rna, "ThreePrimeUntargeted", min_cells_per_gene=10)
    print("after genes:", rna.shape)
    # Optional tighten (ask first if jumping high): filter_genes(..., min_cells_per_gene=100)

# --- Example: imaging when present (lists "imaging" or "CellPaint") ---
if tables & {"imaging", "CellPaint"}:
    check_before(loader.load_tables(table="imaging", concatenated=False))  # memory guard: lazy, raises on FAIL
    img = loader.load_tables(table="imaging", lazy=False)  # alias OK
    ebsd.pp.compute_qc_metrics(img, "Imaging")
    ebsd.pp.filter_cells(img, target_type="Imaging")
    # flags: obsm["filter_flags"]["passes_imaging_intensity"]
    ebsd.pp.generate_filter_report(img)
    # no filter_genes / no normalize

# --- Example: OPS assay when present (DISS identity library; area-focused) ---
run = next(iter(loader._runs.values()))
specialized = run.table_suffixes_for_target_type(("OPS", "SpecializedTargeted"))
# QC_Panel / housekeeping gene panels are also DISS SpecializedTargeted — often labeled OPS
ops_assay = [t for t in specialized if "qc_panel" not in t.lower()]
targeted_diss = [t for t in specialized if t not in ops_assay]
for name in ops_assay:
    check_before(loader.load_tables(table=name, concatenated=False))  # memory guard: lazy, raises on FAIL
    ops = loader.load_tables(table=name, lazy=False, include_target_assignment=True)
    ebsd.pp.compute_qc_metrics(ops, "OPS")
    ebsd.pp.filter_cells(ops, target_type="OPS")  # area defaults; counts off
    # Confidence identity filter → ops-pipeline (OPS assay only)
for name in targeted_diss:
    # Still DISS — depth/detection only; do NOT Confidence-gate; see run-quality-overview
    pass
```


Multi-modality intersection: `ebsd.pp.filter_modalities` — only include **present** modalities; area owner is Transcript or ThreePrimeUntargeted when either is present. After intersection, still run `filter_genes` on RNA tables. See IMPLEMENTATION.md.

## Detail Files

- [IMPLEMENTATION.md](./IMPLEMENTATION.md) — workflows, `filter_modalities`, 3′, imaging, OPS, reports, pitfalls
- [REFERENCE.md](./REFERENCE.md) — defaults, `qc_metrics` columns, filter flags

## Related Skills

- [run-quality-overview](../run-quality-overview/SKILL.md) — prior step (run/well)
- [spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md)
- [normalization-strategies](../normalization-strategies/SKILL.md) — after cell QC (not for imaging / not for OPS)
- [ops-pipeline](../ops-pipeline/SKILL.md) — OPS confidence filtering / merge workflows
- [3prime-untargeted-pipeline](../3prime-untargeted-pipeline/SKILL.md)
- [imaging-pipeline](../imaging-pipeline/SKILL.md) — morphology / IF orchestration after cell QC

## Disclaimer

This content is provided "as-is" without warranty of any kind. AI-generated output may contain
inaccuracies. Review and verify all outputs before relying on them.
