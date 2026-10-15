---
name: ops-pipeline
description: End-to-end workflow for true Optical Pooled Screening (OPS) — DISS SpecializedTargeted identity libraries. Confidence QC, merge cell identity into RNA, filter, UMAP/DE by line×condition. Not for other SpecializedTargeted DISS panels such as QC/housekeeping gene panels (e.g. B01-QC_Panel GAPDH/ACTB) even when metadata says target_type=OPS. Classify first; skip QC_Panel. Typically paired with ThreePrimeUntargeted.
metadata:
  version: 0.2.7
  author: elembio
---

# OPS Pipeline

> **Confirm the SDO first.** Before running, make sure you are using the SDO the user selected in [locate-spatialdata-store](../locate-spatialdata-store/SKILL.md) — a run may have several (instrument-direct, `cells2stats`, `elembio-multiomics-nf`).

## Purpose

End-to-end workflow for **true Optical Pooled Screening (OPS)** — a **DISS** (`DISSBatch` / `SpecializedTargeted`) identity library that assigns cells to a library member (cell line, perturbation, etc.). Almost always paired with `ThreePrimeUntargeted` (also DISS): OPS = identity, RNA = transcriptome. This is **not** legacy `BarcodingBatch` chemistry.

## Prerequisites

- At least one **OPS assay** table (identity / cell-assignment library) — not merely any `SpecializedTargeted` / `target_type=OPS` table
- Python **3.12+** with `elembio_spatialdata_analysis`
- Typically also `ThreePrimeUntargeted`

**Loading:** `load_tables(table=<ops_table_name>, include_target_assignment=True)`. Never `table="OPS"`.

### Guardrail — not every SpecializedTargeted DISS table is OPS

Both OPS and QC/gene panels are **DISS**. The build often labels **targeted QC / housekeeping gene panels** as `target_type=OPS` (e.g. `B01-QC_Panel` with GAPDH+ACTB probes). Those are SpecializedTargeted DISS tables, **not** Optical Pooled Screening identity libraries — do **not** run this skill on them.

| Kind | Examples | Use this skill? | Run-QC metric |
|---|---|---|---|
| **OPS assay** | OPS identity libraries (e.g. `B01-<LibraryName>`) | **Yes**: Confidence filters, merge identity into RNA | `% Confidence>=0.5` |
| **Targeted DISS QC / gene panel** | `QC_Panel`, housekeeping gene panels | **No** — stop; use [run-quality-overview](../run-quality-overview/SKILL.md) depth summaries | depth / detection — **not** Confidence |

Before Step 1, classify discovered suffixes. If **only** QC/gene-panel tables exist, **exit** and do not apply Confidence gates.

## Key Concepts

### OPS table naming — the critical difference from other modalities

OPS table names follow the pattern `{batch_name}-{target_site_name}` where `target_site_name` is a **custom name set in the run manifest** (e.g. `"B01-CustomPanel"`). This means:
- You **cannot** load it with `table="OPS"`
- The name varies between runs depending on the panel configured
- Always discover via `run.table_suffixes_for_target_type(("OPS", "SpecializedTargeted"))`

**How to discover the OPS table name — preferred methods:**

```python
from elembio_spatialdata_analysis.load import Loader

loader = Loader("/path/to/run.zarr")
run = next(iter(loader._runs.values()))

# Preferred: table suffixes (handles OPS + SpecializedTargeted aliases)
specialized = run.table_suffixes_for_target_type(("OPS", "SpecializedTargeted"))
# Drop QC / gene panels — still DISS SpecializedTargeted, not OPS assay
ops_tables = [t for t in specialized if "qc_panel" not in t.lower()]
qc_or_gene_panels = [t for t in specialized if t not in ops_tables]
print("OPS assay tables:", ops_tables)
print("Skip (targeted DISS QC/gene panels):", qc_or_gene_panels)
if not ops_tables:
    raise SystemExit("No OPS assay table — do not run ops-pipeline Confidence workflow")
```

### OPS AnnData structure

| Slot | Contents |
|------|----------|
| `.X` | Sparse int32 DISS polony counts — cells × OPS library members |
| `.obs` | Standard cell columns: `Cell`, `Well`, `WellLabel`, `tile_id`, `X_tile`, `Y_tile`, `cell_id_global`, `well_id` |
| `.var` | One row per library member: `target_id`, `target_gene`, `target_sequence`, `control_type`, `batch_name`, `target_site_name` |
| `obsm["target_assignment"]` | DataFrame with 5 columns — the per-cell identity call; available with `include_target_assignment=True` |
| `obsm["total_counts"]` | DataFrame — per-batch polony counts per cell |
| `obsm["spatial"]` | Region-pixel centroids `X`, `Y` (DataFrame) |
| `obsm["cell_geometry"]` | Cell morphology (`Area`, `NuclearArea`, `Xum_tile`, `Yum_tile`, shape features) |
| `uns["metadata_attrs"]` | `batch_class=DISSBatch`, `target_type=OPS`, `target_site_name=<custom_name>` |

### Assignment fields (`obsm["target_assignment"]` / `obs` after loading)

| Column | Type | Meaning |
|--------|------|---------|
| `MostAbundantTarget` | str | Winning barcode identity (e.g. cell line name), or `"Undetermined"` |
| `MostAbundantTargetCount` | int | Read count for the winning barcode |
| `TotalCount` | int | Total polony reads across all barcodes for this cell |
| `Confidence` | float [0–1] | `MostAbundantTargetCount / (TotalCount + 1)` (+1 pseudocount; 0 when no reads) |
| `target_gene` | str | Gene-level annotation for the winning target |

### Assignment rate QC

| Metric | Typical value | Concern threshold |
|--------|--------------|-------------------|
| `MostAbundantTarget != "Undetermined"` | 55–70% of cells | < 40% warrants investigation |
| `Confidence >= 0.5` | 50–65% of cells | < 35% warrants investigation |
| Median `TotalCount` | 1–5 reads/cell | < 1 suggests loading or sequencing problem |

**Standard filter:** `Confidence >= 0.5` (winning count over total+1; roughly “half or more” of reads once depth is non-trivial).

### Merge OPS identity into RNA **before** any filtering

OPS and RNA tables share identical cell ordering (same `obs_names`, same `cell_id_global`). The correct sequence is:

1. Load OPS table
2. Load RNA table
3. Assert alignment (obs_names, cell_id_global, `obsm["spatial"]` match)
4. Copy OPS fields into `rna.obs` — including the `ops_keep` flag
5. Run RNA QC (`filter_cells`)
6. Apply the `ops_keep` confidence filter as a second pass

This order preserves the ability to report per-well assignment rates before any RNA filtering changes the denominator.

### Pipeline Order

```
Step 0  Setup (imports, thresholds)
Step 1  Loader + discover OPS table name (via run.available_target_types or run.batches)
Step 2  Load OPS table (include_target_assignment=True)
Step 3  Per-well assignment statistics
Step 4  Load ThreePrimeUntargeted RNA (legacy: Transcript)
Step 5  Alignment asserts (obs_names, cell_id_global, spatial)
Step 6  Merge OPS → RNA obs (OPS_line, OPS_Confidence, OPS_TotalCount, ops_keep)
Step 7  RNA QC: compute_qc_metrics + filter_cells
Step 8  Intersect with ops_keep (Confidence >= 0.5)
Step 9  Normalize + PCA + UMAP
Step 10 Visualize by condition and OPS_line
Step 11 Pseudobulk DE by cell line × condition
Step 12 Save results
```

### Batch metadata for OPS

```python
# OPS barcode panel metadata (sequences, CellAssignment flag, primer)
diss_meta = run.diss_targets_metadata
print(diss_meta.columns.tolist())
# e.g. ['TargetSiteName', 'TargetMask', 'CellAssignment', 'RTOligos', 'SequencingPrimer', 'batch_name']

# OPS barcodes (CellAssignment=True)
ops_targets = diss_meta[diss_meta["CellAssignment"] == True]

# OPS barcode rows (by CellAssignment) or filter by SpecializedTargeted batch names
ops_batch_names = run.batches.loc[
    run.batches["batch_type"] == "SpecializedTargeted", "batch_name"
].tolist()
ops_diss = diss_meta[diss_meta["batch_name"].isin(ops_batch_names)] if "batch_name" in diss_meta.columns else diss_meta
```

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
import scipy.sparse as sp
import numpy as np

loader = Loader("/path/to/run.zarr")
print(loader.summary())
print("Available tables:", sorted(loader.available_tables))

# Discover OPS table name
run = next(iter(loader._runs.values()))
print("Target types:", sorted(run.available_target_types))
specialized = run.table_suffixes_for_target_type(("OPS", "SpecializedTargeted"))
ops_tables = [t for t in specialized if "qc_panel" not in t.lower()]
if not ops_tables:
    raise ValueError(
        "No OPS assay table found "
        f"(specialized={specialized}). Targeted DISS QC/gene panels are not ops-pipeline."
    )
ops_table_name = ops_tables[0]
print(f"OPS assay table: {ops_table_name!r}")

# Load OPS with assignment columns
check_before(loader.load_tables(table=ops_table_name, concatenated=False))  # memory guard: lazy, raises on FAIL
ops = loader.load_tables(table=ops_table_name, include_target_assignment=True, lazy=False)
print(ops.obsm["target_assignment"].head())

# Per-well assignment rate
conf = ops.obs["Confidence"] >= 0.5
print(f"Assigned (Conf>=0.5): {conf.sum():,} / {ops.n_obs:,} ({100*conf.mean():.1f}%)")

# Load RNA and merge
check_before(loader.load_tables(table="ThreePrimeUntargeted", concatenated=False))  # memory guard: lazy, raises on FAIL
rna = loader.load_tables(table="ThreePrimeUntargeted", lazy=False)
assert ops.obs_names.equals(rna.obs_names), "obs_names mismatch!"
rna.obs["OPS_line"]       = ops.obs["MostAbundantTarget"].values
rna.obs["OPS_Confidence"] = ops.obs["Confidence"].values
rna.obs["OPS_TotalCount"] = ops.obs["TotalCount"].values
rna.obs["ops_keep"]       = (ops.obs["Confidence"] >= 0.5).values

# RNA QC then OPS filter
ebsd.pp.compute_qc_metrics(rna, "ThreePrimeUntargeted")
ebsd.pp.filter_cells(rna, "ThreePrimeUntargeted")
rna = rna[rna.obs["ops_keep"]].copy()
print(f"After filters: {rna.n_obs:,} cells")
```

## Detail Files

- [IMPLEMENTATION.md](./IMPLEMENTATION.md) — Full step-by-step code for all 12 steps, per-well stats, normalization, UMAP, pseudobulk DE
- [REFERENCE.md](./REFERENCE.md) — AnnData slot inventory, diss_targets_metadata schema, table-name discovery patterns, confidence thresholds, batches DataFrame columns

## Related Skills

- [spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md) — `Loader`, backends, `load_tables`, OPS table loading
- [metadata-and-target-mapping](../metadata-and-target-mapping/SKILL.md) — `RunInfo`, `diss_targets_metadata`, OPS batch discovery
- [3prime-untargeted-pipeline](../3prime-untargeted-pipeline/SKILL.md) — RNA QC, HVG, normalization for the `ThreePrimeUntargeted` table
- [cell-quality-control](../cell-quality-control/SKILL.md) — `compute_qc_metrics`, `filter_cells` API
- [normalization-strategies](../normalization-strategies/SKILL.md) — `pp.normalize` for 3′ data
- [embedding-visualization](../embedding-visualization/SKILL.md) — UMAP coloring by `OPS_line` and `condition`
- [experiment-design](../experiment-design/SKILL.md) — parse `obs["WellLabel"]` (sample-sheet name); `Well` is the plate position
- [pseudobulk-differential-expression](../pseudobulk-differential-expression/SKILL.md) — DESeq2 pseudobulk by cell line × condition

## Disclaimer

This content is provided "as-is" without warranty of any kind. AI-generated output may contain
inaccuracies. Review and verify all outputs before relying on them.
