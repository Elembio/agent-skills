---
name: run-quality-overview
description: Discovery-driven run-level QC for AVITI24. Summarize only present modalities. Classify SpecializedTargeted DISS tables into OPS assay vs targeted QC/gene panels (e.g. QC_Panel) — Confidence ≥ 0.5 only for OPS assay; QC/gene panels get depth/detection only. Soft-interpret; does not filter cells. Use for run quality, cells/well, segmentation, assignment — first QC after loading. Never invent absent modalities.
metadata:
  version: 0.4.11
  author: elembio
---

# Run Quality Overview

## Purpose

**First QC stop after loading:** discover what the run carries, summarize run/well health for **only those modalities**, soft-interpret outliers, then hand off to the matching specialist skills.

```text
Load → discover present modalities → run-quality-overview (present only)
  → cell-quality-control → pipeline(s) for what is present
```

This skill **never** calls `filter_cells` or drops cells. It summarizes and flags.

**Hard rule — do not invent modalities.** After discovery, state the present set once (e.g. `OPS + imaging`). Summarize **only** that set. Do **not** mention, apologize for, speculate about, or give “would look like” guidance for tables/types that are absent (no 3′ depth talk on an OPS+imaging-only run, no OPS talk on a 3′-only run, etc.).

## Prerequisites

- [spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md)
- A local or S3 AVITI24 `.zarr` / `.zarr.zip`

## Key Concepts

### 1. Discover first — modalities are not universal

```python
from elembio_spatialdata_analysis.load import Loader

loader = Loader("/path/to/run.zarr")  # local path or s3:// URI
run = next(iter(loader._runs.values()))
print(sorted(loader.available_tables))
print(sorted(loader.available_target_types))
```

Build an explicit present set, then run **only** matching branches.

**Critical — `SpecializedTargeted` / table `target_type=OPS` is not always Optical Pooled Screening.** These tables are **DISS** (`DISSBatch`). Table metadata uses `target_type=OPS` for several kinds of specialized DISS panel, including **targeted QC / housekeeping gene panels** (e.g. `B01-QC_Panel` = GAPDH+ACTB probes). Those are still DISS — just **not** OPS identity libraries. Classify before any Confidence QC:

```python
tables = set(loader.available_tables)
types = set(loader.available_target_types)
specialized = run.table_suffixes_for_target_type(("OPS", "SpecializedTargeted"))

def classify_specialized_table(name: str, adata=None) -> str:
    """Return 'ops_assay' (OPS identity library) or 'targeted_diss_qc' (QC/gene panel)."""
    n = name.lower()
    if "qc_panel" in n or n.endswith("-qc") or "control_panel" in n:
        return "targeted_diss_qc"
    if adata is not None and "target_gene" in adata.var.columns:
        genes = {str(g).upper() for g in adata.var["target_gene"].dropna()}
        hk = {"GAPDH", "ACTB", "ACTIN", "B2M", "HPRT1"}
        if genes and (genes <= hk or (adata.n_vars <= 20 and genes & hk)):
            return "targeted_diss_qc"
    return "ops_assay"

ops_assay_tables, targeted_diss_tables = [], []
for name in specialized:
    kind = classify_specialized_table(name)
    if kind == "ops_assay":
        # Confirm with var when the name alone is ambiguous
        ad = loader.load_tables(table=name, lazy=True)
        kind = classify_specialized_table(name, ad)
    (targeted_diss_tables if kind == "targeted_diss_qc" else ops_assay_tables).append(name)

present = []
if "ThreePrimeUntargeted" in tables:
    present.append("ThreePrimeUntargeted")
if ops_assay_tables:
    present.append("OPS")  # OPS identity library only
if targeted_diss_tables:
    present.append("targeted_DISS")  # QC / gene panels — still DISS, not ops-pipeline
if tables & {"imaging", "CellPaint"} or "imaging" in types:
    present.append("imaging")
if "Transcript" in tables:
    present.append("Transcript")  # legacy BarcodingBatch
if "Protein" in tables:
    present.append("Protein")  # legacy BarcodingBatch
print("Present:", present)
print("OPS assay tables:", ops_assay_tables)
print("Targeted DISS QC/gene panels:", targeted_diss_tables)
```

Imaging may list `"imaging"` or legacy `"CellPaint"` — both mean imaging (`load_tables(table="imaging")` aliases either way).

### 2. Universal layer (always — not a modality)

Per-well cytometry from `global-table-metadata.uns` (`c2s_runstats_cyto_wells` → fallback `runstats_cyto_wells`):

| Column | Look for |
|---|---|
| `CellCount` | Dead / empty / overloaded wells (compare wells, not a fixed floor) |
| `PercentNucleatedCells` | Low nucleation → debris / nuclear seg issues |
| `PercentConfluency` | Very low sparse; very high over-confluent |
| `MedianCellDiameter` | Drift across wells → focus / seg gradient |

Cyto wells have `WellLocation` (`A1`) and no sample-sheet name. Backfill `WellLabel` from `run.wells` (column `well_label`). Tile-level dropouts: `qc.per_tile_cell_count(adata)` on any cell-shaped table.

**Spatial artifacts:** Check `obsm["spatial_artifacts"]` on any cell-aligned table. If present and cells are flagged in any rollup (`any-imaging`, `any-ThreePrimeUntargeted`, `any-OPS`), report the counts per modality and recommend `filter_spatial_artifacts=True` in [cell-quality-control](../cell-quality-control/SKILL.md). If absent, say nothing.

### 3. Modality branches — skip entirely if not in `present`

| Only if in `present` | Summarize |
|---|---|
| **ThreePrimeUntargeted** | Depth per well (`obsm["total_counts"]` or `X` row sums); genes detected; optional `f2p_3p_runstats_*` uns |
| **OPS** (`ops_assay_tables` only) | `include_target_assignment=True`; **`% Confidence>=0.5`**, undetermined, median `TotalCount`, **by well**. Hand off to [ops-pipeline](../ops-pipeline/SKILL.md) |
| **targeted_DISS** (`targeted_diss_tables`) | Still DISS — **do not use Confidence as a health metric.** Report per-well depth (`TotalCount` or `X` row sums), genes/probes detected, optional gene totals (e.g. GAPDH vs ACTB). Low Confidence on QC_Panel is expected and **not** a failure |
| **Imaging** | Feature inventory; area per well; optional soft S/B — [imaging-pipeline](../imaging-pipeline/SKILL.md) |
| **Legacy Transcript / Protein** | Depth rollups (`BarcodingBatch`, not DISS) |

**Hard guardrail:** Never headline run QC with `Confidence>=0.5` for `QC_Panel` / housekeeping gene panels. That metric is for **OPS assay** (identity library) tables only — both kinds are DISS.

**OPS health (assay only):** cyto `MeanAssignedBarcodingCountPerCell` is a barcoding-assignment metric, not an OPS health metric. Use AnnData `Confidence` / `TotalCount` for assay OPS.

### 4. Soft interpret (present signals only)

Compare wells. Call out outliers **for metrics you actually computed**. Do not cite OPS or 3′ gate bands unless that modality is in `present`.

Universal bands (always OK to use when cyto exists): nucleation ≪ ~85%, confluency outside ~20–70, cell counts orders of magnitude below peer wells.

### 5. Hand off (only for present modalities)

- Always next for filtering: [cell-quality-control](../cell-quality-control/SKILL.md)
- If **OPS** in `present` (assay tables only) → [ops-pipeline](../ops-pipeline/SKILL.md)
- If **targeted_DISS** only (no assay OPS) → **do not** open ops-pipeline; keep as run-QC / depth summary alongside 3′ / imaging
- If imaging in `present` → [cell-quality-control](../cell-quality-control/SKILL.md) then [imaging-pipeline](../imaging-pipeline/SKILL.md)
- If ThreePrimeUntargeted in `present` → [3prime-untargeted-pipeline](../3prime-untargeted-pipeline/SKILL.md)
- Do **not** recommend pipelines for absent modalities.

## Quick Start

<!-- canonical: spatialdata-loading-and-access memory guard. Copy verbatim; edit only in that skill's SKILL.md -->
**Memory guard — before any `lazy=False`, `.to_memory()` or `.compute()`:**

1. Load lazily first and run `check_before` from `scripts/check_load.py` in `spatialdata-loading-and-access`. Run the script unchanged in the Python session that holds the data (`exec` its contents); don't rewrite it. It compares each table's estimated in-memory size with the memory this session can use, lists large session objects you could `del`, and raises on a FAIL.
2. Materialize only the tables and runs the plan needs. With several runs, load `concatenated=False` so it checks each run and their total. A WARN says what to `del` or to load one run at a time. A FAIL means subset before materializing: propose a subset to the user (runs, wells, or a set of features), then slice the lazy table and call `.to_memory()` on the slice, for example `adata[adata.obs["Well"].isin(wells)].to_memory()`. `.copy()` keeps a lazy table lazy. Feature subsets by modality: [ACCESS.md — Subset](../spatialdata-loading-and-access/ACCESS.md#subset-by-well-or-feature).
3. Keep `.X` in the form it was stored (`scipy.sparse.issparse(adata.X)`). Sparse tables (3′, OPS) are mostly zeros: never densify a whole one (`.toarray()`, `.todense()`, `np.asarray(adata.X)`); densify only a subset of cells or features. Dense tables (imaging, legacy barcoding) are float32 with `NaN` where a batch has no data: do not convert them to sparse, which turns missing into zero.
<!-- /canonical -->

```python
import pandas as pd
from elembio_spatialdata_analysis.load import Loader
from elembio_spatialdata_analysis.qc import per_tile_cell_count

STORE = "/path/to/run.zarr"  # local path or s3:// URI
loader = Loader(STORE)
run = next(iter(loader._runs.values()))
tables, types = set(loader.available_tables), set(loader.available_target_types)
specialized = run.table_suffixes_for_target_type(("OPS", "SpecializedTargeted"))
targeted_diss_tables = [t for t in specialized if "qc_panel" in t.lower() or t.lower().endswith("-qc")]
ops_assay_tables = [t for t in specialized if t not in targeted_diss_tables]
present = [m for m, ok in [
    ("ThreePrimeUntargeted", "ThreePrimeUntargeted" in tables),
    ("OPS", bool(ops_assay_tables)),
    ("targeted_DISS", bool(targeted_diss_tables)),
    ("imaging", bool(tables & {"imaging", "CellPaint"}) or "imaging" in types),
    ("Transcript", "Transcript" in tables),
    ("Protein", "Protein" in tables),
] if ok]
print("Present:", present, "OPS assay:", ops_assay_tables, "targeted DISS:", targeted_diss_tables)

_CYTO = ("c2s_runstats_cyto_wells", "runstats_cyto_wells")

def uns_frame(run, keys):
    uns = run._metadata
    for k in keys:
        v = uns.get(k)
        if isinstance(v, pd.DataFrame) and not v.empty:
            return v.reset_index() if v.index.name else v.copy()
    return None

# --- Universal (always) ---
seg = uns_frame(run, _CYTO)
if seg is not None:
    src = run.wells
    label_col = "well_label" if "well_label" in src.columns else "WellLabel"
    if label_col in src.columns:
        mapped = seg["WellLocation"].astype(str).map(src[label_col].astype(str).to_dict())
        if "WellLabel" in seg.columns:
            existing = seg["WellLabel"]
            seg = seg.copy()
            seg["WellLabel"] = existing.where(
                existing.notna() & (existing.astype(str).str.strip() != ""),
                mapped,
            )
        else:
            seg = seg.copy()
            seg.insert(seg.columns.get_loc("WellLocation") + 1, "WellLabel", mapped)
print(seg)

# --- Branches gated on present (skip the rest entirely) ---
if "OPS" in present:  # assay barcode panels only
    for name in ops_assay_tables:
        check_before(loader.load_tables(table=name, concatenated=False))  # memory guard: lazy, raises on FAIL
        ops = loader.load_tables(table=name, lazy=False, include_target_assignment=True)
        conf = ops.obs["Confidence"]
        print(f"OPS assay {name}: Conf>=0.5 = {(conf >= 0.5).mean():.3f}")
        print(ops.obs.assign(a=conf >= 0.5).groupby("Well", observed=True)["a"].mean())

if "targeted_DISS" in present:  # QC / housekeeping — NO Confidence health gate
    for name in targeted_diss_tables:
        check_before(loader.load_tables(table=name, concatenated=False))  # memory guard: lazy, raises on FAIL
        panel = loader.load_tables(table=name, lazy=False, include_target_assignment=True)
        depth = panel.obs["TotalCount"] if "TotalCount" in panel.obs else panel.X.sum(axis=1)
        print(f"Targeted DISS/QC {name}: n_vars={panel.n_vars}", panel.var_names.tolist()[:12])
        print(panel.obs.assign(depth=depth).groupby("Well", observed=True)["depth"].median())
        # Do NOT print or gate on Confidence>=0.5 for these panels

if "ThreePrimeUntargeted" in present:
    check_before(loader.load_tables(table="ThreePrimeUntargeted", concatenated=False))  # memory guard: lazy, raises on FAIL
    rna = loader.load_tables(table="ThreePrimeUntargeted", lazy=False)
    tc = rna.obsm["total_counts"]
    depth = tc.sum(axis=1) if hasattr(tc, "sum") else tc
    print(rna.obs.assign(depth=depth).groupby("Well", observed=True)["depth"].median())

if "imaging" in present:
    check_before(loader.load_tables(table="imaging", concatenated=False))  # memory guard: lazy, raises on FAIL
    img = loader.load_tables(table="imaging", lazy=False)
    if "feature_class" in img.var.columns:
        print(img.var["feature_class"].value_counts(dropna=False))
    print(img.var.get("batch_name", pd.Series(dtype=object)).value_counts(dropna=False))
    # Soft S/B on morphology stains only — full recipe in imaging-pipeline (exclude Bkg/IgG/Cleave)
```

Full recipes: [IMPLEMENTATION.md](./IMPLEMENTATION.md). Reference: [REFERENCE.md](./REFERENCE.md).

## Related Skills

- [spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md)
- [cell-quality-control](../cell-quality-control/SKILL.md) — next: filter cells
- Specialist pipelines — only if that modality is in `present`: [ops-pipeline](../ops-pipeline/SKILL.md), [3prime-untargeted-pipeline](../3prime-untargeted-pipeline/SKILL.md), [imaging-pipeline](../imaging-pipeline/SKILL.md)
- [metadata-and-target-mapping](../metadata-and-target-mapping/SKILL.md), [experiment-design](../experiment-design/SKILL.md)

## Disclaimer

This content is provided "as-is" without warranty of any kind. AI-generated output may contain
inaccuracies. Review and verify all outputs before relying on them.
