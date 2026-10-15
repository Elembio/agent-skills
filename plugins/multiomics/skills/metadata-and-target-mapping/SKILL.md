---
name: metadata-and-target-mapping
description: Wells, tiles, batches, and targets from Loader RunInfo / global-table-metadata. Lists Well (plate position) and WellLabel (sample-sheet name as recorded). DISS-first — diss_targets_metadata (3′ + OPS), targets_imaging (CP-group + optional IF), legacy targets_barcoding. OPS table discovery via table_suffixes_for_target_type(("OPS", "SpecializedTargeted")). Use for run layout, panel inventory, or target mapping — not for inferring condition / replicate groupings.
metadata:
  version: 0.2.4
  author: elembio
---

# Metadata and Target Mapping

## Purpose

Parse **run-level** metadata for AVITI24 stores: wells, tiles, batches, and per-target records via **`Loader` → `RunInfo`** (preferred) or `global-table-metadata.uns`. Current product metadata is **DISS + imaging**; barcoding tables are legacy.

## Prerequisites

- [spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md)
- Store contains **`global-table-metadata`** (recommended)

## Key Concepts

### Where metadata lives

| Source | Access |
|--------|--------|
| **`Loader` / `RunInfo`** (preferred) | `run.wells`, `run.tiles`, `run.batches`, `run.targets_imaging`, `run.diss_targets_metadata`, `run.targets_barcoding` |
| **`global-table-metadata.uns`** | Same keys under `sdata.tables["global-table-metadata"].uns[...]` |

**`loader._runs`** is a private attribute (`dict[str, RunInfo]`) with no public accessor. Use `run = next(iter(loader._runs.values()))` for per-run metadata.

There is **no** single `uns["targets"]` — targets are split by assay family.

### Target metadata tables (current vs legacy)

| Key | Role | Product |
|-----|------|---------|
| **`diss_targets_metadata`** | DISS targets — OPS barcodes (`CellAssignment=True`) and 3′ RT oligos (`CellAssignment=False`) | **Current** when 3′ / OPS present |
| **`targets_imaging`** | Imaging **targets** (channels) for the imaging table — primary morphology group (`CP*`) and additional IF (`P*` / ProteinIF) when present | **Current** when imaging present |
| **`targets_barcoding`** | Legacy Transcript / Protein barcoding panel definitions | **Legacy SDOs only** |
| **`wells` / `tiles` / `batches`** | Layout and batch inventory | Always useful |
| **`build_provenance`**, **`runstats_*`** | Build info / optional stats | Optional |

`targets_imaging` lists **targets**, not a separate “CellPaint modality.” CellPaint vs ProteinIF is target/batch role inside imaging (see [imaging-pipeline](../imaging-pipeline/SKILL.md)).

### `available_target_types` vs `available_tables`

| API | Returns |
|-----|---------|
| `loader.available_tables` | Table suffixes to load (`ThreePrimeUntargeted`, `imaging` or `CellPaint`, `B01-…`, …) |
| `loader.available_target_types` / `run.available_target_types` | `target_type` values from table `metadata_attrs` (e.g. `imaging`, `OPS` / `SpecializedTargeted`, `ThreePrimeUntargeted`) |

Never invent modalities. Discover first. Example current set: `{'ThreePrimeUntargeted', 'OPS', 'imaging'}` — not universal.

### Batches (typical)

| `batch_class` | Typical `batch_type` | Contents |
|---------------|----------------------|----------|
| **DISSBatch** | `ThreePrimeUntargeted`, `SpecializedTargeted` | 3′ / OPS |
| **ImagingBatch** | PreAmp, PostAmp, StainOnly, … | Imaging targets (`CP*` morphology group; `P*` IF when present) |
| **BarcodingBatch** | — | Legacy Transcript / Protein |

## Quick Start

```python
from elembio_spatialdata_analysis.load import Loader

loader = Loader("/path/to/run.zarr")
run = next(iter(loader._runs.values()))

print("Tables:", sorted(loader.available_tables))
print("Types:", sorted(run.available_target_types))

wells_df = run.wells
tiles_df = run.tiles
batches_df = run.batches
diss_meta = run.diss_targets_metadata
targets_imaging = run.targets_imaging
targets_barcoding = run.targets_barcoding  # often empty on current runs

# Labels as recorded — Well is A1; WellLabel is the sample-sheet string
print(adata.obs[["Well", "WellLabel"]].astype(str).drop_duplicates().sort_values("Well").to_string(index=False))
# No table yet: print(wells_df[["well_label"]])
```

### SpecializedTargeted / “OPS” tables — classify before analysis

Prefer **`table_suffixes_for_target_type`** (handles `"OPS"` and `"SpecializedTargeted"`). Do **not** load `table="OPS"`.

**Guardrail:** metadata `target_type=OPS` also appears on **targeted DISS QC/gene panels** (e.g. `B01-QC_Panel` with GAPDH/ACTB). Those are still DISS SpecializedTargeted, but **not** Optical Pooled Screening identity libraries — do not send them to [ops-pipeline](../ops-pipeline/SKILL.md) Confidence workflows. See [run-quality-overview](../run-quality-overview/SKILL.md).

```python
specialized = run.table_suffixes_for_target_type(("OPS", "SpecializedTargeted"))
# e.g. ["B01-<OPSPanel>", "B01-QC_Panel"]: both DISS; split OPS assay vs QC/gene panel
ops_assay = [t for t in specialized if "qc_panel" not in t.lower()]
targeted_qc = [t for t in specialized if t not in ops_assay]

# Optional: batch names only
ops_batches = run.batches_for_target_type("OPS")  # may be incomplete vs SpecializedTargeted alias
# or: run.batches[run.batches["batch_type"] == "SpecializedTargeted"]

if len(diss_meta):
    ops_targets = diss_meta[diss_meta["CellAssignment"] == True]
    rna_oligos = diss_meta[diss_meta["CellAssignment"] == False]
    print("OPS barcodes:", len(ops_targets), "3′ oligos:", len(rna_oligos))
```

### Imaging targets (CP-group vs additional IF)

```python
ti = run.targets_imaging
if len(ti):
    print(ti.columns.tolist())
    for col in ("batch_name", "target_type", "target_name"):
        if col in ti.columns:
            print(ti[col].value_counts(dropna=False).head(20))
# CP* ≈ primary morphology group (may include IF reagents, e.g. histone)
# P* / ProteinIF ≈ additional IF markers when present
```

### Well labels (as recorded)

Print `Well` + `WellLabel` and stop. Do not split those strings into condition / replicate here — that is [experiment-design](../experiment-design/SKILL.md) (needs user confirmation). `run.wells["well_label"]` is the same strings when no cell table is loaded.

## Detail Files

- [IMPLEMENTATION.md](./IMPLEMENTATION.md) — wells/tiles/batches, target tables, DISS/OPS, label listing, pitfalls

## Related Skills

- [spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md) — discovery + `load_tables`
- [experiment-design](../experiment-design/SKILL.md) — condition / replicate from `WellLabel`
- [ops-pipeline](../ops-pipeline/SKILL.md) — OPS end-to-end
- [imaging-pipeline](../imaging-pipeline/SKILL.md) — imaging targets → features / roles
- [3prime-untargeted-pipeline](../3prime-untargeted-pipeline/SKILL.md)
- [anndata-structure-reference](../anndata-structure-reference/SKILL.md)

## Disclaimer

This content is provided "as-is" without warranty of any kind. AI-generated output may contain
inaccuracies. Review and verify all outputs before relying on them.
