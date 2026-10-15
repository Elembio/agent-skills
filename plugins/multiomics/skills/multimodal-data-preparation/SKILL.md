---
name: multimodal-data-preparation
description: Aligning present AVITI24 modalities for joint analysis — 3′ / OPS / imaging (CellPaint alias) and legacy Transcript/Protein. Intersection vs union, per-modality QC and transforms (count normalize vs imaging MAD-z), MuData. Use when combining modalities. Defers joint embedding stats to integrated-analysis; imaging details to imaging-pipeline.
metadata:
  version: 0.2.3
  author: elembio
---

# Multimodal Data Preparation

## Purpose

Align **modalities that exist in the store** for joint analysis. Typical tables: **ThreePrimeUntargeted**, **OPS** (custom table names), **imaging** (`CellPaint` alias). Legacy stores: `Transcript` / `Protein` tables. Discover first; never invent absent modalities.

**IO:** `Loader.load_tables` per modality. Images/labels for viz → `spatialdata.read_zarr` ([spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md)).

## Prerequisites

- [spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md)
- [run-quality-overview](../run-quality-overview/SKILL.md) / [cell-quality-control](../cell-quality-control/SKILL.md)
- Count transforms: [normalization-strategies](../normalization-strategies/SKILL.md)
- Imaging / IF: [imaging-pipeline](../imaging-pipeline/SKILL.md)

## Key Concepts

### Modality overview (present-only)

| Modality | Table | Transform before joint DR |
|---|---|---|
| **ThreePrimeUntargeted** | `ThreePrimeUntargeted` | `filter_genes` → `ebsd.pp.normalize` |
| **Imaging** | `imaging` / `CellPaint` | Cell QC → MAD-z feature matrix (**not** `pp.normalize`); exclude Bkg/IgG/Cleave from embeddings |
| **OPS** | custom `{batch}-{panel}` | Area/confidence filters; **do not** count-normalize; often merge identity into RNA `obs` |
| **Transcript / Protein** (legacy) | `Transcript` / `Protein` | `ebsd.pp.normalize` |

### Intersection vs union

| Strategy | When |
|---|---|
| **Intersection** | Default for joint embedding — cells in all selected modalities |
| **Union** | Keep cells missing a modality (NaNs) — only if methods tolerate it |

### Key rule: transform per modality after alignment

Normalize / MAD-z **each modality independently after** aligning cell IDs. Do not run count `normalize` on imaging. Do not put imaging Bkg/IgG/Cleave into a joint feature matrix used for UMAP.

---

## Quick Start

```python
from elembio_spatialdata_analysis.load import Loader

loader = Loader("/path/to/run.zarr")
tables = set(loader.available_tables)
print("Tables:", sorted(tables))

modalities = {}
if "ThreePrimeUntargeted" in tables:
    modalities["ThreePrimeUntargeted"] = loader.load_tables(table="ThreePrimeUntargeted", lazy=False)
if tables & {"imaging", "CellPaint"}:
    modalities["Imaging"] = loader.load_tables(table="imaging", lazy=False)
run = next(iter(loader._runs.values()))
specialized = run.table_suffixes_for_target_type(("OPS", "SpecializedTargeted"))
# OPS assay only (DISS identity library) — skip QC_Panel / housekeeping gene panels
ops_assay = [t for t in specialized if "qc_panel" not in t.lower()]
for name in ops_assay:
    modalities[name] = loader.load_tables(table=name, lazy=False, include_target_assignment=True)

cell_sets = {n: set(a.obs_names.astype(str)) for n, a in modalities.items()}
common = set.intersection(*cell_sets.values()) if cell_sets else set()
print(f"Common cells: {len(common)}")
```

Then QC each present modality; `ebsd.pp.normalize` only counts; imaging → [imaging-pipeline](../imaging-pipeline/SKILL.md). Joint embedding → [integrated-analysis](../integrated-analysis/SKILL.md).

## Detail Files

- [IMPLEMENTATION.md](./IMPLEMENTATION.md) — alignment helpers, per-modality transforms (counts vs imaging MAD-z), MuData, pitfalls
- [REFERENCE.md](./REFERENCE.md) — checklists

## Related Skills

- [spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md)
- [cell-quality-control](../cell-quality-control/SKILL.md)
- [normalization-strategies](../normalization-strategies/SKILL.md)
- [imaging-pipeline](../imaging-pipeline/SKILL.md)
- [ops-pipeline](../ops-pipeline/SKILL.md)
- [integrated-analysis](../integrated-analysis/SKILL.md)

## Disclaimer

This content is provided "as-is" without warranty of any kind. AI-generated output may contain
inaccuracies. Review and verify all outputs before relying on them.
