# Metadata and Target Mapping — Implementation

Detailed implementation steps for the [metadata-and-target-mapping](./SKILL.md) skill.

## Table of Contents

1. [Load Metadata](#1-load-metadata)
2. [Parse Well Metadata](#2-parse-well-metadata)
3. [Parse Tile Metadata](#3-parse-tile-metadata)
4. [Parse Batch Metadata](#4-parse-batch-metadata)
5. [Target Metadata and ID → Name Maps](#5-target-metadata-and-id--name-maps)
6. [Control Targets](#6-control-targets)
7. [DISS / OPS Metadata + table discovery](#7-diss--ops-metadata--table-discovery)
8. [Well labels as recorded](#8-well-labels-as-recorded)
9. [Target–Batch Lookup](#9-targetbatch-lookup)
10. [Parameters and Customization](#10-parameters-and-customization)
11. [Expected Outputs](#11-expected-outputs)
12. [Common Pitfalls](#12-common-pitfalls)

---

## 1. Load Metadata

### Preferred: `Loader` / `RunInfo`

```python
import pandas as pd
from pathlib import Path
from elembio_spatialdata_analysis.load import Loader

store_path = Path("/path/to/run.zarr")
loader = Loader(store_path)
run_key, run = next(iter(loader._runs.items()))  # one run_key per store

wells_df = run.wells
tiles_df = run.tiles
batches_df = run.batches
targets_barcoding = run.targets_barcoding
targets_imaging = run.targets_imaging
diss_targets_metadata = run.diss_targets_metadata
```

### Alternative: `SpatialData` + `global-table-metadata.uns`

```python
import spatialdata as sd
import pandas as pd

sdata = sd.read_zarr("/path/to/run.zarr")
metadata_table = sdata.tables["global-table-metadata"]
uns = metadata_table.uns

wells_df = uns["wells"]
tiles_df = uns["tiles"]
batches_df = uns["batches"]
targets_barcoding = uns["targets_barcoding"]
targets_imaging = uns.get("targets_imaging", pd.DataFrame())
diss_targets_metadata = uns.get("diss_targets_metadata", pd.DataFrame())
```

---

## 2. Parse Well Metadata

```python
print(wells_df.head())

well_a1 = wells_df.loc["A1"]  # index is well name (A1)
label_col = "well_label" if "well_label" in wells_df.columns else "WellLabel"
print(f"Well A1 sample-sheet name: {well_a1[label_col]}")
well_names = wells_df.index.tolist()
# On a cell table the same string is obs["WellLabel"]; Well is the plate position.
```

---

## 3. Parse Tile Metadata

```python
well_a1_tiles = tiles_df[tiles_df["well_name"] == "A1"]
tile_positions = tiles_df[["tile_name", "x_mm", "y_mm"]].set_index("tile_name")
```

---

## 4. Parse Batch Metadata

```python
barcoding_batches = batches_df[batches_df["batch_class"] == "BarcodingBatch"]
imaging_batches   = batches_df[batches_df["batch_class"] == "ImagingBatch"]
diss_batches      = batches_df[batches_df["batch_class"] == "DISSBatch"]
transcript_batch_names = barcoding_batches["batch_name"].tolist()
```

### Batch classes (typical)

| Batch class | Typical `batch_type` | Contents |
|-------------|----------------------|----------|
| **DISSBatch** | ThreePrimeUntargeted, SpecializedTargeted | 3′ / OPS (**current**) |
| **ImagingBatch** | PreAmp, PostAmp, StainOnly, … | Imaging targets — CP-group + optional IF |
| **BarcodingBatch** | — | Legacy Transcript / Protein |

---

## 5. Target Metadata and ID → Name Maps

**Current:** `diss_targets_metadata` + `targets_imaging`. **Legacy:** `targets_barcoding`. Column names vary — inspect first.

```python
ti = targets_imaging
diss = diss_targets_metadata
tb = targets_barcoding

# Imaging targets by role (when columns exist)
if len(ti) and "batch_name" in ti.columns:
    print(ti["batch_name"].astype(str).str.slice(0, 2).value_counts())  # CP* vs P* prefix peek
if len(ti) and "target_type" in ti.columns:
    print(ti["target_type"].value_counts(dropna=False))

# DISS: OPS barcodes vs 3′ oligos
if len(diss) and "CellAssignment" in diss.columns:
    ops_targets = diss[diss["CellAssignment"] == True]
    rna_oligos = diss[diss["CellAssignment"] == False]

# Legacy barcoding (often empty on current runs)
if len(tb) and "target_type" in tb.columns:
    transcript_rows = tb[tb["target_type"] == "Transcript"]
    protein_rows = tb[tb["target_type"] == "Protein"]
```

---

## 6. Control Targets

In **metadata** DataFrames, **`control_type`** describes panel design — not per-cell counts. On loaded AnnData, NSB / positive-control layout is modality-specific (legacy Transcript often has `obsm["control_variables"]`; see [anndata-structure-reference](../anndata-structure-reference/SKILL.md)).

```python
if len(tb) and "control_type" in tb.columns:
    control_rows = tb[tb["control_type"].astype(str).str.len() > 0]
```

---

## 7. DISS / OPS Metadata + table discovery

```python
# Canonical SpecializedTargeted DISS suffixes (never table="OPS")
# Split OPS assay vs QC/gene panels — both are DISS; only assay goes to ops-pipeline
specialized = run.table_suffixes_for_target_type(("OPS", "SpecializedTargeted"))
ops_tables = [t for t in specialized if "qc_panel" not in t.lower()]
qc_or_gene_panels = [t for t in specialized if t not in ops_tables]

if len(diss_targets_metadata):
    print(diss_targets_metadata.head())
    print(diss_targets_metadata.columns.tolist())
```

---

## 8. Well labels as recorded

List the sample-sheet names. Do not split them into condition / replicate here.

```python
# Preferred — cell table already loaded
print(
    adata.obs[["Well", "WellLabel"]]
    .astype(str)
    .drop_duplicates()
    .sort_values("Well")
    .to_string(index=False)
)

# Metadata-only fallback (no AnnData yet)
# col = "well_label" if "well_label" in wells_df.columns else "WellLabel"
# print(wells_df[[col]].to_string())
```

To infer `condition` / `replicate` and write those columns, open [experiment-design](../experiment-design/SKILL.md) and confirm with the user.

### `metadata_attrs` on a loaded table

Each `global-table-*` carries `adata.uns["metadata_attrs"]` with `target_type`, `batch_name`, `batch_class`, `merge_strategy`.

```python
meta_attrs = adata.uns.get("metadata_attrs", {})
```

---

## 9. Target–Batch Lookup

```python
batch_info = batches_df.set_index("batch_name")[["batch_class", "batch_type"]]

# Prefer the metadata frame that matches the run (imaging / DISS / legacy barcoding)
panel = targets_imaging if len(targets_imaging) else (
    diss_targets_metadata if len(diss_targets_metadata) else targets_barcoding
)
name_col = "target_name" if "target_name" in panel.columns else panel.columns[0]
target_lookup = panel.merge(batch_info, left_on="batch_name", right_index=True, how="left")

def targets_in_batch(batch_name: str) -> list[str]:
    return target_lookup[target_lookup["batch_name"] == batch_name][name_col].astype(str).tolist()
```

---

## 10. Parameters and Customization

### Modality → metadata table

```python
def targets_for_modality(modality: str, tb, ti, diss) -> pd.DataFrame:
    m = modality.lower().strip()
    if m in ("imaging", "cellpaint", "cell_paint", "proteinif"):
        return ti
    if m in ("ops", "diss", "specializedtargeted", "threeprimeuntargeted", "threeprime"):
        return diss
    if m in ("transcript", "protein"):  # legacy barcoding
        return tb
    raise ValueError(modality)
```

---

## 11. Expected Outputs

After running the recipes above you should have:

1. **`wells_df`**, **`tiles_df`**, **`batches_df`** from `RunInfo` or `uns`
2. **`targets_barcoding`**, **`targets_imaging`**, **`diss_targets_metadata`**
3. **Well / WellLabel listing** (as recorded). Condition / replicate columns come from [experiment-design](../experiment-design/SKILL.md)
4. **Lookups** for target IDs and batch membership

### Validation

```python
assert len(wells_df) > 0, "No wells"
assert (
    len(targets_barcoding) > 0
    or len(targets_imaging) > 0
    or len(diss_targets_metadata) > 0
), "No target metadata"

# Validate against a table that exists
t = "ThreePrimeUntargeted" if "ThreePrimeUntargeted" in loader.available_tables else sorted(loader.available_tables)[0]
adata = loader.load_tables(table=t, lazy=False)
data_wells = set(adata.obs["Well"].unique())
metadata_wells = set(wells_df.index)
missing = data_wells - metadata_wells
assert len(missing) == 0, f"Wells in table but not in metadata: {missing}"
```

---

## 12. Common Pitfalls

1. **Missing `global-table-metadata`**: `Loader` warns; `RunInfo` metadata properties return empty defaults.

2. **Old key names**: `uns["targets"]` / `uns["diss_targets"]` are gone — use `targets_barcoding`, `targets_imaging`, `diss_targets_metadata`.

3. **OPS discovery**: use `table_suffixes_for_target_type(("OPS", "SpecializedTargeted"))`, not `table="OPS"` and not batch-name `startswith` alone.

4. **`targets_imaging` ≠ CellPaint modality** — it is the imaging target list; CP-group vs IF is role inside imaging.

5. **Empty `targets_barcoding`** on current DISS runs is normal — not a broken store.

6. **Well naming**: Match `"A1"` vs `"A01"` consistently across `obs` and `wells_df`.

7. **Duplicate target names across batches**: Disambiguate with `batch_name` in `var` and in metadata.

8. **Parsing `WellLabel` into condition here** — list the labels only. Inferring groupings belongs in [experiment-design](../experiment-design/SKILL.md).
