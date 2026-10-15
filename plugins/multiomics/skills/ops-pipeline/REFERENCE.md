# OPS Pipeline — Reference

Quick-reference tables for OPS AnnData layout, metadata schemas, table naming, and confidence thresholds.

---

## Table of Contents

- [OPS AnnData Slot Inventory](#ops-anndata-slot-inventory)
- [target_assignment Column Reference](#target_assignment-column-reference)
- [uns["metadata_attrs"] for OPS](#unsmetadata_attrs-for-ops)
- [Table Name Discovery Pattern](#table-name-discovery-pattern)
- [batches DataFrame Columns Relevant to OPS](#batches-dataframe-columns-relevant-to-ops)
- [diss_targets_metadata Schema](#diss_targets_metadata-schema)
- [Confidence Threshold Reference](#confidence-threshold-reference)
- [OPS vs Other Modalities](#ops-vs-other-modalities)
- [load_tables() Parameters for OPS](#load_tables-parameters-for-ops)

---

## OPS AnnData Slot Inventory

| Slot | Key / Column | Type | Contents |
|------|-------------|------|---------|
| `.X` | — | `scipy.sparse.csr_matrix` (int32) | Raw polony read counts: cells × OPS barcodes |
| `.obs` | `Cell` | int | Cell identifier |
| `.obs` | `Well` | str | Plate position (e.g. `"A1"`) |
| `.obs` | `WellLabel` | str | Sample-sheet name on that well. Color/groupby uses this, not `well_label`. |
| `.obs` | `tile_id` | str | Tile identifier |
| `.obs` | `X_tile`, `Y_tile` | float | Tile-local centroid coordinates (pixels) |
| `.obs` | `cell_id_local` | int | Local cell ID within tile |
| `.obs` | `cell_id_global` | int | Global cell ID (`tile_index << 16 \| local_id`) |
| `.obs` | `region_key` | str | SpatialData region reference |
| `.obs` | `well_id` | str | Well identifier string |
| `.obs` | `MostAbundantTarget` | str | **Added by `include_target_assignment=True`** — OPS identity call or `"Undetermined"` |
| `.obs` | `MostAbundantTargetCount` | int | **Added by `include_target_assignment=True`** — winning library-member polony count |
| `.obs` | `TotalCount` | int | **Added by `include_target_assignment=True`** — total OPS polony reads for this cell |
| `.obs` | `Confidence` | float [0–1] | **Added by `include_target_assignment=True`** — `MostAbundantTargetCount / (TotalCount + 1)` |
| `.obs` | `target_gene` | str | **Added by `include_target_assignment=True`** — gene-level annotation for the winning target |
| `.var` | `target_id` | int | Numeric target identifier |
| `.var` | `target_gene` | str | Gene / barcode identity label |
| `.var` | `target_sequence` | str | Barcode DNA sequence |
| `.var` | `control_type` | str | Control classification (empty for analysis targets) |
| `.var` | `batch_name` | str | Which batch this barcode belongs to (e.g. `"B01"`) |
| `.var` | `target_site_name` | str | Custom panel name from the run manifest |
| `obsm["target_assignment"]` | — | `pd.DataFrame` (cells × 5) | Full assignment DataFrame — present before `include_target_assignment` merging into obs |
| `obsm["total_counts"]` | batch column | `pd.DataFrame` | Per-batch total polony counts |
| `obsm["data_missing"]` | batch column | `pd.DataFrame` | Missing data flags per batch |
| `obsm["spatial"]` | `X`, `Y` | `pd.DataFrame` | Region-pixel centroids (convert to ndarray before squidpy) |
| `obsm["cell_geometry"]` | `Area`, `NuclearArea`, `Xum_tile`, `Yum_tile`, etc. | `pd.DataFrame` | Cell morphology features |
| `uns["metadata_attrs"]` | — | dict | Batch/target metadata (see below) |

---

## target_assignment Column Reference

`obsm["target_assignment"]` is a DataFrame with 5 columns, available after `load_tables(..., include_target_assignment=True)`. These columns are also merged into `obs` when that flag is set.

| Column | Type | Meaning |
|--------|------|---------|
| `MostAbundantTarget` | str | Winning barcode name (e.g. cell line or perturbation identifier), or `"Undetermined"` if no reads |
| `MostAbundantTargetCount` | int | Read count for the winning barcode |
| `TotalCount` | int | Total polony reads summed across all barcodes for this cell |
| `Confidence` | float [0–1] | `MostAbundantTargetCount / (TotalCount + 1)` (+1 pseudocount from build; equals 0 when `TotalCount == 0`) |
| `target_gene` | str | Gene-level annotation for the winning target |

**Undetermined** cells have `TotalCount == 0` (no polony reads at all) or `Confidence < some threshold`. A cell can have `TotalCount > 0` but still be `Undetermined` if reads are split across many barcodes below the confidence cutoff.

---

## uns["metadata_attrs"] for OPS

```python
ops.uns["metadata_attrs"]
# {
#     "batch_class":      "DISSBatch",
#     "batch_name":       "B01",         # or np.array(["B01"], dtype=StringDType())
#     "target_type":      "OPS",         # always "OPS" for this batch type
#     "target_site_name": "<custom>",    # custom name set in the run manifest
# }
```

Key fields:

| Field | Value | Meaning |
|-------|-------|---------|
| `batch_class` | `"DISSBatch"` | OPS uses the DISS (in-situ sequencing) chemistry |
| `target_type` | `"OPS"` | Identifies this as an OPS batch — use this for programmatic detection |
| `batch_name` | e.g. `"B01"` | Batch identifier |
| `target_site_name` | custom string | Panel / OPS library name from the run manifest — varies between runs |

---

## Table Name Discovery Pattern

OPS table names do **not** follow the standard `"Transcript"` / `"Protein"` pattern. They encode both the batch name and the target site name:

```
global-table-{batch_name}-{target_site_name}
```

As seen by `loader.available_tables` (which strips the `global-table-` prefix):

```
{batch_name}-{target_site_name}
e.g. "B01-<CustomPanelName>"
```

**Discovery recipes (in order of preference):**

```python
run = next(iter(loader._runs.values()))

# Preferred — then drop targeted DISS QC/gene panels (often labeled target_type=OPS)
specialized = run.table_suffixes_for_target_type(("OPS", "SpecializedTargeted"))
ops_tables = [t for t in specialized if "qc_panel" not in t.lower()]
# e.g. keep B01-<LibraryName>; drop B01-QC_Panel (GAPDH/ACTB — DISS, not OPS assay)

# Optional: batches DataFrame
ops_batches = run.batches[run.batches["batch_type"] == "SpecializedTargeted"]

# Optional: metadata_attrs scan (still classify — QC panels also say OPS)
ops_tables_via_attrs = []
for t in loader.available_tables:
    attrs = run._read_table_metadata_attrs(f"global-table-{t}")
    if attrs.get("target_type") in ("OPS", "SpecializedTargeted") and "qc_panel" not in t.lower():
        ops_tables_via_attrs.append(t)
```

---

## batches DataFrame Columns Relevant to OPS

`run.batches` is a DataFrame with one row per batch. Key columns for OPS:

| Column | OPS value | Description |
|--------|-----------|-------------|
| `batch_name` | e.g. `"B01"` | Batch identifier |
| `batch_class` | `"DISSBatch"` | Chemistry class — always `DISSBatch` for OPS |
| `batch_type` | `"SpecializedTargeted"` | Identifies OPS batches (vs `ThreePrimeUntargeted`) |
| `min_cycles` | e.g. `25` | Sequencing cycles |
| `num_diss_targets` | e.g. `25` | Number of DISS targets for this batch |

**Finding OPS batches:**

```python
# Filter by batch_type
ops_batches = run.batches[run.batches["batch_type"] == "SpecializedTargeted"]

# Preferred for **table** names (not batch names alone); exclude QC/gene panels
specialized = run.table_suffixes_for_target_type(("OPS", "SpecializedTargeted"))
ops_tables = [t for t in specialized if "qc_panel" not in t.lower()]
```

---

## diss_targets_metadata Schema

`run.diss_targets_metadata` is a DataFrame from the `global-table-metadata` table. It contains per-barcode metadata for all DISS (OPS + 3′ untargeted) batches. The exact columns depend on the run manifest schema.

| Column | Description |
|--------|-------------|
| `TargetSiteName` | Name of the barcode panel / target site (custom, from run manifest) |
| `TargetMask` | Sequencing cycle mask |
| `CellAssignment` | `True` for OPS barcodes (used for cell identity); `False` for 3′ RT oligos |
| `RTOligos` | Reverse transcription oligo sequences |
| `SequencingPrimer` | Sequencing primer sequence |
| `batch_name` | Which batch this target belongs to |

**Filtering to OPS-only barcodes:**

```python
diss_meta = run.diss_targets_metadata

# OPS barcodes (used for cell identity assignment)
ops_targets = diss_meta[diss_meta["CellAssignment"] == True]

# 3′ untargeted RT oligos (not used for cell assignment)
rna_targets = diss_meta[diss_meta["CellAssignment"] == False]

# Filter DISS metadata to SpecializedTargeted batches
ops_batch_names = run.batches.loc[
    run.batches["batch_type"] == "SpecializedTargeted", "batch_name"
].tolist()
ops_batch_targets = diss_meta[diss_meta["batch_name"].isin(ops_batch_names)]
```

---

## Confidence Threshold Reference

| Threshold | Meaning | Recommended use |
|-----------|---------|----------------|
| `Confidence > 0` | Cell received at least one OPS read | Loose — includes noisy calls |
| `Confidence >= 0.3` | Minimal signal; 30% of reads support the call | Use only for exploratory analysis |
| `Confidence >= 0.5` | **Standard filter** — `max / (total + 1) ≥ 0.5` | Default for analysis and DE |
| `Confidence >= 0.7` | High-confidence calls only | Stricter — reduces cell count |
| `Confidence >= 0.9` | Stringent — nearly unambiguous calls | Very low yield; rarely used |

**Run-level interpretation:**

| Pct confident (>= 0.5) | Assessment |
|------------------------|-----------|
| > 65% | Excellent OPS data quality |
| 50–65% | Good; typical for OPS runs |
| 35–50% | Acceptable; check polony density |
| < 35% | Poor assignment; investigate sequencing or panel |

**Why `Undetermined`?**

A cell is `Undetermined` when:
1. No polony reads mapped to the cell (`TotalCount == 0`) — no OPS DISS signal in that cell
2. Reads are split across many barcodes and no single barcode dominates (low `Confidence`)

---

## OPS vs Other Modalities

| Aspect | OPS (`SpecializedTargeted`) | `ThreePrimeUntargeted` | `Transcript` |
|--------|---------------------------|------------------------|-------------|
| Table name | Custom: `{batch}-{site_name}` | `"ThreePrimeUntargeted"` | `"Transcript"` |
| `.X` content | OPS library-member counts (DISS) | Gene counts (DISS) | Gene/protein counts |
| Assignment in obs | `MostAbundantTarget`, `Confidence` | N/A | N/A |
| Purpose | Cell identity (which barcode/perturbation) | Transcriptome profiling | Targeted gene panel |
| `batch_type` | `SpecializedTargeted` | `ThreePrimeUntargeted` | `BarcodingBatch` |
| `batch_class` | `DISSBatch` | `DISSBatch` | `BarcodingBatch` |
| `include_target_assignment` | **Must be `True`** | Ignored | Ignored |
| `diss_targets_metadata` | Used (`CellAssignment=True`) | Used (`CellAssignment=False`) | Not used |
| Typical var count | 10–100 barcodes | 10,000–25,000 genes | 100–500 genes |
| Discovery via | `table_suffixes_for_target_type(("OPS", "SpecializedTargeted"))` then drop `QC_Panel` / gene panels | `"ThreePrimeUntargeted" in loader.available_tables` | `"Transcript" in loader.available_tables` |

---

## load_tables() Parameters for OPS

```python
# Discover the OPS assay table name first (custom, not "OPS"; skip QC/gene panels)
specialized = run.table_suffixes_for_target_type(("OPS", "SpecializedTargeted"))
ops_tables = [t for t in specialized if "qc_panel" not in t.lower()]
ops_table_name = ops_tables[0]

ops = loader.load_tables(
    table=ops_table_name,            # Required: use the discovered custom table name, not "OPS"
    lazy=False,                       # Required for downstream Scanpy / pp operations
    include_target_assignment=True,   # Required: merges assignment columns into obs
    include_cell_geometry=False,      # Optional: True appends geometry columns to X/var
    concatenated=True,                # Default: single AnnData across runs
)
```

| Parameter | OPS recommendation | Reason |
|-----------|-------------------|--------|
| `table` | Discovered custom name (e.g. `"{batch}-{site}"`) | Table never stored as `"OPS"` |
| `lazy` | `False` | Assignment columns and sparse X need materialisation |
| `include_target_assignment` | `True` | Essential — merges `MostAbundantTarget`, `Confidence`, etc. into `obs` |
| `include_cell_geometry` | Usually `False` | Geometry available via `obsm["cell_geometry"]`; only merge if needed for morphology analysis |
