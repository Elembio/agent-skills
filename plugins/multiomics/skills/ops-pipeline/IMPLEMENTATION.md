# OPS Pipeline — Implementation

Full step-by-step code for the AVITI24 OPS + RNA workflow. Steps 0–12 correspond to the pipeline order in [SKILL.md](./SKILL.md).

---

## Table of Contents

- [Step 0: Environment Setup](#step-0-environment-setup)
- [Step 1: Open Store and Discover OPS Table](#step-1-open-store-and-discover-ops-table)
- [Step 2: Load OPS Table](#step-2-load-ops-table)
- [Step 3: Per-Well Assignment Statistics](#step-3-per-well-assignment-statistics)
- [Step 4: Load RNA Table](#step-4-load-rna-table)
- [Step 5: Alignment Asserts](#step-5-alignment-asserts)
- [Step 6: Merge OPS into RNA obs](#step-6-merge-ops-into-rna-obs)
- [Step 7: RNA QC](#step-7-rna-qc)
- [Step 8: Intersect with OPS Confidence Filter](#step-8-intersect-with-ops-confidence-filter)
- [Step 9: Normalize + UMAP](#step-9-normalize--umap)
- [Step 10: Visualize by Condition and OPS Line](#step-10-visualize-by-condition-and-ops-line)
- [Step 11: Pseudobulk DE by Cell Line × Condition](#step-11-pseudobulk-de-by-cell-line--condition)
- [Step 12: Save Results](#step-12-save-results)
- [Common Pitfalls](#common-pitfalls)

---

## Step 0: Environment Setup

```python
from __future__ import annotations

import os
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import seaborn as sns

import elembio_spatialdata_analysis as ebsd
from elembio_spatialdata_analysis.load import Loader

warnings.filterwarnings("ignore")

STORE_PATH = Path("/path/to/run.zarr.zip")
os.makedirs("figures", exist_ok=True)
os.makedirs("results", exist_ok=True)

# Thresholds
CONF_THR          = 0.5    # Confidence >= this to call a cell as assigned
MIN_CELLS_PER_LINE = 20    # Minimum cells per cell line for pseudobulk
MIN_REP_DESEQ     = 2      # Minimum wells per condition for DESeq2
PADJ_THR          = 0.05
MIN_LFC_ABS       = 1.0

sns.set_theme(style="whitegrid", context="notebook")
plt.rcParams["figure.dpi"] = 110
sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=110, facecolor="white")

print("Setup complete")
```

---

## Step 1: Open Store and Discover OPS Table

OPS table names follow the pattern `{batch_name}-{target_site_name}` where `target_site_name` is a custom name set in the run manifest. They are **not** accessible as `table="OPS"`. Always discover the name at runtime from the loader metadata.

```python
loader = Loader(STORE_PATH)
print(loader.summary())
print("Available tables:", sorted(loader.available_tables))

run = next(iter(loader._runs.values()))
print("Target types:", run.available_target_types)

# Inspect batch metadata
print("\nBatches:")
print(run.batches[["batch_name", "batch_class", "batch_type"]].to_string())

# Find OPS batches using available_target_types (preferred) or batches DataFrame
if "OPS" not in run.available_target_types:
    raise ValueError("No OPS batch found. Check run.available_target_types.")

# Preferred: table suffixes (handles OPS + SpecializedTargeted aliases)
specialized = run.table_suffixes_for_target_type(("OPS", "SpecializedTargeted"))
# Guardrail: both OPS and QC/gene panels are DISS SpecializedTargeted. Build metadata
# often labels QC/housekeeping panels as target_type=OPS (e.g. B01-QC_Panel = GAPDH+ACTB).
# Those are NOT Optical Pooled Screening identity libraries — do not Confidence-gate them.
ops_tables = [t for t in specialized if "qc_panel" not in t.lower() and not t.lower().endswith("-qc")]
qc_or_gene_panels = [t for t in specialized if t not in ops_tables]
print(f"\nSpecialized DISS tables: {specialized}")
print(f"OPS assay tables: {ops_tables}")
print(f"Skip (targeted DISS QC/gene panels): {qc_or_gene_panels}")
if not ops_tables:
    raise ValueError(
        "No OPS assay table found. "
        f"SpecializedTargeted/OPS suffixes={specialized}. "
        "Targeted DISS QC/gene panels are not ops-pipeline — use run-quality-overview depth summaries."
    )

ops_table_name = ops_tables[0]
print(f"OPS assay table name: {ops_table_name!r}")

# Inspect OPS barcode panel via diss_targets_metadata
diss_meta = run.diss_targets_metadata
ops_diss = diss_meta[diss_meta["batch_name"].isin(ops_batch_names)]
print(f"\nOPS targets in diss_targets_metadata: {len(ops_diss)}")
print(ops_diss.head(3).to_string())
```

---

## Step 2: Load OPS Table

<!-- canonical: spatialdata-loading-and-access memory guard. Copy verbatim; edit only in that skill's SKILL.md -->
**Memory guard — before any `lazy=False`, `.to_memory()` or `.compute()`:**

1. Load lazily first and run `check_before` from `scripts/check_load.py` in `spatialdata-loading-and-access`. Run the script unchanged in the Python session that holds the data (`exec` its contents); don't rewrite it. It compares each table's estimated in-memory size with the memory this session can use, lists large session objects you could `del`, and raises on a FAIL.
2. Materialize only the tables and runs the plan needs. With several runs, load `concatenated=False` so it checks each run and their total. A WARN says what to `del` or to load one run at a time. A FAIL means subset before materializing: propose a subset to the user (runs, wells, or a set of features), then slice the lazy table and call `.to_memory()` on the slice, for example `adata[adata.obs["Well"].isin(wells)].to_memory()`. `.copy()` keeps a lazy table lazy. Feature subsets by modality: [ACCESS.md — Subset](../spatialdata-loading-and-access/ACCESS.md#subset-by-well-or-feature).
3. Keep `.X` in the form it was stored (`scipy.sparse.issparse(adata.X)`). Sparse tables (3′, OPS) are mostly zeros: never densify a whole one (`.toarray()`, `.todense()`, `np.asarray(adata.X)`); densify only a subset of cells or features. Dense tables (imaging, legacy barcoding) are float32 with `NaN` where a batch has no data: do not convert them to sparse, which turns missing into zero.
<!-- /canonical -->

Always pass `include_target_assignment=True` to merge the per-cell assignment into `obs`.

```python
check_before(loader.load_tables(table=ops_table_name, concatenated=False))  # memory guard: lazy, raises on FAIL
ops = loader.load_tables(
    table=ops_table_name,
    include_target_assignment=True,
    lazy=False,
)

print(ops)
print(f"\nShape: {ops.n_obs:,} cells × {ops.n_vars} barcodes")
print(f"Sparse X: {sp.issparse(ops.X)}")

# Verify assignment columns arrived in obs
assignment_cols = ["MostAbundantTarget", "MostAbundantTargetCount", "TotalCount", "Confidence"]
print("\nOPS obs assignment columns:")
print(ops.obs[assignment_cols].head())

# Quick barcode panel view
print("\nBarcodes (var):")
print(ops.var[["target_id", "target_gene", "control_type", "batch_name"]].head(10))

# uns metadata
print("\nmetadata_attrs:", ops.uns["metadata_attrs"])
```

---

## Step 3: Per-Well Assignment Statistics

Compute and display assignment rates per well before any filtering. This is the primary OPS QC step.

```python
def compute_well_assignment_stats(obs: pd.DataFrame, conf_thr: float = 0.5) -> pd.DataFrame:
    """Per-well assignment statistics for an OPS obs DataFrame."""
    g = obs.copy()
    g["_assigned"] = g["MostAbundantTarget"] != "Undetermined"
    g["_confident"] = g["Confidence"] >= conf_thr
    stats = g.groupby("Well").agg(
        n_cells=("_assigned", "count"),
        n_assigned=("_assigned", "sum"),
        n_confident=("_confident", "sum"),
        median_total_count=("TotalCount", "median"),
        median_confidence=("Confidence", "median"),
    ).assign(
        pct_assigned=lambda d: (100 * d["n_assigned"] / d["n_cells"]).round(1),
        pct_confident=lambda d: (100 * d["n_confident"] / d["n_cells"]).round(1),
    )
    return stats

stats = compute_well_assignment_stats(ops.obs, conf_thr=CONF_THR)
print("=== OPS Assignment Statistics (per well) ===")
print(stats.to_string())

# Run-level summary
n_total = ops.n_obs
n_assigned = int((ops.obs["MostAbundantTarget"] != "Undetermined").sum())
n_confident = int((ops.obs["Confidence"] >= CONF_THR).sum())
print(f"\nRun total: {n_total:,} cells")
print(f"  Assigned (TotalCount >= 1):  {n_assigned:,}  ({100*n_assigned/n_total:.1f}%)")
print(f"  Confident (Conf >= {CONF_THR}):    {n_confident:,}  ({100*n_confident/n_total:.1f}%)")

# Cell line composition (excluding Undetermined)
line_counts = ops.obs.loc[ops.obs["MostAbundantTarget"] != "Undetermined", "MostAbundantTarget"].value_counts()
print(f"\nTop cell lines (assigned cells):\n{line_counts.head(10).to_string()}")

# Confidence distribution plot
fig, axes = plt.subplots(1, 2, figsize=(12, 4))
ops.obs["Confidence"].hist(bins=50, ax=axes[0])
axes[0].axvline(CONF_THR, color="red", linestyle="--", label=f"threshold={CONF_THR}")
axes[0].set_xlabel("Confidence")
axes[0].set_title("OPS Confidence distribution")
axes[0].legend()

ops.obs["TotalCount"].clip(upper=20).hist(bins=20, ax=axes[1])
axes[1].set_xlabel("TotalCount (clipped at 20)")
axes[1].set_title("Polony reads per cell")
plt.tight_layout()
plt.savefig("figures/ops_qc_distributions.png", dpi=150)
plt.show()
```

---

## Step 4: Load RNA Table

Load the paired RNA table (usually `ThreePrimeUntargeted`; may be `Transcript` for targeted panels).

```python
# Determine RNA table
rna_table = "ThreePrimeUntargeted" if "ThreePrimeUntargeted" in loader.available_tables else "Transcript"
print(f"Loading RNA table: {rna_table!r}")

check_before(loader.load_tables(table=rna_table, concatenated=False))  # memory guard: lazy, raises on FAIL
rna = loader.load_tables(table=rna_table, lazy=False)
print(rna)
print(f"\nRNA: {rna.n_obs:,} cells × {rna.n_vars:,} genes")
print("RNA obs columns:", rna.obs.columns.tolist())
```

---

## Step 5: Alignment Asserts

OPS and RNA tables from the same store always have identical cell ordering. Verify this before merging.

```python
assert ops.n_obs == rna.n_obs, \
    f"Cell count mismatch: ops={ops.n_obs:,}, rna={rna.n_obs:,}"

assert ops.obs_names.equals(rna.obs_names), \
    "obs_names do not match — tables are not aligned!"

cig_match = (ops.obs["cell_id_global"].values == rna.obs["cell_id_global"].values).all()
assert cig_match, "cell_id_global arrays differ between OPS and RNA tables!"

def _spatial_xy(adata):
    sp = adata.obsm["spatial"]
    return (
        sp[["X", "Y"]].to_numpy(dtype=float)
        if hasattr(sp, "columns")
        else np.asarray(sp, dtype=float)
    )

assert np.allclose(_spatial_xy(ops), _spatial_xy(rna), equal_nan=True), \
    "obsm['spatial'] differs between OPS and RNA tables!"

print("Alignment checks passed:")
print(f"  n_obs:          {ops.n_obs:,} ✓")
print(f"  obs_names:      ✓")
print(f"  cell_id_global: ✓")
print(f"  spatial:        ✓")
```

---

## Step 6: Merge OPS into RNA obs

Copy OPS assignment fields into `rna.obs` **before** any filtering. This preserves full-denominator assignment rates and allows the confidence filter to be applied after RNA QC.

```python
# Number of distinct barcodes with reads per cell (signal diversity)
X_ops = ops.X.toarray() if sp.issparse(ops.X) else ops.X
n_barcodes_with_reads = (X_ops > 0).sum(axis=1)

rna.obs["OPS_line"]             = ops.obs["MostAbundantTarget"].astype(str).values
rna.obs["OPS_Confidence"]       = ops.obs["Confidence"].astype(float).values
rna.obs["OPS_TotalCount"]       = ops.obs["TotalCount"].astype(float).values
rna.obs["OPS_MostAbundantCount"] = ops.obs["MostAbundantTargetCount"].astype(float).values
rna.obs["OPS_n_barcodes"]       = n_barcodes_with_reads
rna.obs["ops_keep"]             = (ops.obs["Confidence"] >= CONF_THR).values

# Prefer experiment-design. If labeling from cells, parse WellLabel (not well_label):
def condition_from_well_label(label: str) -> str:
    """Extract condition from obs['WellLabel'] (sample-sheet name)."""
    s = str(label)
    # Example: 'Sample1_DMSO_1' -> 'DMSO', 'Sample1_DrugA_1' -> 'DrugA'
    for keyword in ["DMSO", "DrugA", "DrugB", "Treated", "Control", "Vehicle"]:
        if keyword.lower() in s.lower():
            return keyword
    return s  # return full label if no keyword matched

rna.obs["condition"] = rna.obs["WellLabel"].astype(str).map(condition_from_well_label)

print(f"OPS fields merged into rna.obs. ops_keep={rna.obs['ops_keep'].sum():,} cells")
print("\nCondition breakdown:")
print(rna.obs["condition"].value_counts().to_string())
print("\nCell line × condition (top 10 lines):")
ct = pd.crosstab(rna.obs["OPS_line"], rna.obs["condition"])
print(ct.head(10).to_string())
```

---

## Step 7: RNA QC

Run standard RNA quality control. For `ThreePrimeUntargeted`, this includes `compute_qc_metrics` and `filter_cells`. See [cell-quality-control](../cell-quality-control/SKILL.md) and [3prime-untargeted-pipeline](../3prime-untargeted-pipeline/SKILL.md) for parameter details.

```python
ebsd.pp.compute_qc_metrics(rna, rna_table)

# Optional: visualize QC before filtering
qc = rna.obsm["qc_metrics"]
fig, ax = plt.subplots(figsize=(10, 4))
plot_df = pd.DataFrame({
    "Well": rna.obs["Well"].astype(str).values,
    "total_counts": qc["total_counts"].values,
})
sns.boxplot(data=plot_df, x="Well", y="total_counts", ax=ax, color="lavender")
ax.set_title(f"{rna_table} total_counts by well (pre-filter)")
ax.tick_params(axis="x", rotation=45)
plt.tight_layout()
plt.savefig("figures/rna_qc_boxplot_prefilter.png", dpi=150)
plt.show()

# Apply RNA QC filter
n_before_rna = rna.n_obs
ebsd.pp.filter_cells(rna, rna_table)
print(f"After RNA QC: {rna.n_obs:,} (was {n_before_rna:,})")
ebsd.pp.generate_filter_report(rna)
```

---

## Step 8: Intersect with OPS Confidence Filter

After RNA QC, apply the OPS confidence filter as a second pass. The `ops_keep` flag was set before filtering so it still references the original cell population.

```python
n_before_ops = rna.n_obs
rna = rna[rna.obs["ops_keep"]].copy()
print(f"After OPS confidence filter (Conf >= {CONF_THR}): {rna.n_obs:,} (was {n_before_ops:,})")

# Final cell line counts (must have enough cells for DE)
line_counts_final = rna.obs["OPS_line"].value_counts()
line_counts_final = line_counts_final[line_counts_final >= MIN_CELLS_PER_LINE]
usable_lines = line_counts_final.index.tolist()
print(f"\nUsable cell lines (>= {MIN_CELLS_PER_LINE} cells): {len(usable_lines)}")
print(line_counts_final.to_string())

# Filter to usable lines
rna = rna[rna.obs["OPS_line"].isin(usable_lines)].copy()
print(f"After cell-line filter: {rna.n_obs:,} cells, {len(usable_lines)} lines")
```

---

## Step 9: Normalize + UMAP

For `ThreePrimeUntargeted`: use `ebsd.pp.normalize` then HVG selection. For targeted panels: use `ebsd.pp.normalize` directly without HVG.

```python
ebsd.pp.normalize(rna, rna_table)
rna = rna.to_memory()

if rna_table == "ThreePrimeUntargeted":
    # HVG selection for 3′ untargeted data
    sc.pp.highly_variable_genes(rna, n_top_genes=2000, subset=True)
    print(f"HVG selected: {rna.n_vars:,} genes")

sc.pp.pca(rna, svd_solver="arpack")
sc.pp.neighbors(rna, n_neighbors=15, n_pcs=min(40, rna.obsm["X_pca"].shape[1]))
sc.tl.umap(rna)
print("UMAP complete")
```

---

## Step 10: Visualize by Condition and OPS Line

```python
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

sc.pl.umap(
    rna, color="condition", ax=axes[0],
    show=False, title="UMAP by Condition",
    palette="tab10",
)
sc.pl.umap(
    rna, color="OPS_line", ax=axes[1],
    show=False, title="UMAP by Cell Line (OPS)",
    legend_fontsize=6, legend_loc="on data",
)
plt.tight_layout()
plt.savefig("figures/umap_condition_and_line.png", dpi=150, bbox_inches="tight")
plt.show()

# Optionally color by OPS confidence
fig, ax = plt.subplots(figsize=(7, 6))
sc.pl.umap(rna, color="OPS_Confidence", ax=ax, show=False, title="OPS Confidence", vmin=0, vmax=1)
plt.tight_layout()
plt.savefig("figures/umap_ops_confidence.png", dpi=150)
plt.show()
```

---

## Step 11: Pseudobulk DE by Cell Line × Condition

Pseudobulk DE compares conditions within each cell line. Requires `>= MIN_REP_DESEQ` wells per condition after filtering. See [pseudobulk-differential-expression](../pseudobulk-differential-expression/SKILL.md) for full API details.

```python
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import anndata as ad

# Pseudobulk key: condition × well (one pseudo-replicate per well)
rna.obs["pseudo_key"] = rna.obs["condition"] + "__" + rna.obs["Well"].astype(str)

de_results = {}

for cell_line in usable_lines:
    sub = rna[rna.obs["OPS_line"] == cell_line].copy()

    # Check replicates per condition
    rep_counts = sub.obs.groupby("condition")["Well"].nunique()
    valid_conditions = rep_counts[rep_counts >= MIN_REP_DESEQ].index.tolist()
    if len(valid_conditions) < 2:
        print(f"  Skipping {cell_line}: not enough replicates in conditions {valid_conditions}")
        continue

    sub = sub[sub.obs["condition"].isin(valid_conditions)].copy()

    # Pseudobulk aggregation — sum counts per (condition, well)
    raw_layer = sub.layers.get("raw", sub.X)
    raw_X = raw_layer.toarray() if sp.issparse(raw_layer) else raw_layer.copy()

    agg_rows = []
    for key in sub.obs["pseudo_key"].unique():
        mask = sub.obs["pseudo_key"] == key
        if mask.sum() < 3:
            continue
        agg_rows.append({
            "pseudo_key": key,
            "condition": sub.obs.loc[mask, "condition"].iloc[0],
            "Well": sub.obs.loc[mask, "Well"].iloc[0],
            "n_cells": int(mask.sum()),
            "counts": raw_X[mask].sum(axis=0),
        })

    if len(agg_rows) < 4:
        print(f"  Skipping {cell_line}: too few pseudobulk replicates ({len(agg_rows)})")
        continue

    pb_obs = pd.DataFrame([{k: v for k, v in r.items() if k != "counts"} for r in agg_rows])
    pb_obs.index = pb_obs["pseudo_key"].values
    pb_X = np.vstack([r["counts"] for r in agg_rows])

    pb_adata = ad.AnnData(
        X=pb_X.astype(int),
        obs=pb_obs,
        var=sub.var[[]].copy(),
    )

    try:
        dds = DeseqDataSet(
            adata=pb_adata,
            design_factors="condition",
            refit_cooks=True,
        )
        dds.deseq2()

        reference = "DMSO" if "DMSO" in valid_conditions else valid_conditions[0]
        treatment = [c for c in valid_conditions if c != reference][0]

        stat_res = DeseqStats(dds, contrast=["condition", treatment, reference])
        stat_res.summary()
        results_df = stat_res.results_df.copy()
        results_df["cell_line"] = cell_line
        results_df["treatment"] = treatment
        results_df["reference"] = reference
        de_results[cell_line] = results_df

        sig = results_df.dropna(subset=["padj"])
        sig = sig[(sig["padj"] < PADJ_THR) & (sig["log2FoldChange"].abs() >= MIN_LFC_ABS)]
        print(f"  {cell_line}: {len(sig)} DE genes (|LFC|>={MIN_LFC_ABS}, padj<{PADJ_THR})")

    except Exception as e:
        print(f"  {cell_line}: DESeq2 failed — {e}")

# Combine all results
if de_results:
    all_de = pd.concat(de_results.values(), axis=0)
    all_de.to_csv("results/ops_de_by_cell_line.csv")
    print(f"\nDE results saved: {len(all_de):,} rows across {len(de_results)} cell lines")
```

---

## Step 12: Save Results

```python
import os

# Save the filtered, normalized AnnData
out_path = "results/rna_ops_filtered.h5ad"
rna.write_h5ad(out_path)
print(f"Saved: {out_path}  ({rna.n_obs:,} cells × {rna.n_vars:,} features)")

# Save assignment stats
stats.to_csv("results/ops_assignment_stats_per_well.csv")
print("Saved: results/ops_assignment_stats_per_well.csv")

# Save cell line counts per condition
ct_final = pd.crosstab(rna.obs["OPS_line"], rna.obs["condition"])
ct_final.to_csv("results/ops_cell_line_condition_counts.csv")
print("Saved: results/ops_cell_line_condition_counts.csv")
```

---

## Common Pitfalls

| Pitfall | Cause | Fix |
|---------|-------|-----|
| `table="OPS"` not found | OPS table has a custom name like `B01-CustomPanel` | Discover via `run.table_suffixes_for_target_type(("OPS", "SpecializedTargeted"))` |
| All cells `Confidence=0` | OPS alignment pipeline was not run; `target_name` absent from points | Confirm `target_name` is present in the points table; if it is missing, contact Element support |
| `obs_names` mismatch between OPS and RNA | Using different Loader calls with different filters | Always load both tables from the same `Loader` instance without any filter |
| Low assignment rate (< 40%) | Low polony density, poor sequencing, or panel mismatch | Inspect `median_total_count` in stats; check `diss_targets_metadata` for panel info |
| No cells pass `MIN_CELLS_PER_LINE` | Rare cell line in pool, or over-filtering | Lower `CONF_THR` to 0.3 tentatively and label results as exploratory; if assignment stays low, contact Element support |
| DESeq2 convergence warning | Very few pseudobulk replicates | Require `MIN_REP_DESEQ >= 2`; skip lines with < 2 wells per condition |
| Merging OPS after RNA filter | OPS flag applied to already-filtered cells | Always copy `ops_keep` into `rna.obs` before calling `filter_cells` |
| `include_target_assignment=False` | Default; assignment columns not in `obs` | Always pass `include_target_assignment=True` for OPS tables |
