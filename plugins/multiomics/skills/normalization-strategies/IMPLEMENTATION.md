# Normalization Strategies — Implementation

## Table of Contents

1. [3′ Untargeted / RNA Normalization Workflow](#1-3-untargeted--rna-normalization-workflow)
2. [Protein Normalization](#2-protein-normalization)
3. [Imaging / IF — use imaging-pipeline](#3-imaging--if--use-imaging-pipeline)
4. [Fixed Target Count Normalization](#4-fixed-target-count-normalization)
5. [Multi-Run Normalization](#5-multi-run-normalization)
6. [Before/After Distribution Plots](#6-beforeafter-distribution-plots)
7. [Access Intermediate Layers](#7-access-intermediate-layers)
8. [Per-Gene Normalization Comparison Plot](#8-per-gene-normalization-comparison-plot)
9. [Validation Asserts](#9-validation-asserts)
10. [Parameters Reference and Expected Outputs](#10-parameters-reference-and-expected-outputs) → [REFERENCE.md](./REFERENCE.md)
11. [Common Pitfalls](#11-common-pitfalls)

---

## 1. 3′ Untargeted / RNA Normalization Workflow

The recommended entry point is `ebsd.pp.normalize`. It dispatches by `target_type` and runs the full three-step pipeline (agg-count scaling → median-of-ratios → log1p). Intermediate layers are opt-in.

```python
import elembio_spatialdata_analysis as ebsd

# Start from a filtered AnnData
adata = filtered_adata.copy()

# Full pipeline: agg-count scaling → median-of-ratios → log1p
ebsd.pp.normalize(
    adata,
    "ThreePrimeUntargeted",        # preferred; use "Transcript" for legacy targeted panels
    batch_colname="batch_name",
    well_colname="Well",
    run_colname=None,              # auto-detects run_name / run_id from Loader output
    layer=None,                    # operate on .X
    log_base=None,                 # natural log (recommended)
    raw_layer_name="raw",          # snapshot before any step
    agg_norm_layer_name="norm_agg_counts",      # opt-in intermediate
    mor_norm_layer_name="norm_med_of_ratios",   # opt-in intermediate
    inplace=True,
    target_counts=None,            # data-driven median (recommended)
)

print("Layers after normalization:", list(adata.layers.keys()))
print(f"Raw snapshot range : {adata.layers['raw'].min():.0f} – {adata.layers['raw'].max():.0f}")
print(f"Final .X range     : {adata.X.min():.2f} – {adata.X.max():.2f}")
```

`normalize_transcripts` and `normalize_proteins` are also public for cases where you need a kwarg not forwarded by `ebsd.pp.normalize`. Low-level helpers (`_normalize_agg_counts`, etc.) are private.

---

## 2. Protein Normalization

Protein uses aggregated-count scaling (optional) followed by log1p. There is no CLR step.

```python
from elembio_spatialdata_analysis.load import Loader

loader = Loader("/path/to/run.zarr")
protein_adata = loader.load_tables(table="Protein", lazy=False)

# Default: agg-count scaling → log1p
ebsd.pp.normalize(protein_adata, "Protein")

# Skip agg-count scaling (log1p only)
ebsd.pp.normalize(protein_adata, "Protein", apply_agg_norm=False)

print(f"Protein .X range: {protein_adata.X.min():.2f} – {protein_adata.X.max():.2f}")
```

---

## 3. Imaging / IF — use imaging-pipeline

`ebsd.pp.normalize(adata, "Imaging")` (and `"CellPaint"`) raises `NotImplementedError`. Fluorescence and CellProfiler features are not counts.

**Canonical path:** [imaging-pipeline](../imaging-pipeline/SKILL.md) — build a feature matrix (base + CP + AreaShape; **exclude Bkg / IgG / Cleave**), then **MAD-robust z against the reference the design implies**: control wells when biology is constant within wells and shared controls exist, global when it is constant with no control, per-well only when wells are internally heterogeneous (pooled / OPS). Per-well centering on a one-condition-per-well plate removes the condition. Do not overwrite imaging `X` with count-style normalize. Keep raw intensities for S/B and localization plots.

---

## 4. Fixed Target Count Normalization

Use `target_counts` to scale all cells to the same total rather than the per-batch median. Useful when comparing datasets where the median differs.

```python

adata = filtered_adata.copy()

ebsd.pp.normalize(
    adata,
    "ThreePrimeUntargeted",        # or "Transcript" for legacy targeted panels
    batch_colname="batch_name",
    well_colname="Well",
    target_counts=10_000,
    inplace=True,
)
```

Omit `target_counts` (or pass `None`) to use the data-driven per-batch median — the default and recommended approach for single experiments.

---

## 5. Multi-Run Normalization

When data is loaded via `Loader` with `concatenated=True`, the Loader adds `run_name` (and `source_key`) to `.obs`. `ebsd.pp.normalize` auto-detects this column:

```python

# run_colname=None → auto-detects "run_name" or "run_id" from obs when present
ebsd.pp.normalize(
    adata,
    "ThreePrimeUntargeted",        # or "Transcript" for legacy targeted panels
    batch_colname="batch_name",
    well_colname="Well",
    run_colname=None,
    inplace=True,
)
```

Pass an explicit string (e.g. `run_colname="run_id"`) if your column has a different name.

---

## 6. Before/After Distribution Plots

Requires the intermediate `norm_agg_counts` layer to have been written (pass `agg_norm_layer_name="norm_agg_counts"` to `ebsd.pp.normalize`).

```python
import numpy as np
import scipy.sparse as sp
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

raw = adata.layers["raw"]
agg = adata.layers["norm_agg_counts"]

total_raw = np.asarray(raw.sum(axis=1)).ravel() if sp.issparse(raw) else raw.sum(axis=1)
total_agg = np.asarray(agg.sum(axis=1)).ravel() if sp.issparse(agg) else agg.sum(axis=1)

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

sns.histplot(total_raw, ax=axes[0, 0], bins=50)
axes[0, 0].set_title("Raw: Total Counts per Cell")
axes[0, 0].set_xlabel("Total Counts")

sns.histplot(total_agg, ax=axes[0, 1], bins=50)
axes[0, 1].set_title("After Agg Norm: Total Counts per Cell")
axes[0, 1].set_xlabel("Total Counts")

plot_df = pd.DataFrame({
    "Well": adata.obs["Well"].values,
    "Total_Raw": total_raw,
    "Total_Agg": total_agg,
})

sns.boxplot(data=plot_df, x="Well", y="Total_Raw", ax=axes[1, 0])
axes[1, 0].set_xticklabels(axes[1, 0].get_xticklabels(), rotation=45)
axes[1, 0].set_title("Raw Counts by Well")

sns.boxplot(data=plot_df, x="Well", y="Total_Agg", ax=axes[1, 1])
axes[1, 1].set_xticklabels(axes[1, 1].get_xticklabels(), rotation=45)
axes[1, 1].set_title("After Agg Norm by Well")

plt.tight_layout()
plt.show()
```

---

## 7. Access Intermediate Layers

```python
import scipy.sparse as sp

# Always present
raw_counts = adata.layers["raw"]

# Present only if you set agg_norm_layer_name / mor_norm_layer_name
agg_normalized = adata.layers.get("norm_agg_counts")
mor_normalized = adata.layers.get("norm_med_of_ratios")

# Final normalized result
log_normalized = adata.X

# --- Prepare a copy for pseudobulk DE (needs raw integer counts) ---
adata_for_de = adata.copy()
adata_for_de.X = adata_for_de.layers["raw"]

# --- Use .X directly for clustering / visualization ---
adata_for_viz = adata  # .X is already log-normalized
```

---

## 8. Per-Gene Normalization Comparison Plot

Visualize how normalization transforms the distribution of a single gene. Requires `norm_med_of_ratios` to be present for a three-panel view; falls back to two panels otherwise.

```python
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt

gene = "GAPDH"
if gene in adata.var_names:
    gene_idx = adata.var_names.get_loc(gene)
    mor_normalized = adata.layers.get("norm_med_of_ratios")
    ncols = 3 if mor_normalized is not None else 2

    fig, axes = plt.subplots(1, ncols, figsize=(5 * ncols, 4))

    def _col(m):
        return m[:, gene_idx].toarray().flatten() if sp.issparse(m) else m[:, gene_idx]

    axes[0].hist(_col(adata.layers["raw"]), bins=50)
    axes[0].set_title(f"{gene} — Raw")

    if mor_normalized is not None:
        axes[1].hist(_col(mor_normalized), bins=50)
        axes[1].set_title(f"{gene} — After MoR")
        axes[2].hist(_col(adata.X), bins=50)
        axes[2].set_title(f"{gene} — log-normalized (.X)")
    else:
        axes[1].hist(_col(adata.X), bins=50)
        axes[1].set_title(f"{gene} — log-normalized (.X)")

    plt.tight_layout()
    plt.show()
```

---

## 9. Validation Asserts

```python
import numpy as np
import scipy.sparse as sp

# Layer presence
assert "raw" in adata.layers, "Raw snapshot not stored — check raw_layer_name parameter"

# Sanity-check final values
assert adata.X.max() < 20, "Values seem too high for log-normalized data; expected < 20"

# Total-count CV after aggregated scaling (requires norm_agg_counts layer)
X_agg = adata.layers.get("norm_agg_counts")
if X_agg is None:
    raise AssertionError(
        "norm_agg_counts layer not found. "
        "Pass agg_norm_layer_name='norm_agg_counts' to ebsd.pp.normalize() to enable this check."
    )
totals = np.asarray(X_agg.sum(axis=1)).ravel() if sp.issparse(X_agg) else X_agg.sum(axis=1)
cv = np.std(totals) / np.mean(totals)
assert cv < 0.5, f"High variation in total counts after agg normalization: CV={cv:.2f}"
```

---

## 10. Parameters Reference and Expected Outputs

See [REFERENCE.md](./REFERENCE.md) — full `ebsd.pp.normalize()` parameter table, expected outputs by modality, and the log1p primitive.

---

## 11. Common Pitfalls

1. **Normalizing before filtering** — Outlier cells shift the median used for aggregated-count scaling. Always run QC filtering first.

2. **Using `adata.X` for DE after normalization** — PyDESeq2 and DESeq2 require raw integer counts. Once `ebsd.pp.normalize` has run, `.X` holds log-normalized values; pass `adata.layers["raw"]` instead. If you have not called `ebsd.pp.normalize`, `.X` still contains raw counts and is the correct input.

3. **Missing batch column** — If `batch_name` is absent from `.obs`, normalization is applied globally across all cells. Verify the column exists before calling `ebsd.pp.normalize`.

4. **Memory spike from sparse → dense conversion** — `log1p` may densify a sparse matrix. For large datasets (>500k cells), monitor memory and consider chunking or using a machine with more RAM.

5. **Accessing opt-in intermediate layers before enabling them** — `adata.layers.get("norm_agg_counts")` returns `None` unless you passed `agg_norm_layer_name="norm_agg_counts"` to `ebsd.pp.normalize`. Use `.get()` defensively and check for `None` before plotting.

---

