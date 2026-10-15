# Pseudobulk DE — Implementation

Full implementation reference for pseudobulk differential expression with PyDESeq2.

## Table of Contents

1. [Imports](#1-imports)
2. [Raw Layer Handling](#2-raw-layer-handling)
3. [Define Experimental Design](#3-define-experimental-design)
4. [create_pseudobulk — Full Implementation](#4-create_pseudobulk--full-implementation)
5. [Run PyDESeq2](#5-run-pydeseq2)
6. [Extract Results and Contrasts](#6-extract-results-and-contrasts)
7. [Meso-Level Loop with run_meso_level_de](#7-meso-level-loop-with-run_meso_level_de)
8. [Covariates (Batch Correction)](#8-covariates-batch-correction)
9. [Save Results](#9-save-results)
10. [Validation Plots](#10-validation-plots)
11. [Parameters Table](#11-parameters-table)
12. [Expected Outputs](#12-expected-outputs)
13. [Common Pitfalls](#13-common-pitfalls)

---

## 1. Imports

```python
import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
```

---

## 2. Raw Layer Handling

DESeq2 requires raw (un-normalized) integer counts. Always extract from `layers["raw"]`.

```python
# Start from your filtered, annotated AnnData
adata = filtered_normalized_adata.copy()

# Prefer layers["raw"]; fall back only if the layer is genuinely absent
if "raw" in adata.layers:
    counts = adata.layers["raw"]
elif "raw_counts" in adata.layers:
    counts = adata.layers["raw_counts"]
else:
    counts = adata.X  # only if .X was never normalized

# Densify sparse matrix
if sp.issparse(counts):
    counts = counts.toarray()

# DESeq2 requires integers — round if fractional
if not np.allclose(counts, counts.astype(int)):
    print("Warning: non-integer counts detected — rounding.")
    counts = np.round(counts).astype(int)

adata.layers["raw"] = counts
print(f"Data shape: {adata.shape}")
print(f"Wells: {adata.obs['Well'].nunique()}")
```

---

## 3. Define Experimental Design

```python
# Option A: condition already in obs
# adata.obs["experimental_condition"] must exist

# Option B: map from well labels
if "experimental_condition" not in adata.obs.columns:
    well_to_condition = {
        "A1": "Control",
        "A2": "Control",
        "A3": "Control",
        "B1": "Treatment",
        "B2": "Treatment",
        "B3": "Treatment",
    }
    adata.obs["experimental_condition"] = adata.obs["Well"].map(well_to_condition)

# Verify design
design_df = (
    adata.obs
    .groupby(["experimental_condition", "Well"])
    .size()
    .reset_index(name="n_cells")
)
print("Experimental design:")
print(design_df)

replicates_per_condition = design_df.groupby("experimental_condition")["Well"].nunique()
print(f"\nReplicates per condition:\n{replicates_per_condition}")

if (replicates_per_condition < 2).any():
    print("WARNING: some conditions have < 2 replicates. DE will not be reliable.")
```

---

## 4. `create_pseudobulk` — Full Implementation

Aggregates cells to pseudobulk samples (samples × genes).

```python
def create_pseudobulk(adata, sample_col, condition_col,
                      min_cells=10, layer="raw"):
    """
    Aggregate single-cell data to pseudobulk.

    Parameters
    ----------
    adata : AnnData
        Single-cell data. Raw counts must be in `adata.layers[layer]`.
    sample_col : str
        obs column identifying biological samples (e.g. "Well").
    condition_col : str
        obs column identifying experimental conditions.
    min_cells : int
        Minimum cells a sample must contribute; samples below threshold are dropped.
    layer : str
        Layer containing raw counts. Defaults to "raw".

    Returns
    -------
    AnnData
        Pseudobulk object: one row per sample (samples × genes).
    """
    X = adata.layers[layer] if layer in adata.layers else adata.X
    if sp.issparse(X):
        X = X.toarray()

    samples = adata.obs[sample_col].unique()
    pb_counts, pb_samples, pb_conditions, pb_n_cells = [], [], [], []

    for sample in samples:
        mask = adata.obs[sample_col] == sample
        n_cells = mask.sum()

        if n_cells < min_cells:
            print(f"  Dropping {sample}: {n_cells} cells < min_cells={min_cells}")
            continue

        pb_counts.append(X[mask].sum(axis=0))
        pb_samples.append(sample)
        pb_conditions.append(adata.obs.loc[mask, condition_col].iloc[0])
        pb_n_cells.append(n_cells)

    pb_adata = ad.AnnData(
        X=np.vstack(pb_counts),
        obs=pd.DataFrame({
            sample_col: pb_samples,
            condition_col: pb_conditions,
            "n_cells": pb_n_cells,
        }),
        var=adata.var.copy(),
    )
    pb_adata.obs.index = pb_adata.obs[sample_col].values
    return pb_adata


# Macro-level pseudobulk
pseudobulk = create_pseudobulk(
    adata,
    sample_col="Well",
    condition_col="experimental_condition",
    min_cells=10,
    layer="raw",
)

print(f"Pseudobulk shape: {pseudobulk.shape}")
print(pseudobulk.obs)
```

---

## 5. Run PyDESeq2

```python
counts_df = pd.DataFrame(
    pseudobulk.X.astype(int).clip(min=0),
    index=pseudobulk.obs_names,
    columns=pseudobulk.var_names,
)

metadata = pseudobulk.obs.copy()

dds = DeseqDataSet(
    counts=counts_df,
    metadata=metadata,
    design_factors="experimental_condition",
    refit_cooks=True,   # refit outlier genes (recommended)
    min_mu=0.5,         # minimum mean for shrinkage
)

dds.deseq2()

print("DESeq2 complete")
print(f"Size factors: {dds.obsm['size_factors']}")
```

---

## 6. Extract Results and Contrasts

### Single contrast

```python
stat_res = DeseqStats(
    dds,
    contrast=["experimental_condition", "Treatment", "Control"],
)
stat_res.summary()

# Optional: apply LFC shrinkage (recommended for noisy data)
stat_res.lfc_shrink()

results_df = stat_res.results_df.sort_values("padj")
print(results_df.head(20))

# Significant DEGs
significant = results_df[
    (results_df["padj"] < 0.05) &
    (np.abs(results_df["log2FoldChange"]) > 1)
].copy()

print(f"Significant DEGs (padj<0.05, |log2FC|>1): {len(significant)}")
print(f"  Upregulated:   {(significant['log2FoldChange'] > 0).sum()}")
print(f"  Downregulated: {(significant['log2FoldChange'] < 0).sum()}")
```

### Multiple pairwise contrasts

```python
conditions = pseudobulk.obs["experimental_condition"].unique()
all_results = {}

for i, cond1 in enumerate(conditions):
    for cond2 in conditions[i + 1:]:
        key = f"{cond1}_vs_{cond2}"
        sr = DeseqStats(dds, contrast=["experimental_condition", cond1, cond2])
        sr.summary()
        all_results[key] = sr.results_df.copy()
        n_sig = (sr.results_df["padj"] < 0.05).sum()
        print(f"{key}: {n_sig} significant genes")
```

---

## 7. Meso-Level Loop with `run_meso_level_de`

Runs DE independently for each cluster. Requires `create_pseudobulk` from step 4.

```python
def run_meso_level_de(adata, cluster_col, sample_col, condition_col,
                      min_cells=10, layer="raw"):
    """
    Run pseudobulk DE per cluster (meso level).

    Parameters
    ----------
    adata : AnnData
        Full single-cell dataset with cluster labels.
    cluster_col : str
        obs column containing cluster/cell-type labels (e.g. "leiden").
    sample_col : str
        obs column identifying biological samples (e.g. "Well").
    condition_col : str
        obs column identifying experimental conditions.
    min_cells : int
        Minimum cells per sample within a cluster to include that sample.
    layer : str
        Layer containing raw counts.

    Returns
    -------
    dict
        Keys are cluster labels; values are dicts with keys:
        "results" (DataFrame), "n_cells_total", "n_samples", "n_sig".
    """
    meso_results = {}

    for cluster in adata.obs[cluster_col].unique():
        print(f"\n--- Cluster: {cluster} ---")
        cl_adata = adata[adata.obs[cluster_col] == cluster].copy()

        # Filter samples with too few cells
        cells_per_sample = cl_adata.obs.groupby(sample_col).size()
        valid_samples = cells_per_sample[cells_per_sample >= min_cells].index

        if len(valid_samples) < 4:
            print(f"  Skip: only {len(valid_samples)} valid samples (need ≥4)")
            continue

        cl_adata = cl_adata[cl_adata.obs[sample_col].isin(valid_samples)]

        pb = create_pseudobulk(
            cl_adata,
            sample_col=sample_col,
            condition_col=condition_col,
            min_cells=min_cells,
            layer=layer,
        )

        if pb.n_obs < 4:
            print(f"  Skip: {pb.n_obs} pseudobulk samples (need ≥4)")
            continue

        cond_counts = pb.obs[condition_col].value_counts()
        if (cond_counts < 2).any():
            print(f"  Skip: fewer than 2 replicates in at least one condition")
            continue

        try:
            counts_df = pd.DataFrame(
                pb.X.astype(int).clip(min=0),
                index=pb.obs_names,
                columns=pb.var_names,
            )
            dds = DeseqDataSet(
                counts=counts_df,
                metadata=pb.obs.copy(),
                design_factors=condition_col,
            )
            dds.deseq2()

            conds = pb.obs[condition_col].unique()
            sr = DeseqStats(dds, contrast=[condition_col, conds[0], conds[1]])
            sr.summary()

            meso_results[cluster] = {
                "results": sr.results_df.copy(),
                "n_cells_total": len(cl_adata),
                "n_samples": pb.n_obs,
                "n_sig": (sr.results_df["padj"] < 0.05).sum(),
            }
            print(f"  Significant DEGs: {meso_results[cluster]['n_sig']}")

        except Exception as exc:
            print(f"  Error: {exc}")

    return meso_results


# Run meso-level DE (requires leiden or equivalent cluster column)
if "leiden" in adata.obs.columns:
    meso_results = run_meso_level_de(
        adata,
        cluster_col="leiden",
        sample_col="Well",
        condition_col="experimental_condition",
        min_cells=10,
        layer="raw",
    )
```

---

## 8. Covariates (Batch Correction)

Add covariates to the `design_factors` list to control for batch effects or other nuisance variables. The primary contrast of interest is still `experimental_condition`.

```python
# Derive batch from well row (example only — adapt to your metadata)
pseudobulk.obs["batch"] = pseudobulk.obs["Well"].str[0]

counts_df = pd.DataFrame(
    pseudobulk.X.astype(int).clip(min=0),
    index=pseudobulk.obs_names,
    columns=pseudobulk.var_names,
)

dds_with_batch = DeseqDataSet(
    counts=counts_df,
    metadata=pseudobulk.obs.copy(),
    design_factors=["batch", "experimental_condition"],  # batch modeled first
)
dds_with_batch.deseq2()

# Contrast is still on experimental_condition
stat_res_batch = DeseqStats(
    dds_with_batch,
    contrast=["experimental_condition", "Treatment", "Control"],
)
stat_res_batch.summary()
```

> **Note:** Batch and condition must not be perfectly confounded (e.g., all "batch A" wells are "Control"). If they are, the model is not identifiable and will fail.

---

## 9. Save Results

```python
os.makedirs("results/de", exist_ok=True)
os.makedirs("figures", exist_ok=True)

# Macro-level
results_df.to_csv("results/de/macro_de_results.csv")
significant.to_csv("results/de/macro_de_significant.csv")

# Meso-level
for cluster, data in meso_results.items():
    safe_name = str(cluster).replace(" ", "_").replace("/", "-")
    data["results"].to_csv(f"results/de/meso_{safe_name}_de_results.csv")

# Summary table
summary_rows = [
    {"level": "macro", "cluster": "all",
     "n_tested": len(results_df), "n_significant": len(significant)}
]
for cluster, data in meso_results.items():
    summary_rows.append({
        "level": "meso",
        "cluster": cluster,
        "n_tested": len(data["results"]),
        "n_significant": data["n_sig"],
    })

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv("results/de/de_summary.csv", index=False)
print(summary_df)
```

---

## 10. Validation Plots

Visualize top DE genes as boxplots across conditions using pseudobulk log-normalized counts.

```python
def validate_de_results(adata, results_df, sample_col, condition_col,
                        n_genes=5, layer="raw"):
    """
    Plot top up- and down-regulated genes as per-sample boxplots.

    Parameters
    ----------
    adata : AnnData
        Full single-cell AnnData (used to re-create pseudobulk for plotting).
    results_df : DataFrame
        DE results from DeseqStats.results_df.
    sample_col : str
        obs column for biological samples.
    condition_col : str
        obs column for conditions.
    n_genes : int
        Number of top up- and down-regulated genes to plot.
    layer : str
        Layer containing raw counts for pseudobulk aggregation.
    """
    sig = results_df[results_df["padj"] < 0.05].copy()
    top_up   = sig[sig["log2FoldChange"] > 0].head(n_genes).index.tolist()
    top_down = sig[sig["log2FoldChange"] < 0].head(n_genes).index.tolist()
    top_genes = top_up + top_down

    if not top_genes:
        print("No significant genes to validate.")
        return

    pb = create_pseudobulk(adata, sample_col, condition_col, layer=layer)
    log_counts = np.log1p(pb.X)

    n_cols = min(3, len(top_genes))
    n_rows = int(np.ceil(len(top_genes) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
    axes = np.array(axes).flatten()

    for ax, gene in zip(axes, top_genes):
        gene_idx = pb.var_names.get_loc(gene)
        plot_df = pd.DataFrame({
            "Expression (log1p)": log_counts[:, gene_idx],
            "Condition": pb.obs[condition_col].values,
        })
        sns.boxplot(data=plot_df, x="Condition", y="Expression (log1p)", ax=ax)
        sns.swarmplot(data=plot_df, x="Condition", y="Expression (log1p)",
                      color="black", size=4, ax=ax)
        stats = results_df.loc[gene]
        ax.set_title(
            f"{gene}\nlog2FC={stats['log2FoldChange']:.2f}, "
            f"padj={stats['padj']:.2e}"
        )

    for ax in axes[len(top_genes):]:
        ax.axis("off")

    plt.tight_layout()
    plt.savefig("figures/de_validation.pdf")
    plt.show()


validate_de_results(adata, results_df, "Well", "experimental_condition")
```

---

## 11. Parameters Table

### PyDESeq2 parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `design_factors` | required | Factor(s) for the design matrix. String or list (list = multiple covariates). |
| `refit_cooks` | `True` | Refit genes flagged as Cook's distance outliers. Recommended. |
| `min_mu` | `0.5` | Minimum fitted mean for shrinkage estimation. |
| `contrast` | required | `[factor, numerator, denominator]` — e.g. `["condition", "Treatment", "Control"]`. |

### Significance thresholds

| Stringency | `padj` | `\|log2FC\|` | Typical use |
|------------|--------|------------|-------------|
| Strict | < 0.01 | > 1.5 | Publication, follow-up validation |
| Moderate | < 0.05 | > 1.0 | Standard exploratory analysis |
| Lenient | < 0.10 | > 0.5 | Hypothesis generation, small experiments |

```python
thresholds = {
    "strict":   {"padj": 0.01, "lfc": 1.5},
    "moderate": {"padj": 0.05, "lfc": 1.0},
    "lenient":  {"padj": 0.10, "lfc": 0.5},
}

for name, t in thresholds.items():
    n = ((results_df["padj"] < t["padj"]) &
         (np.abs(results_df["log2FoldChange"]) > t["lfc"])).sum()
    print(f"{name}: {n} DEGs")
```

### `create_pseudobulk` parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `adata` | AnnData | required | Single-cell dataset |
| `sample_col` | str | required | obs column for biological samples |
| `condition_col` | str | required | obs column for conditions |
| `min_cells` | int | `10` | Minimum cells per sample; samples below are excluded |
| `layer` | str | `"raw"` | Layer with raw counts |

---

## 12. Expected Outputs

| Object | Description |
|--------|-------------|
| `pseudobulk` | AnnData (samples × genes) with summed raw counts |
| `results_df` | DataFrame: `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj` |
| `significant` | Filtered `results_df` by chosen padj + log2FC thresholds |
| `meso_results` | Dict keyed by cluster; each value contains `results`, `n_cells_total`, `n_samples`, `n_sig` |

Saved files:

| Path | Contents |
|------|----------|
| `results/de/macro_de_results.csv` | All genes, macro-level |
| `results/de/macro_de_significant.csv` | Significant DEGs only |
| `results/de/meso_{cluster}_de_results.csv` | Per-cluster results |
| `results/de/de_summary.csv` | Summary counts across levels |
| `figures/de_validation.pdf` | Boxplots of top DE genes |

---

## 13. Common Pitfalls

1. **Using normalized counts.** DESeq2 performs its own normalization internally. Feeding log-normalized or scaled values produces incorrect size factor estimates and inflated significance. Always use `layers["raw"]`.

2. **Insufficient replicates.** Fewer than 2 replicates per condition makes dispersion estimation impossible; fewer than 3 makes it unreliable. The meso-level loop skips clusters that fail this check.

3. **Too few cells per pseudobulk sample.** Pseudobulk aggregates with < 10 cells are dominated by noise. Increase `min_cells` if you have very small clusters.

4. **Batch–condition confounding.** If all wells from batch A are "Control" and batch B are "Treatment", the batch covariate cannot be separated from condition. The model will fail or produce misleading results. Check your design matrix before running.

5. **Using raw p-values.** Always filter on `padj` (Benjamini–Hochberg adjusted). Raw `pvalue` columns lead to a high false discovery rate.

6. **Cluster labels with special characters.** Slashes or spaces in cluster names break file paths. The save-results recipe uses `str.replace` to sanitize names — ensure cluster labels are clean upstream.
