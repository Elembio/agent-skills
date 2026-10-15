# Batch Effect Assessment — Implementation

## Table of Contents

1. [Prepare Data for Assessment](#1-prepare-data-for-assessment)
2. [PCA Coloring by Technical Variables](#2-pca-coloring-by-technical-variables)
3. [UMAP Coloring by Technical and Biological Variables](#3-umap-coloring-by-technical-and-biological-variables)
4. [Variance Explained Calculation](#4-variance-explained-calculation)
5. [Silhouette Scores by Well and Condition](#5-silhouette-scores-by-well-and-condition)
6. [Replicate Pseudobulk Correlation Heatmap](#6-replicate-pseudobulk-correlation-heatmap)
7. [Harmony Integration](#7-harmony-integration)
8. [ComBat Expression Correction](#8-combat-expression-correction)
9. [Batch Report Generator](#9-batch-report-generator)
10. [Parameters and Thresholds](#10-parameters-and-thresholds)
11. [Expected Outputs](#11-expected-outputs)
12. [Common Pitfalls](#12-common-pitfalls)

---

## 1. Prepare Data for Assessment

```python
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Start with normalized data
adata = normalized_adata.copy()

# Resolve tile column name (tile_id or Tile, depending on the export)
_tile = "tile_id" if "tile_id" in adata.obs.columns else "Tile"

print("Technical variables available:")
tech_vars = ["Well", _tile, "run_id", "run_name", "batch_name"]
for var in tech_vars:
    if var in adata.obs.columns:
        print(f"  {var}: {adata.obs[var].nunique()} unique values")

if "experimental_condition" not in adata.obs.columns:
    print("\nNote: Add experimental_condition to .obs for proper assessment")
```

---

## 2. PCA Coloring by Technical Variables

```python
_tile = "tile_id" if "tile_id" in adata.obs.columns else "Tile"

sc.pp.pca(adata, n_comps=50, svd_solver="arpack")

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

sc.pl.pca(adata, color="Well", ax=axes[0, 0], show=False, title="PCA by Well")

if _tile in adata.obs.columns and adata.obs[_tile].nunique() <= 20:
    sc.pl.pca(adata, color=_tile, ax=axes[0, 1], show=False, title=f"PCA by {_tile}")
else:
    n_tiles = adata.obs[_tile].nunique() if _tile in adata.obs.columns else 0
    axes[0, 1].text(0.5, 0.5, f"Too many tiles ({n_tiles})",
                    ha="center", va="center", transform=axes[0, 1].transAxes)
    axes[0, 1].set_title("PCA by tile")

if "experimental_condition" in adata.obs.columns:
    sc.pl.pca(adata, color="experimental_condition", ax=axes[1, 0],
              show=False, title="PCA by Condition")
else:
    axes[1, 0].text(0.5, 0.5, "No condition metadata",
                    ha="center", va="center", transform=axes[1, 0].transAxes)

# Variance scree plot
ax = axes[1, 1]
variance_ratio = adata.uns["pca"]["variance_ratio"]
ax.bar(range(1, 21), variance_ratio[:20] * 100)
ax.set_xlabel("PC")
ax.set_ylabel("Variance Explained (%)")
ax.set_title("Variance Explained by PC")

plt.tight_layout()
plt.show()
```

---

## 3. UMAP Coloring by Technical and Biological Variables

```python
_tile = "tile_id" if "tile_id" in adata.obs.columns else "Tile"

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)

fig, axes = plt.subplots(2, 3, figsize=(18, 12))

sc.pl.umap(adata, color="Well", ax=axes[0, 0], show=False, title="UMAP by Well")

if "run_id" in adata.obs.columns:
    sc.pl.umap(adata, color="run_id", ax=axes[0, 1], show=False, title="UMAP by Run")
else:
    axes[0, 1].axis("off")

# Log total counts from the raw-count layer
raw_layer = "raw" if "raw" in adata.layers else "raw_counts"
totals = np.asarray(adata.layers[raw_layer].sum(axis=1)).ravel()
adata.obs["log_total_counts"] = np.log1p(totals)
sc.pl.umap(adata, color="log_total_counts", ax=axes[0, 2], show=False,
           title="UMAP by Total Counts (raw)")

if "experimental_condition" in adata.obs.columns:
    sc.pl.umap(adata, color="experimental_condition", ax=axes[1, 0],
               show=False, title="UMAP by Condition")

if "cell_type" in adata.obs.columns:
    sc.pl.umap(adata, color="cell_type", ax=axes[1, 1], show=False,
               title="UMAP by Cell Type")

if "Area" in adata.obs.columns:
    sc.pl.umap(adata, color="Area", ax=axes[1, 2], show=False, title="UMAP by Cell Area")

plt.tight_layout()
plt.show()
```

---

## 4. Variance Explained Calculation

```python
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import LabelEncoder, OneHotEncoder


def variance_explained_by_factor(adata, factor_col, n_pcs=30):
    """Calculate mean R² of a categorical factor across the top PCs."""
    if factor_col not in adata.obs.columns:
        return None

    X = adata.obsm["X_pca"][:, :n_pcs]

    le = LabelEncoder()
    y = le.fit_transform(adata.obs[factor_col].astype(str))

    ohe = OneHotEncoder(sparse_output=False, drop="first")
    y_onehot = ohe.fit_transform(y.reshape(-1, 1))

    r2_values = []
    for i in range(X.shape[1]):
        lr = LinearRegression()
        lr.fit(y_onehot, X[:, i])
        r2_values.append(lr.score(y_onehot, X[:, i]))

    return np.mean(r2_values)


_tile = "tile_id" if "tile_id" in adata.obs.columns else "Tile"
factors = ["Well", _tile, "experimental_condition"]
for extra in ("run_id", "run_name"):
    if extra in adata.obs.columns:
        factors.append(extra)

print("Variance explained by technical/biological factors (top 30 PCs):")
results = {}
for factor in factors:
    if factor in adata.obs.columns:
        r2 = variance_explained_by_factor(adata, factor)
        results[factor] = r2
        print(f"  {factor}: {r2:.1%}")

# Bar chart
fig, ax = plt.subplots(figsize=(8, 5))
factors_present = [f for f in factors if f in results]
values = [results[f] for f in factors_present]
tech_keys = {"Well", "Tile", "tile_id", "run_id", "run_name"}
colors = ["red" if f in tech_keys else "blue" for f in factors_present]
ax.barh(factors_present, values, color=colors)
ax.set_xlabel("Variance Explained (R²)")
ax.set_title("Variance Explained by Factors\n(Red=Technical, Blue=Biological)")
plt.tight_layout()
plt.show()
```

---

## 5. Silhouette Scores by Well and Condition

```python
from sklearn.metrics import silhouette_score, adjusted_rand_score
from sklearn.preprocessing import LabelEncoder


def assess_technical_clustering(adata, tech_var, cluster_key="leiden"):
    """
    Silhouette score on UMAP coordinates measures separation of technical groups.
    High silhouette (> 0.3) suggests cells cluster by the technical variable.
    High ARI with leiden clusters suggests clusters align with batch, which is problematic.
    """
    if tech_var not in adata.obs.columns:
        return None

    umap_coords = adata.obsm["X_umap"]
    labels = LabelEncoder().fit_transform(adata.obs[tech_var].astype(str))
    n_groups = len(np.unique(labels))

    if n_groups < 2 or n_groups > 50:
        return {"n_groups": n_groups, "silhouette": None, "ari": None}

    silhouette = silhouette_score(umap_coords, labels)

    ari = None
    if cluster_key in adata.obs.columns:
        cluster_labels = LabelEncoder().fit_transform(adata.obs[cluster_key].astype(str))
        ari = adjusted_rand_score(labels, cluster_labels)

    return {"n_groups": n_groups, "silhouette": silhouette, "ari": ari}


print("Technical clustering assessment:")
print("High silhouette = cells cluster by technical variable = potential batch effect")
print("High ARI = clusters align with technical groups = problematic\n")

_tile = "tile_id" if "tile_id" in adata.obs.columns else "Tile"
for tech_var in ["Well", _tile, "run_id", "run_name"]:
    if tech_var not in adata.obs.columns:
        continue
    result = assess_technical_clustering(adata, tech_var)
    if result:
        print(f"{tech_var}:")
        print(f"  Groups: {result['n_groups']}")
        if result["silhouette"] is not None:
            print(f"  Silhouette: {result['silhouette']:.3f}")
        if result["ari"] is not None:
            print(f"  ARI with clusters: {result['ari']:.3f}")
```

---

## 6. Replicate Pseudobulk Correlation Heatmap

```python
def compare_replicates(adata, condition_col, replicate_col="Well"):
    """
    Compute pseudobulk profiles per (condition, replicate) pair and compare
    within-condition vs between-condition Pearson correlations.
    """
    if condition_col not in adata.obs.columns:
        print(f"Column '{condition_col}' not found")
        return None

    pseudobulk = {}
    for (cond, rep), group in adata.obs.groupby([condition_col, replicate_col]):
        cells = group.index
        X_sub = adata[cells].X
        if hasattr(X_sub, "toarray"):
            X_sub = X_sub.toarray()
        pseudobulk[(cond, rep)] = X_sub.mean(axis=0)

    pb_df = pd.DataFrame(pseudobulk).T
    pb_df.index = pd.MultiIndex.from_tuples(pb_df.index, names=[condition_col, replicate_col])

    conditions = adata.obs[condition_col].unique()
    within_corrs, between_corrs = [], []

    for cond in conditions:
        cond_samples = [k for k in pb_df.index if k[0] == cond]
        for i, s1 in enumerate(cond_samples):
            for s2 in cond_samples[i + 1:]:
                within_corrs.append(np.corrcoef(pb_df.loc[s1], pb_df.loc[s2])[0, 1])

    for i, cond1 in enumerate(conditions):
        for cond2 in conditions[i + 1:]:
            for s1 in [k for k in pb_df.index if k[0] == cond1]:
                for s2 in [k for k in pb_df.index if k[0] == cond2]:
                    between_corrs.append(np.corrcoef(pb_df.loc[s1], pb_df.loc[s2])[0, 1])

    return {
        "within_condition_corr": np.mean(within_corrs) if within_corrs else None,
        "between_condition_corr": np.mean(between_corrs) if between_corrs else None,
        "pseudobulk_df": pb_df,
    }


if "experimental_condition" in adata.obs.columns:
    rep_comparison = compare_replicates(adata, "experimental_condition", "Well")
    if rep_comparison and rep_comparison["within_condition_corr"] is not None:
        print(f"Within-condition correlation:  {rep_comparison['within_condition_corr']:.3f}")
        print(f"Between-condition correlation: {rep_comparison['between_condition_corr']:.3f}")
        print("(Higher within than between indicates good replicate consistency)")

        # Heatmap of all pseudobulk correlations
        pb_df = rep_comparison["pseudobulk_df"]
        corr_matrix = pb_df.T.corr()
        fig, ax = plt.subplots(figsize=(10, 8))
        sns.heatmap(corr_matrix, cmap="vlag", center=0, vmin=-1, vmax=1,
                    annot=True, fmt=".2f", ax=ax)
        ax.set_title("Pseudobulk Pearson Correlation (condition × well)")
        plt.tight_layout()
        plt.show()
```

---

## 7. Harmony Integration

Use Harmony when batch effects are present but you want to preserve expression values (correction lives in the embedding, not `.X`).

```python
try:
    import harmonypy as hm

    ho = hm.run_harmony(
        adata.obsm["X_pca"],
        adata.obs,
        vars_use=["Well"],   # variables to correct for
        max_iter_harmony=20,
    )

    adata.obsm["X_pca_harmony"] = ho.Z_corr.T

    # Rebuild neighbor graph and UMAP on the corrected embedding
    sc.pp.neighbors(adata, use_rep="X_pca_harmony")
    sc.tl.umap(adata)

    print("Harmony correction applied. Use X_pca_harmony for downstream steps.")
except ImportError:
    print("harmonypy not installed. Install with: pip install harmonypy")
```

---

## 8. ComBat Expression Correction

Use ComBat when downstream steps require corrected expression values in `.X` (e.g., differential expression without a batch covariate). Only apply when effects are severe.

```python
adata_combat = adata.copy()

# Modifies .X in place; requires batch variable in .obs
sc.pp.combat(adata_combat, key="Well")

print("ComBat correction applied. Corrected values are in adata_combat.X")

# Visual comparison
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
sc.pl.umap(adata, color="Well", ax=axes[0], show=False, title="Before ComBat")

sc.pp.pca(adata_combat, n_comps=50, svd_solver="arpack")
sc.pp.neighbors(adata_combat, n_pcs=30)
sc.tl.umap(adata_combat)
sc.pl.umap(adata_combat, color="Well", ax=axes[1], show=False, title="After ComBat")

plt.tight_layout()
plt.show()
```

---

## 9. Batch Report Generator

```python
def generate_batch_report(adata, tech_vars=None, bio_vars=None):
    """Generate a console summary report of batch effect metrics."""
    if tech_vars is None:
        _tile = "tile_id" if "tile_id" in adata.obs.columns else "Tile"
        tech_vars = [v for v in ["Well", _tile, "run_id", "run_name"]
                     if v in adata.obs.columns]
    if bio_vars is None:
        bio_vars = [v for v in ["experimental_condition"] if v in adata.obs.columns]

    report = {"n_cells": adata.n_obs, "n_genes": adata.n_vars}

    all_vars = tech_vars + bio_vars
    for var in all_vars:
        if var in adata.obs.columns:
            report[f"n_{var}"] = adata.obs[var].nunique()
            r2 = variance_explained_by_factor(adata, var)
            report[f"var_explained_{var}"] = r2

    for var in tech_vars:
        result = assess_technical_clustering(adata, var)
        if result and result["silhouette"] is not None:
            report[f"silhouette_{var}"] = result["silhouette"]

    print("=" * 50)
    print("BATCH EFFECT ASSESSMENT REPORT")
    print("=" * 50)
    print(f"Dataset: {report['n_cells']} cells x {report['n_genes']} genes\n")

    print("Groups:")
    for var in all_vars:
        if f"n_{var}" in report:
            print(f"  {var}: {report[f'n_{var}']} groups")

    print("\nVariance Explained (on PCA):")
    for var in all_vars:
        key = f"var_explained_{var}"
        if key in report:
            val = report[key]
            flag = " [!]" if var in tech_vars and val > 0.1 else ""
            print(f"  {var}: {val:.1%}{flag}")

    print("\nTechnical Clustering (silhouette; higher = more batch effect):")
    for var in tech_vars:
        key = f"silhouette_{var}"
        if key in report:
            val = report[key]
            flag = " [!]" if val > 0.3 else ""
            print(f"  {var}: {val:.3f}{flag}")

    print("=" * 50)
    return report


report = generate_batch_report(adata)
```

---

## 10. Parameters and Thresholds

| Metric | Low concern | Moderate | High concern |
|--------|-------------|----------|--------------|
| Technical factor R² (variance explained) | < 5% | 5–15% | > 15% |
| Technical silhouette score | < 0.1 | 0.1–0.3 | > 0.3 |
| Within/between replicate correlation ratio | > 1.2 | 1.0–1.2 | < 1.0 |

**PCA components used for variance calculation:** 30 (adjust `n_pcs` for datasets with fewer cells).

**Harmony `max_iter_harmony`:** 20 is sufficient for most datasets; increase to 50 for large, heterogeneous experiments.

**ComBat:** Avoid when conditions are fully confounded with wells (one condition per well). The model will fail or produce meaningless results.

---

## 11. Expected Outputs

After running these recipes you should have:

1. **Visualizations** — PCA and UMAP plots colored by `Well`, `tile_id`/`Tile`, `run_id`, total counts, cell type, and experimental condition.
2. **Variance metrics** — Mean R² per technical and biological factor across the top 30 PCs.
3. **Silhouette scores** — Per-well and per-run separation on UMAP coordinates.
4. **Replicate correlations** — Within- vs. between-condition pseudobulk Pearson correlations with a heatmap.
5. **A correction decision** — Whether to apply Harmony (embedding-level) or ComBat (expression-level).
6. **Corrected embeddings or expression** (if applicable) — `adata.obsm["X_pca_harmony"]` or `adata_combat.X`.
7. **Batch report** — Console summary with flagged metrics exceeding concern thresholds.

---

## 12. Common Pitfalls

**Confounded design** — If each biological condition occupies exactly one well, you cannot distinguish batch from biology. No correction method can recover this; redesign the experiment.

**Over-correction** — Aggressive correction (especially ComBat) can remove genuine biological signal. Always compare corrected and uncorrected UMAP/clustering side by side before proceeding.

**Wrong correction variable** — Correcting for `Well` when the true batch is `run_id` (or vice versa) may be ineffective. Check which variable has the highest R² and silhouette score first.

**Post-hoc correction on confounded data** — Batch correction cannot rescue fundamentally confounded experiments. Check design balance before applying any correction.

**`sparse_output` deprecation** — `OneHotEncoder(sparse=False, ...)` is deprecated in recent scikit-learn; use `sparse_output=False` as shown in Recipe 4.
