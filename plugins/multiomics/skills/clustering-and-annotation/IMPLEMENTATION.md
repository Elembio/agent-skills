# Clustering and Annotation — Implementation

Detailed implementation recipes for the [Clustering and Annotation skill](./SKILL.md).

## Table of Contents

1. [Run Leiden Clustering](#1-run-leiden-clustering)
2. [Multi-Resolution Sweep](#2-multi-resolution-sweep)
3. [Silhouette-Based Resolution Selection](#3-silhouette-based-resolution-selection)
4. [Full Marker Gene Identification](#4-full-marker-gene-identification)
5. [Marker Gene Visualization](#5-marker-gene-visualization)
6. [Manual Cluster Annotation](#6-manual-cluster-annotation)
7. [Cell Cycle Scoring](#7-cell-cycle-scoring)
8. [Gene Set Scoring](#8-gene-set-scoring)
9. [Cluster Composition Analysis](#9-cluster-composition-analysis)
10. [Export](#10-export)
11. [Parameters Reference and Expected Outputs](#11-parameters-reference-and-expected-outputs) → [REFERENCE.md](./REFERENCE.md)
12. [Common Pitfalls](#12-common-pitfalls)

---

## 1. Run Leiden Clustering

```python
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Ensure neighbor graph exists
if "neighbors" not in adata.uns:
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)

sc.tl.leiden(
    adata,
    resolution=0.5,        # Start with balanced default
    key_added="leiden",
    flavor="igraph", n_iterations=2, directed=False,  # fast, memory-lean backend — important at scale
    random_state=42,        # Set for reproducibility
)

print(f"Number of clusters: {adata.obs['leiden'].nunique()}")
print(f"Cluster sizes:\n{adata.obs['leiden'].value_counts()}")

sc.pl.umap(adata, color="leiden", legend_loc="on data")
```

---

## 2. Multi-Resolution Sweep

Test a range of resolutions to understand the clustering landscape before committing to one.

```python
resolutions = [0.1, 0.3, 0.5, 0.8, 1.0, 1.5]

for res in resolutions:
    key = f"leiden_{res}"
    sc.tl.leiden(adata, resolution=res, key_added=key,
                 flavor="igraph", n_iterations=2, directed=False, random_state=42)
    n_clusters = adata.obs[key].nunique()
    print(f"Resolution {res}: {n_clusters} clusters")

# Panel UMAP for all resolutions
fig, axes = plt.subplots(2, 3, figsize=(18, 12))
for ax, res in zip(axes.flat, resolutions):
    key = f"leiden_{res}"
    sc.pl.umap(
        adata, color=key, ax=ax, show=False,
        title=f"Resolution={res} ({adata.obs[key].nunique()} clusters)",
        legend_loc="on data",
    )
plt.tight_layout()
plt.savefig("figures/leiden_resolutions.pdf")
plt.show()
```

---

## 3. Silhouette-Based Resolution Selection

Compute silhouette score and cluster size statistics across resolutions to guide the final choice.

```python
from sklearn.metrics import silhouette_score

def evaluate_clustering(adata, key):
    """Return quality metrics for a given cluster column."""
    labels = adata.obs[key].astype("category").cat.codes
    coords = adata.obsm["X_umap"]

    sil = silhouette_score(coords, labels) if len(np.unique(labels)) > 1 else 0

    sizes = adata.obs[key].value_counts()
    return {
        "n_clusters": len(np.unique(labels)),
        "silhouette": sil,
        "size_cv": sizes.std() / sizes.mean(),
        "min_cluster_size": sizes.min(),
    }

resolutions = [0.1, 0.3, 0.5, 0.8, 1.0, 1.5]
results = []
for res in resolutions:
    key = f"leiden_{res}"
    metrics = evaluate_clustering(adata, key)
    metrics["resolution"] = res
    results.append(metrics)

results_df = pd.DataFrame(results)
print(results_df)

fig, axes = plt.subplots(1, 3, figsize=(15, 4))

axes[0].plot(results_df["resolution"], results_df["n_clusters"], "o-")
axes[0].set(xlabel="Resolution", ylabel="Number of Clusters", title="Cluster Count")

axes[1].plot(results_df["resolution"], results_df["silhouette"], "o-")
axes[1].set(xlabel="Resolution", ylabel="Silhouette Score", title="Clustering Quality")

axes[2].plot(results_df["resolution"], results_df["min_cluster_size"], "o-")
axes[2].set(xlabel="Resolution", ylabel="Min Cluster Size", title="Smallest Cluster")

plt.tight_layout()
plt.show()

# Choose the resolution where silhouette is high, min cluster size is reasonable,
# and the number of clusters aligns with biological expectations.
```

---

## 4. Full Marker Gene Identification

Run on `adata`, not `adata_hvg`. `adata.X` must be log-normalized here — not z-scaled.
For 3′ data this is guaranteed because only `adata_hvg` was scaled (see [dimensionality-reduction §2](../dimensionality-reduction/IMPLEMENTATION.md#2-scale-data)).
For targeted panels, restore from `layers["normalized"]` if `adata` was scaled in-place.

```python
sc.tl.rank_genes_groups(
    adata,
    groupby="leiden",
    method="wilcoxon",             # Recommended; 't-test' and 'logreg' also available
    use_raw=False,                 # adata.raw is not set by normalization-strategies; use adata.X
    pts=True,                      # Include percent-expressing columns (pts, pts_rest)
    min_fold_change=1.5,
    min_in_group_fraction=0.25,
    key_added="rank_genes_leiden",
)

# Overview panel
sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False, key="rank_genes_leiden")

# Full DataFrame — pts = fraction expressing in group, pts_rest = fraction in rest
# Use these directly; do not recompute percent-expressing from adata.X
marker_df = sc.get.rank_genes_groups_df(adata, group=None, key="rank_genes_leiden")
print(marker_df.head(20))

# Per-cluster top markers
for cluster in sorted(adata.obs["leiden"].unique()):
    top = marker_df[marker_df["group"] == cluster].head(5)
    print(f"\nCluster {cluster} top markers:")
    print(top[["names", "scores", "pvals_adj", "logfoldchanges", "pts", "pts_rest"]].to_string(index=False))
```

---

## 5. Marker Gene Visualization

### Dot Plot

```python
sc.pl.rank_genes_groups_dotplot(
    adata,
    n_genes=5,
    key="rank_genes_leiden",
    groupby="leiden",
)
```

### Heatmap

```python
sc.pl.rank_genes_groups_heatmap(
    adata,
    n_genes=5,
    key="rank_genes_leiden",
    groupby="leiden",
    show_gene_labels=True,
)
```

### Stacked Violin

```python
top_genes = marker_df.groupby("group").head(2)["names"].unique()[:12]
sc.pl.stacked_violin(adata, top_genes, groupby="leiden", rotation=90)
```

---

## 6. Manual Cluster Annotation

Annotation is experiment-specific and requires domain knowledge of expected cell types.

```python
cluster_annotation = {
    "0": "Cycling_cells",
    "1": "Stressed_cells",
    "2": "Quiescent",
    "3": "Apoptotic",
    # Extend with the clusters found in your data
}

adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_annotation).fillna("Unknown")

sc.pl.umap(adata, color="cell_type", legend_loc="right margin")
```

Tip: leave clusters labeled `"Unknown"` rather than forcing a label when marker evidence is ambiguous.

---

## 7. Cell Cycle Scoring

```python
# S/G2M gene lists from Tirosh et al. 2015 — see cell-cycle-analysis/REFERENCE.md for full lists
S_GENES = [...]   # 43 S-phase genes
G2M_GENES = [...]  # 54 G2M genes

s_genes   = [g for g in S_GENES   if g in adata.var_names]
g2m_genes = [g for g in G2M_GENES if g in adata.var_names]

print(f"S phase genes found:  {len(s_genes)}/{len(S_GENES)}")
print(f"G2M genes found:      {len(g2m_genes)}/{len(G2M_GENES)}")

if s_genes and g2m_genes:
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
    # Adds: adata.obs["S_score"], adata.obs["G2M_score"], adata.obs["phase"]
    print(adata.obs["phase"].value_counts())

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    sc.pl.umap(adata, color="S_score",   ax=axes[0], show=False, cmap="RdBu_r")
    sc.pl.umap(adata, color="G2M_score", ax=axes[1], show=False, cmap="RdBu_r")
    sc.pl.umap(adata, color="phase",     ax=axes[2], show=False)
    plt.tight_layout()
    plt.show()
```

See also: [Cell Cycle Analysis skill](../cell-cycle-analysis/SKILL.md) for deeper cell cycle workflows.

---

## 8. Gene Set Scoring

Score cells for arbitrary gene sets such as pathway activity signatures.

```python
gene_sets = {
    "apoptosis":     ["BAX", "BCL2", "CASP3", "CASP9"],
    "proliferation": ["MKI67", "PCNA", "CDK1", "CCNB1"],
}

for name, genes in gene_sets.items():
    genes_present = [g for g in genes if g in adata.var_names]
    if genes_present:
        sc.tl.score_genes(adata, genes_present, score_name=f"{name}_score")
```

Use at least 10-20 genes per set; smaller sets produce noisy scores.

---

## 9. Cluster Composition Analysis

Check whether clusters are driven by biology or by batch/well effects.

```python
import seaborn as sns

def analyze_cluster_composition(adata, cluster_col, condition_col):
    """Return and plot cluster proportions across conditions."""
    crosstab = pd.crosstab(
        adata.obs[cluster_col],
        adata.obs[condition_col],
        normalize="columns",
    )

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    crosstab.T.plot(kind="bar", stacked=True, ax=axes[0], colormap="tab20")
    axes[0].set(xlabel=condition_col, ylabel="Proportion",
                title="Cluster Composition by Condition")
    axes[0].legend(title=cluster_col, bbox_to_anchor=(1.05, 1))

    sns.heatmap(crosstab, annot=True, fmt=".2f", cmap="YlOrRd", ax=axes[1])
    axes[1].set_title("Cluster Proportions")

    plt.tight_layout()
    return crosstab

if "experimental_condition" in adata.obs.columns:
    composition = analyze_cluster_composition(adata, "leiden", "experimental_condition")
    plt.show()
```

If one cluster is dominated by a single well, consider [Batch Effect Assessment](../batch-effect-assessment/SKILL.md).

---

## 10. Export

```python
# Cluster assignments
cols = ["leiden", "cell_type"]
if "phase" in adata.obs.columns:
    cols.append("phase")
adata.obs[cols].to_csv("results/cluster_assignments.csv")

# Marker genes table
marker_df = sc.get.rank_genes_groups_df(adata, group=None, key="rank_genes_leiden")
marker_df.to_csv("results/cluster_markers.csv", index=False)

# Per-cluster summary statistics
cluster_summary = adata.obs.groupby("leiden").agg(
    n_cells=("leiden", "count"),
    mean_area=("Area", "mean"),
    std_area=("Area", "std"),
).reset_index()
print(cluster_summary)
cluster_summary.to_csv("results/cluster_summary.csv", index=False)
```

---

## 11. Parameters Reference and Expected Outputs

See [REFERENCE.md](./REFERENCE.md) — `sc.tl.leiden` and `sc.tl.rank_genes_groups` parameter tables, expected outputs, and validation checks.

---

## 12. Common Pitfalls

1. **Resolution too high** — Too many clusters dilutes biological signal and produces singleton clusters. Start at 0.5, then adjust.

2. **Ignoring batch effects** — Clusters may reflect well-of-origin rather than biology. Always run cluster composition analysis (Recipe 9) and compare against known batch variables.

3. **Over-annotating** — Not every cluster needs a distinct biological label. Using `"Unknown"` is scientifically valid and preferable to a speculative label.

4. **Small gene sets** — Sets with fewer than 10-20 genes produce noisy `score_genes` output. Prefer curated pathway databases over ad-hoc small lists.

5. **Forgetting `random_state`** — Leiden is non-deterministic by default. Always set `random_state=42` (or any fixed value) for reproducible cluster labels.

6. **Using UMAP for silhouette** — UMAP coordinates are a non-linear embedding; silhouette computed on them is an approximation. Use it for relative comparison across resolutions, not as an absolute quality measure.
