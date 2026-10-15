# 3′ Untargeted Pipeline — Implementation

Step-by-step implementation for the full 3′ untargeted analysis pipeline.

## Table of Contents

1. [Setup](#1-setup)
2. [Load 3′ Untargeted Data](#2-load-3-untargeted-data)
3. [Compute QC Metrics](#3-compute-qc-metrics)
4. [Visualize QC Metrics](#4-visualize-qc-metrics)
5. [Filter Cells and Genes](#5-filter-cells-and-genes)
6. [Normalization](#6-normalization)
7. [Highly Variable Gene Selection](#7-highly-variable-gene-selection)
8. [Scale and PCA](#8-scale-and-pca)
9. [Batch Effect Assessment](#9-batch-effect-assessment)
10. [Compute Embeddings](#10-compute-embeddings)
11. [Clustering](#11-clustering)
12. [Marker Gene Identification](#12-marker-gene-identification)
13. [Cell Type Annotation](#13-cell-type-annotation)
14. [Cell Cycle Scoring (Optional)](#14-cell-cycle-scoring-optional)
15. [Pseudobulk Differential Expression (Optional)](#15-pseudobulk-differential-expression-optional)
16. [Save Results](#16-save-results)
- [Expected Outputs](#expected-outputs)
- [Common Pitfalls](#common-pitfalls)

---

## 1. Setup

```python
import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

os.makedirs("figures", exist_ok=True)
os.makedirs("results", exist_ok=True)
os.makedirs("results/de", exist_ok=True)

plt.rcParams["figure.dpi"] = 100
sc.settings.verbosity = 2

print("Setup complete")
```

---

## 2. Load 3′ Untargeted Data

<!-- canonical: spatialdata-loading-and-access memory guard. Copy verbatim; edit only in that skill's SKILL.md -->
**Memory guard — before any `lazy=False`, `.to_memory()` or `.compute()`:**

1. Load lazily first and run `check_before` from `scripts/check_load.py` in `spatialdata-loading-and-access`. Run the script unchanged in the Python session that holds the data (`exec` its contents); don't rewrite it. It compares each table's estimated in-memory size with the memory this session can use, lists large session objects you could `del`, and raises on a FAIL.
2. Materialize only the tables and runs the plan needs. With several runs, load `concatenated=False` so it checks each run and their total. A WARN says what to `del` or to load one run at a time. A FAIL means subset before materializing: propose a subset to the user (runs, wells, or a set of features), then slice the lazy table and call `.to_memory()` on the slice, for example `adata[adata.obs["Well"].isin(wells)].to_memory()`. `.copy()` keeps a lazy table lazy. Feature subsets by modality: [ACCESS.md — Subset](../spatialdata-loading-and-access/ACCESS.md#subset-by-well-or-feature).
3. Keep `.X` in the form it was stored (`scipy.sparse.issparse(adata.X)`). Sparse tables (3′, OPS) are mostly zeros: never densify a whole one (`.toarray()`, `.todense()`, `np.asarray(adata.X)`); densify only a subset of cells or features. Dense tables (imaging, legacy barcoding) are float32 with `NaN` where a batch has no data: do not convert them to sparse, which turns missing into zero.
<!-- /canonical -->

```python
# See: spatialdata-loading-and-access skill
from elembio_spatialdata_analysis.load import Loader

STORE_PATH = "/path/to/run.zarr"
loader = Loader(STORE_PATH)

if "ThreePrimeUntargeted" not in loader.available_tables:
    raise ValueError("No ThreePrimeUntargeted table — check loader.available_tables")

check_before(loader.load_tables(table="ThreePrimeUntargeted", concatenated=False))  # memory guard: lazy, raises on FAIL
adata = loader.load_tables(table="ThreePrimeUntargeted", lazy=False)
print(f"3′ Untargeted data: {adata.shape}")
print(f"Genes detected: {adata.n_vars}")

# Full SpatialData (only if you need non-table elements):
# import spatialdata as sd
# sdata = sd.read_zarr(STORE_PATH)
# adata = sdata.tables["global-table-ThreePrimeUntargeted"].copy()
```

---

## 3. Compute QC Metrics

```python
# See: cell-quality-control skill
import elembio_spatialdata_analysis as ebsd

ebsd.pp.compute_qc_metrics(adata, "ThreePrimeUntargeted")
# Writes adata.obsm["qc_metrics"] — required before ebsd.pp.filter_cells / filter plots.
# For 3′, MT/ribo/HB metrics are filled via Scanpy if missing.

qc = adata.obsm["qc_metrics"]
print("\nQC metrics summary (obsm['qc_metrics']):")
print(qc.describe())

# Check for spatial artifacts
if "spatial_artifacts" in adata.obsm:
    art = adata.obsm["spatial_artifacts"]
    col = "any-ThreePrimeUntargeted"
    if col in art.columns:
        n_flagged = art[col].sum()
        print(f"Spatial artifacts: {n_flagged} / {adata.n_obs} cells flagged in {col}")
        # If flagged > 0, ask the user whether to opt in during filter_cells:
        #   ebsd.pp.filter_cells(adata, "ThreePrimeUntargeted", filter_spatial_artifacts=True)
```

---

## 4. Visualize QC Metrics

```python
qc = adata.obsm["qc_metrics"]

fig, axes = plt.subplots(1, 4, figsize=(20, 5))

# total_counts and n_genes_by_counts are always written for ThreePrimeUntargeted, so they are safe to index unguarded here.
# The qc_metrics column set is modality-specific, though — when reusing this across modalities, guard with `if col in qc.columns` (see cell-quality-control).
sns.violinplot(y=qc["total_counts"], ax=axes[0], inner="box")
axes[0].set_ylabel("Total Counts")

sns.violinplot(y=qc["n_genes_by_counts"], ax=axes[1], inner="box")
axes[1].set_ylabel("Genes Detected")

if "pct_counts_mt" in qc.columns:
    sns.violinplot(y=qc["pct_counts_mt"], ax=axes[2], inner="box")
    axes[2].set_ylabel("MT%")
else:
    axes[2].set_visible(False)

ax = axes[3]
scat = ax.scatter(
    qc["total_counts"],
    qc["n_genes_by_counts"],
    c=qc["pct_counts_mt"] if "pct_counts_mt" in qc.columns else "C0",
    cmap="RdYlGn_r",
    s=1,
    alpha=0.5,
)
ax.set_xlabel("Total Counts")
ax.set_ylabel("Genes Detected")
ax.set_title("Counts vs Genes (colored by MT%)" if "pct_counts_mt" in qc.columns else "Counts vs Genes")
if "pct_counts_mt" in qc.columns:
    plt.colorbar(scat, ax=ax, label="MT%")

plt.tight_layout()
plt.savefig("figures/qc_metrics.pdf")
plt.show()
```

---

## 5. Filter Cells and Genes

```python

n_before = adata.n_obs

# max_pct_mt is on by default at 0.10 for 3′. Inspect the MT% distribution before overriding it.
adata_filtered, diagnostics = ebsd.pp.filter_cells(
    adata,
    "ThreePrimeUntargeted",
    diagnostics=True,
    inplace=False,
)

print(f"Before: {n_before}")
print(f"After:  {adata_filtered.n_obs}")
print(f"Removed: {n_before - adata_filtered.n_obs} ({(n_before - adata_filtered.n_obs) / n_before:.1%})")

filter_summary = (~diagnostics.flags).sum()
print("\nCells failing each filter:")
print(filter_summary.sort_values(ascending=False))

adata = adata_filtered

# Gene-level QC after cell filter (per-well logic)
adata = ebsd.pp.filter_genes(
    adata,
    "ThreePrimeUntargeted",
    min_cells_per_gene=10,
    inplace=False,
)
print(f"Genes remaining: {adata.n_vars}")
```

---

## 6. Normalization

```python
# See: normalization-strategies skill
# ebsd.pp.normalize is well-aware: it scales each cell to the per-well+batch
# median total counts, then applies median-of-ratios (DESeq2-style composition
# correction), then log1p. layers["raw"] is written automatically before any transformation.

ebsd.pp.normalize(adata, "ThreePrimeUntargeted")  # inplace=True by default

# Verify
assert "raw" in adata.layers, "raw layer missing — normalization may have failed"
print("Normalization complete")
print(f"Raw range  : {adata.layers['raw'].max():.0f}")
print(f".X range   : {adata.X.max():.2f}")
```

---

## 7. Highly Variable Gene Selection

```python
# See: dimensionality-reduction skill; REFERENCE.md for method comparison

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    flavor="seurat_v3",   # use "seurat" if working on log-normalized data
    batch_key="Well",     # account for per-well batch effects
    subset=False,
)

print(f"Highly variable genes selected: {adata.var['highly_variable'].sum()}")

sc.pl.highly_variable_genes(adata)
plt.savefig("figures/hvg_selection.pdf")
plt.show()

# Subset to HVGs for downstream dimensionality reduction
adata_hvg = adata[:, adata.var["highly_variable"]].copy()
print(f"HVG subset shape: {adata_hvg.shape}")
```

---

## 8. Scale and PCA

```python
sc.pp.scale(adata_hvg, max_value=10)
sc.pp.pca(adata_hvg, n_comps=50)

# Transfer PCA results to the full object
adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
adata.uns["pca"] = adata_hvg.uns["pca"]

# Scree plot
fig, ax = plt.subplots(figsize=(10, 5))
variance = adata.uns["pca"]["variance_ratio"]
ax.bar(range(1, 31), variance[:30] * 100)
ax.set_xlabel("PC")
ax.set_ylabel("Variance Explained (%)")
ax.set_title("Scree Plot")
plt.savefig("figures/pca_variance.pdf")
plt.show()

print(f"Variance explained by top 30 PCs: {sum(variance[:30]):.1%}")
```

---

## 9. Batch Effect Assessment

```python
# See: batch-effect-assessment skill

fig, axes = plt.subplots(1, 3, figsize=(18, 5))

sc.pl.pca(adata, color="Well", ax=axes[0], show=False, title="PCA by Well")
sc.pl.pca(adata, color="pct_counts_mt", ax=axes[1], show=False,
          title="PCA by MT%", cmap="RdYlGn_r")
sc.pl.pca(adata, color="total_counts", ax=axes[2], show=False,
          title="PCA by Total Counts", cmap="viridis")

plt.tight_layout()
plt.savefig("figures/pca_batch_assessment.pdf")
plt.show()

# Quantify batch mixing with NMI
from sklearn.metrics import normalized_mutual_info_score
from sklearn.cluster import KMeans

kmeans = KMeans(n_clusters=5, random_state=42)
adata.obs["_temp_cluster"] = kmeans.fit_predict(adata.obsm["X_pca"][:, :20])

nmi_well = normalized_mutual_info_score(
    adata.obs["_temp_cluster"].astype(str),
    adata.obs["Well"].astype(str),
)
print(f"NMI(cluster, Well): {nmi_well:.3f}  (>0.3 → consider batch correction)")

del adata.obs["_temp_cluster"]
```

---

## 10. Compute Embeddings

```python
# See: dimensionality-reduction skill

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata, min_dist=0.3)

print("Neighbor graph and UMAP computed")
```

---

## 11. Clustering

```python
# See: clustering-and-annotation skill

for res in [0.3, 0.5, 0.8, 1.0]:
    sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{res}",
                 flavor="igraph", n_iterations=2, directed=False)
    print(f"Resolution {res}: {adata.obs[f'leiden_{res}'].nunique()} clusters")

chosen_res = 0.5
adata.obs["leiden"] = adata.obs[f"leiden_{chosen_res}"]
print(f"\nUsing resolution {chosen_res}: {adata.obs['leiden'].nunique()} clusters")

fig, axes = plt.subplots(1, 3, figsize=(18, 5))
sc.pl.umap(adata, color="leiden", ax=axes[0], show=False,
           legend_loc="on data", title="Leiden Clusters")
sc.pl.umap(adata, color="Well", ax=axes[1], show=False, title="Well")
sc.pl.umap(adata, color="pct_counts_mt", ax=axes[2], show=False,
           cmap="RdYlGn_r", title="MT%")
plt.tight_layout()
plt.savefig("figures/clustering.pdf")
plt.show()
```

---

## 12. Marker Gene Identification

```python
sc.tl.rank_genes_groups(
    adata, groupby="leiden", method="wilcoxon",
    use_raw=False,   # adata.raw is not set; use adata.X (log-normalized)
    pts=True,        # fraction expressing in group (pts) and rest (pts_rest)
)

sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False)
plt.savefig("figures/marker_genes.pdf")
plt.show()

sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, groupby="leiden")
plt.savefig("figures/marker_dotplot.pdf")
plt.show()

markers = sc.get.rank_genes_groups_df(adata, group=None)
# pts = fraction expressing in cluster; pts_rest = fraction expressing in all other clusters
# Use these directly — do not recompute percent-expressing from adata.X
print(markers[["group", "names", "logfoldchanges", "pvals_adj", "pts", "pts_rest"]].head(20))
markers.to_csv("results/cluster_markers.csv", index=False)
print(f"Saved {len(markers)} marker gene results")
```

---

## 13. Cell Type Annotation

```python
# See: clustering-and-annotation skill; REFERENCE.md for reference database options.
# Annotation is experiment-specific — customize based on your marker genes.

marker_to_type = {
    "KRT19": "Epithelial",
    "VIM": "Mesenchymal",
    "PTPRC": "Immune",   # CD45
    "CDH1": "Epithelial",
    "MKI67": "Proliferating",
}

known_markers = [m for m in marker_to_type if m in adata.var_names]
if known_markers:
    sc.pl.dotplot(adata, known_markers, groupby="leiden")
    plt.savefig("figures/known_markers.pdf")
    plt.show()

# Map cluster IDs to cell type labels
cluster_annotation = {
    "0": "Type_A",
    "1": "Type_B",
    # extend per your analysis...
}

adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_annotation).fillna("Unknown")
print(f"Cell type counts:\n{adata.obs['cell_type'].value_counts()}")
```

---

## 14. Cell Cycle Scoring (Optional)

```python
# S/G2M gene lists from Tirosh et al. 2015 — see cell-cycle-analysis/REFERENCE.md for full lists
S_GENES = [...]   # 43 S-phase genes
G2M_GENES = [...]  # 54 G2M genes

s_genes = [g for g in S_GENES if g in adata.var_names]
g2m_genes = [g for g in G2M_GENES if g in adata.var_names]

print(f"S-phase genes found: {len(s_genes)}/{len(S_GENES)}")
print(f"G2M genes found:     {len(g2m_genes)}/{len(G2M_GENES)}")

if s_genes and g2m_genes:
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    sc.pl.umap(adata, color="phase", ax=axes[0], show=False, title="Cell Cycle Phase")

    phase_by_cluster = pd.crosstab(adata.obs["leiden"], adata.obs["phase"], normalize="index")
    phase_by_cluster.plot(kind="bar", stacked=True, ax=axes[1])
    axes[1].set_ylabel("Proportion")
    axes[1].set_title("Cell Cycle by Cluster")

    plt.tight_layout()
    plt.savefig("figures/cell_cycle.pdf")
    plt.show()
```

---

## 15. Pseudobulk Differential Expression (Optional)

```python
# See: pseudobulk-differential-expression skill
# Requires an "experimental_condition" column in adata.obs.

if "experimental_condition" in adata.obs.columns:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    import scipy.sparse as sp

    def create_pseudobulk(adata, sample_col, condition_col):
        raw_layer = "raw" if "raw" in adata.layers else "raw_counts"
        X = adata.layers[raw_layer]
        if sp.issparse(X):
            X = X.toarray()

        samples = adata.obs[sample_col].unique()
        pb_counts, pb_meta = [], []

        for sample in samples:
            mask = adata.obs[sample_col] == sample
            if mask.sum() < 10:
                continue
            pb_counts.append(X[mask].sum(axis=0))
            pb_meta.append({
                sample_col: sample,
                condition_col: adata.obs.loc[mask, condition_col].iloc[0],
            })

        return (
            pd.DataFrame(
                np.vstack(pb_counts),
                index=[m[sample_col] for m in pb_meta],
                columns=adata.var_names,
            ).astype(int),
            pd.DataFrame(pb_meta),
        )

    counts, meta = create_pseudobulk(adata, "Well", "experimental_condition")
    meta.index = meta["Well"]

    conditions = meta["experimental_condition"].unique()
    if len(conditions) >= 2:
        dds = DeseqDataSet(
            counts=counts.clip(lower=0),
            metadata=meta,
            design_factors="experimental_condition",
        )
        dds.deseq2()

        stat_res = DeseqStats(
            dds,
            contrast=["experimental_condition", conditions[0], conditions[1]],
        )
        stat_res.summary()

        de_results = stat_res.results_df
        sig = de_results[
            (de_results["padj"] < 0.05) & (np.abs(de_results["log2FoldChange"]) > 1)
        ]
        print(f"Significant DE genes: {len(sig)}")

        de_results.to_csv("results/de/3prime_de_results.csv")
```

---

## 16. Save Results

```python
adata.write("results/3prime_processed.h5ad")

summary = f"""
# 3′ Untargeted Analysis Summary

## Data
- Cells: {adata.n_obs}
- Genes (total): {adata.n_vars}
- HVGs: {adata.var['highly_variable'].sum()}
- Wells: {adata.obs['Well'].nunique()}

## Clustering
- Resolution: {chosen_res}
- Clusters: {adata.obs['leiden'].nunique()}

## Files
- results/3prime_processed.h5ad
- results/cluster_markers.csv
- figures/*.pdf
"""

print(summary)
with open("results/3prime_summary.md", "w") as f:
    f.write(summary)
```

---

## Expected Outputs

| File | Description |
|------|-------------|
| `results/3prime_processed.h5ad` | Fully processed AnnData with all layers, embeddings, and cluster labels |
| `results/cluster_markers.csv` | Wilcoxon marker genes per cluster |
| `results/de/3prime_de_results.csv` | DESeq2 results (if DE step run) |
| `results/3prime_summary.md` | Run summary with cell/gene counts and key parameters |
| `figures/qc_metrics.pdf` | Total counts, genes detected, MT% violin + scatter |
| `figures/hvg_selection.pdf` | HVG dispersion plot |
| `figures/pca_variance.pdf` | Scree plot |
| `figures/pca_batch_assessment.pdf` | PCA colored by Well / MT% / counts |
| `figures/clustering.pdf` | UMAP by cluster, well, and MT% |
| `figures/marker_genes.pdf` | Top 10 marker genes per cluster |
| `figures/marker_dotplot.pdf` | Dot plot of top 5 markers |
| `figures/cell_cycle.pdf` | Cell cycle phase UMAP and per-cluster bar chart (optional) |

---

## Common Pitfalls

| Pitfall | Resolution |
|---------|-----------|
| `obsm["qc_metrics"]` missing | Always run `ebsd.pp.compute_qc_metrics` before `ebsd.pp.filter_cells` or QC plots |
| Assuming 3′ MT filtering is off | It is on at `max_pct_mt=0.10` by default; pass a value only to override |
| HVG flavor mismatch | Use `"seurat_v3"` on raw counts; use `"seurat"` on log-normalized data |
| `batch_key` in `highly_variable_genes` fails | Ensure the column exists in `adata.obs` (commonly `"Well"`) |
| PCA computed on wrong object | Scale and PCA run on `adata_hvg`; transfer `X_pca` back to `adata` before `neighbors` |
| DE step fails with sparse matrix | `toarray()` before summing into pseudobulk counts |
| `chosen_res` undefined in save step | Assign `chosen_res` before the clustering loop or explicitly reassign before saving |
