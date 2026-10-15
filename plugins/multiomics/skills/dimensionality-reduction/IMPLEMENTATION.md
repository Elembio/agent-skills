# Dimensionality Reduction — Implementation

Detailed implementation steps for the dimensionality-reduction skill. For the concise overview and quick start, see [SKILL.md](./SKILL.md).

## Table of Contents

1. [HVG Selection (3′ untargeted only)](#1-hvg-selection-3-untargeted-only)
2. [Scale Data](#2-scale-data)
3. [Run PCA](#3-run-pca)
4. [Variance / Elbow Analysis](#4-variance--elbow-analysis)
5. [Compute Neighborhood Graph](#5-compute-neighborhood-graph)
6. [Run UMAP](#6-run-umap)
7. [UMAP Parameter Grid](#7-umap-parameter-grid) → [REFERENCE.md](./REFERENCE.md)
8. [t-SNE (when to use)](#8-t-sne-when-to-use)
9. [Gene Loading Analysis](#9-gene-loading-analysis)
10. [Saving Results](#10-saving-results)
11. [Expected Outputs / Validation](#11-expected-outputs--validation) → [REFERENCE.md](./REFERENCE.md)
12. [Common Pitfalls](#12-common-pitfalls)

---

## 1. HVG Selection (3′ untargeted only)

Skip this step for targeted panels — all genes are signal. Apply only when `adata.n_vars > 500`.

```python
import scanpy as sc

# flavor="seurat_v3" requires raw counts in adata.X
# Use flavor="seurat" if adata.X already holds log-normalized values
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    flavor="seurat_v3",   # or "seurat" for log-normalized
    batch_key="Well",     # compute per well, then combine — prevents well-specific genes dominating
    subset=False,         # keep all genes; HVG flag stored in adata.var["highly_variable"]
)

print(f"Highly variable genes: {adata.var['highly_variable'].sum()}")
sc.pl.highly_variable_genes(adata)

# Subset to HVGs before scaling — sc.pp.scale densifies whatever it is given
adata_hvg = adata[:, adata.var["highly_variable"]].copy()
print(f"HVG subset shape: {adata_hvg.shape}")
```

Keep `subset=False` so the HVG flag and the full gene set stay on `adata` for marker
detection and DE later. `adata_hvg` is used only for scale → PCA; `X_pca` is transferred
back to `adata` so neighbors / UMAP / Leiden find it there. Those steps read `obsm["X_pca"]`,
not `adata.X`. See [Run PCA](#3-run-pca).

**Parameter options for `highly_variable_genes`:**

| Parameter | Default | Notes |
|-----------|---------|-------|
| `n_top_genes` | 2000 | Increase to 3000–5000 for very heterogeneous datasets |
| `flavor` | `"seurat_v3"` | Use `"seurat"` for log-normalized; `"cell_ranger"` for Cell Ranger output |
| `batch_key` | `None` | Recommended when multiple wells/samples exist |
| `min_mean` / `max_mean` | auto | Fine-tune to exclude lowly / extremely expressed genes |
| `subset` | `False` | Set `True` to immediately drop non-HVGs from `adata` |

---

## 2. Scale Data

Centering and scaling gives each gene equal weight in PCA. Run **after** normalization (and HVG selection if applicable).

`sc.pp.scale` zero-centers, which destroys sparsity — the result is a dense
`float32` array of `n_cells × n_genes × 4` bytes regardless of how sparse the input
was. Scale the **HVG subset**, never the full gene set, on 3′ data:

```python
# 3′ untargeted: scale the HVG subset only — adata.X stays log-normalized
# This is intentional: adata.X must remain log-normalized for marker identification
# and percent-expressing calculations downstream (rank_genes_groups, score_genes).
sc.pp.scale(adata_hvg, max_value=10)

# Targeted panel (< ~500 genes, no HVG step): scale a copy, not adata itself
# adata.X must stay log-normalized for rank_genes_groups downstream
adata_scaled = adata.copy()
sc.pp.scale(adata_scaled, max_value=10)

# Verify
print(f"Mean of first gene: {adata_hvg.X[:, 0].mean():.4f}")   # ≈ 0
print(f"Std  of first gene: {adata_hvg.X[:, 0].std():.4f}")    # ≈ 1
```

**Check the dense footprint before scaling** when the object is large:

```python
gb = adata.n_obs * adata.n_vars * 4 / 1e9
print(f"Dense footprint if scaled as-is: {gb:.1f} GB")
```

A typical 3′ run of ~100k cells × ~13k post-QC genes is ~5 GB for the array alone,
and `sc.pp.scale` needs room for the copy on top of that — enough to kill the kernel
on a machine with limited memory. The same cells restricted to 2000 HVGs are under 1 GB. PCA only
uses the HVGs anyway (`use_highly_variable=True`), so nothing is lost.

> **Note:** Always scale a copy (`adata_hvg` for 3′, `adata_scaled` for targeted panels) — never scale `adata` in-place. `adata.X` must remain log-normalized for `rank_genes_groups` and `score_genes` downstream. Transfer `X_pca` back to `adata` after PCA (see [Run PCA](#3-run-pca)).

---

## 3. Run PCA

For 3′ data, run PCA on `adata_hvg` and transfer `X_pca` back to `adata`. Neighbors, UMAP, and Leiden all read from `obsm["X_pca"]` — not from `adata.X` — so they work correctly on `adata`. Marker identification then has access to all genes, not just the 2000 HVGs in `adata_hvg`:

```python
sc.pp.pca(adata_hvg, n_comps=50, svd_solver="arpack", random_state=42)
adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
adata.uns["pca"] = adata_hvg.uns["pca"]

print(f"PCA shape: {adata.obsm['X_pca'].shape}")  # (n_cells, 50)
```

For targeted panels / imaging MAD-z matrix, run PCA on the scaled copy and transfer back:

```python
sc.pp.pca(adata_scaled, n_comps=50, svd_solver="arpack", random_state=42)
adata.obsm["X_pca"] = adata_scaled.obsm["X_pca"]
adata.uns["pca"] = adata_scaled.uns["pca"]
```

**Key outputs:**

| Location | Content |
|----------|---------|
| `adata.obsm["X_pca"]` | Cell embeddings |
| `adata.varm["PCs"]` | Gene loadings |
| `adata.uns["pca"]["variance_ratio"]` | Fraction of variance per PC |

> **Memory (large runs):** once PCA is computed, downstream steps (neighbors, UMAP, Leiden) read `obsm["X_pca"]`, not `X`, so `adata_hvg` / `adata_scaled` are the largest dense objects you may no longer need — release them with `del adata_hvg` / `del adata_scaled` and `gc.collect()` if memory is tight. `adata.X` is log-normalized and stays in memory for marker calling.

**`svd_solver` options:**

| Value | When to use |
|-------|-------------|
| `"arpack"` | Default; good for sparse matrices |
| `"randomized"` | Faster for large dense matrices |
| `"auto"` | Let scanpy choose |

---

## 4. Variance / Elbow Analysis

Use the scree plot to decide how many PCs to pass to `sc.pp.neighbors`.

```python
import numpy as np
import matplotlib.pyplot as plt

variance_ratio = adata.uns["pca"]["variance_ratio"]
cumulative = np.cumsum(variance_ratio) * 100

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

# Scree plot
axes[0].bar(range(1, len(variance_ratio) + 1), variance_ratio * 100)
axes[0].set_xlabel("Principal Component")
axes[0].set_ylabel("Variance Explained (%)")
axes[0].set_title("Scree Plot")

# Cumulative variance
axes[1].plot(range(1, len(cumulative) + 1), cumulative, "b-o", markersize=3)
axes[1].axhline(80, color="red",    linestyle="--", label="80 %")
axes[1].axhline(90, color="orange", linestyle="--", label="90 %")
axes[1].set_xlabel("Principal Component")
axes[1].set_ylabel("Cumulative Variance (%)")
axes[1].set_title("Cumulative Variance")
axes[1].legend()

plt.tight_layout()
plt.show()

# Quick read-outs
n_pcs_80 = int(np.argmax(cumulative > 80)) + 1
n_pcs_90 = int(np.argmax(cumulative > 90)) + 1
print(f"PCs for 80 % variance: {n_pcs_80}")
print(f"PCs for 90 % variance: {n_pcs_90}")
```

**Automatic elbow detection (optional — requires `kneed`):**

```python
from kneed import KneeLocator

kl = KneeLocator(
    range(1, len(variance_ratio) + 1),
    variance_ratio,
    curve="convex",
    direction="decreasing",
)
print(f"Elbow at PC{kl.elbow}")
```

**Rule of thumb:** Use 20–50 PCs. For targeted panels, `n_pcs=30` is almost always sufficient; for 3′ untargeted data, check the elbow but rarely exceed 50.

---

## 5. Compute Neighborhood Graph

```python
sc.pp.neighbors(
    adata,
    n_neighbors=15,       # k in kNN; higher → smoother, more global structure
    n_pcs=30,             # PCs to use; adjust after elbow analysis
    metric="euclidean",   # distance metric
    method="umap",        # graph construction algorithm
)

print(f"n_neighbors: {adata.uns['neighbors']['params']['n_neighbors']}")
print(f"n_pcs:       {adata.uns['neighbors']['params']['n_pcs']}")
```

**Outputs:**

| Key | Content |
|-----|---------|
| `adata.obsp["distances"]` | Sparse kNN distance matrix |
| `adata.obsp["connectivities"]` | Sparse connectivity matrix |
| `adata.uns["neighbors"]` | Parameters |

**`n_neighbors` guidance:**

| Value | Effect |
|-------|--------|
| 5–10 | Fine-grained local structure; noisier |
| 15 | Recommended default |
| 30–50 | Smoother; better global topology |

---

## 6. Run UMAP

```python
sc.tl.umap(
    adata,
    min_dist=0.5,       # minimum distance between embedded points (higher = more spread)
    spread=1.0,         # scale of the embedding
    n_components=2,     # 2 for visualization; 3 for 3-D plots
    random_state=42,    # set for reproducibility
)

print(f"UMAP shape: {adata.obsm['X_umap'].shape}")  # (n_cells, 2)

sc.pl.umap(adata, color="Well")
```

---

## 7. UMAP Parameter Grid

See [REFERENCE.md](./REFERENCE.md) for the parameter cheat-sheet (`min_dist`, `n_neighbors`, `spread`).

Explore `min_dist` and `n_neighbors` to find the best layout for your data:

```python
# --- Effect of min_dist ---
fig, axes = plt.subplots(2, 2, figsize=(12, 12))
for ax, min_dist in zip(axes.flat, [0.1, 0.3, 0.5, 1.0]):
    tmp = adata.copy()
    sc.tl.umap(tmp, min_dist=min_dist, random_state=42)
    sc.pl.umap(tmp, color="Well", ax=ax, show=False, title=f"min_dist={min_dist}")
plt.tight_layout()
plt.show()

# --- Effect of n_neighbors ---
fig, axes = plt.subplots(2, 2, figsize=(12, 12))
for ax, k in zip(axes.flat, [5, 15, 30, 50]):
    tmp = adata.copy()
    sc.pp.neighbors(tmp, n_neighbors=k, n_pcs=30)
    sc.tl.umap(tmp, random_state=42)
    sc.pl.umap(tmp, color="Well", ax=ax, show=False, title=f"n_neighbors={k}")
plt.tight_layout()
plt.show()
```

---

## 8. t-SNE (when to use)

Prefer UMAP for most workflows. Use t-SNE when:
- Comparing with legacy analyses that used t-SNE
- You need maximum separation of dense, well-defined clusters

```python
sc.tl.tsne(
    adata,
    n_pcs=30,
    perplexity=30,   # analogous to n_neighbors; typical range 5–50
    random_state=42,
)

# Side-by-side comparison
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
sc.pl.umap(adata, color="Well", ax=axes[0], show=False, title="UMAP")
sc.pl.tsne(adata, color="Well", ax=axes[1], show=False, title="t-SNE")
plt.tight_layout()
plt.show()
```

**t-SNE caveats:** Slower than UMAP; does not preserve global distances; cluster sizes and inter-cluster distances are not meaningful.

---

## 9. Gene Loading Analysis

Identify which genes drive each principal component.

```python
import numpy as np

def get_top_genes_per_pc(adata, n_genes=10, n_pcs=5):
    """Return top genes by absolute loading for the first n_pcs PCs."""
    loadings = adata.varm["PCs"]  # (n_genes_in_pca, n_pcs)

    # When HVG subsetting was used, varm["PCs"] rows correspond to HVGs only
    if "highly_variable" in adata.var.columns and adata.var["highly_variable"].sum() < adata.n_vars:
        gene_names = adata.var_names[adata.var["highly_variable"]]
    else:
        gene_names = adata.var_names

    results = {}
    for pc in range(n_pcs):
        vals = loadings[:, pc]
        top_pos = np.argsort(vals)[-n_genes:][::-1]
        top_neg = np.argsort(vals)[:n_genes]
        results[f"PC{pc+1}_positive"] = list(zip(gene_names[top_pos], vals[top_pos]))
        results[f"PC{pc+1}_negative"] = list(zip(gene_names[top_neg], vals[top_neg]))
    return results


top_genes = get_top_genes_per_pc(adata)
print("Top genes driving PC1:")
print("  Positive:", [g for g, _ in top_genes["PC1_positive"][:5]])
print("  Negative:", [g for g, _ in top_genes["PC1_negative"][:5]])

# Built-in scanpy plot
sc.pl.pca_loadings(adata, components=[1, 2, 3])
```

---

## 10. Saving Results

All embeddings are stored inside `adata` automatically. Persist with:

```python
# Save the full object
adata.write_h5ad("processed_with_dr.h5ad")

# Quick inventory
print("Dimensionality reduction complete")
print(f"  PCA:        {adata.obsm['X_pca'].shape}")
print(f"  UMAP:       {adata.obsm['X_umap'].shape}")
print(f"  Neighbors:  {adata.obsp['connectivities'].shape}")
```

**What is stored and where:**

| Location | Content |
|----------|---------|
| `adata.obsm["X_pca"]` | PCA cell embeddings |
| `adata.obsm["X_umap"]` | UMAP cell coordinates |
| `adata.varm["PCs"]` | Gene loadings |
| `adata.uns["pca"]` | PCA metadata (variance ratio, etc.) |
| `adata.uns["neighbors"]` | Neighbor-graph parameters |
| `adata.obsp["distances"]` | Sparse kNN distance matrix |
| `adata.obsp["connectivities"]` | Sparse connectivity matrix |

---

## 11. Expected Outputs / Validation

See [REFERENCE.md](./REFERENCE.md) — key parameter summary, expected output shapes, and validation assert checks.

---

## 12. Common Pitfalls

| Pitfall | Fix |
|---------|-----|
| **Scaling all genes on 3′ data** | `sc.pp.scale` densifies: `n_cells × n_genes × 4` bytes. ~100k cells × ~13k genes is ~5 GB plus a copy, which can exhaust memory on a small machine. Subset to HVGs first (`adata_hvg = adata[:, adata.var["highly_variable"]].copy()`) and scale that. |
| **Scaling raw counts** | Always normalize first (`ebsd.pp.normalize`) before `sc.pp.scale`. Raw counts are dominated by highly expressed genes and will distort PCA. |
| **HVG on normalized data with `seurat_v3`** | Use `flavor="seurat_v3"` for raw counts; switch to `flavor="seurat"` if `adata.X` is already log-normalized. |
| **Too few PCs** | Using < 15 PCs risks losing biological signal. Start with 30–50 and check the elbow plot. |
| **Lazy AnnData passed to `scale` / PCA** | Reload with `lazy=False` (`Loader.load_tables` defaults to `lazy=True` and keeps `X` on disk) or call `.to_memory()`. Materialize at **load**, not right before `sc.pp.scale` — by then you pull in every gene and then densify it. |
| **`TypeError: scale() got an unexpected keyword argument 'return_mean_std'`** | Not a version problem. `sc.pp.scale` dispatches on `type(adata.X)` and only implements `np.ndarray`, scipy `csr`/`csc`, and dask `Array`; anything else hits the catch-all, which lacks that internal argument. Print `type(adata.X)` — if it is a backed / lazy type (`_CSRDataset`, zarr-backed), materialize it. **Do not hand-roll the standardization**: it densifies all genes, and scanpy's std uses `correction=1` while a naive `X.std(0)` uses `ddof=0`, so the values no longer match `sc.pp.scale` output. |
| **Missing `random_state`** | UMAP and t-SNE are stochastic. Always set `random_state=42` (or any fixed seed) for reproducibility. |
| **Over-interpreting UMAP** | UMAP distorts global distances and cluster sizes. Do not draw quantitative conclusions from proximity or cluster area alone. |
| **HVG selection on targeted panels** | Targeted panels are pre-selected for signal; applying HVG selection can inadvertently discard informative genes. Skip it. |
| **Imaging via count normalize / HVG / second scale** | Imaging uses MAD-z from [imaging-pipeline](../imaging-pipeline/SKILL.md); no `pp.normalize`, no HVG, no extra `sc.pp.scale`. Exclude Bkg/IgG/Cleave from the embedding matrix. |
