---
name: dimensionality-reduction
description: Canonical PCA, UMAP, t-SNE, HVG (3′ only), and neighbors. Count modalities after ebsd.pp.normalize; imaging/IF after the scaled profiling matrix from imaging-pipeline. Use when computing embeddings — not for coloring plots (embedding-visualization) or building the imaging feature matrix (imaging-pipeline).
metadata:
  version: 0.4.2
  author: elembio
---

# Dimensionality Reduction

## Purpose

**Canonical skill for computing** PCA / neighbors / UMAP (and HVG for 3′). Does not build imaging feature matrices or color plots.

## Prerequisites

- [cell-quality-control](../cell-quality-control/SKILL.md)
- **Counts:** [normalization-strategies](../normalization-strategies/SKILL.md)
- **Imaging / IF:** profiling AnnData already MAD-z–scaled via [imaging-pipeline](../imaging-pipeline/SKILL.md) (exclude Bkg/IgG/Cleave upstream). That skill owns **which reference** the scaling used — per-well centering is invalid when condition is constant within a well
- `lazy=False` before PCA / neighbors

## Key Concepts

### Which pipeline?

| Data | Before PCA | HVG? | Extra `sc.pp.scale`? |
|---|---|---|---|
| 3′ | `ebsd.pp.normalize` | Yes (~2000) | Yes — `sc.pp.scale` on the **HVG subset**, then PCA |
| Transcript / Protein (targeted) | `ebsd.pp.normalize` | No | Yes (`sc.pp.scale`) then PCA |
| **Imaging / IF / morphology** | **MAD-z** matrix from imaging-pipeline, scaled against the reference the design implies | **No** | **No** — already scaled; then PCA → neighbors → UMAP |

Imaging matrix must already exclude Bkg / IgG / Cleave upstream ([imaging-pipeline](../imaging-pipeline/SKILL.md)). Prefer a separate profiling `AnnData` so raw imaging stays for S/B plots.

**Before running PCA on an imaging matrix, check `prof.uns["scaling_reference"]`.** A per-well reference on a plate where condition is constant within wells has already removed the effect — the embedding will show no separation no matter what you do here. The reference choice belongs to [imaging-pipeline](../imaging-pipeline/SKILL.md#scaling); carry the label onto the figures. On lean imaging feature sets, show PC1/PC2 alongside the UMAP — a clean linear split in early PCs is often softened by UMAP.

### HVG selection: only for 3′/untargeted data

For **targeted panels** (< ~500 genes) skip HVG. For **3′**, select ~2000 HVGs and subset to them before scaling and PCA. Never HVG on imaging features.

### Pipeline order

```
3′ counts:      HVG → subset to HVGs → sc.pp.scale → PCA → transfer X_pca → neighbors → UMAP
Targeted panel: sc.pp.scale → PCA → neighbors → UMAP
Imaging:        MAD-z feature matrix (no Bkg/IgG/Cleave; reference set by design) → PCA → neighbors → UMAP
```

**Scale the HVG subset, not the full object.** `sc.pp.scale` zero-centers and therefore
densifies: `n_cells × n_genes × 4` bytes. A ~100k-cell 3′ run with ~13k post-QC genes is
~5 GB before the copy scanpy makes, which will kill a small kernel; the 2000-HVG subset
is under 1 GB and is all PCA uses anyway.

### Where results live

| Key | Content |
|-----|---------|
| `adata.obsm["X_pca"]` | Cell PCA coordinates *(n_cells × n_pcs)* |
| `adata.obsm["X_umap"]` | Cell UMAP coordinates *(n_cells × 2)* |
| `adata.uns["neighbors"]` | Neighbor-graph parameters |
| `adata.uns["pca"]["variance_ratio"]` | Variance explained per PC |
| `adata.varm["PCs"]` | Gene loadings *(n_genes × n_pcs)* |

### n_comps vs n_pcs

These are different parameters with different effects:

| Parameter | Where | What it controls |
|---|---|---|
| `n_comps` | `sc.pp.pca(n_comps=50)` | How many PCs are **computed** — compute more than you need so you can inspect the elbow |
| `n_pcs` | `sc.pp.neighbors(n_pcs=30)` | How many PCs are **used for the graph** — this is what drives UMAP structure |

**Targeted Transcript panels** — few features, so most variance is captured early. Start with `n_pcs=30` for neighbors; rarely need more.

**3′ untargeted** — more genes, more variance in higher PCs. Use 30–50 PCs for neighbors; check the elbow plot and look for where cumulative variance flattens (~70–90%).

To find the elbow:
```python
import matplotlib.pyplot as plt
plt.plot(adata.uns["pca"]["variance_ratio"].cumsum())
plt.xlabel("PC"); plt.ylabel("Cumulative variance explained")
```

**Rule**: cap `n_pcs` at `adata.obsm["X_pca"].shape[1]`; if late PCs are driven by batch noise, use fewer.

### Recommended defaults

- `n_comps=min(50, adata.n_vars - 1)` — compute up to 50 PCs
- `n_pcs=30` — use 30 for the neighbor graph (targeted panels)
- `n_neighbors=15` — good balance of local and global structure

## Quick Start — targeted panel

```python
import numpy as np
import scanpy as sc
from scipy import sparse

# adata is normalized. Materialize at LOAD — Loader.load_tables defaults to lazy=True
# and leaves X on disk:  loader.load_tables(table=..., lazy=False)
if not isinstance(adata.X, (np.ndarray, sparse.csr_matrix, sparse.csc_matrix)):
    print("X is", type(adata.X).__name__)  # backed / lazy → adata = adata.to_memory()

sc.pp.scale(adata, max_value=10)
sc.pp.pca(adata, n_comps=min(50, adata.n_vars - 1), svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata, min_dist=0.5, random_state=42)
# Color / plot → embedding-visualization
```

**If `sc.pp.scale` raises `TypeError: ... unexpected keyword argument 'return_mean_std'`,** `X` is not a type `scale` implements — usually because it was never materialized. `scale` dispatches on `type(adata.X)` and only implements ndarray / scipy csr-csc / dask; everything else falls to a catch-all that lacks that internal argument, so the message names the wrong problem. Fix the type — never reimplement the scaling by hand ([IMPLEMENTATION.md](./IMPLEMENTATION.md#12-common-pitfalls)).

For 3′ HVG, elbow, t-SNE, loadings — see **[IMPLEMENTATION.md](./IMPLEMENTATION.md)**. For imaging MAD-z matrix construction — [imaging-pipeline](../imaging-pipeline/SKILL.md). Plotting — [embedding-visualization](../embedding-visualization/SKILL.md).

## Detail Files

- [IMPLEMENTATION.md](./IMPLEMENTATION.md) — HVG, elbow, parameter grids, pitfalls
- [REFERENCE.md](./REFERENCE.md) — parameter tables and validation

## Related Skills

- [normalization-strategies](../normalization-strategies/SKILL.md) — count normalize before DR
- [imaging-pipeline](../imaging-pipeline/SKILL.md) — build imaging MAD-z matrix and choose its scaling reference (then come here for PCA/UMAP)
- [embedding-visualization](../embedding-visualization/SKILL.md) — color / plot embeddings (after this skill)
- [cell-quality-control](../cell-quality-control/SKILL.md)
- [spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md)

## Disclaimer

This content is provided "as-is" without warranty of any kind. AI-generated output may contain
inaccuracies. Review and verify all outputs before relying on them.
