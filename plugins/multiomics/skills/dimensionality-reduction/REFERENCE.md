# Dimensionality Reduction — Reference

## Contents

1. [Key parameter summary](#1-key-parameter-summary)
2. [UMAP parameter cheat-sheet](#2-umap-parameter-cheat-sheet)
3. [Expected outputs and validation](#3-expected-outputs-and-validation)

---

## 1. Key parameter summary

| Parameter | Function | Typical value | Notes |
|---|---|---|---|
| `n_comps` | `sc.pp.pca` | `min(50, n_vars - 1)` | How many PCs to **compute** — compute more than you need |
| `n_pcs` | `sc.pp.neighbors` | 30 (targeted), 30–50 (3′) | How many PCs to **use for the graph** — drives UMAP structure |
| `n_neighbors` | `sc.pp.neighbors` | 15 | Graph connectivity; lower = finer local detail |
| `svd_solver` | `sc.pp.pca` | `"arpack"` | Deterministic; preferred for reproducibility |
| `min_dist` | `sc.tl.umap` | 0.5 | Minimum distance between embedded points |
| `spread` | `sc.tl.umap` | 1.0 | Scale of the embedding |
| `random_state` | `sc.tl.umap`, `sc.tl.tsne` | 42 | Always set for reproducibility |
| `max_value` | `sc.pp.scale` | 10 | Clips scaled values to prevent outlier genes dominating PCA |

---

## 2. UMAP parameter cheat-sheet

| Parameter | Low value effect | High value effect |
|---|---|---|
| `min_dist` | Tight, compact clusters | Loose, uniform distribution |
| `n_neighbors` | Fine local detail | Coarser, more global topology |
| `spread` | Compact embedding | Expanded embedding |

---

## 3. Expected outputs and validation

After a successful run, `adata` should contain:

| Check | Expected |
|---|---|
| `adata.obsm["X_pca"].shape[1]` | 30–50 |
| Cumulative variance (30 PCs) | > 60% for targeted; varies for 3′ |
| `adata.obsm["X_umap"].shape[1]` | 2 |
| NaN count in UMAP | 0 |

Validation checks:

```python
import numpy as np

assert "X_pca" in adata.obsm, "PCA not computed"
assert "X_umap" in adata.obsm, "UMAP not computed"
assert not np.isnan(adata.obsm["X_umap"]).any(), "NaN values in UMAP coordinates"
assert "connectivities" in adata.obsp, "Neighbor graph not computed"

total_var = adata.uns["pca"]["variance_ratio"].sum()
assert total_var > 0.3, f"Total variance explained too low: {total_var:.1%}"
```
