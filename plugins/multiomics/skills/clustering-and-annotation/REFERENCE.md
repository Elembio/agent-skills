# Clustering and Annotation — Reference

## Contents

1. [sc.tl.leiden parameters](#1-sctlleiden-parameters)
2. [sc.tl.rank_genes_groups parameters](#2-sctlrank_genes_groups-parameters)
3. [Expected outputs and validation](#3-expected-outputs-and-validation)

---

## 1. sc.tl.leiden parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| `resolution` | `1.0` | 0.1 (few) to 2.0 (many) clusters |
| `random_state` | `0` | Set for deterministic results |
| `key_added` | `"leiden"` | Column written to `adata.obs` |
| `flavor` | `"igraph"` | Algorithm backend |
| `n_iterations` | `2` | `-1` runs until convergence |

---

## 2. sc.tl.rank_genes_groups parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| `method` | — | `"wilcoxon"` recommended; also `"t-test"`, `"logreg"` |
| `pts` | `False` | Compute percent-expressing columns |
| `min_fold_change` | `0` | Pre-filter by fold change |
| `min_in_group_fraction` | `0` | Min fraction expressing in cluster |
| `key_added` | `"rank_genes_groups"` | Key in `adata.uns` |

---

## 3. Expected outputs and validation

After running the full workflow, `adata` should contain:

| Location | Key | Description |
|----------|-----|-------------|
| `adata.obs` | `"leiden"` | Cluster integer labels |
| `adata.obs` | `"cell_type"` | Annotated biological labels |
| `adata.uns` | `"rank_genes_leiden"` | Marker gene results |
| `adata.obs` | `"S_score"`, `"G2M_score"` | Cell cycle scores (if computed) |
| `adata.obs` | `"<name>_score"` | Custom gene set scores |

Validation checks:

```python
assert "leiden" in adata.obs.columns, "Clustering not performed"
assert adata.obs["leiden"].nunique() > 1, "Only one cluster found — lower resolution or check neighbor graph"
assert "rank_genes_leiden" in adata.uns, "Marker genes not computed"

# Sanity: no cluster smaller than ~10 cells
min_size = adata.obs["leiden"].value_counts().min()
if min_size < 10:
    print(f"Warning: smallest cluster has {min_size} cells — consider lower resolution")
```
