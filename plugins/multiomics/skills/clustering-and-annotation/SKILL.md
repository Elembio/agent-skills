---
name: clustering-and-annotation
description: Unsupervised cell clustering using Leiden algorithm, marker gene identification with Wilcoxon tests, and cell type annotation strategies. Resolution selection, cluster quality assessment, and cell state scoring. Use when the user asks about clustering, finding cell types, annotating clusters, identifying marker genes, or visualizing clusters on UMAP. Works on any normalized AnnData with a computed neighbor graph — AVITI24, scRNA-seq, or AnnData built from a CSV/parquet count matrix.
metadata:
  version: 0.2.4
  author: elembio
---

# Clustering and Annotation

## Purpose

This skill covers unsupervised cell clustering using the Leiden algorithm, marker gene identification, and cell type annotation strategies. It explains resolution selection, cluster quality assessment, and cell state scoring for single-cell expression data — AVITI24 spatial transcriptomics, scRNA-seq, or any AnnData built from a generic CSV/parquet count matrix.

## Prerequisites

- Completed: [Cell Quality Control](../cell-quality-control/SKILL.md) and [Normalization Strategies](../normalization-strategies/SKILL.md)
- Completed: [Dimensionality Reduction](../dimensionality-reduction/SKILL.md)
- An AnnData object with a computed neighbor graph (`adata.uns["neighbors"]`)

## Key Concepts

### Leiden vs Louvain

| Algorithm | Advantages | Use When |
|-----------|------------|----------|
| **Leiden** | Guaranteed connected clusters, faster | Default choice |
| **Louvain** | Legacy, widely published | Reproducibility with older studies |

### Resolution Parameter

- **Lower (0.1–0.5)**: Fewer, broader clusters — suitable for major cell types
- **Medium (0.5–1.0)**: Balanced clustering — good starting point
- **Higher (1.0–2.0)**: More, finer clusters — suitable for cell states

### Large datasets (100K+ cells)

On 100K+ cells, **use `flavor="igraph", n_iterations=2`** (shown in Quick Start). The legacy `leidenalg` backend is markedly slower and heavier on large graphs; the `igraph` backend is the intended default. It can assign different cluster labels than `leidenalg`, so results shift versus prior legacy-backend runs (still deterministic under `random_state`).

### Workflow Overview

1. **Neighborhood Graph** — Built in PCA space (already computed upstream)
2. **Community Detection** — Leiden algorithm partitions cells into clusters
3. **Marker Genes** — Wilcoxon rank-sum test identifies per-cluster enriched genes
4. **Annotation** — Biological labels assigned to clusters from marker evidence

## Quick Start

```python
import scanpy as sc

# Ensure neighbor graph exists
if "neighbors" not in adata.uns:
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)

# Cluster at balanced resolution.
# flavor="igraph" + n_iterations=2 is the fast, memory-lean backend — important at scale.
sc.tl.leiden(adata, resolution=0.5, key_added="leiden",
             flavor="igraph", n_iterations=2, directed=False, random_state=42)
print(f"Clusters: {adata.obs['leiden'].nunique()}")
sc.pl.umap(adata, color="leiden", legend_loc="on data")

# Identify marker genes
sc.tl.rank_genes_groups(
    adata, groupby="leiden", method="wilcoxon",
    use_raw=False,                   # adata.raw is not set; use adata.X (log-normalized)
    pts=True, key_added="rank_genes_leiden",
)
sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False, key="rank_genes_leiden")
```

## Detail Files

- [IMPLEMENTATION.md](./IMPLEMENTATION.md) — Multi-resolution sweep, silhouette selection, full marker visualization (dotplot/heatmap/violin), manual annotation mapping, cell cycle scoring, gene set scoring, cluster composition analysis, export, and common pitfalls
- [REFERENCE.md](./REFERENCE.md) — `sc.tl.leiden` and `sc.tl.rank_genes_groups` parameter tables, expected outputs, and validation checks

## Related Skills

- [Dimensionality Reduction](../dimensionality-reduction/SKILL.md) — Compute neighbor graph
- [Cluster Validation](../cluster-validation/SKILL.md) — Validate cluster robustness
- [Pseudobulk Differential Expression](../pseudobulk-differential-expression/SKILL.md) — DE analysis by cluster

## Disclaimer

This content is provided "as-is" without warranty of any kind. AI-generated output may contain
inaccuracies. Review and verify all outputs before relying on them.
