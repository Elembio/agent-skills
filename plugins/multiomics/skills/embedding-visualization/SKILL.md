---
name: embedding-visualization
description: Visualizing UMAP/PCA embeddings — color by genes, clusters, Well (plate position), WellLabel or condition (sample name), QC, or imaging metrics (S/B, nuclear fraction). Always use when plotting UMAP or PCA. Canonical skill for embedding plots; compute embeddings via dimensionality-reduction first.
metadata:
  version: 0.4.3
  author: elembio
---

# Embedding Visualization

## Purpose

**Canonical skill for coloring/plotting** PCA/UMAP embeddings (counts or imaging). Does not build the embedding — that is [dimensionality-reduction](../dimensionality-reduction/SKILL.md).

## Prerequisites

- Embeddings already in `obsm` from [dimensionality-reduction](../dimensionality-reduction/SKILL.md)
- Upstream prep: counts → [normalization-strategies](../normalization-strategies/SKILL.md); imaging feature matrix → [imaging-pipeline](../imaging-pipeline/SKILL.md) then DR

## Key Concepts

### Metadata Vocabulary

| Generic Term | Example Mappings |
|---|---|
| `Well` | Plate position (`A1`) — technique / replicate unit |
| `WellLabel` | Sample-sheet string — use this (or `condition`) for biology |
| `tile_id` | field of view, tile, region |
| `leiden` | cluster label from leiden/louvain |
| `cell_type` | annotated cell type string |
| `scores` / imaging metrics | `sb_*`, `nucfrac_*` (from imaging-pipeline prep), QC scores |
| `condition` | from [experiment-design](../experiment-design/SKILL.md) |

### Colormap conventions — red is hotter/higher

| Data type | Colormap | Notes |
|---|---|---|
| Gene expression, counts, nonneg metrics | `viridis` | Sequential; high = yellow |
| Diverging — z-scores, fold changes, batch scores (centered at 0) | `RdBu_r` + `vcenter=0` | `_r` reverses so red = high; without `_r`, blue is high |
| Nuclear fraction, similarity (0–1) | `coolwarm`, `vmin=0, vmax=1` | Red = fully nuclear / maximally similar |
| Pass / fail rates | `RdYlGn` | Semantic exception — green = passing |
| Categorical | scanpy default palette | Do not set `cmap` |

> Never use `RdBu` without `_r`. Plain `RdBu` puts **blue at the high end**.

### Imaging / IF coloring tips

- Color **`Well` and `WellLabel` (or `condition`) together, always** — `Well` is plate position; `WellLabel` is the sample name. Either alone hides whether a separation is biology or well technique. Then `sb_*` / `nucfrac_*`. There is no `obs["well_label"]`.
- **Label the figure with `prof.uns["scaling_reference"]`.** An imaging embedding is uninterpretable without knowing whether it was scaled per-well, against controls, or globally; a per-well embedding on a plate where condition is constant within wells shows no separation by construction. See [imaging-pipeline](../imaging-pipeline/SKILL.md#scaling).
- Reading well-vs-condition structure → [batch-effect-assessment](../batch-effect-assessment/SKILL.md).
- Prefer `obs` metrics over assuming stain names are in `var_names`.
- UMAP is overview only — for “what differs by treatment?” use well-level contrasts in [imaging-pipeline](../imaging-pipeline/SKILL.md), not the plot alone.

### Visualization Type Selection

| Type | Use Case | Example Variables |
|---|---|---|
| **Categorical** | Discrete groups | `leiden`, `cell_type`, `Well`, `WellLabel` |
| **Continuous** | Quantitative values | gene expression, `total_counts` |
| **Multi-panel** | Compare multiple variables | condition + well + gene |
| **Split / facet** | Separate plot per group | facet by condition or leiden |

## Quick Start

```python
import scanpy as sc

# Categorical: cluster labels or metadata
sc.pl.umap(adata, color="leiden", legend_loc="on data", frameon=False)
sc.pl.umap(adata, color="Well", frameon=False)          # plate position
label = "condition" if "condition" in adata.obs.columns else "WellLabel"
sc.pl.umap(adata, color=label, frameon=False)           # sample-sheet name

# Gene expression — uses adata.X (normalized) by default
gene = "GAPDH"
if gene not in adata.var_names:
    hits = adata.var_names[adata.var_names.str.contains(gene, case=False)]
    print(hits[:10])
else:
    sc.pl.umap(adata, color=gene, cmap="viridis", vmin=0, frameon=False)

# Clip contrast to 99th percentile (useful when bright outlier cells dominate)
sc.pl.umap(adata, color=gene, cmap="viridis", vmax="p99", frameon=False)

# Plot from a specific layer (e.g. raw counts)
sc.pl.umap(adata, color=gene, layer="raw", cmap="viridis", frameon=False)

# Multi-variable in one call — mix genes and obs columns
sc.pl.umap(adata, color=["leiden", "Well", label, gene], ncols=4, frameon=False)
```

## Detail Files

- [IMPLEMENTATION.md](./IMPLEMENTATION.md) — metadata mapping patterns, multi-panel UMAPs, gene expression overlays, faceting by condition, density/contour plots, publication-quality rcParams, interactive plotly, condition comparison, and common pitfalls
- [REFERENCE.md](./REFERENCE.md) — `sc.pl.umap` key options, recommended colormaps, expected figure outputs

## Related Skills

- [dimensionality-reduction](../dimensionality-reduction/SKILL.md) — compute PCA/UMAP
- [clustering-and-annotation](../clustering-and-annotation/SKILL.md) — visualize clusters
- [batch-effect-assessment](../batch-effect-assessment/SKILL.md) — assess technical effects

## Disclaimer

This content is provided "as-is" without warranty of any kind. AI-generated output may contain
inaccuracies. Review and verify all outputs before relying on them.
