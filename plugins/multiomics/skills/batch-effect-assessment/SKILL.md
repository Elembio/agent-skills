---
name: batch-effect-assessment
description: Detecting and assessing batch effects in single-cell or expression-matrix data — across wells, runs, molecular batches, or any user-defined batch variable. Visualizing technical variation, quantifying batch contribution with ANOVA/silhouette scores, and strategies for batch correction (ComBat, Harmony). Use when the user asks about batch effects, technical variation, samples/wells clustering together, cells separating by run or batch instead of condition, or whether to apply batch correction. Works for AVITI24 AnnData as well as generic AnnData or tabular data (CSV/parquet) that carries a batch label.
metadata:
  version: 0.3.4
  author: elembio
---

# Batch Effect Assessment

## Purpose

This skill covers detecting and assessing batch effects in single-cell or expression-matrix data. It applies to AVITI24 AnnData as well as generic AnnData or tabular data (CSV/parquet) that carries one or more batch labels. It explains how to identify technical variation across wells, runs, and batches, and strategies for batch correction when needed. Apply it after QC and normalization to decide whether correction is necessary before downstream analysis.

## Prerequisites

- Completed: [experiment-design](../experiment-design/SKILL.md) — `condition` column confirmed with user (needed to separate biological from technical variation)
- Completed: [cell-quality-control](../cell-quality-control/SKILL.md)
- Completed: [normalization-strategies](../normalization-strategies/SKILL.md)
- A normalized AnnData object with `.obs` columns identifying one or more batch variables (for AVITI24 data: `Well`, `tile_id` (or `Tile`), and optionally `run_id` / `run_name`; for generic data: any user-defined batch label such as `batch`, `sample_id`, `plate`, or `donor`)

## Key Concepts

### Sources of Technical Variation

| Batch class | What it can introduce |
|-------------|----------------------|
| **Well** | Edge effects, reagent distribution gradients (moderate impact) |
| **Run** | Library prep differences, instrument-to-instrument variation (high impact) |
| **Batch (molecular)** | Probe efficiency differences across lot numbers (variable impact) |
| **Tile** | Focus and illumination non-uniformity within a well (low–moderate impact) |

### Batch vs. Biology

Batch effects are **problematic** when:
- Technical groupings cluster more strongly than biological groupings in PCA/UMAP
- The same biological condition splits by well or run rather than clustering together
- Control samples from different wells fail to co-cluster

Good experimental design that keeps the two separable includes biological replicates across wells/runs, technical replicates within wells, and balanced conditions across batches.

### When to Apply Correction

**Correct** when technical variables explain more variance than biological variables, or when cells cluster primarily by batch rather than condition.

**Do not correct** when the batch effect is minor (< 5% variance explained), when correction removes biological signal, or when the design is confounded (condition constant within each well — correction cannot fix this).

### Confounded designs: reading well structure you cannot correct

When `condition` is constant within every well — the usual AVITI24 plate layout, flagged by [experiment-design](../experiment-design/SKILL.md) — well technique and biology are aliased. Any per-well correction (per-well centering, ComBat on `Well`, Harmony on `Well`) removes the condition along with the batch. **Replicate wells do not change this**; the criterion is heterogeneity *within* a well, not the number of wells per condition.

What you can still do is read the residual structure, given ≥2 wells per condition and an embedding built **without** per-well centering (for imaging: the global or control-referenced matrix from [imaging-pipeline](../imaging-pipeline/SKILL.md)):

| Pattern on the embedding | Reading |
|---|---|
| Conditions separate, same-condition wells overlap | Biology dominates; well effect is small relative to it |
| Conditions separate, wells form islands **within** each condition | Real biology **plus** residual well technique — global scaling does not remove it, and per-well would delete the biology |
| Wells separate but conditions do not | Well technique dominates, or there is no condition effect — not distinguishable here |
| One well per condition | **Cannot attribute.** Report the confounding; do not present the split as a phenotype |

Corroborate rather than correct: check that the driving features are biologically coherent, and look for the same split in an orthogonal modality on the same cells (3′ RNA / OPS identity). Cross-modality agreement is far stronger evidence than any embedding geometry.

Always color the embedding by **both** `Well` (plate position) and `condition` or `WellLabel` (sample name) — one alone cannot show this. There is no `obs["well_label"]`.

## Quick Start

Color UMAP by Well and compute per-well cluster composition (assumes UMAP already computed):

```python
import scanpy as sc

# If UMAP not yet computed:
# sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
# sc.tl.umap(adata)

cond = next(
    (c for c in ("condition", "experimental_condition", "WellLabel") if c in adata.obs.columns),
    "WellLabel",
)
sc.pl.umap(adata, color=["Well", cond],
           ncols=2, title=["UMAP by Well", f"UMAP by {cond}"])

# Per-well cluster composition (requires leiden clustering)
sc.tl.leiden(adata, resolution=0.5, flavor="igraph", n_iterations=2, directed=False)
composition = (adata.obs
               .groupby(["Well", "leiden"])
               .size()
               .unstack(fill_value=0)
               .apply(lambda r: r / r.sum(), axis=1))
print(composition)
```

## Detail Files

- [IMPLEMENTATION.md](./IMPLEMENTATION.md) — Full implementations: PCA/UMAP batch coloring, variance explained, silhouette scores, replicate pseudobulk correlation heatmap, Harmony integration, ComBat expression correction, batch report generator.
- [REFERENCE.md](./REFERENCE.md) — Parameters/thresholds table, expected outputs, common pitfalls.

## Related Skills

- [Normalization Strategies](../normalization-strategies/SKILL.md) — normalize before assessment
- [Dimensionality Reduction](../dimensionality-reduction/SKILL.md) — PCA/UMAP for visualization
- [Experiment Design](../experiment-design/SKILL.md) — reports whether condition is confounded with Well
- [Imaging Pipeline](../imaging-pipeline/SKILL.md) — imaging scaling reference (per-well centering is itself a well correction)
- [Pseudobulk Differential Expression](../pseudobulk-differential-expression/SKILL.md) — model batch effects in DE

## Disclaimer

This content is provided "as-is" without warranty of any kind. AI-generated output may contain
inaccuracies. Review and verify all outputs before relying on them.
