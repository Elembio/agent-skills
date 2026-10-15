---
name: pseudobulk-differential-expression
description: Pseudobulk differential expression analysis using PyDESeq2. Aggregating single-cell data to pseudobulk by sample / replicate / well, running DE with biological replicates, multi-level analysis (macro/meso), and per-cell-type DE. Works on any AnnData with raw integer counts in a layer — AVITI24, scRNA-seq, or AnnData built from a CSV count matrix plus sample metadata. Use when the user asks about differential expression, which genes are up or down, comparing gene expression between conditions or treatments, or DE analysis.
metadata:
  version: 0.2.2
  author: elembio
---

# Pseudobulk Differential Expression Analysis

## Purpose

This skill covers pseudobulk differential expression (DE) analysis using PyDESeq2. It handles both macro-level (condition vs condition across all cells) and meso-level (per-cluster, cell-type-specific) analysis. Raw counts from `layers["raw"]` are always used — never normalized values.

## Prerequisites

- `condition`, `replicate`, `group` columns present in `adata.obs` with `groupby` and `reference_group` established. For AVITI24 data use [experiment-design](../experiment-design/SKILL.md) to derive these from well labels; for generic data the user supplies them directly as a metadata CSV or DataFrame.
- QC and filtering applied ([cell-quality-control](../cell-quality-control/SKILL.md) for AVITI24, or equivalent upstream filtering).
- Raw integer counts preserved in a layer (AVITI24: `.layers["raw"]` after [normalization-strategies](../normalization-strategies/SKILL.md); generic: whatever layer you populated from the source counts — pass it via the `layer` argument).
- Cluster labels in `adata.obs` for meso-level analysis ([clustering-and-annotation](../clustering-and-annotation/SKILL.md) or any external clustering).

## Key Concepts

### Macro vs Meso Analysis

| Level | Aggregation | What it answers | When to use |
|-------|-------------|-----------------|-------------|
| **Macro** | All cells per well/sample | Overall condition effects across the full dataset | First-pass DE, treatment vs control |
| **Meso** | Cells per cluster per well | Cell-type-specific condition effects | After clustering/annotation; reveals heterogeneous responses |

### The Raw Counts Rule

DESeq2 requires raw (un-normalized) integer counts. Always use `layers["raw"]`, never `adata.X` after normalization. The `create_pseudobulk` helper accepts a `layer` argument that defaults to `"raw"`.

### Minimum Replicates Guideline

At least **3 biological replicates per condition** are recommended (hard minimum: 2). Pseudobulk samples with fewer than 10 cells are dropped before aggregation.

---

## Quick Start

`create_pseudobulk` is a helper defined in [IMPLEMENTATION.md §4](./IMPLEMENTATION.md#4-create_pseudobulk--full-implementation) — copy it into your notebook before use.

```python
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# create_pseudobulk — see IMPLEMENTATION.md §4 for full definition
# Always uses layers["raw"] (raw integer counts) — never normalized .X

# 1. Aggregate to pseudobulk (macro level)
pseudobulk = create_pseudobulk(
    adata,
    sample_col="Well",
    condition_col="experimental_condition",
    min_cells=10,
    layer="raw",
)

# 2. Build DESeq2 dataset and run
counts_df = pd.DataFrame(
    pseudobulk.X.astype(int),
    index=pseudobulk.obs_names,
    columns=pseudobulk.var_names,
)
dds = DeseqDataSet(counts=counts_df, metadata=pseudobulk.obs.copy(),
                   design_factors="experimental_condition")
dds.deseq2()

# 3. Extract results
stat_res = DeseqStats(dds, contrast=["experimental_condition", "Treatment", "Control"])
stat_res.summary()
results_df = stat_res.results_df.sort_values("padj")
print(f"Significant DEGs (padj < 0.05): {(results_df['padj'] < 0.05).sum()}")
```

## Detail Files

- [IMPLEMENTATION.md](./IMPLEMENTATION.md) — raw layer handling, design setup, `create_pseudobulk` full definition, multiple contrasts, meso-level loop, covariates, saving results, validation plots.
- [REFERENCE.md](./REFERENCE.md) — PyDESeq2 parameters, significance thresholds, `create_pseudobulk` parameters, expected outputs, common pitfalls.

---

## Related Skills

- [normalization-strategies](../normalization-strategies/SKILL.md) — preserve raw counts in `layers["raw"]`
- [clustering-and-annotation](../clustering-and-annotation/SKILL.md) — define clusters for meso-level analysis
- [embedding-visualization](../embedding-visualization/SKILL.md) — visualize DE results on UMAP

## Disclaimer

This content is provided "as-is" without warranty of any kind. AI-generated output may contain
inaccuracies. Review and verify all outputs before relying on them.
