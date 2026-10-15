---
name: de-visualization-and-interpretation
description: Visualizing and interpreting differential expression results — volcano plots, MA plots, heatmaps of top DEGs, gene set enrichment analysis (GSEA), and pathway analysis. Accepts any DE results table (CSV, parquet, or DataFrame) with `log2FoldChange` / `pvalue` / `padj` columns, whether produced by pseudobulk-differential-expression, DESeq2, edgeR, limma, or another tool. Always use this skill when the user asks to make a volcano plot, plot DE results, show significant or differentially expressed genes, interpret fold changes, run pathway enrichment, or understand what DE results mean biologically.
metadata:
  version: 0.2.1
  author: elembio
---

# DE Visualization and Interpretation

## Purpose

This skill covers visualizing differential expression results, including volcano plots, MA plots, heatmaps, and gene set enrichment analysis (GSEA). It explains how to interpret and communicate DE findings using standard statistical thresholds and effect-size reasoning.

## Prerequisites

- A DE results DataFrame, CSV, or parquet with `log2FoldChange`, `pvalue`, `padj` columns — produced by [pseudobulk-differential-expression](../pseudobulk-differential-expression/SKILL.md) or any equivalent tool (DESeq2, edgeR, limma, etc.).  If upstream pseudobulk is used, it assumes [cell-quality-control](../cell-quality-control/SKILL.md) and [normalization-strategies](../normalization-strategies/SKILL.md) in turn.

## Key Concepts

### Plot Type Selection

| Plot Type | Shows | Best For |
|-----------|-------|----------|
| **Volcano plot** | log2FC vs -log10(padj) or -log10(pvalue) | Overall DE landscape, identifying outlier genes. Prefer `padj` on Y; use `pvalue` only if most genes have `padj=NaN` |
| **MA plot** | Mean expression vs log2FC | Detecting expression-dependent bias or LFC inflation |
| **Heatmap** | Expression of top DEGs across samples | Pattern discovery and sample clustering |
| **GSEA** | Ranked enrichment of gene sets | Biological pathway interpretation using all genes |
| **Enrichment barplot** | Top significant pathways by -log10(padj) | Summarizing ORA results at a glance |

### Interpretation Guidelines

- **Adjusted p-value threshold**: Use `padj < 0.05` (Benjamini-Hochberg); raw p-values are not corrected for multiple testing.
- **LFC magnitude matters**: `|log2FC| > 1` is a common effect-size filter (2× change), but domain context determines relevance. Statistical significance alone does not imply biological importance.
- **Watch for LFC inflation at low expression**: Genes with very low baseMean often show extreme fold changes. The MA plot exposes this pattern — use a minimum expression filter (e.g., `baseMean > 10`) if needed.

## Quick Start

Use the `plot_volcano()` helper from [IMPLEMENTATION.md §2](./IMPLEMENTATION.md#2-volcano-plot) — copy it into your notebook before calling.

```python
import pandas as pd

results_df = pd.read_csv("results/de/macro_de_results.csv", index_col=0)
print(f"Total genes tested: {len(results_df)}")
print(f"Significant (padj < 0.05): {(results_df['padj'] < 0.05).sum()}")

# plot_volcano — defined in IMPLEMENTATION.md §2
fig = plot_volcano(
    results_df,
    lfc_thresh=0.5,       # default ~1.4× change; adjust to 1.0 (2×) or 1.5 (2.8×) if needed
    pval_thresh=0.05,
    top_genes=10,         # number of genes to label by name
    title="Treatment vs Control\nTargeted transcript panel · PyDESeq2 (pseudobulk by well)",
)
fig.savefig("figures/volcano_plot.pdf", bbox_inches="tight")
```

Style: gray = not significant, red = up, blue = down. Horizontal dashed line at `padj` threshold. Gene labels via `adjustText` (no overlap). Requires `pip install adjustText`.

## Detail Files

- [IMPLEMENTATION.md](./IMPLEMENTATION.md) — parameterized `plot_volcano()` with gene labeling, `plot_ma()` for bias detection, `plot_de_heatmap()` with z-score clustering, GSEA preranked and ORA via `gseapy`, `plot_enrichment_barplot()`, cross-cluster DE comparison heatmap, Markdown report generator.
- [REFERENCE.md](./REFERENCE.md) — plot threshold parameters, GSEA library names, significance thresholds, expected output files, common pitfalls.

## Related Skills

- [pseudobulk-differential-expression](../pseudobulk-differential-expression/SKILL.md) — Run DE analysis
- [clustering-and-annotation](../clustering-and-annotation/SKILL.md) — Cell type annotation
- [embedding-visualization](../embedding-visualization/SKILL.md) — General visualization

## Disclaimer

This content is provided "as-is" without warranty of any kind. AI-generated output may contain
inaccuracies. Review and verify all outputs before relying on them.
