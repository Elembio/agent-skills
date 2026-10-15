# DE Visualization and Interpretation — Reference

## Interpretation Guidelines

| Criterion | Threshold | Notes |
|---|---|---|
| Adjusted p-value | `padj < 0.05` | Benjamini-Hochberg corrected; never use raw `pvalue` |
| Effect size | `\|log2FC\| > 1` | 2× change; domain context determines relevance |
| Base expression filter | `baseMean > 10` | Remove LFC-inflated low-expression genes (visible on MA plot) |

## Plot Threshold Parameters

Default is `lfc_thresh=0.5` (~1.4× change). Tell the user the threshold being used and the resulting gene counts, then ask if they want to adjust.

```python
# Default (show broadly, good starting point for targeted panels)
plot_volcano(results_df, lfc_thresh=0.5, pval_thresh=0.05)

# Standard / publication
plot_volcano(results_df, lfc_thresh=1.0, pval_thresh=0.05)

# Strict publication
plot_volcano(results_df, lfc_thresh=1.5, pval_thresh=0.01)

# Color by padj only (no LFC filter)
plot_volcano(results_df, lfc_thresh=0, pval_thresh=0.05)
```

## GSEA Gene Set Libraries (via `gseapy`)

```python
import gseapy as gp
gp.get_library_name()  # list all available
```

Common choices:

| Library | Description |
|---|---|
| `KEGG_2021_Human` | Metabolic and signaling pathways |
| `GO_Biological_Process_2021` | Gene Ontology biological processes |
| `MSigDB_Hallmark_2020` | 50 curated hallmark gene sets |
| `Reactome_2022` | Reaction-level pathway database |
| `WikiPathways_2021_Human` | Community-curated pathways |

## Expected Output Files

| File | Description |
|---|---|
| `figures/volcano_plot.pdf` | Volcano plot |
| `figures/ma_plot.pdf` | MA plot |
| `figures/de_heatmap.pdf` | Heatmap of top DEGs |
| `figures/enrichment_up.pdf` | Upregulated pathway barplot |
| `figures/enrichment_down.pdf` | Downregulated pathway barplot |
| `figures/de_cluster_comparison.pdf` | Cross-cluster log2FC heatmap |
| `results/gsea/` | GSEA output directory |
| `results/de/de_report.md` | Markdown summary report |

## Common Pitfalls

1. **Ignoring effect size** — Statistical significance (`padj`) alone does not imply biological importance. Always consider `|log2FC|` magnitude alongside `padj`.

2. **LFC inflation at low expression** — Genes with `baseMean < 10` often have extreme, unreliable fold changes. Inspect the MA plot and apply a minimum expression filter when necessary.

3. **Over-interpreting pathways** — Pathway databases have significant term overlap. Treat enrichment as hypothesis generation; validate key findings independently.

4. **Gene name mismatches** — Ensure gene identifiers use HGNC symbols. ENSEMBL IDs or aliases will fail to match most pathway databases.

5. **Too few permutations in GSEA** — Use `permutation_num=100` only for exploration; increase to 1000 for publication-quality FDR estimates.
