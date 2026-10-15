# Pseudobulk Differential Expression — Reference

## PyDESeq2 Parameters

| Parameter | Default | Description |
|---|---|---|
| `design_factors` | required | Factor(s) for the design matrix. String or list (list = multiple covariates). |
| `refit_cooks` | `True` | Refit genes flagged as Cook's distance outliers. Recommended. |
| `min_mu` | `0.5` | Minimum fitted mean for shrinkage estimation. |
| `contrast` | required | `[factor, numerator, denominator]` — e.g. `["condition", "Treatment", "Control"]`. |

## `create_pseudobulk` Parameters

Defined in [IMPLEMENTATION.md §4](./IMPLEMENTATION.md#4-create_pseudobulk--full-implementation).

| Parameter | Type | Default | Description |
|---|---|---|---|
| `adata` | AnnData | required | Single-cell dataset |
| `sample_col` | str | required | `.obs` column for biological samples (e.g. `"Well"`) |
| `condition_col` | str | required | `.obs` column for conditions |
| `min_cells` | int | `10` | Minimum cells per sample; samples below are excluded |
| `layer` | str | `"raw"` | Layer with raw integer counts — always use `"raw"` |

## Significance Thresholds

| Stringency | `padj` | `\|log2FC\|` | Typical use |
|---|---|---|---|
| Strict | < 0.01 | > 1.5 | Publication, follow-up validation |
| Moderate | < 0.05 | > 1.0 | Standard exploratory analysis |
| Lenient | < 0.10 | > 0.5 | Hypothesis generation, small experiments |

```python
thresholds = {
    "strict":   {"padj": 0.01, "lfc": 1.5},
    "moderate": {"padj": 0.05, "lfc": 1.0},
    "lenient":  {"padj": 0.10, "lfc": 0.5},
}
for name, t in thresholds.items():
    n = ((results_df["padj"] < t["padj"]) &
         (results_df["log2FoldChange"].abs() > t["lfc"])).sum()
    print(f"{name}: {n} DEGs")
```

## Macro vs Meso Analysis

| Level | Aggregation | What it answers | When to use |
|---|---|---|---|
| **Macro** | All cells per well/sample | Overall condition effects across the full dataset | First-pass DE, treatment vs control |
| **Meso** | Cells per cluster per well | Cell-type-specific condition effects | After clustering; reveals heterogeneous responses |

## Expected Outputs

| Object | Description |
|---|---|
| `pseudobulk` | AnnData (samples × genes) with summed raw counts |
| `results_df` | DataFrame: `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj` |
| `significant` | Filtered `results_df` by chosen padj + log2FC thresholds |
| `meso_results` | Dict keyed by cluster; each value contains `results`, `n_cells_total`, `n_samples`, `n_sig` |

Saved files:

| Path | Contents |
|---|---|
| `results/de/macro_de_results.csv` | All genes, macro-level |
| `results/de/macro_de_significant.csv` | Significant DEGs only |
| `results/de/meso_{cluster}_de_results.csv` | Per-cluster results |
| `results/de/de_summary.csv` | Summary counts across levels |
| `figures/de_validation.pdf` | Boxplots of top DE genes |

## Common Pitfalls

1. **Using normalized counts** — DESeq2 performs its own normalization internally. Feeding log-normalized or scaled values produces incorrect size factor estimates and inflated significance. Always use `layers["raw"]`.

2. **Insufficient replicates** — Fewer than 2 replicates per condition makes dispersion estimation impossible; fewer than 3 makes it unreliable. The meso-level loop skips clusters that fail this check.

3. **Too few cells per pseudobulk sample** — Aggregates with < 10 cells are dominated by noise. Increase `min_cells` if you have very small clusters.

4. **Batch–condition confounding** — If all wells from batch A are "Control" and batch B are "Treatment", the batch covariate cannot be separated from condition. Check your design matrix before running.

5. **Using raw p-values** — Always filter on `padj` (Benjamini–Hochberg adjusted). Raw `pvalue` columns lead to a high false discovery rate.

6. **Cluster labels with special characters** — Slashes or spaces in cluster names break file paths. Sanitize cluster labels upstream.
