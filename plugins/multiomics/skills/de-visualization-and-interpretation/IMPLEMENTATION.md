# DE Visualization and Interpretation — Implementation

## Table of Contents

1. [Load DE Results](#1-load-de-results)
2. [Volcano Plot](#2-volcano-plot)
3. [MA Plot](#3-ma-plot)
4. [Heatmap of Top DEGs](#4-heatmap-of-top-degs)
5. [GSEA Preranked](#5-gsea-preranked)
6. [Over-Representation Analysis (ORA)](#6-over-representation-analysis-ora)
7. [Enrichment Barplot](#7-enrichment-barplot)
8. [Cross-Cluster DE Comparison](#8-cross-cluster-de-comparison)
9. [Markdown Report Generator](#9-markdown-report-generator)
10. [Parameters](#10-parameters)
11. [Expected Outputs](#11-expected-outputs)
12. [Common Pitfalls](#12-common-pitfalls)

---

## 1. Load DE Results

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load DE results produced by pseudobulk-differential-expression skill
results_df = pd.read_csv("results/de/macro_de_results.csv", index_col=0)

# Inspect summary counts
print(f"Total genes tested:                      {len(results_df)}")
print(f"Significant (padj < 0.05):               {(results_df['padj'] < 0.05).sum()}")
print(f"Upregulated   (log2FC > 1, padj < 0.05): "
      f"{((results_df['padj'] < 0.05) & (results_df['log2FoldChange'] > 1)).sum()}")
print(f"Downregulated (log2FC < -1, padj < 0.05): "
      f"{((results_df['padj'] < 0.05) & (results_df['log2FoldChange'] < -1)).sum()}")
```

---

## 2. Volcano Plot

Style conventions (canonical — always follow this):
- **Y-axis: `-log10(padj)`** (adjusted p-value). Use `-log10(pvalue)` only if the majority of genes have `padj=NaN` (e.g. too few replicates for BH correction).
- Colors: gray `#BBBBBB` = not significant, red `#E41A1C` = up, blue `#377EB8` = down
- Horizontal dashed line at `padj` threshold; no vertical lines
- No vertical LFC threshold lines by default (set `lfc_thresh=0` to color by padj only)
- Seaborn whitegrid background
- Legend shows counts: `Up (padj < 0.05, n=N)`
- Gene labels via `adjustText` to avoid overlap

```python
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text

def plot_volcano(results_df, lfc_thresh=0.5, pval_thresh=0.05,
                 top_genes=10, title="Volcano Plot"):
    """
    Volcano plot of DE results.

    Parameters
    ----------
    results_df : pd.DataFrame
        DE results with columns log2FoldChange, padj.
    lfc_thresh : float
        Absolute log2FC threshold for significance coloring. Default 0.5 (~1.4× change).
        Set to 0 to color by padj only. Common alternatives: 1.0 (2×), 1.5 (2.8×).
    pval_thresh : float
        Adjusted p-value threshold.
    top_genes : int
        Number of top significant genes to label by name.
    title : str
        Plot title (multi-line supported with \\n).
    """
    sns.set_style("whitegrid")

    df = results_df.copy()
    df["-log10(padj)"] = -np.log10(df["padj"].clip(lower=1e-300))

    up_mask   = (df["padj"] < pval_thresh) & (df["log2FoldChange"] >  lfc_thresh)
    down_mask = (df["padj"] < pval_thresh) & (df["log2FoldChange"] < -lfc_thresh)
    df["significant"] = "Not sig."
    df.loc[up_mask,   "significant"] = "Up"
    df.loc[down_mask, "significant"] = "Down"

    colors = {"Not sig.": "#BBBBBB", "Up": "#E41A1C", "Down": "#377EB8"}
    labels = {
        "Not sig.": f"Not sig. (padj ≥ {pval_thresh})",
        "Up":       f"Up (padj < {pval_thresh}, n={up_mask.sum()})",
        "Down":     f"Down (padj < {pval_thresh}, n={down_mask.sum()})",
    }

    fig, ax = plt.subplots(figsize=(10, 8))
    for sig_type in ["Not sig.", "Down", "Up"]:
        mask = df["significant"] == sig_type
        ax.scatter(
            df.loc[mask, "log2FoldChange"],
            df.loc[mask, "-log10(padj)"],
            c=colors[sig_type],
            label=labels[sig_type],
            alpha=0.6,
            s=20,
            linewidths=0,
        )

    ax.axhline(-np.log10(pval_thresh), color="gray", linestyle="--", alpha=0.6,
               label=f"p = {pval_thresh}")

    # Gene labels with adjustText to avoid overlap
    sig = df[up_mask | down_mask]
    if len(sig) > 0:
        top = sig.nlargest(top_genes, "-log10(padj)")
        texts = [
            ax.text(row["log2FoldChange"], row["-log10(padj)"], gene, fontsize=8)
            for gene, row in top.iterrows()
        ]
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="gray", lw=0.5))

    ax.set_xlabel("log₂ fold change (treatment / control)")
    ax.set_ylabel("-log₁₀(Adjusted P-value)")
    ax.set_title(title, fontsize=11)
    ax.legend(loc="upper right", framealpha=0.9, fontsize=9)
    plt.tight_layout()
    return fig


lfc_thresh = 0.5   # default: ~1.4× change; adjust to 1.0 (2×) or 1.5 (2.8×) if needed
pval_thresh = 0.05

print(f"Volcano thresholds: |log2FC| > {lfc_thresh}, padj < {pval_thresh}")
print(f"  Up:   {((results_df['padj'] < pval_thresh) & (results_df['log2FoldChange'] >  lfc_thresh)).sum()} genes")
print(f"  Down: {((results_df['padj'] < pval_thresh) & (results_df['log2FoldChange'] < -lfc_thresh)).sum()} genes")
print("Adjust lfc_thresh or pval_thresh if the number of highlighted genes is too high or too low.")

fig = plot_volcano(
    results_df,
    lfc_thresh=lfc_thresh,
    pval_thresh=pval_thresh,
    top_genes=10,
    title="Treatment vs Control\nTargeted transcript panel · PyDESeq2 (pseudobulk by well)",
)
fig.savefig("figures/volcano_plot.pdf", bbox_inches="tight")
plt.show()
```

---

## 3. MA Plot

```python
def plot_ma(results_df, lfc_thresh=1.0, pval_thresh=0.05, title="MA Plot"):
    """
    MA plot: mean expression vs log fold change.

    Useful for detecting LFC inflation at low baseMean values.

    Parameters
    ----------
    results_df : pd.DataFrame
        DE results with columns baseMean, log2FoldChange, padj.
    lfc_thresh : float
        Absolute log2FC threshold for significance coloring.
    pval_thresh : float
        Adjusted p-value threshold.
    title : str
        Plot title.
    """
    df = results_df.copy()
    df["log10(baseMean)"] = np.log10(df["baseMean"].clip(lower=1))

    df["significant"] = "Not significant"
    df.loc[(df["padj"] < pval_thresh) & (df["log2FoldChange"] > lfc_thresh), "significant"] = "Up"
    df.loc[(df["padj"] < pval_thresh) & (df["log2FoldChange"] < -lfc_thresh), "significant"] = "Down"

    colors = {"Not significant": "#BBBBBB", "Up": "#E41A1C", "Down": "#377EB8"}

    fig, ax = plt.subplots(figsize=(10, 8))
    for sig_type in ["Not significant", "Down", "Up"]:
        mask = df["significant"] == sig_type
        ax.scatter(
            df.loc[mask, "log10(baseMean)"],
            df.loc[mask, "log2FoldChange"],
            c=colors[sig_type],
            label=f"{sig_type} ({mask.sum()})",
            alpha=0.6,
            s=10,
        )

    ax.axhline(0, color="black", linestyle="-", alpha=0.3)
    ax.axhline(lfc_thresh, color="gray", linestyle="--", alpha=0.5)
    ax.axhline(-lfc_thresh, color="gray", linestyle="--", alpha=0.5)

    ax.set_xlabel("log₁₀(Mean Expression)")
    ax.set_ylabel("log₂(Fold Change)")
    ax.set_title(title)
    ax.legend(loc="upper right")
    plt.tight_layout()
    return fig


fig = plot_ma(results_df)
plt.savefig("figures/ma_plot.pdf")
plt.show()
```

---

## 4. Heatmap of Top DEGs

```python
def plot_de_heatmap(adata, results_df, sample_col, condition_col,
                    n_genes=50, layer="raw"):
    """
    Clustered heatmap of top DE genes using pseudobulk z-scores.

    Parameters
    ----------
    adata : AnnData
        Annotated data object.
    results_df : pd.DataFrame
        DE results with log2FoldChange and padj columns.
    sample_col : str
        adata.obs column identifying biological samples (e.g. "Well").
    condition_col : str
        adata.obs column identifying condition groups.
    n_genes : int
        Number of top significant genes (by |log2FC|) to display.
    layer : str
        adata.layers key for expression matrix, or "X" to use adata.X.
    """
    sig = results_df[results_df["padj"] < 0.05].copy()
    sig["abs_lfc"] = sig["log2FoldChange"].abs()
    top_genes = sig.nlargest(n_genes, "abs_lfc").index.tolist()
    top_genes = [g for g in top_genes if g in adata.var_names]

    if not top_genes:
        print("No significant genes found.")
        return

    samples = adata.obs[sample_col].unique()
    X = adata.layers.get(layer, adata.X) if layer != "X" else adata.X
    if hasattr(X, "toarray"):
        X = X.toarray()

    gene_idx = [adata.var_names.get_loc(g) for g in top_genes]
    pb_data = {}
    for sample in samples:
        mask = adata.obs[sample_col] == sample
        pb_data[sample] = X[mask][:, gene_idx].sum(axis=0)

    pb_df = pd.DataFrame(pb_data, index=top_genes).T
    pb_log = np.log1p(pb_df)
    pb_scaled = (pb_log - pb_log.mean()) / pb_log.std()

    sample_conditions = adata.obs.groupby(sample_col)[condition_col].first()
    conditions = sample_conditions[pb_scaled.index].values
    unique_conditions = pd.unique(conditions)
    condition_colors = dict(zip(unique_conditions,
                                sns.color_palette("Set1", len(unique_conditions))))
    col_colors = [condition_colors[c] for c in conditions]

    g = sns.clustermap(
        pb_scaled.T,
        cmap="RdBu_r",
        center=0,
        col_colors=col_colors,
        figsize=(12, 10),
        dendrogram_ratio=0.15,
        cbar_pos=(0.02, 0.8, 0.03, 0.15),
    )
    g.ax_heatmap.set_ylabel("Gene")
    g.ax_heatmap.set_xlabel("Sample")
    plt.title("Top DE Genes Heatmap", pad=100)

    handles = [plt.Rectangle((0, 0), 1, 1, color=color)
               for color in condition_colors.values()]
    plt.legend(handles, condition_colors.keys(),
               loc="upper left", bbox_to_anchor=(1.05, 1))

    plt.savefig("figures/de_heatmap.pdf", bbox_inches="tight")
    plt.show()


# plot_de_heatmap(adata, results_df, sample_col="Well", condition_col="experimental_condition")
```

---

## 5. GSEA Preranked

```python
# pip install gseapy
import gseapy as gp

# Rank all genes by signed -log10(p-value)
results_df["rank_score"] = (
    -np.log10(results_df["pvalue"]) * np.sign(results_df["log2FoldChange"])
)
ranked_genes = results_df["rank_score"].sort_values(ascending=False)

print(f"Top upregulated:   {ranked_genes.head(5).index.tolist()}")
print(f"Top downregulated: {ranked_genes.tail(5).index.tolist()}")

gene_sets = [
    "KEGG_2021_Human",
    "GO_Biological_Process_2021",
    "MSigDB_Hallmark_2020",
]

try:
    gsea_results = gp.prerank(
        rnk=ranked_genes,
        gene_sets=gene_sets[0],
        processes=4,
        permutation_num=1000,   # Use 100 for exploration, 1000 for publication
        outdir="results/gsea",
        format="pdf",
        seed=42,
    )

    gsea_df = gsea_results.res2d.sort_values("NES", ascending=False)
    print("Top enriched (positive NES = upregulated):")
    print(gsea_df[["Term", "NES", "FDR q-val"]].head(10))
    print("\nTop depleted (negative NES = downregulated):")
    print(gsea_df[["Term", "NES", "FDR q-val"]].tail(10))

except Exception as e:
    print(f"GSEA error: {e}")
    print("Install with: pip install gseapy")
```

---

## 6. Over-Representation Analysis (ORA)

```python
significant_up = results_df[
    (results_df["padj"] < 0.05) & (results_df["log2FoldChange"] > 1)
].index.tolist()

significant_down = results_df[
    (results_df["padj"] < 0.05) & (results_df["log2FoldChange"] < -1)
].index.tolist()

print(f"Genes for ORA — Up: {len(significant_up)}, Down: {len(significant_down)}")

if significant_up:
    try:
        enr_up = gp.enrichr(
            gene_list=significant_up,
            gene_sets=["KEGG_2021_Human", "GO_Biological_Process_2021"],
            organism="Human",
            outdir="results/enrichr_up",
        )
        print("\nTop pathways (upregulated genes):")
        print(enr_up.results[["Term", "Adjusted P-value", "Overlap"]].head(10))
    except Exception as e:
        print(f"Enrichr error: {e}")

if significant_down:
    try:
        enr_down = gp.enrichr(
            gene_list=significant_down,
            gene_sets=["KEGG_2021_Human", "GO_Biological_Process_2021"],
            organism="Human",
            outdir="results/enrichr_down",
        )
        print("\nTop pathways (downregulated genes):")
        print(enr_down.results[["Term", "Adjusted P-value", "Overlap"]].head(10))
    except Exception as e:
        print(f"Enrichr error: {e}")
```

---

## 7. Enrichment Barplot

```python
def plot_enrichment_barplot(enr_results, n_terms=15, title="Pathway Enrichment"):
    """
    Horizontal barplot of top enrichment results.

    Parameters
    ----------
    enr_results : gseapy Enrichr result object
        Returned by gp.enrichr().
    n_terms : int
        Number of top terms to display (sorted by adjusted p-value).
    title : str
        Plot title.
    """
    df = enr_results.results.copy()
    df = df.nsmallest(n_terms, "Adjusted P-value")
    df["-log10(padj)"] = -np.log10(df["Adjusted P-value"].clip(lower=1e-50))
    df = df.sort_values("-log10(padj)")
    df["Term_short"] = df["Term"].str[:50]

    fig, ax = plt.subplots(figsize=(10, 8))
    colors = plt.cm.viridis(df["-log10(padj)"] / df["-log10(padj)"].max())
    ax.barh(df["Term_short"], df["-log10(padj)"], color=colors)
    ax.axvline(-np.log10(0.05), color="red", linestyle="--", alpha=0.5, label="padj = 0.05")
    ax.set_xlabel("-log₁₀(Adjusted P-value)")
    ax.set_title(title)
    ax.legend()
    plt.tight_layout()
    return fig


# fig = plot_enrichment_barplot(enr_up, title="Upregulated Pathways")
# plt.savefig("figures/enrichment_up.pdf"); plt.show()
```

---

## 8. Cross-Cluster DE Comparison

```python
def plot_de_comparison_across_clusters(meso_results, n_genes=20):
    """
    Heatmap comparing log2FC values across cell types/clusters.

    Parameters
    ----------
    meso_results : dict
        Mapping of cluster name → dict with key "results" (DE DataFrame).
        Produced by running DE per cluster in the pseudobulk skill.
    n_genes : int
        Number of top genes to show (by maximum absolute log2FC across clusters).
    """
    all_genes = set()
    for cluster, data in meso_results.items():
        sig_genes = data["results"][data["results"]["padj"] < 0.05].index
        all_genes.update(sig_genes)

    if not all_genes:
        print("No significant genes across any cluster.")
        return

    lfc_matrix = pd.DataFrame(index=list(all_genes))
    for cluster, data in meso_results.items():
        lfc_matrix[cluster] = data["results"].reindex(lfc_matrix.index)["log2FoldChange"].values

    lfc_matrix["max_abs"] = lfc_matrix.abs().max(axis=1)
    top_genes = lfc_matrix.nlargest(n_genes, "max_abs").index
    lfc_plot = lfc_matrix.loc[top_genes].drop(columns=["max_abs"])

    fig, ax = plt.subplots(figsize=(12, 10))
    sns.heatmap(
        lfc_plot,
        cmap="RdBu_r",
        center=0,
        ax=ax,
        annot=True,
        fmt=".1f",
        cbar_kws={"label": "log₂(Fold Change)"},
    )
    ax.set_xlabel("Cluster")
    ax.set_ylabel("Gene")
    ax.set_title("Differential Expression Across Cell Types")
    plt.tight_layout()
    plt.savefig("figures/de_cluster_comparison.pdf")
    plt.show()


# plot_de_comparison_across_clusters(meso_results)
```

---

## 9. Markdown Report Generator

```python
def generate_de_report(results_df, condition_contrast,
                       output_file="results/de/de_report.md"):
    """
    Write a markdown summary report of DE results.

    Parameters
    ----------
    results_df : pd.DataFrame
        DE results DataFrame.
    condition_contrast : str
        Human-readable label for the contrast (e.g. "Treatment vs Control").
    output_file : str
        Path to write the report.
    """
    sig = results_df[results_df["padj"] < 0.05]
    up = sig[sig["log2FoldChange"] > 1]
    down = sig[sig["log2FoldChange"] < -1]

    report = f"""# Differential Expression Report

## Analysis: {condition_contrast}

### Summary Statistics

- **Total genes tested**: {len(results_df)}
- **Significant genes (padj < 0.05)**: {len(sig)}
- **Upregulated (log2FC > 1)**: {len(up)}
- **Downregulated (log2FC < -1)**: {len(down)}

### Top Upregulated Genes

| Gene | log2FC | padj |
|------|--------|------|
"""
    for gene, row in up.nlargest(10, "log2FoldChange").iterrows():
        report += f"| {gene} | {row['log2FoldChange']:.2f} | {row['padj']:.2e} |\n"

    report += """
### Top Downregulated Genes

| Gene | log2FC | padj |
|------|--------|------|
"""
    for gene, row in down.nsmallest(10, "log2FoldChange").iterrows():
        report += f"| {gene} | {row['log2FoldChange']:.2f} | {row['padj']:.2e} |\n"

    report += """
### Figures

- `figures/volcano_plot.pdf` — Volcano plot
- `figures/ma_plot.pdf` — MA plot
- `figures/de_heatmap.pdf` — Heatmap of top DEGs
- `figures/enrichment_up.pdf` — Upregulated pathway enrichment

### Methods

Differential expression analysis was performed using pseudobulk aggregation
and PyDESeq2. Cells were aggregated per sample and a negative binomial model
was fit. Adjusted p-values were calculated using the Benjamini-Hochberg method.
"""
    import os
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, "w") as f:
        f.write(report)
    print(f"Report saved to {output_file}")


generate_de_report(results_df, "Treatment vs Control")
```

---

## 10. Parameters

### Volcano and MA Plot Thresholds

```python
# Loosen thresholds for exploratory work
plot_volcano(results_df, lfc_thresh=0.5, pval_thresh=0.1)

# Stricter thresholds for publication
plot_volcano(results_df, lfc_thresh=1.5, pval_thresh=0.01)
```

### GSEA Gene Set Libraries

```python
# List all available libraries
gp.get_library_name()

# Common choices:
# - KEGG_2021_Human
# - GO_Biological_Process_2021
# - MSigDB_Hallmark_2020
# - Reactome_2022
# - WikiPathways_2021_Human
```

### Expression Filter Before Heatmap

```python
# Remove low-expression genes before downstream analysis
results_df_filtered = results_df[results_df["baseMean"] > 10]
```

---

## 11. Expected Outputs

| File | Description |
|------|-------------|
| `figures/volcano_plot.pdf` | Volcano plot |
| `figures/ma_plot.pdf` | MA plot |
| `figures/de_heatmap.pdf` | Heatmap of top DEGs |
| `figures/enrichment_up.pdf` | Upregulated pathway barplot |
| `figures/enrichment_down.pdf` | Downregulated pathway barplot |
| `figures/de_cluster_comparison.pdf` | Cross-cluster log2FC heatmap |
| `results/gsea/` | GSEA output directory |
| `results/de/de_report.md` | Markdown summary report |

---

## 12. Common Pitfalls

1. **Ignoring effect size**: Statistical significance (padj) alone does not imply biological importance. Always consider `|log2FC|` magnitude alongside padj.

2. **LFC inflation at low expression**: Genes with `baseMean < 10` often have extreme, unreliable fold changes. Inspect the MA plot and apply a minimum expression filter when necessary.

3. **Over-interpreting pathways**: Pathway databases have significant term overlap. Treat enrichment as hypothesis generation; validate key findings independently.

4. **Gene name mismatches**: Ensure gene identifiers use HGNC symbols. ENSEMBL IDs or aliases will fail to match most pathway databases.

5. **Too few permutations in GSEA**: Use `permutation_num=100` only for exploration; increase to 1000 for publication-quality FDR estimates.
