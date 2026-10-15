# Embedding Visualization — Implementation

## Table of Contents

1. [Setup Generic Metadata](#1-setup-generic-metadata)
2. [Basic UMAP Visualization](#2-basic-umap-visualization)
3. [Imaging / IF embedding colors](#3-imaging--if-embedding-colors)
4. [Multi-Panel Comparison](#4-multi-panel-comparison)
5. [Gene Expression Overlay](#5-gene-expression-overlay)
6. [Faceting by Condition](#6-faceting-by-condition)
7. [Continuous Variable Visualization](#7-continuous-variable-visualization)
8. [Custom Color Palettes](#8-custom-color-palettes)
9. [Density and Contour Plots](#9-density-and-contour-plots)
10. [Publication-Quality Figures](#10-publication-quality-figures)
11. [Interactive Plotly Visualization](#11-interactive-plotly-visualization)
12. [Condition Comparison Side-by-Side](#12-condition-comparison-side-by-side)
13. [Parameters Reference and Expected Outputs](#13-parameters-reference-and-expected-outputs) → [REFERENCE.md](./REFERENCE.md)
14. [Common Pitfalls](#14-common-pitfalls)

---

## 1. Setup Generic Metadata

Map experiment-specific column names to generic terms so the rest of the recipes work without modification.

```python
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

adata = adata_with_umap.copy()

def add_generic_metadata(adata, mappings):
    """Add generic metadata columns from specific columns."""
    for generic_name, specific_name in mappings.items():
        if specific_name in adata.obs.columns:
            adata.obs[generic_name] = adata.obs[specific_name]
        else:
            print(f"Warning: {specific_name} not found, skipping {generic_name}")
    return adata

# Customize these mappings for your experiment
metadata_mappings = {
    "experimental_condition": "WellLabel",  # sample-sheet name; or "condition" after experiment-design
    "sample_id": "Well",
    "technical_batch": "run_id",            # or "batch_name"
}

adata = add_generic_metadata(adata, metadata_mappings)
```

---

## 2. Basic UMAP Visualization

```python
# Single categorical variable
sc.pl.umap(
    adata,
    color="experimental_condition",
    title="Cells by Experimental Condition",
    legend_loc="on data",     # 'right margin', 'on data', None
    legend_fontsize=8,
    size=10,                  # point size
    frameon=False,
    save="_by_condition.pdf", # saves to figures/
)
```


---

## 3. Imaging / IF embedding colors

Use when `adata` is a profiling AnnData from [imaging-pipeline](../imaging-pipeline/SKILL.md) with `X_umap` already computed via [dimensionality-reduction](../dimensionality-reduction/SKILL.md). Prefer `obs` metrics over assuming stain names are gene-like `var_names`.

```python
import scanpy as sc

# Technical structure (plate position) and sample name together
sc.pl.umap(adata, color="Well", frameon=False)
label = "condition" if "condition" in adata.obs.columns else "WellLabel"
sc.pl.umap(adata, color=label, frameon=False)

# Derived localization / S/B written by imaging-pipeline prep
sb = [c for c in adata.obs.columns if c.startswith("sb_")]
nuc = [c for c in adata.obs.columns if c.startswith("nucfrac_")]
if sb:
    sc.pl.umap(adata, color=sb[:4], ncols=2, cmap="viridis", frameon=False)
if nuc:
    sc.pl.umap(adata, color=nuc[:4], ncols=2, cmap="coolwarm", vmin=0, vmax=1, frameon=False)

# Condition / treatment (from experiment-design) — overview only;
# treatment hits → well-level contrasts in imaging-pipeline, not this plot alone
if "condition" in adata.obs.columns:
    sc.pl.umap(adata, color="condition", frameon=False)
```

---

## 4. Multi-Panel Comparison

```python
import scipy.sparse as sp

fig, axes = plt.subplots(2, 3, figsize=(18, 12))

# Row 1: sample name / derived factors (WellLabel, not well_label)
cond = next(
    (c for c in ("experimental_condition", "condition", "WellLabel") if c in adata.obs.columns),
    None,
)
for col_key, ax, title in [
    (cond,                     axes[0, 0], "Sample / condition"),
    ("condition_level",        axes[0, 1], "Level/Dose"),
    ("temporal_variable",      axes[0, 2], "Time"),
]:
    if col_key and col_key in adata.obs.columns:
        sc.pl.umap(adata, color=col_key, ax=ax, show=False, title=title)
    else:
        ax.axis("off")

# Row 2: plate position (Well) vs other technical variables
sc.pl.umap(adata, color="Well", ax=axes[1, 0], show=False, title="Well")

if "technical_batch" in adata.obs.columns:
    sc.pl.umap(adata, color="technical_batch", ax=axes[1, 1],
               show=False, title="Technical Batch")
else:
    axes[1, 1].axis("off")

# Total counts quality metric
if "total_counts" not in adata.obs.columns:
    X = adata.layers.get("raw", adata.layers.get("raw_counts", adata.X))
    if sp.issparse(X):
        adata.obs["total_counts"] = np.asarray(X.sum(axis=1)).ravel()
    else:
        adata.obs["total_counts"] = X.sum(axis=1)

sc.pl.umap(adata, color="total_counts", ax=axes[1, 2],
           show=False, title="Total Counts", cmap="viridis")

plt.tight_layout()
plt.savefig("figures/umap_overview.pdf", dpi=150, bbox_inches="tight")
plt.show()
```

---

## 5. Gene Expression Overlay

```python
genes_of_interest = ["GAPDH", "MYC", "TP53", "ACTB"]
genes_present = [g for g in genes_of_interest if g in adata.var_names]

if genes_present:
    sc.pl.umap(
        adata,
        color=genes_present,
        ncols=2,
        cmap="viridis",
        vmin=0,
        use_raw=False,  # use normalized values
        title=[f"{g} Expression" for g in genes_present],
    )
else:
    print("No genes of interest found in data")

# Top variable genes
sc.pl.umap(adata, color=adata.var_names[:4].tolist(), ncols=2, cmap="viridis")
```

---

## 6. Faceting by Condition

```python
def plot_umap_split_by(adata, split_by, color_by, ncols=3, **kwargs):
    """Plot UMAP faceted by a categorical variable."""
    categories = adata.obs[split_by].unique()
    nrows = int(np.ceil(len(categories) / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 5 * nrows))
    axes = axes.flatten() if len(categories) > 1 else [axes]

    for ax, cat in zip(axes, categories):
        mask = adata.obs[split_by] == cat
        sc.pl.umap(
            adata[mask],
            color=color_by,
            ax=ax,
            show=False,
            title=f"{split_by}={cat}",
            **kwargs,
        )

    for ax in axes[len(categories):]:
        ax.axis("off")

    plt.tight_layout()
    return fig


if "experimental_condition" in adata.obs.columns:
    fig = plot_umap_split_by(
        adata,
        split_by="experimental_condition",
        color_by="Well",
        ncols=3,
    )
    plt.savefig("figures/umap_by_condition.pdf")
    plt.show()
```

---

## 7. Continuous Variable Visualization

```python
# Cell area — check obsm["cell_geometry"] first, then obs["Area"]
def _area_for_plot(adata):
    if "cell_geometry" in adata.obsm and "Area" in adata.obsm["cell_geometry"].columns:
        return adata.obsm["cell_geometry"]["Area"]
    if "Area" in adata.obs.columns:
        return adata.obs["Area"]
    return None


area_s = _area_for_plot(adata)
if area_s is not None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    adata.obs["_plot_Area"] = area_s.values
    sc.pl.umap(adata, color="_plot_Area", ax=axes[0], show=False,
               cmap="YlOrRd", title="Cell Area")
    adata.obs["log_Area"] = np.log1p(adata.obs["_plot_Area"])
    sc.pl.umap(adata, color="log_Area", ax=axes[1], show=False,
               cmap="YlOrRd", title="Cell Area (log)")
    del adata.obs["_plot_Area"]
    plt.tight_layout()
    plt.show()

# Assignment rate
if "assignment_rate" in adata.obs.columns:
    sc.pl.umap(adata, color="assignment_rate", cmap="RdYlGn",
               title="Assignment Rate", vmin=0, vmax=1)
```

### Ordered numeric levels (timepoints, doses)

A numeric column renders as a continuous colorbar, which misreads unevenly spaced
levels — 0, 5, 10, 30, 60, 120 min squeezes everything under 30 min into the bottom
quarter of the bar. Cast to an ordered categorical for equal visual weight and a
discrete legend:

```python
levels = sorted(adata.obs["timepoint_min"].unique())
adata.obs["timepoint_cat"] = pd.Categorical(adata.obs["timepoint_min"], categories=levels)
sc.pl.umap(adata, color="timepoint_cat", palette="plasma")
```

Derive the column before subsetting, or it will be missing on the subset.

---

## 8. Custom Color Palettes

```python
condition_colors = {
    "Control":     "#1f77b4",
    "Treatment_A": "#ff7f0e",
    "Treatment_B": "#2ca02c",
    "Treatment_C": "#d62728",
}

if "experimental_condition" in adata.obs.columns:
    conditions = adata.obs["experimental_condition"].unique()
    palette = [condition_colors.get(c, "#7f7f7f") for c in conditions]
    sc.pl.umap(
        adata,
        color="experimental_condition",
        palette=palette,
        title="Experimental Conditions",
    )

# Seaborn palettes work directly
sc.pl.umap(adata, color="Well", palette="tab20")
```

---

## 9. Density and Contour Plots

```python
from scipy.stats import gaussian_kde

# Built-in scanpy density by group
sc.pl.embedding_density(
    adata,
    basis="umap",
    groupby="experimental_condition" if "experimental_condition" in adata.obs.columns else "Well",
)


def plot_umap_with_density(adata, color_by=None, levels=10):
    """Scatter UMAP overlaid with Gaussian KDE contour lines."""
    fig, ax = plt.subplots(figsize=(10, 8))
    x = adata.obsm["X_umap"][:, 0]
    y = adata.obsm["X_umap"][:, 1]

    if color_by and color_by in adata.obs.columns:
        cats = adata.obs[color_by].astype("category")
        for cat in cats.cat.categories:
            mask = cats == cat
            ax.scatter(x[mask], y[mask], label=cat, alpha=0.5, s=5)
        ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    else:
        ax.scatter(x, y, alpha=0.5, s=5, c="steelblue")

    z = gaussian_kde(np.vstack([x, y]))(np.vstack([x, y]))
    ax.tricontour(x, y, z, levels=levels, colors="black", linewidths=0.5, alpha=0.5)

    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title("UMAP with Density Contours")
    plt.tight_layout()
    return fig


plot_umap_with_density(adata, color_by="Well")
plt.show()
```

---

## 10. Publication-Quality Figures

```python
plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.sans-serif":   ["Arial", "DejaVu Sans"],
    "font.size":         10,
    "axes.titlesize":    12,
    "axes.labelsize":    10,
    "xtick.labelsize":   8,
    "ytick.labelsize":   8,
    "legend.fontsize":   8,
    "figure.dpi":        150,
    "savefig.dpi":       300,
    "savefig.format":    "pdf",
})

fig, axes = plt.subplots(1, 3, figsize=(12, 4))

# Panel A: condition or well
cond_col = "experimental_condition" if "experimental_condition" in adata.obs.columns else "Well"
sc.pl.umap(adata, color=cond_col, ax=axes[0], show=False,
           title="", frameon=False, legend_loc="right margin")
axes[0].set_title("A. Experimental Condition", loc="left", fontweight="bold")

# Panel B: key gene
if "MYC" in adata.var_names:
    sc.pl.umap(adata, color="MYC", ax=axes[1], show=False,
               cmap="viridis", title="", frameon=False)
    axes[1].set_title("B. MYC Expression", loc="left", fontweight="bold")
else:
    axes[1].axis("off")

# Panel C: well / technical
sc.pl.umap(adata, color="Well", ax=axes[2], show=False,
           title="", frameon=False, legend_loc="none")
axes[2].set_title("C. Well", loc="left", fontweight="bold")

plt.tight_layout()
plt.savefig("figures/figure1_umap.pdf", bbox_inches="tight")
plt.savefig("figures/figure1_umap.png", bbox_inches="tight")
plt.show()
```

---

## 11. Interactive Plotly Visualization

```python
try:
    import plotly.express as px

    plot_df = pd.DataFrame({
        "UMAP1": adata.obsm["X_umap"][:, 0],
        "UMAP2": adata.obsm["X_umap"][:, 1],
        "Well":  adata.obs["Well"].values,
    })

    if "experimental_condition" in adata.obs.columns:
        plot_df["Condition"] = adata.obs["experimental_condition"].values
        color_col = "Condition"
    else:
        color_col = "Well"

    fig = px.scatter(
        plot_df,
        x="UMAP1", y="UMAP2",
        color=color_col,
        hover_data=["Well"],
        opacity=0.7,
        title="Interactive UMAP",
    )
    fig.update_layout(width=800, height=600)
    fig.show()

except ImportError:
    print("plotly not installed. Install with: pip install plotly")
```

---

## 12. Condition Comparison Side-by-Side

```python
def plot_condition_comparison(adata, condition_col, gene, ncols=2):
    """Plot gene expression per condition with a shared colorbar range."""
    conditions = adata.obs[condition_col].unique()
    nrows = int(np.ceil(len(conditions) / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 5 * nrows))
    axes = axes.flatten() if len(conditions) > 1 else [axes]

    if gene in adata.var_names:
        gene_expr = adata[:, gene].X
        if hasattr(gene_expr, "toarray"):
            gene_expr = gene_expr.toarray().flatten()
        vmin, vmax = np.percentile(gene_expr, [5, 95])
    else:
        vmin, vmax = None, None

    for ax, cond in zip(axes, conditions):
        adata_sub = adata[adata.obs[condition_col] == cond]
        if gene in adata.var_names:
            sc.pl.umap(adata_sub, color=gene, ax=ax, show=False,
                       vmin=vmin, vmax=vmax, cmap="viridis",
                       title=f"{cond}\n{gene} expression")
        else:
            sc.pl.umap(adata_sub, ax=ax, show=False, title=cond)

    for ax in axes[len(conditions):]:
        ax.axis("off")

    plt.tight_layout()
    return fig


if "experimental_condition" in adata.obs.columns and "MYC" in adata.var_names:
    fig = plot_condition_comparison(adata, "experimental_condition", "MYC")
    plt.savefig("figures/MYC_by_condition.pdf")
    plt.show()
```

---

## 13. Parameters Reference and Expected Outputs

See [REFERENCE.md](./REFERENCE.md) — `sc.pl.umap` key options, recommended colormaps, and expected figure outputs.

---


## 14. Common Pitfalls

1. **Too many categories** — if > 20 categories, use `legend_loc="none"` or aggregate groups before plotting.

2. **Overplotting** — for datasets > 100k cells, reduce `size` and `alpha`, or switch to `sc.pl.embedding_density`.

3. **Inconsistent color ranges** — always set matching `vmin`/`vmax` when comparing the same continuous variable across conditions.

4. **NaN values in color column** — scanpy will silently drop cells with NaN; check `adata.obs[col].isna().sum()` before plotting.

5. **`use_raw` mismatch** — if `adata.raw` exists, `sc.pl.umap` defaults to `use_raw=True`; pass `use_raw=False` to plot from the normalized layer instead.

6. **`AttributeError: module 'matplotlib.cm' has no attribute 'get_cmap'`** — `matplotlib.cm.get_cmap` and `plt.get_cmap` were both removed in matplotlib 3.9. Pass colormaps as strings everywhere possible (`cmap="viridis"` works in `sc.pl.umap`, `plt.scatter`, `plt.imshow`). If you genuinely need a colormap object: `matplotlib.colormaps["viridis"]`. For a discrete N-color palette: `matplotlib.colormaps["tab10"].resampled(N)`.

7. **Detaching the embedding** — `X_umap` lives in `obsm` and survives slicing, so subset the object (`adata[mask].copy()`) and keep calling `sc.pl.umap`. Pulling coordinates into a bare array beside a separate `obs` frame forces hand-rolled matplotlib for no gain.

8. **`color="well_label"`** — that column is only on `run.wells`. On AnnData use `WellLabel` (sample-sheet name) or `condition`. `Well` is the plate position (`A1`), not the sample name.
