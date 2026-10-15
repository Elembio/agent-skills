# Embedding Visualization — Reference

## Contents

1. [sc.pl.umap key options](#1-scplumap-key-options)
2. [Recommended colormaps](#2-recommended-colormaps)
3. [Expected outputs](#3-expected-outputs)

---

## 1. sc.pl.umap key options

```python
sc.pl.umap(
    adata,
    color="variable",           # categorical or continuous obs column, or gene name
    palette="tab20",            # color palette for categorical variables
    cmap="viridis",             # colormap for continuous variables
    vmin=0, vmax=10,            # clamp range for continuous colormap
    size=10,                    # point size
    alpha=0.7,                  # point transparency
    frameon=False,              # remove axes frame
    legend_loc="right margin",  # 'on data', 'right margin', None
    legend_fontsize=8,
    title="Custom Title",
    save="_filename.pdf",       # saves under figures/
    show=True,
)
```

---

## 2. Recommended colormaps

| Data Type | Colormap | Notes |
|---|---|---|
| Gene expression | `viridis`, `magma` | Perceptually uniform |
| Diverging | `RdBu_r`, `coolwarm` | For +/− centered values |
| Categorical | `tab20`, `Set2` | Visually distinct groups |
| Quality metric | `RdYlGn` | Red = bad, green = good |

---

## 3. Expected outputs

| Output | Description |
|---|---|
| `figures/umap_overview.pdf` | 2×3 multi-panel overview across all variables |
| `figures/umap_by_condition.pdf` | Faceted panels, one per condition |
| `figures/MYC_by_condition.pdf` | Gene expression split by condition |
| `figures/figure1_umap.pdf/.png` | Publication-ready 3-panel figure |
