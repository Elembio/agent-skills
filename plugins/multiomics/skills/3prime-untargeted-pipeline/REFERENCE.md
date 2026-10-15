# 3′ Untargeted Pipeline — Reference

Quick-lookup tables for HVG selection and cell type reference databases.

## Table of Contents

1. [HVG Selection Method Comparison](#hvg-selection-method-comparison)
2. [Reference Databases for Cell Type Annotation](#reference-databases-for-cell-type-annotation)

---

## HVG Selection Method Comparison

Used with `sc.pp.highly_variable_genes(adata, flavor=...)`.

| Method | Input Data | When to Use | Notes |
|--------|-----------|-------------|-------|
| `seurat_v3` | Raw counts (integers) | Default for 3′ untargeted | Most stringent; models count variance directly; requires `layer="raw"` or raw `.X` |
| `seurat` | Log-normalized counts | When `.X` already log-normalized | Classic dispersion-based method; less sensitive to library size differences |
| `cell_ranger` | Log-normalized counts | Rarely needed | Selects top N by normalized dispersion; typically not relevant for DISS data |

**Recommended default:** `seurat_v3` on raw counts, with `batch_key="Well"` to account for per-well library size variation.

```python
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    flavor="seurat_v3",
    batch_key="Well",
    subset=False,        # keeps all genes; use adata.var["highly_variable"] to subset later
)
```

**Choosing `n_top_genes`:** 2000 is a common default. For low-complexity samples (fewer total genes detected), 1000–1500 may suffice. For very deep datasets, up to 3000 is reasonable.

---

## Reference Databases for Cell Type Annotation

Use after marker gene identification (Step 12 in [IMPLEMENTATION.md](./IMPLEMENTATION.md#12-marker-gene-identification)) to map cluster markers to known cell types.

| Database | URL | Coverage | Best For |
|----------|-----|----------|----------|
| **CellMarker 2.0** | [cellmarker.com](http://cellmarker.com) | Human + mouse, 100+ cell types | Curated cell-type marker genes; downloadable CSV |
| **PanglaoDB** | [panglaodb.se](https://panglaodb.se) | Human + mouse, broad | Quick lookup; includes sensitivity/specificity scores per marker |
| **Human Cell Atlas (HCA)** | [humancellatlas.org](https://www.humancellatlas.org) | Human, single-cell resolution | Gold-standard human single-cell reference; use `scArches` or `scANVI` for mapping |
| **Tabula Muris** | [tabula-muris.ds.czbiohub.org](https://tabula-muris.ds.czbiohub.org) | Mouse, broad multi-organ atlas | Best mouse reference; FACS + droplet datasets |
| **Custom reference** | — | Experiment-matched | Highest accuracy when matched to the same biological context and protocol; use `scVI` transfer or `ingest` |

**Programmatic lookup example (PanglaoDB):**

```python
import pandas as pd

# Download from: https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz
panglao = pd.read_csv("PanglaoDB_markers_27_Mar_2020.tsv.gz", sep="\t")
human = panglao[panglao["species"].str.contains("Hs")]

def lookup_cell_type(markers, db=human, top_n=5):
    """Return cell types enriched for the provided marker gene list."""
    hits = db[db["official gene symbol"].isin(markers)]
    return hits.groupby("cell type")["official gene symbol"].count().sort_values(ascending=False).head(top_n)

# Example: look up top markers for cluster 0
top_markers = markers[markers["group"] == "0"].head(20)["names"].tolist()
print(lookup_cell_type(top_markers))
```

**Automated annotation tools:**

| Tool | Approach | Notes |
|------|----------|-------|
| `celltypist` | Logistic regression on log-normalized data | Fast; pre-trained models for 45+ cell types |
| `scVI` / `scANVI` | Deep generative model transfer | Best accuracy; requires reference AnnData |
| `SingleR` (R) | Correlation-based | Widely used; wrap via `rpy2` if needed |
