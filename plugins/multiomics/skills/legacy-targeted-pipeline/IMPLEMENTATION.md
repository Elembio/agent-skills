# Legacy Targeted Pipeline — Implementation

Legacy Transcript / Protein barcoding path only. For imaging / CellPaint use [imaging-pipeline](../imaging-pipeline/SKILL.md). For 3′ / OPS use those pipelines.

Full end-to-end code for the AVITI24 standard analysis pipeline (Steps 0–13).
Copy blocks into a notebook or script, or run sequentially as a standalone pipeline.

---

## Table of Contents

- [Step 0: Environment Setup](#step-0-environment-setup)
- [Step 1: Open Store with Loader](#step-1-open-store-with-loader)
- [Step 2: Parse RunInfo Metadata](#step-2-parse-runinfo-metadata)
- [Step 3: Load Transcript / Protein Tables](#step-3-load-transcript--protein-tables)
- [Step 4: Add Experimental Metadata](#step-4-add-experimental-metadata)
- [Step 5: Quality Control Filtering](#step-5-quality-control-filtering)
- [Step 6: Normalize Data](#step-6-normalize-data)
- [Step 7: Assess Batch Effects](#step-7-assess-batch-effects)
- [Step 8: Dimensionality Reduction](#step-8-dimensionality-reduction)
- [Step 9: Embedding Visualization](#step-9-embedding-visualization)
- [Step 10: Clustering](#step-10-clustering)
- [Step 11: Differential Expression Analysis](#step-11-differential-expression-analysis)
- [Step 12: Save Results and Final Report](#step-12-save-results-and-final-report)
- *(Optional)* [Cluster Validation](../cluster-validation/SKILL.md) — invoke on demand only

---

## Memory budget for large datasets (50K+ cells)

In a long-running session memory **accumulates across these steps**, so a large run can exhaust available memory even when no single step is large. **When memory is the binding constraint**, run the pipeline memory-lean:

- **`float32` for large count / expression matrices** — cast `adata.X` after load (Step 3) and re-assert after normalize/scale; it halves every dense matrix. Not for CellPaint morphology / precision-sensitive stats.
- **Keep one working object** — at QC filtering (Step 5) the `inplace=True` default avoids a duplicate; keep the unfiltered object only when you will compare or re-filter at new thresholds.
- **Free intermediates you will not revisit** — `gc.collect()` between heavy steps. The scaled `X` is the largest such intermediate after PCA; releasing it is opt-in (see the Dimensionality Reduction note) — do not overwrite `.X` mid-session in interactive work.
- **Leiden (Step 10) is the spike** — on 100K+ cells use `flavor="igraph", n_iterations=2`, free memory first, and if your runtime caps execution time, raise it or split the work.

Any memory-forced deletion that could lose work the user wants (the unfiltered object, a layer, precision) should be surfaced to the user as a choice, not done silently; default to keeping it when you cannot ask.

---

## Step 0: Environment Setup

```python
import elembio_spatialdata_analysis as ebsd
from elembio_spatialdata_analysis.load import Loader
import spatialdata as sd
import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

os.makedirs("figures", exist_ok=True)
os.makedirs("results", exist_ok=True)
os.makedirs("results/de", exist_ok=True)

plt.rcParams["figure.dpi"] = 100
sc.settings.verbosity = 2

print("Setup complete")
```

---

## Step 1: Open Store with Loader

> Skill: [spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md)

```python
STORE_PATH = "/path/to/run.zarr"
loader = Loader(STORE_PATH)

print(loader.summary())
print("Tables:", sorted(loader.available_tables))
print("Target types (points):", sorted(loader.available_target_types))
```

**Alternative — full `SpatialData` (images / labels / spatial elements):**

```python
# Use only when the workflow needs sdata.images, sdata.labels, sdata.points as linked elements.
sdata = sd.read_zarr(STORE_PATH)
print(f"Images: {len(sdata.images)}, Labels: {len(sdata.labels)}")
loader = Loader(sdata)  # wraps in-memory store; same load_tables API
```

---

## Step 2: Parse RunInfo Metadata

> Skill: [metadata-and-target-mapping](../metadata-and-target-mapping/SKILL.md)

`loader._runs` is a `dict[str, RunInfo]`. For a single store:

```python
run_key, run = next(iter(loader._runs.items()))

wells_df = run.wells
tiles_df = run.tiles
batches_df = run.batches
targets_barcoding = run.targets_barcoding
targets_imaging = run.targets_imaging

print(f"Wells: {len(wells_df)}")
print(f"Tiles: {len(tiles_df)}")
print(f"Batches: {len(batches_df)}")
print(f"Targets (barcoding): {len(targets_barcoding)}")
print(f"Targets (imaging): {len(targets_imaging)}")

if "batch_class" in batches_df.columns:
    barcoding_batches = batches_df[batches_df["batch_class"] == "BarcodingBatch"]
    print(f"\nBarcoding batches: {barcoding_batches['batch_name'].tolist()}")
```

For multi-run stores, iterate with `for run_key, run in loader._runs.items()`.

---

## Step 3: Load Transcript / Protein Tables

> Skill: [spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md)

Use `lazy=False` for this pipeline — the default `lazy=True` (zarr-backed `X`) conflicts
with `sc.pp.pca`, neighbors, UMAP, and most `pp` steps.

```python
# lazy=False: required for pp / Scanpy steps (lazy=True is default, fine for exploration)
adata = loader.load_tables(table="Transcript", lazy=False)
print(f"Transcript data: {adata.shape}")

if "Protein" in loader.available_tables:
    protein_adata = loader.load_tables(table="Protein", lazy=False)
    print(f"Protein data: {protein_adata.shape}")
else:
    protein_adata = None
    print("No protein table in store")

print(f"\n.obs columns:  {adata.obs.columns.tolist()}")
print(f".var columns:  {adata.var.columns.tolist()}")
print(f".obsm keys:    {list(adata.obsm.keys())}")
print(f".layers keys:  {list(adata.layers.keys())}")
```

---

## Step 4: Add Experimental Metadata

> Skill: [metadata-and-target-mapping](../metadata-and-target-mapping/SKILL.md)

Map well labels to experimental conditions. Customize the parsing logic for your naming convention.

```python
def parse_well_conditions(wells_df):
    """Parse experimental conditions from sample-sheet labels."""
    label_col = next((c for c in ("WellLabel", "well_label") if c in wells_df.columns), None)
    conditions = {}
    for well_name, row in wells_df.iterrows():
        label = str(row[label_col]) if label_col else str(well_name)
        conditions[well_name] = {
            "experimental_condition": label.split("_")[0] if "_" in label else label,
            "condition_level": label.split("_")[1] if len(label.split("_")) > 1 else None,
        }
    return pd.DataFrame(conditions).T

well_conditions = parse_well_conditions(wells_df)

adata.obs["experimental_condition"] = adata.obs["Well"].map(
    well_conditions["experimental_condition"]
)
adata.obs["condition_level"] = adata.obs["Well"].map(
    well_conditions["condition_level"]
)

print(f"Experimental conditions: {adata.obs['experimental_condition'].unique()}")
print(f"Condition levels:        {adata.obs['condition_level'].unique()}")
```

---

## Step 5: Quality Control Filtering

> Skill: [cell-quality-control](../cell-quality-control/SKILL.md)

Run `ebsd.pp.compute_qc_metrics` first — required before filter plots and reused by `ebsd.pp.filter_cells`
to avoid recomputation.

```python
ebsd.pp.compute_qc_metrics(adata, "Transcript")
ebsd.pp.filter_cells(adata, target_type="Transcript", inplace=True)

print(f"Cells after filtering: {adata.n_obs}")

# Per-filter breakdown is always in adata.uns["filter_flags_full"] — no diagnostics=True needed
filter_summary = (~adata.uns["filter_flags_full"]).sum()
print(f"\nCells failing each filter:\n{filter_summary.sort_values(ascending=False)}")
```

> **Decision:** If rejection rate exceeds ~30 %, inspect `adata.uns["filter_flags_full"]` per filter and
> adjust thresholds — see the `cell-quality-control` skill for threshold parameters.

---

## Step 6: Normalize Data

> Skill: [normalization-strategies](../normalization-strategies/SKILL.md)

`ebsd.pp.normalize` is a dispatcher: it snapshots raw counts to `layers["raw"]` and runs the
modality-appropriate pipeline (agg-count scaling → median-of-ratios → log1p for Transcript; agg-count scaling → log1p for Protein).

```python

ebsd.pp.normalize(
    adata,
    "Transcript",
    batch_colname="batch_name",
    well_colname="Well",
    run_colname=None,  # set if multi-run (Loader adds run_name to obs for multi-run tables)
    inplace=True,
)

print(f"Layers: {list(adata.layers.keys())}")
print(f"Expression range: {adata.X.min():.2f} to {adata.X.max():.2f}")
```

For Protein (if loaded):

```python
if protein_adata is not None:
    ebsd.pp.normalize(
        protein_adata,
        "Protein",
        batch_colname="batch_name",
        well_colname="Well",
        inplace=True,
    )
    print(f"Protein layers: {list(protein_adata.layers.keys())}")
```

---

## Step 7: Assess Batch Effects

> Skill: [batch-effect-assessment](../batch-effect-assessment/SKILL.md)

```python
sc.pp.pca(adata, n_comps=50)

fig, axes = plt.subplots(1, 3, figsize=(18, 5))

sc.pl.pca(adata, color="Well", ax=axes[0], show=False, title="PCA by Well")

if "experimental_condition" in adata.obs.columns:
    sc.pl.pca(
        adata, color="experimental_condition",
        ax=axes[1], show=False, title="PCA by Condition"
    )
else:
    axes[1].axis("off")

axes[2].bar(range(1, 21), adata.uns["pca"]["variance_ratio"][:20] * 100)
axes[2].set_xlabel("PC")
axes[2].set_ylabel("Variance Explained (%)")
axes[2].set_title("Scree Plot")

plt.tight_layout()
plt.savefig("figures/pca_batch_assessment.pdf")
plt.show()
```

> **Decision:** If wells cluster strongly in PCA space, apply Harmony or scVI before computing
> neighbors — see the `batch-effect-assessment` skill.

---

## Step 8: Dimensionality Reduction

> Skill: [dimensionality-reduction](../dimensionality-reduction/SKILL.md)

```python
sc.pp.scale(adata, max_value=10)
sc.pp.pca(adata, n_comps=min(50, adata.n_vars - 1), svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata, min_dist=0.5, random_state=42)

print(f"PCA:  {adata.obsm['X_pca'].shape}")
print(f"UMAP: {adata.obsm['X_umap'].shape}")
```

> For 3′ untargeted data: run `sc.pp.highly_variable_genes` before `sc.pp.scale` — see [dimensionality-reduction](../dimensionality-reduction/SKILL.md).

---

## Step 9: Embedding Visualization

> Skill: [embedding-visualization](../embedding-visualization/SKILL.md)

```python
fig, axes = plt.subplots(2, 2, figsize=(14, 14))

sc.pl.umap(adata, color="Well", ax=axes[0, 0], show=False, title="Well")

if "experimental_condition" in adata.obs.columns:
    sc.pl.umap(
        adata, color="experimental_condition",
        ax=axes[0, 1], show=False, title="Condition"
    )
else:
    axes[0, 1].axis("off")

# Total counts from raw snapshot (ebsd.pp.normalize stores raw as layers["raw"])
if "total_counts" not in adata.obs.columns:
    raw_layer = adata.layers.get("raw", adata.layers.get("raw_counts", adata.X))
    adata.obs["total_counts"] = np.asarray(raw_layer.sum(axis=1)).ravel()

sc.pl.umap(
    adata, color="total_counts",
    ax=axes[1, 0], show=False, title="Total Counts", cmap="viridis"
)

if adata.n_vars > 0:
    gene = adata.var_names[0]
    sc.pl.umap(
        adata, color=gene,
        ax=axes[1, 1], show=False, title=f"{gene} Expression", cmap="viridis"
    )
else:
    axes[1, 1].axis("off")

plt.tight_layout()
plt.savefig("figures/umap_overview.pdf")
plt.show()
```

---

## Step 10: Clustering

> Skill: [clustering-and-annotation](../clustering-and-annotation/SKILL.md)

```python
sc.tl.leiden(adata, resolution=0.5, key_added="leiden",
             flavor="igraph", n_iterations=2, directed=False)

print(f"Clusters: {adata.obs['leiden'].nunique()}")
print(f"Cluster sizes:\n{adata.obs['leiden'].value_counts()}")

sc.tl.rank_genes_groups(
    adata, groupby="leiden", method="wilcoxon",
    use_raw=False,   # adata.raw is not set; use adata.X (log-normalized)
    pts=True,        # fraction expressing in group (pts) and rest (pts_rest)
)

fig, axes = plt.subplots(1, 2, figsize=(16, 6))

sc.pl.umap(
    adata, color="leiden",
    ax=axes[0], show=False, legend_loc="on data", title="Leiden Clusters"
)
sc.pl.rank_genes_groups(adata, n_genes=5, sharey=False, ax=axes[1], show=False)

plt.tight_layout()
plt.savefig("figures/clustering.pdf")
plt.show()
```

> **Decision:** Too many or too few clusters? Rerun with a different `resolution` (lower → fewer
> clusters, higher → more granular). See `clustering-and-annotation` skill.

---

## Step 11: Differential Expression Analysis

> Skill: [pseudobulk-differential-expression](../pseudobulk-differential-expression/SKILL.md)

Uses pseudobulk aggregation + DESeq2. Requires ≥ 2 biological replicates per condition.
Aggregate from `layers["raw"]` (integer counts stored by `ebsd.pp.normalize`).

```python
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import scipy.sparse as sp

def create_pseudobulk(adata, sample_col, condition_col, layer="raw"):
    """Aggregate cells to pseudobulk samples from integer count layer."""
    X = adata.layers[layer]
    if sp.issparse(X):
        X = X.toarray()

    pb_counts, pb_meta = [], []
    for sample in adata.obs[sample_col].unique():
        mask = adata.obs[sample_col] == sample
        if mask.sum() < 10:
            continue
        pb_counts.append(X[mask].sum(axis=0))
        pb_meta.append({
            sample_col: sample,
            condition_col: adata.obs.loc[mask, condition_col].iloc[0],
            "n_cells": int(mask.sum()),
        })

    pb = ad.AnnData(
        X=np.vstack(pb_counts).astype(int),
        obs=pd.DataFrame(pb_meta),
        var=adata.var.copy(),
    )
    pb.obs.index = pb.obs[sample_col].values
    return pb

sig = pd.DataFrame()  # default (used in Step 13 report)

if "experimental_condition" in adata.obs.columns:
    conditions = adata.obs["experimental_condition"].dropna().unique()

    if len(conditions) >= 2:
        pb = create_pseudobulk(adata, "Well", "experimental_condition")
        print(f"Pseudobulk samples: {pb.n_obs}")

        counts_df = pd.DataFrame(
            pb.X.clip(min=0).astype(int),
            index=pb.obs_names,
            columns=pb.var_names,
        )

        dds = DeseqDataSet(
            counts=counts_df,
            metadata=pb.obs.copy(),
            design_factors="experimental_condition",
        )
        dds.deseq2()

        stat_res = DeseqStats(
            dds,
            contrast=["experimental_condition", conditions[0], conditions[1]],
        )
        stat_res.summary()

        de_results = stat_res.results_df.sort_values("padj")

        sig = de_results[
            (de_results["padj"] < 0.05) &
            (np.abs(de_results["log2FoldChange"]) > 1)
        ]
        print(f"\nSignificant DEGs: {len(sig)}")
        print(f"  Upregulated:   {(sig['log2FoldChange'] > 0).sum()}")
        print(f"  Downregulated: {(sig['log2FoldChange'] < 0).sum()}")

        de_results.to_csv("results/de/differential_expression.csv")
        sig.to_csv("results/de/significant_degs.csv")
    else:
        print("Need at least 2 conditions for DE analysis")
else:
    print("No experimental_condition column — skipping DE")
```

---

## Step 12: Save Results and Final Report

```python
adata.write("results/processed_adata.h5ad")

n_conditions = adata.obs.get("experimental_condition", pd.Series()).nunique()
n_clusters = adata.obs["leiden"].nunique() if "leiden" in adata.obs.columns else "N/A"
classifier_acc = f"{scores.mean():.3f}" if "scores" in dir() else "N/A"
n_sig = len(sig) if not sig.empty else "N/A"

report = f"""# Analysis Summary

## Data
- Cells: {adata.n_obs}
- Genes: {adata.n_vars}
- Wells: {adata.obs['Well'].nunique()}
- Conditions: {n_conditions}

## Clustering
- Clusters: {n_clusters}
- Classifier accuracy: {classifier_acc}

## Differential Expression
- Significant genes (padj < 0.05, |LFC| > 1): {n_sig}

## Files Generated
- results/processed_adata.h5ad
- results/de/differential_expression.csv
- results/de/significant_degs.csv
- figures/pca_batch_assessment.pdf
- figures/umap_overview.pdf
- figures/clustering.pdf
"""

print(report)

with open("results/analysis_summary.md", "w") as f:
    f.write(report)
```
