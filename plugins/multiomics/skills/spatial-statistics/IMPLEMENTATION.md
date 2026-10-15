# Spatial Statistics — Implementation

## Table of Contents

1. [Prepare Spatial Coordinates](#1-prepare-spatial-coordinates)
2. [Build Spatial Neighborhood Graph](#2-build-spatial-neighborhood-graph)
3. [Neighborhood Enrichment Analysis](#3-neighborhood-enrichment-analysis)
4. [Interpret Neighborhood Enrichment](#4-interpret-neighborhood-enrichment)
5. [Spatial Autocorrelation — Moran's I](#5-spatial-autocorrelation--morans-i)
6. [Geary's C — Alternative Autocorrelation](#6-gearys-c--alternative-autocorrelation)
7. [Co-occurrence Analysis](#7-co-occurrence-analysis)
8. [Ripley's Statistics](#8-ripleys-statistics)
9. [Centrality Scores](#9-centrality-scores)
10. [Visualize Spatial Clusters](#10-visualize-spatial-clusters)
11. [Per-Well Spatial Analysis](#11-per-well-spatial-analysis)
12. [Export Spatial Statistics Summary](#12-export-spatial-statistics-summary)
13. [Parameters](#13-parameters)
14. [Expected Outputs](#14-expected-outputs)
15. [Common Pitfalls](#15-common-pitfalls)

---

## 1. Prepare Spatial Coordinates

The store writes `obsm["spatial"]` as a DataFrame (`X`, `Y` in the region
pixel frame). Squidpy and `[:, 0]` indexing need an `(n, 2)` ndarray.
Discover keys; do not invent `obs["X"]` / `obs["Y"]`.

```python
import squidpy as sq
import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from elembio_spatialdata_analysis.load import Loader

loader = Loader("/path/to/run.zarr")  # zip is fine — tables only
print(loader.available_tables)
adata = loader.load_tables(table="ThreePrimeUntargeted", lazy=False)
print(list(adata.obsm.keys()))


def spatial_xy(adata):
    if "spatial" in adata.obsm:
        sp = adata.obsm["spatial"]
        if hasattr(sp, "columns"):
            cols = [c for c in ("X", "Y") if c in sp.columns]
            if len(cols) < 2:
                cols = list(sp.columns[:2])
            return sp[cols].to_numpy(dtype=float)
        return np.asarray(sp, dtype=float)
    cg = adata.obsm.get("cell_geometry")
    if cg is not None and {"Xum_tile", "Yum_tile"} <= set(cg.columns):
        return cg[["Xum_tile", "Yum_tile"]].to_numpy(dtype=float)
    if {"X_tile", "Y_tile"} <= set(adata.obs.columns):
        return adata.obs[["X_tile", "Y_tile"]].to_numpy(dtype=float)
    raise KeyError("No spatial coordinates — print(list(adata.obsm.keys())) / obs.columns")


adata.obsm["spatial"] = spatial_xy(adata)

# Wells are separate physical spaces — never build one graph for the plate
adata = adata[adata.obs["Well"].astype(str) == "A1"].copy()

xy = adata.obsm["spatial"]
print(f"Spatial extent: X=[{xy[:,0].min():.0f}, {xy[:,0].max():.0f}]")
print(f"                Y=[{xy[:,1].min():.0f}, {xy[:,1].max():.0f}]")
```

If `read_zarr` is already open, use
`sdata.tables["global-table-ThreePrimeUntargeted"]` (or `global-table-imaging`)
instead of `Loader`. The same `spatial_xy` conversion applies.

`obsm["spatial"]` units are pixels. For tile-local microns, pass
`cell_geometry[["Xum_tile", "Yum_tile"]].to_numpy()` as `obsm["spatial"]`
(still per well).

---

## 2. Build Spatial Neighborhood Graph

```python
# K-nearest neighbors graph (recommended for variable cell density)
sq.gr.spatial_neighbors(
    adata,
    coord_type='generic',  # generic = arbitrary 2-D coordinates, not a regular grid
    n_neighs=10,
    spatial_key='spatial'
)

# Check the graph
n_edges = adata.obsp['spatial_connectivities'].nnz
avg_neighbors = n_edges / adata.n_obs
print(f"Spatial graph: {n_edges} edges, {avg_neighbors:.1f} avg neighbors/cell")

# Radius is in the same units as obsm["spatial"] (pixels if you used the build table)
sq.gr.spatial_neighbors(
    adata,
    coord_type='generic',
    radius=100,
    spatial_key='spatial'
)
```

---

## 3. Neighborhood Enrichment Analysis

```python
# Run clustering upstream if needed
if 'leiden' not in adata.obs.columns:
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.leiden(adata, resolution=0.5, flavor="igraph", n_iterations=2, directed=False)

# Compute neighborhood enrichment
sq.gr.nhood_enrichment(adata, cluster_key='leiden')

# Results stored in adata.uns['leiden_nhood_enrichment']
#   'zscore': z-score matrix (clusters × clusters)
#   'count':  observed count matrix

# Visualize
sq.pl.nhood_enrichment(adata, cluster_key='leiden', mode='zscore', figsize=(8, 6))
plt.savefig('nhood_enrichment.png', dpi=150, bbox_inches='tight')
plt.close()
```

---

## 4. Interpret Neighborhood Enrichment

```python
# Extract enrichment matrix as a labelled DataFrame
enrichment = adata.uns['leiden_nhood_enrichment']['zscore']
clusters = adata.obs['leiden'].cat.categories

enrichment_df = pd.DataFrame(enrichment, index=clusters, columns=clusters)

fig, ax = plt.subplots(figsize=(10, 8))
sns.heatmap(
    enrichment_df,
    cmap='RdBu_r',
    center=0,
    annot=True,
    fmt='.1f',
    ax=ax
)
ax.set_title('Neighborhood Enrichment (Z-scores)')
ax.set_xlabel('Neighbor Cluster')
ax.set_ylabel('Cluster')
plt.tight_layout()
plt.savefig('nhood_enrichment_annotated.png', dpi=150, bbox_inches='tight')
plt.close()

# Interpretation:
#   Positive values → clusters co-localize more than expected by chance
#   Negative values → clusters avoid each other
#   Near zero       → random spatial distribution
```

---

## 5. Spatial Autocorrelation — Moran's I

```python
# Moran's I identifies genes whose expression is spatially coherent
sq.gr.spatial_autocorr(
    adata,
    mode='moran',
    genes=adata.var_names[:100],  # or a specific gene list
    n_perms=100,
    n_jobs=4
)

# Results in adata.uns['moranI'] — columns: I, pval_norm, pval_z_sim, ...
moran_df = adata.uns['moranI']
print(moran_df.sort_values('I', ascending=False).head(20))

# Plot top spatially variable genes
top_spatial_genes = moran_df.sort_values('I', ascending=False).head(6).index.tolist()

fig, axes = plt.subplots(2, 3, figsize=(15, 10))
for ax, gene in zip(axes.flat, top_spatial_genes):
    sc.pl.spatial(adata, color=gene, spot_size=10, ax=ax, show=False)
    ax.set_title(f'{gene}\nMoran I = {moran_df.loc[gene, "I"]:.3f}')

plt.tight_layout()
plt.savefig('spatially_variable_genes.png', dpi=150, bbox_inches='tight')
plt.close()
```

---

## 6. Geary's C — Alternative Autocorrelation

```python
# Geary's C is more sensitive to local variation than Moran's I
sq.gr.spatial_autocorr(
    adata,
    mode='geary',
    genes=['GAPDH', 'ACTB', 'TFRC'],
    n_perms=100
)

# Results in adata.uns['gearyC']
geary_df = adata.uns['gearyC']
print(geary_df)

# Interpretation:
#   C < 1 → positive autocorrelation (clustering)
#   C = 1 → no autocorrelation (random)
#   C > 1 → negative autocorrelation (dispersion)
```

---

## 7. Co-occurrence Analysis

```python
# Co-occurrence probability as a function of distance between cell types
sq.gr.co_occurrence(
    adata,
    cluster_key='leiden',
    spatial_key='spatial',
    n_splits=10  # number of distance bins
)

# Results in adata.uns['leiden_co_occurrence']
#   'occ':      co-occurrence ratios at each distance
#   'interval': distance interval boundaries

sq.pl.co_occurrence(
    adata,
    cluster_key='leiden',
    clusters=['0', '1', '2']  # specific clusters to display
)
plt.savefig('co_occurrence.png', dpi=150, bbox_inches='tight')
plt.close()
```

---

## 8. Ripley's Statistics

```python
# Ripley's L: normalised version of Ripley's K for point pattern analysis
# Tests whether a cell type is clustered, randomly distributed, or dispersed
sq.gr.ripley(
    adata,
    cluster_key='leiden',
    mode='L',
    spatial_key='spatial',
    max_dist=500  # maximum distance to analyse (µm)
)

# Results in adata.uns['leiden_ripley_L']

sq.pl.ripley(adata, cluster_key='leiden', mode='L')
plt.savefig('ripley_L.png', dpi=150, bbox_inches='tight')
plt.close()

# Interpretation:
#   L(r) > r → clustering at distance r
#   L(r) = r → complete spatial randomness (CSR)
#   L(r) < r → dispersion / inhibition at distance r
```

---

## 9. Centrality Scores

```python
# Graph-theoretic centrality metrics per cluster
sq.gr.centrality_scores(adata, cluster_key='leiden')

# Results in adata.uns['leiden_centrality_scores']
centrality = adata.uns['leiden_centrality_scores']
print(centrality)

# Metrics:
#   degree_centrality     — normalised number of graph neighbours
#   closeness_centrality  — inverse average path length to all other nodes
#   clustering_coefficient — local connectivity / triangles in the graph
```

---

## 10. Visualize Spatial Clusters

```python
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# Left: continuous colour-coded scatter
ax = axes[0]
scatter = ax.scatter(
    adata.obsm['spatial'][:, 0],
    adata.obsm['spatial'][:, 1],
    c=adata.obs['leiden'].astype(int),
    s=2,
    alpha=0.7,
    cmap='tab20'
)
ax.set_xlabel('X (px)')
ax.set_ylabel('Y (px)')
ax.set_title('Spatial Distribution by Cluster')
ax.set_aspect('equal')
ax.invert_yaxis()
plt.colorbar(scatter, ax=ax, label='Cluster')

# Right: per-cluster overlay
ax = axes[1]
for cluster in adata.obs['leiden'].cat.categories:
    mask = adata.obs['leiden'] == cluster
    ax.scatter(
        adata.obsm['spatial'][mask, 0],
        adata.obsm['spatial'][mask, 1],
        s=1, alpha=0.5, label=cluster
    )
ax.set_xlabel('X (px)')
ax.set_ylabel('Y (px)')
ax.set_title('Clusters Overlay')
ax.set_aspect('equal')
ax.invert_yaxis()
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', markerscale=5)

plt.tight_layout()
plt.savefig('spatial_clusters.png', dpi=150, bbox_inches='tight')
plt.close()
```

---

## 11. Per-Well Spatial Analysis

Section 1 already subsets to one well. To compare wells, start from the
full table (after `spatial_xy`), then loop. Do not build neighbors on the
unsubsetted plate.

```python
# `full` = table after spatial_xy, all wells still present
wells = full.obs['Well'].astype(str).unique()
enrichment_by_well = {}

for well in wells:
    well_adata = full[full.obs['Well'].astype(str) == well].copy()

    # Rebuild spatial graph for this well's coordinate extent
    sq.gr.spatial_neighbors(well_adata, coord_type='generic', n_neighs=10)
    sq.gr.nhood_enrichment(well_adata, cluster_key='leiden')

    enrichment_by_well[well] = well_adata.uns['leiden_nhood_enrichment']['zscore']

# Compare across wells
fig, axes = plt.subplots(1, len(wells), figsize=(5 * len(wells), 4))
if len(wells) == 1:
    axes = [axes]

for ax, well in zip(axes, wells):
    sns.heatmap(enrichment_by_well[well], cmap='RdBu_r', center=0, ax=ax, cbar=False)
    ax.set_title(f'Well {well}')

plt.tight_layout()
plt.savefig('nhood_enrichment_by_well.png', dpi=150, bbox_inches='tight')
plt.close()
```

---

## 12. Export Spatial Statistics Summary

```python
import json

spatial_summary = {
    'n_cells': int(adata.n_obs),
    'n_clusters': int(adata.obs['leiden'].nunique()),
    'spatial_extent': {
        'x_range': [float(adata.obsm['spatial'][:, 0].min()),
                    float(adata.obsm['spatial'][:, 0].max())],
        'y_range': [float(adata.obsm['spatial'][:, 1].min()),
                    float(adata.obsm['spatial'][:, 1].max())]
    },
    'avg_neighbors': float(adata.obsp['spatial_connectivities'].sum() / adata.n_obs),
    'top_spatially_variable_genes': (
        moran_df.sort_values('I', ascending=False).head(10).to_dict('index')
    )
}

with open('spatial_summary.json', 'w') as f:
    json.dump(spatial_summary, f, indent=2)

print("Spatial analysis summary saved to spatial_summary.json")
```

---

## 13. Parameters

### Neighborhood graph construction

```python
# K-nearest neighbors — recommended when cell density varies across the sample
sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=10)

# Radius-based — radius is in obsm["spatial"] units (pixels for the build table)
sq.gr.spatial_neighbors(adata, coord_type='generic', radius=50)

# Delaunay triangulation — captures natural spatial topology
sq.gr.spatial_neighbors(adata, coord_type='generic', delaunay=True)
```

### Large-dataset subsampling

```python
# For datasets with >50 k cells, subsample before computationally heavy analyses
if adata.n_obs > 50_000:
    np.random.seed(42)
    idx = np.random.choice(adata.n_obs, size=50_000, replace=False)
    adata_subset = adata[idx].copy()
    sq.gr.spatial_autocorr(adata_subset, mode='moran', genes=genes)
```

---

## 14. Expected Outputs

After running this workflow the AnnData object will contain:

| Location | Content |
|----------|---------|
| `.obsp['spatial_connectivities']` | Sparse spatial neighborhood graph |
| `.uns['{cluster_key}_nhood_enrichment']` | Z-scores and counts for neighborhood enrichment |
| `.uns['moranI']` | Moran's I scores and p-values per gene |
| `.uns['gearyC']` | Geary's C scores per gene |
| `.uns['{cluster_key}_co_occurrence']` | Co-occurrence ratios across distance bins |
| `.uns['{cluster_key}_ripley_L']` | Ripley's L function values |
| `.uns['{cluster_key}_centrality_scores']` | Degree, closeness, and clustering coefficient per cluster |

Saved files: `nhood_enrichment.png`, `nhood_enrichment_annotated.png`, `spatially_variable_genes.png`, `co_occurrence.png`, `ripley_L.png`, `spatial_clusters.png`, `nhood_enrichment_by_well.png`, `spatial_summary.json`.

---

## 15. Common Pitfalls

1. **Coordinates and units** — Use `obsm["spatial"]` (`X`, `Y` pixels). Convert a DataFrame to an ndarray. A radius is in those units. Tile-local microns are `cell_geometry` `Xum_tile` / `Yum_tile`.
2. **Mixing wells** — Do not build one graph for the plate.
3. **DataFrame in `obsm["spatial"]`** — squidpy and `xy[:, 0]` fail until you assign `spatial_xy(adata)`.
4. **Edge effects** — Cells at sample boundaries have fewer neighbours. Interpret statistics for boundary clusters with caution.
5. **Tile-local vs global** — `obs["X_tile"]` / `obs["Y_tile"]` are within-tile. Neighborhoods use `obsm["spatial"]`.
6. **Multiple testing** — When scanning all genes for spatial autocorrelation, apply FDR correction (`pval_norm_fdr_bh` column in `moranI` results).
7. **Compute time** — Spatial autocorrelation scales with `n_genes × n_perms`. Parallelize with `n_jobs` or test on a gene subset before running the full panel.
