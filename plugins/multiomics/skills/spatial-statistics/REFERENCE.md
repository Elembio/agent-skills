# Spatial Statistics — Reference

## Neighborhood Graph Parameters

```python
# K-nearest neighbors — recommended when cell density varies across the sample
sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=10)

# Radius-based — radius is in obsm["spatial"] units (pixels for the build table)
sq.gr.spatial_neighbors(adata, coord_type="generic", radius=50)

# Delaunay triangulation — captures natural spatial topology
sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)
```

## Key squidpy Function Parameters

| Function | Key parameters | Notes |
|---|---|---|
| `sq.gr.spatial_neighbors` | `coord_type`, `n_neighs`, `radius`, `delaunay` | First step; required for all analyses |
| `sq.gr.spatial_autocorr` | `mode` (`"moran"` / `"geary"`), `genes`, `n_perms`, `n_jobs` | Scale `n_jobs` to CPU count |
| `sq.gr.nhood_enrichment` | `cluster_key` | Requires leiden cluster labels |
| `sq.gr.co_occurrence` | `cluster_key`, `interval` | Computationally expensive for large N |
| `sq.gr.ripley` | `cluster_key`, `mode` (`"L"`) | Per-cluster point pattern analysis |
| `sq.gr.centrality_scores` | `cluster_key` | Graph-theoretic hub detection |

## Large-Dataset Subsampling

For datasets with > 50k cells, subsample before computationally heavy analyses:

```python
if adata.n_obs > 50_000:
    import numpy as np
    np.random.seed(42)
    idx = np.random.choice(adata.n_obs, size=50_000, replace=False)
    adata_subset = adata[idx].copy()
    sq.gr.spatial_autocorr(adata_subset, mode="moran", genes=adata.var_names)
```

## Expected Output Slots

| Location | Content |
|---|---|
| `.obsp["spatial_connectivities"]` | Sparse spatial neighborhood graph |
| `.uns["{cluster_key}_nhood_enrichment"]` | Z-scores and counts for neighborhood enrichment |
| `.uns["moranI"]` | Moran's I scores and p-values per gene |
| `.uns["gearyC"]` | Geary's C scores per gene |
| `.uns["{cluster_key}_co_occurrence"]` | Co-occurrence ratios across distance bins |
| `.uns["{cluster_key}_ripley_L"]` | Ripley's L function values |
| `.uns["{cluster_key}_centrality_scores"]` | Degree, closeness, and clustering coefficient per cluster |

## Common Pitfalls

1. **Coordinates and units** — `obsm["spatial"]` is region-pixel `X`/`Y` (often a DataFrame). Convert to an `(n, 2)` array before squidpy. Radius is in those units. Tile-local microns: `obsm["cell_geometry"]` `Xum_tile`/`Yum_tile`.
2. **Mixing wells** — one graph on the whole plate is invalid; wells do not share a coordinate space for neighbors.
3. **Edge effects** — Cells at sample boundaries have fewer neighbours. Interpret statistics for boundary clusters with caution.
4. **Tile-local vs global** — `obs["X_tile"]` / `obs["Y_tile"]` are within-tile. Neighborhoods use `obsm["spatial"]` (region pixel frame).
5. **Multiple testing** — When scanning all genes for spatial autocorrelation, apply FDR correction (`pval_norm_fdr_bh` column in `moranI` results).
6. **Compute time** — Spatial autocorrelation scales with `n_genes × n_perms`. Parallelize with `n_jobs` or test on a gene subset before running the full panel.
