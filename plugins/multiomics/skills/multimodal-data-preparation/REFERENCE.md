# Multimodal Data Preparation — Reference

## Alignment Strategy Comparison

| Strategy | Behaviour | When to use |
|---|---|---|
| **Intersection** | Keep only cells measured in **all** modalities | Strictest; no missing values; best for joint embedding |
| **Union** | Keep all cells; absent modalities become **NaN / zeros** | Permissive; retains rare cells; requires modality-aware methods downstream |

Default: **intersection**. Union-padded zeros look like genuine zero expression — use the `_has_data` flag to mask those cells before statistics.

## Cell ID Alignment Checklist

1. Cell IDs are `adata.obs_names` (the AnnData index). If a `Cell` column exists in `.obs`, it should match.
2. Cast IDs to `str` before set operations — IDs may be integers in one modality and strings in another.
3. After alignment, verify `n_obs` is identical across all modalities and `obs_names` arrays match exactly.

## Modality Transform Order

Always transform **after** alignment, never before:

| Step | Action |
|---|---|
| 1 | Load present modalities with `lazy=False` |
| 2 | QC + filter each modality independently |
| 3 | Align cells (intersection or union) |
| 4 | **Counts:** `ebsd.pp.normalize` — [normalization-strategies](../normalization-strategies/SKILL.md) |
| 4b | **Imaging:** MAD-z feature matrix — [imaging-pipeline](../imaging-pipeline/SKILL.md) (not `pp.normalize`; not global StandardScaler) |
| 5 | Create MuData or concatenate for joint analysis |

Imaging / CellPaint: `ebsd.pp.normalize(..., "Imaging")` raises `NotImplementedError` by design.
## Expected Outputs

| Object | Type | Description |
|---|---|---|
| `aligned` | `dict[str, AnnData]` | Modality AnnDatas with identical cell order and `_cell_id` column |
| `normalized` | `dict[str, AnnData]` | Each modality normalized by its own pipeline |
| `combined` | `AnnData` | Wide matrix with all features (optional; useful for joint PCA) |
| `mdata` | `MuData` | Recommended container; write to `.h5mu` for persistence |

Validation checks:
- All modalities report the same `n_obs`
- `obs_names` arrays are identical across modalities
- Metadata on `.obs` (`Well`, tile, `X_tile`, `Y_tile`) is present in every modality
- Region centroids stay in `obsm["spatial"]` (`X`, `Y`) — already on each table; not an `.obs` transfer

## Common Pitfalls

1. **Cell ID type mismatch** — Cell IDs are sometimes integers in one modality and strings in another. Cast with `.astype(str)` before building sets.

2. **Index vs column alignment** — pandas aligns on index by default. When building a new AnnData, set `obs.index` explicitly rather than relying on positional order.

3. **Union + zeros vs NaN** — The union helper zero-pads absent cells, which looks like genuine zero expression. Use `_has_data` to mask those cells before statistics.

4. **Normalize after alignment, not before** — Normalizing separate populations then aligning can produce different scale factors per modality.

5. **Memory for dense concatenation** — `concatenate_modalities` densifies sparse matrices. For large experiments (> 500k cells × 5k genes) consider keeping modalities sparse and using MuData instead.
