# Multimodal Data Preparation — Implementation

Implementation details for [SKILL.md](./SKILL.md).

## Table of Contents

1. [Step 3 — Align modalities](#step-3--align-modalities)
2. [Step 4 — Union strategy and missing-data handling](#step-4--union-strategy-and-missing-data-handling)
3. [Step 5 — Transfer metadata across modalities](#step-5--transfer-metadata-across-modalities)
4. [Step 6 — Transform each modality independently](#step-6--transform-each-modality-independently)
5. [Step 7 — Concatenate into a combined feature matrix](#step-7--concatenate-into-a-combined-feature-matrix)
6. [Step 8 — Create MuData object](#step-8--create-mudata-object)
7. [Step 9 — Handle batch-specific features](#step-9--handle-batch-specific-features)
8. [Step 10 — Validate multimodal alignment](#step-10--validate-multimodal-alignment)
9. [Expected Outputs](#expected-outputs)
10. [Common Pitfalls](#common-pitfalls)

---

> **Steps 1–2** (loading modalities and identifying common cells) live in the Quick Start section of [SKILL.md](./SKILL.md).

---

## Step 3 — Align modalities

`align_modalities` reorders (and optionally pads) every AnnData in a dict so that all objects share the same cell order.

```python
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

def align_modalities(modalities_dict, cell_col="Cell", strategy="intersection"):
    """
    Align multiple AnnData objects to have the same cells.

    Parameters
    ----------
    modalities_dict : dict
        Mapping of modality name -> AnnData.
    cell_col : str
        Column in .obs containing cell identifiers.
    strategy : {"intersection", "union"}
        "intersection" keeps only cells present in every modality.
        "union" keeps all cells and zero-pads missing entries.

    Returns
    -------
    dict
        Aligned AnnData objects with matching cell order.
    """
    # Attach a canonical _cell_id column
    cell_ids = {}
    for name, adata in modalities_dict.items():
        if cell_col in adata.obs.columns:
            cell_ids[name] = set(adata.obs[cell_col].astype(str).values)
            adata.obs["_cell_id"] = adata.obs[cell_col].astype(str)
        else:
            cell_ids[name] = set(adata.obs_names.astype(str))
            adata.obs["_cell_id"] = adata.obs_names.astype(str)

    if strategy == "intersection":
        target_cells = sorted(set.intersection(*cell_ids.values()))
    elif strategy == "union":
        target_cells = sorted(set.union(*cell_ids.values()))
    else:
        raise ValueError(f"Unknown strategy: {strategy!r}")

    print(f"Target cells ({strategy}): {len(target_cells)}")

    aligned = {}
    for name, adata in modalities_dict.items():
        present_cells = adata.obs["_cell_id"].values

        if strategy == "intersection":
            mask = adata.obs["_cell_id"].isin(target_cells)
            adata_aligned = adata[mask].copy()
            order = adata_aligned.obs["_cell_id"].apply(target_cells.index)
            adata_aligned = adata_aligned[order.argsort()].copy()

        else:  # union
            X_full = np.zeros((len(target_cells), adata.n_vars))
            for i, cell in enumerate(target_cells):
                if cell in present_cells:
                    idx = np.where(present_cells == cell)[0][0]
                    row = adata.X[idx]
                    X_full[i] = row.toarray().flatten() if sp.issparse(row) else row

            obs_full = pd.DataFrame(index=target_cells)
            adata_aligned = ad.AnnData(X=X_full, obs=obs_full, var=adata.var.copy())
            adata_aligned.obs["_cell_id"] = target_cells
            adata_aligned.obs["_has_data"] = [c in present_cells for c in target_cells]

        adata_aligned.obs.drop(columns=["_order"], errors="ignore", inplace=True)
        aligned[name] = adata_aligned

    return aligned

# Align with intersection (recommended default)
aligned = align_modalities(available_modalities, cell_col="Cell", strategy="intersection")

for name, adata in aligned.items():
    print(f"{name}: {adata.shape}")

# Confirm cell order is identical across modalities
first = aligned[next(iter(aligned))].obs["_cell_id"].values
all_match = all(
    (adata.obs["_cell_id"].values == first).all() for adata in aligned.values()
)
print(f"Cells aligned: {all_match}")
```

---

## Step 4 — Union strategy and missing-data handling

```python
aligned_union = align_modalities(available_modalities, cell_col="Cell", strategy="union")

for name, adata in aligned_union.items():
    if "_has_data" in adata.obs.columns:
        print(f"{name}: {adata.obs['_has_data'].mean():.1%} coverage")

# Downstream options for union output:
#   1. Subset to cells with all modalities  -> effectively intersection
#   2. Use modality-aware methods that tolerate NaN
#   3. Impute (use with caution — introduces signal)
```

---

## Step 5 — Transfer metadata across modalities

After alignment, propagate spatial / well metadata from a **donor** modality (prefer 3′ when present, else first aligned table) so every AnnData carries the same obs columns.

```python
def transfer_metadata(source_adata, target_adata, columns, cell_col="_cell_id"):
    """Transfer obs columns from source to target using _cell_id as key."""
    source_meta = source_adata.obs[columns + [cell_col]].set_index(cell_col)
    for col in columns:
        target_adata.obs[col] = target_adata.obs[cell_col].map(source_meta[col])
    return target_adata

# Prefer ThreePrimeUntargeted as donor when present
donor_name = next(
    (k for k in ("ThreePrimeUntargeted", "Transcript", "imaging", "CellPaint") if k in aligned),
    next(iter(aligned)),
)
donor = aligned[donor_name]
_tile = "tile_id" if "tile_id" in donor.obs.columns else "Tile"
metadata_cols = [
    c for c in ["Well", _tile, "X_tile", "Y_tile"] if c in donor.obs.columns
]
# Area lives in obsm["cell_geometry"], not obs — do not invent obs["Area"] / obs["X"]

for name, adata in aligned.items():
    if name != donor_name:
        aligned[name] = transfer_metadata(donor, adata, columns=metadata_cols)

print("Metadata donor:", donor_name)
for name, adata in aligned.items():
    print(f"  {name}: {adata.obs.columns.tolist()}")
```

---

## Step 6 — Transform each modality independently

**QC first:** `compute_qc_metrics` + `filter_cells` per modality — see [cell-quality-control](../cell-quality-control/SKILL.md).

| Modality | Transform | Skill |
|---|---|---|
| `ThreePrimeUntargeted` / `Transcript` / `Protein` | `ebsd.pp.normalize` | [normalization-strategies](../normalization-strategies/SKILL.md) |
| `imaging` / `CellPaint` | **MAD-z** feature matrix (exclude Bkg/IgG/Cleave), reference chosen from the design | [imaging-pipeline](../imaging-pipeline/SKILL.md) |
| OPS | Identity filters only — **do not** count-normalize | [ops-pipeline](../ops-pipeline/SKILL.md) |

Do **not** call `ebsd.pp.normalize` on imaging (`NotImplementedError`). Do **not** `StandardScaler` the full imaging `X` (pulls in Bkg/IgG/Cleave and skips the curated matrix). Do **not** assume a per-well reference: when condition is constant within each well, per-well centering removes it — scale against control wells or globally instead ([imaging-pipeline → Scaling](../imaging-pipeline/SKILL.md#scaling)). This matters doubly here, since a joint embedding will otherwise show the RNA modality carrying condition structure the imaging block has had removed.

```python
import elembio_spatialdata_analysis as ebsd

prepared = {}

# --- Count modalities: ebsd.pp.normalize ---
for name, target_type in (
    ("ThreePrimeUntargeted", "ThreePrimeUntargeted"),
    ("Transcript", "Transcript"),
    ("Protein", "Protein"),
):
    if name not in aligned:
        continue
    ad = aligned[name].copy()
    ebsd.pp.normalize(
        ad, target_type, batch_colname="batch_name", well_colname="Well", inplace=True
    )
    prepared[name] = ad

# --- Imaging: MAD-z profiling AnnData via imaging-pipeline (not pp.normalize) ---
# Keys may be "imaging" or "CellPaint" depending on alignment naming.
img_key = next((k for k in ("imaging", "CellPaint") if k in aligned), None)
if img_key is not None:
    # Follow imaging-pipeline IMPLEMENTATION.md §§1–6 on aligned[img_key]:
    # inventory → S/B → feature matrix (exclude Bkg/IgG/Cleave) → MAD-z vs design reference → `prof`
    # Then: prepared[img_key] = prof
    # Do not StandardScaler(aligned[img_key].X) here.
    pass  # replace with `prof` from imaging-pipeline

# OPS tables (custom names): keep filtered counts as-is — do not normalize
for name, ad in aligned.items():
    if name in prepared or name in ("imaging", "CellPaint"):
        continue
    prepared[name] = ad.copy()

normalized = prepared  # name kept for downstream steps that expect this variable
```

---

## Step 7 — Concatenate into a combined feature matrix

Produces a single wide AnnData where every feature is prefixed with its modality name. Useful for PCA / UMAP over all features jointly.

```python
def concatenate_modalities(modalities_dict, modality_key="_modality"):
    """
    Horizontally concatenate modality AnnData objects.

    Cells must already be aligned (same order, same _cell_id values).
    Feature names are prefixed as "<Modality>_<feature>" to avoid collisions.
    uns["modality_boundaries"] records the column slice for each modality.
    """
    first_cells = modalities_dict[next(iter(modalities_dict))].obs["_cell_id"].values
    for name, adata in modalities_dict.items():
        if not np.array_equal(adata.obs["_cell_id"].values, first_cells):
            raise ValueError(f"Cells not aligned for modality '{name}'")

    X_list, var_list = [], []
    for name, adata in modalities_dict.items():
        X = adata.X.toarray() if sp.issparse(adata.X) else adata.X
        X_list.append(X)

        var = adata.var.copy()
        var[modality_key] = name
        var["feature_name"] = adata.var_names
        var.index = [f"{name}_{fn}" for fn in adata.var_names]
        var_list.append(var)

    combined = ad.AnnData(
        X=np.hstack(X_list),
        obs=modalities_dict[next(iter(modalities_dict))].obs.copy(),
        var=pd.concat(var_list),
    )

    boundaries, start = {}, 0
    for name, adata in modalities_dict.items():
        boundaries[name] = (start, start + adata.n_vars)
        start += adata.n_vars
    combined.uns["modality_boundaries"] = boundaries

    return combined

combined = concatenate_modalities(normalized)
print(f"Combined shape: {combined.shape}")
print(f"Modality boundaries: {combined.uns['modality_boundaries']}")
print(combined.var["_modality"].value_counts())
```

---

## Step 8 — Create MuData object

MuData is the recommended container for multimodal data: each modality lives in `mdata.mod[name]` while shared observations are accessible at `mdata.obs`.

```python
try:
    import mudata as md

    mdata = md.MuData({name: adata for name, adata in normalized.items()})
    print(mdata)

    for mod in mdata.mod:
        print(f"{mod} shape: {mdata[mod].shape}")
    print(f"Shared cells:     {mdata.n_obs}")

    mdata.write("multimodal_data.h5mu")

except ImportError:
    print("mudata not installed. Install with: pip install mudata")
    mdata = None
```

---

## Step 9 — Handle batch-specific features

Some protein panels or gene sets differ across batches. Use this helper to audit which features are present in which batches before concatenating across runs.

```python
def get_batch_feature_matrix(adata, batch_col="batch_name"):
    """Return a boolean DataFrame: rows = features, columns = batches."""
    if batch_col not in adata.var.columns:
        print(f"No '{batch_col}' column in .var")
        return None

    batches = adata.var[batch_col].unique()
    feature_presence = {b: set(adata.var_names[adata.var[batch_col] == b]) for b in batches}

    return pd.DataFrame(
        {b: [f in feature_presence[b] for f in adata.var_names] for b in batches},
        index=adata.var_names,
    )

if "batch_name" in transcript_adata.var.columns:
    presence = get_batch_feature_matrix(transcript_adata)
    print("Features per batch:\n", presence.sum())
```

---

## Step 10 — Validate multimodal alignment

Run this after alignment and normalization to catch mismatches before downstream analysis.

```python
def validate_multimodal_data(modalities_dict):
    """Print and return a validation summary for aligned modalities."""
    checks = {
        "n_modalities": len(modalities_dict),
        "modalities": list(modalities_dict.keys()),
        "cell_counts": {n: a.n_obs for n, a in modalities_dict.items()},
        "feature_counts": {n: a.n_vars for n, a in modalities_dict.items()},
    }

    checks["cells_aligned"] = len(set(checks["cell_counts"].values())) == 1

    if all("_cell_id" in a.obs.columns for a in modalities_dict.values()):
        first = modalities_dict[next(iter(modalities_dict))].obs["_cell_id"].values
        checks["cell_ids_match"] = all(
            np.array_equal(a.obs["_cell_id"].values, first)
            for a in modalities_dict.values()
        )
    else:
        checks["cell_ids_match"] = "N/A (no _cell_id column)"

    sep = "=" * 50
    print(sep)
    print("MULTIMODAL DATA VALIDATION")
    print(sep)
    print(f"Modalities:      {checks['modalities']}")
    print(f"Cells aligned:   {checks['cells_aligned']}")
    print(f"Cell IDs match:  {checks['cell_ids_match']}")
    print(f"Cells/modality:  {checks['cell_counts']}")
    print(f"Feats/modality:  {checks['feature_counts']}")
    print(sep)

    return checks

validation = validate_multimodal_data(normalized)
```

---

## Expected Outputs

After completing steps 3–10 you should have:

| Object | Type | Description |
|---|---|---|
| `aligned` | `dict[str, AnnData]` | Modality AnnDatas with identical cell order and `_cell_id` column |
| `normalized` | `dict[str, AnnData]` | Each modality normalized by its own pipeline |
| `combined` | `AnnData` | Wide matrix with all features (optional; useful for joint PCA) |
| `mdata` | `MuData` | Recommended container; write to `.h5mu` for persistence |

Validation checks:
- All modalities report the same `n_obs`
- `_cell_id` arrays are identical across modalities
- Metadata on `.obs` (`Well`, tile, `X_tile`, `Y_tile`) is present in every modality
- Region centroids stay in `obsm["spatial"]` (`X`, `Y`) — already on each table; not an `.obs` transfer

---

## Common Pitfalls

1. **Cell ID type mismatch** — cell IDs are sometimes integers in one modality and strings in another. Cast with `.astype(str)` before building sets.

2. **Index vs column alignment** — pandas aligns on index by default. When building a new AnnData, set `obs.index` explicitly rather than relying on positional order.

3. **Union + zeros vs NaN** — the union helper zero-pads absent cells, which looks like genuine zero expression. Use `_has_data` to mask those cells before statistics.

4. **Normalize after alignment, not before** — normalizing separate populations then aligning can produce different scale factors per modality. Always align first.

5. **Memory for dense concatenation** — `concatenate_modalities` densifies sparse matrices. For large experiments (> 500k cells × 5k genes) consider keeping modalities sparse and using MuData instead.
