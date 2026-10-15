# Imaging Pipeline — Implementation

Concrete workflows for imaging tables. Prefer `table="imaging"` (alias `CellPaint`). Discover channel names from `var` / `select_morphology_channels` — panels differ by run.

## Table of Contents

1. [Load and inventory](#1-load-and-inventory)
2. [Run QC then cell QC](#2-run-qc-then-cell-qc)
3. [Annotate var with roles](#3-annotate-var-with-roles)
4. [Scope selection](#4-scope-selection)
   - [4a. P01 NaN/zero check](#4a-p01-nanzero-check)
5. [S/B and localization diagnostics](#5-sb-and-localization-diagnostics)
6. [Build analysis matrix and scale](#6-build-analysis-matrix-and-scale)
7. [Hand off to DR + embedding visualization](#7-hand-off-to-dr--embedding-visualization)
8. [Treatment contrasts (well-level)](#8-treatment-contrasts-well-level)
9. [Common pitfalls](#9-common-pitfalls)

## 1. Load and inventory

<!-- canonical: spatialdata-loading-and-access memory guard. Copy verbatim; edit only in that skill's SKILL.md -->
**Memory guard — before any `lazy=False`, `.to_memory()` or `.compute()`:**

1. Load lazily first and run `check_before` from `scripts/check_load.py` in `spatialdata-loading-and-access`. Run the script unchanged in the Python session that holds the data (`exec` its contents); don't rewrite it. It compares each table's estimated in-memory size with the memory this session can use, lists large session objects you could `del`, and raises on a FAIL.
2. Materialize only the tables and runs the plan needs. With several runs, load `concatenated=False` so it checks each run and their total. A WARN says what to `del` or to load one run at a time. A FAIL means subset before materializing: propose a subset to the user (runs, wells, or a set of features), then slice the lazy table and call `.to_memory()` on the slice, for example `adata[adata.obs["Well"].isin(wells)].to_memory()`. `.copy()` keeps a lazy table lazy. Feature subsets by modality: [ACCESS.md — Subset](../spatialdata-loading-and-access/ACCESS.md#subset-by-well-or-feature).
3. Keep `.X` in the form it was stored (`scipy.sparse.issparse(adata.X)`). Sparse tables (3′, OPS) are mostly zeros: never densify a whole one (`.toarray()`, `.todense()`, `np.asarray(adata.X)`); densify only a subset of cells or features. Dense tables (imaging, legacy barcoding) are float32 with `NaN` where a batch has no data: do not convert them to sparse, which turns missing into zero.
<!-- /canonical -->

```python
import numpy as np
import pandas as pd
from elembio_spatialdata_analysis.load import Loader
from elembio_spatialdata_analysis.tl.morphology import select_morphology_channels

loader = Loader("/path/to/run.zarr")
tables = set(loader.available_tables)
assert tables & {"imaging", "CellPaint"}, "No imaging table"

check_before(loader.load_tables(table="imaging", concatenated=False))  # memory guard: lazy, raises on FAIL
img = loader.load_tables(table="imaging", lazy=False)
v = img.var

print("shape:", img.shape)
print("var columns:", v.columns.tolist())
print("geometry:", img.obsm["cell_geometry"].columns.tolist())

for col in ("feature_class", "batch_name", "suffix", "target_type", "is_background"):
    if col in v.columns:
        print(v[col].value_counts(dropna=False))

print("QC morphology channels:", select_morphology_channels(v))
bn = v["batch_name"].dropna().astype(str) if "batch_name" in v.columns else pd.Series(dtype=str)
tt = set(v["target_type"].dropna()) if "target_type" in v.columns else set()
has_protein = "ProteinIF" in tt or bn.str.startswith("P").any()
print("has ProteinIF:", has_protein)
```

### Target inventory

The CP/P batch split is technical — classify targets by what they stain, not which batch they are in. CP01 can contain antibodies; P01 "morphology panels" contain structural/organelle markers.

```python
# Classify targets for the user — base whole-cell features only
if "feature_class" in v.columns:
    base_wc = v[(v["feature_class"] == "base") &
                (v.get("suffix", pd.Series("", index=v.index)).fillna("") == "")]
else:
    base_wc = v[~v.index.to_series().str.contains("_Nuclear|_Cytoplasmic|Fraction|Location_|Granularity|Texture|Intensity_|RadialDistribution")]

_CTRL = re.compile(r"(?i)(Bkg|IgG|Cleave|FinalCleave)")
for _, row in base_wc.iterrows():
    batch = row.get("batch_name", "?")
    ttype = row.get("target_type", "?")
    is_ctrl = bool(_CTRL.search(row.name))
    label = "control" if is_ctrl else f"{batch} / {ttype}"
    print(f"  {row.name:30s}  {label}")
```

Present this inventory to the user. When a target's biology is unclear from the name, ask what it stains — channel names may be abbreviated, suffixed, or non-standard. See [REFERENCE.md](./REFERENCE.md#product-panels-proteinif) for the standard product panels.

## 2. Run QC then cell QC

Prefer [run-quality-overview](../run-quality-overview/SKILL.md) first, then:

```python
import elembio_spatialdata_analysis as ebsd

ebsd.pp.compute_qc_metrics(img, "Imaging")  # total_counts = CP* morphology sum only

# Check for spatial artifacts before filtering
has_artifacts = "spatial_artifacts" in img.obsm
if has_artifacts:
    art = img.obsm["spatial_artifacts"]
    n_flagged = art["any-imaging"].sum() if "any-imaging" in art.columns else 0
    print(f"Spatial artifacts: {n_flagged} / {img.n_obs} cells flagged")
    # If cells are flagged, ask the user whether to opt in:
    #   ebsd.pp.filter_cells(img, "Imaging", filter_spatial_artifacts=True)
    # If zero flagged, the gate is a no-op either way.

ebsd.pp.filter_cells(img, "Imaging")         # NaN/empty/intensity scoped to CP* morphology base
ebsd.pp.generate_filter_report(img)
# Do NOT ebsd.pp.normalize(img, ...)
```

**Spatial artifacts:** If `obsm["spatial_artifacts"]` is present and cells are flagged in `any-imaging`, tell the user how many are affected and ask whether to opt in with `filter_spatial_artifacts=True`. If absent, say nothing.

## 3. Annotate var with roles

Set role flags directly on `img.var`, then use AnnData views.

```python
import re

_CONTROL = re.compile(r"(?i)(^Bkg_|_Bkg_|Bkg\b|IgG|Cleave|FinalCleave)")

def annotate_var_roles(adata):
    """Annotate img.var with role and profiling flags."""
    v = adata.var
    names = v.index.to_series()

    # Feature class
    if "feature_class" in v.columns:
        base = v["feature_class"] == "base"
        cp = v["feature_class"] == "cellprofiler"
    else:
        cp = names.str.startswith(("Granularity_", "Texture_", "Intensity_",
                                   "RadialDistribution_", "Location_"))
        base = ~cp

    # Background
    bkg = v["is_background"].fillna(False) if "is_background" in v.columns \
        else names.str.contains(r"(?i)bkg|background")

    # Controls (Bkg, IgG, Cleave)
    is_control = bkg | names.str.contains(_CONTROL)

    # Absolute Location_* (FOV coordinates — not phenotype)
    location_abs = names.str.startswith("Location_")

    # Batch / scope
    if "batch_name" in v.columns:
        bn = v["batch_name"].fillna("").astype(str).replace("nan", "")
        # CellProfiler features have no batch_name; take it from the stain their name ends with
        stain_batch = bn[base & (bn != "")]
        for stain in sorted(stain_batch.index, key=len, reverse=True):
            bn[cp & (bn == "") & names.str.endswith("_" + stain)] = stain_batch[stain]
        is_cp = bn.str.startswith("CP")
        is_protein = bn.str.startswith("P")
    else:
        is_cp = pd.Series(True, index=v.index)
        is_protein = pd.Series(False, index=v.index)

    # Profiling flag: included in the analysis matrix
    v["_is_base"] = base
    v["_is_cp"] = cp
    v["_is_control"] = is_control
    v["_location_abs"] = location_abs
    v["_scope_morphology"] = is_cp & ~is_control & ~location_abs
    v["_scope_protein"] = is_protein & ~is_control & ~location_abs
    v["_profiling"] = (~is_control) & (~location_abs)

annotate_var_roles(img)
print("profiling features:", img.var["_profiling"].sum())
print("morphology scope:", img.var["_scope_morphology"].sum())
print("protein scope:", img.var["_scope_protein"].sum())
```

The `_scope_morphology` and `_scope_protein` flags follow the technical batch boundary. When the user's question spans both (e.g., "morphology" but the run has a morphology panel in P01), expand the scope by combining flags or by selecting specific targets from the inventory.

## 4. Scope selection

Select features based on the user's question. Always state the scope. After the default selection, offer to expand or narrow based on the target inventory.

```python
def select_scope(adata, scope: str = "morphology"):
    """Return an AnnData view for the requested scope.

    scope: "morphology", "protein", or "all"
    The batch-based flags are the default. For cross-batch scoping (e.g.,
    including morphology-panel P01 targets in a morphology analysis),
    build a custom mask from the target inventory.
    """
    v = adata.var
    if scope == "morphology":
        mask = v["_scope_morphology"]
    elif scope == "protein":
        mask = v["_scope_protein"]
    else:
        mask = v["_profiling"]
    return adata[:, mask]

# Example: morphology analysis (CP01 default)
scoped = select_scope(img, "morphology")
print(f"Scope: morphology (CP01) — {scoped.shape[1]} features")

# Example: expand morphology to include morphology-panel P01 targets
# expanded_mask = img.var["_scope_morphology"] | (img.var["_scope_protein"] & custom_target_mask)
# scoped = img[:, expanded_mask]
```

### 4a. P01 NaN/zero check

`filter_cells("Imaging")` scopes QC to CP01 channels only. When P01 targets are in scope, check for NaN before building the analysis matrix. NaN in a P01 column typically means the antibody was not measured for that cell — it is not zero signal.

```python
if scope in ("protein", "all"):
    p01_cols = img.var_names[img.var["_scope_protein"]].tolist()
    if p01_cols:
        p01_data = np.asarray(img[:, p01_cols].X)
        nan_per_col = np.isnan(p01_data).sum(axis=0)
        for i, col in enumerate(p01_cols):
            if nan_per_col[i] > 0:
                print(f"  {col}: {nan_per_col[i]} cells with NaN ({nan_per_col[i]/img.n_obs*100:.1f}%)")
        nan_any = np.isnan(p01_data).any(axis=1)
        if nan_any.sum() > 0:
            print(f"Total cells with NaN in P01 scope: {nan_any.sum()}")
            # Drop cells with NaN in scoped columns before building the matrix
            img = img[~nan_any].copy()
            print(f"After dropping P01 NaN: {img.n_obs} cells remain")
```

## 5. S/B and localization diagnostics

Compute into `obs` for channel-health plots. **These are diagnostics — never add them to the analysis matrix.**

```python
def _col(adata, name):
    return np.asarray(adata[:, name].X).ravel()

def add_diagnostics(adata, channels: list[str], eps: float = 1e-6):
    """Write obs columns sb_<ch> and nucfrac_<ch> when counterparts exist."""
    for ch in channels:
        bkg = f"Bkg_{ch}"
        if ch in adata.var_names and bkg in adata.var_names:
            adata.obs[f"sb_{ch}"] = _col(adata, ch) / np.maximum(_col(adata, bkg), eps)
        nuc, cyto = f"{ch}_Nuclear", f"{ch}_Cytoplasmic"
        if nuc in adata.var_names and cyto in adata.var_names:
            n, c = _col(adata, nuc), _col(adata, cyto)
            adata.obs[f"nucfrac_{ch}"] = n / np.maximum(n + c, eps)

# Compute for morphology channels
morph_ch = [c for c in select_morphology_channels(img.var)
            if not _CONTROL.search(c)]
add_diagnostics(img, morph_ch)

# If ProteinIF targets exist, compute for those too
if has_protein:
    protein_base = img.var.index[
        img.var["_is_base"] & img.var["_scope_protein"]
        & (img.var.get("suffix", pd.Series("", index=img.var.index)).fillna("") == "")
    ].tolist()
    add_diagnostics(img, protein_base)

# Plot S/B per well — a flat or near-zero channel is dead
sb_cols = [c for c in img.obs.columns if c.startswith("sb_")]
if sb_cols:
    print(img.obs.groupby("Well", observed=True)[sb_cols].median())
```

### IgG validation (any antibody target)

When IgG columns are present, check that antibody signal exceeds the isotype control. This applies to all P01 targets and any antibody-based CP01 targets (use the target inventory to identify them):

```python
igg_cols = [c for c in img.var_names if "IgG" in c and not c.startswith("Bkg_")]
if igg_cols:
    igg_signal = np.asarray(img[:, igg_cols[0]].X).ravel()
    for ch in protein_base:
        ratio = _col(img, ch) / np.maximum(igg_signal, 1e-6)
        print(f"{ch} / IgG median: {np.median(ratio):.1f}x")
        # Ratio >> 1 means the antibody is working; ratio ~ 1 means non-specific
```

## 6. Build analysis matrix and scale

### 6a. Confounding check

```python
def well_confounding(obs, condition_key="condition", well_key="Well"):
    """Is the biology constant within each well?"""
    ct = pd.crosstab(obs[well_key], obs[condition_key])
    conds_per_well = (ct > 0).sum(axis=1)
    confounded = bool((conds_per_well <= 1).all())
    wells_per_cond = (ct > 0).sum(axis=0)
    return {
        "confounded": confounded,
        "min_wells_per_condition": int(wells_per_cond.min()),
        "per_well_centering_allowed": not confounded,
    }

if "condition" in img.obs.columns:
    info = well_confounding(img.obs)
    if info["confounded"]:
        print("Condition is constant within wells — do NOT center per well.")
```

### 6b. Robust-z with guardrails

One core routine. The three scaling options differ only in which rows define the reference.

```python
import anndata as ad

def robust_z(df, ref, global_ref=None, rel_floor=1e-3, clip=10.0):
    """Robust z-score df against median/MAD of ref.

    global_ref: run-wide data for the MAD floor (pass the full matrix when ref
    is a subgroup so near-constant features don't blow up).
    """
    med = ref.median(axis=0)
    mad = (ref - med).abs().median(axis=0) * 1.4826

    g = df if global_ref is None else global_ref
    global_mad = (g - g.median(axis=0)).abs().median(axis=0) * 1.4826
    floor = (global_mad * rel_floor).replace(0.0, np.nan)

    # Degenerate MAD: fall back to std
    degenerate = mad <= floor.fillna(0.0)
    if degenerate.any():
        mad = mad.mask(degenerate, ref.std(axis=0, ddof=0))

    scale = mad.where(mad > 0, other=floor).fillna(1.0)
    out = (df - med) / scale
    out = out.replace([np.inf, -np.inf], np.nan).fillna(0.0)
    return out.clip(-clip, clip) if clip else out
```

### 6c. Assemble and scale

```python
# Select scope
scope = "morphology"  # or "protein" or "all" — based on user intent
scoped = select_scope(img, scope)

# Dense matrix from X
feat = pd.DataFrame(
    np.asarray(scoped.X), index=img.obs_names,
    columns=scoped.var_names.astype(str),
)

# Append geometry (morphology scope default; protein scope only if asked)
if scope in ("morphology", "all"):
    cg = img.obsm["cell_geometry"]
    shape_cols = [c for c in cg.columns if c.startswith("AreaShape_")]
    for c in shape_cols:
        feat[c] = cg[c].to_numpy()

feat = feat.replace([np.inf, -np.inf], np.nan).fillna(0.0)

# Scale — reference from the design
if not info["confounded"]:
    # Wells internally heterogeneous → per-well is a valid technical correction
    scaled = robust_z(feat, feat)
    for well, idx in img.obs.groupby("Well", observed=True).groups.items():
        block = feat.loc[idx]
        scaled.loc[idx] = robust_z(block, block, global_ref=feat)
    scaling_used = "per-well"
elif "control_wells" in dir() and control_wells:
    ref_mask = img.obs["Well"].isin(list(control_wells))
    scaled = robust_z(feat, feat.loc[ref_mask.values], global_ref=feat)
    scaling_used = f"control-referenced ({sorted(control_wells)})"
else:
    scaled = robust_z(feat, feat)
    scaling_used = "global"

print(f"Scaling: {scaling_used} | features: {scaled.shape[1]} | max |z|: {np.abs(scaled.values).max():.1f}")

# AnnData for DR — keep img for raw interpretation
prof = ad.AnnData(X=scaled.values.astype(np.float32), obs=img.obs.copy())
prof.var_names = scaled.columns.astype(str)
prof.uns["scaling_reference"] = scaling_used
prof.uns["scope"] = scope
```

## 7. Hand off to DR + embedding visualization

With `prof` (scaled AnnData) ready:

- Compute embeddings → [dimensionality-reduction](../dimensionality-reduction/SKILL.md) (PCA → neighbors → UMAP; no HVG, no second scale)
- Color / plot → [embedding-visualization](../embedding-visualization/SKILL.md) — **Well and condition on every imaging embedding**, then S/B columns
- Spatial overlays → [spatial-multimodal-visualization](../spatial-multimodal-visualization/SKILL.md)

Label every figure with `prof.uns["scaling_reference"]` and `prof.uns["scope"]`.

On lean feature sets, compare PC1/PC2 alongside UMAP — a clean linear split in early PCs is often clearer than UMAP.

## 8. Treatment contrasts (well-level)

Requires `condition` (and replicates) from [experiment-design](../experiment-design/SKILL.md).

```python
# Scope features by intent — morphology, protein, or both
scope = "protein"  # e.g. "did my protein change with treatment?"
scoped = select_scope(img, scope)

feat = img.obs[["Well", "condition"]].copy()

# Add scoped raw intensities (whole-cell base)
base_mask = img.var["_is_base"] & img.var[f"_scope_{scope}"]
if "suffix" in img.var.columns:
    base_mask &= img.var["suffix"].fillna("").eq("")
base_cols = img.var_names[base_mask].tolist()
for c in base_cols:
    feat[c] = np.asarray(img[:, c].X).ravel()

# Add S/B diagnostics
for c in img.obs.columns:
    if c.startswith("sb_"):
        feat[c] = img.obs[c]

# Add geometry if morphology scope
if scope in ("morphology", "all"):
    cg = img.obsm["cell_geometry"]
    for c in cg.columns:
        if c.startswith("AreaShape_") or c in ("Area", "NuclearArea"):
            feat[c] = cg[c].to_numpy()

# Well-level medians → effect sizes
test_cols = [c for c in feat.columns if c not in ("Well", "condition")]
well_feat = feat.groupby(["Well", "condition"], observed=True)[test_cols].median().reset_index()

groups = well_feat["condition"].unique()
if len(groups) == 2:
    a, b = groups
    A = well_feat.loc[well_feat["condition"] == a, test_cols]
    B = well_feat.loc[well_feat["condition"] == b, test_cols]
    effect = (B.median() - A.median()).sort_values(key=np.abs, ascending=False)
    print(effect.head(20))
```

After the primary scope, offer to check the other: "I can also check whether morphology / protein changed alongside — want me to run that?"

## 9. Common pitfalls

1. **Calling `ebsd.pp.normalize` on imaging** — raises `NotImplementedError`. Use MAD-z.
2. **Ignoring `obsm["cell_geometry"]`** — AreaShape features are first-class morphology data.
3. **Putting Bkg / IgG / Cleave into the analysis matrix** — exclude via var flags. Use Bkg for S/B diagnostics; IgG for antibody validation.
4. **Keeping absolute `Location_*` in UMAP** — FOV pixel coordinates dominate variance via plating layout. Always drop; use `RadialDistribution_*` for within-cell localization.
5. **Centering per well when the well *is* the condition** — per-well robust z removes every between-well effect. Run `well_confounding` first.
6. **Absolute MAD floors** — `mad.clip(lower=1e-6)` inflates near-constant features to ~10^6 and hands PCA to a few columns. Floor relative to global MAD, fall back to std, clip to ±10.
7. **Reading biology off scaled values** — keep raw `X` for S/B and localization plots.
8. **Presenting a confounded plate as a phenotype result** — with one well per condition, well technique and biology are aliased. Report the confounding.
