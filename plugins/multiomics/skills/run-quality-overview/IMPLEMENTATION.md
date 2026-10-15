# Run Quality Overview — Implementation

Discovery-driven recipes for run-level QC.

**Hard rule:** after discovery, compute a `present` list and only execute /
narrate sections for modalities in that list. Do not mention absent modalities.

## Table of Contents

1. [Setup helpers](#0-setup-helpers)
2. [Universal — segmentation / cells per well](#1-universal--segmentation--cells-per-well)
3. [OPS assay (DISS identity library) — Confidence by panel and well](#2a-ops-assay-diss-identity-library--confidence-by-panel-and-well)
4. [Targeted DISS QC / gene panels — depth only (no Confidence gate)](#2b-targeted-diss-qc--gene-panels--depth-only-no-confidence-gate)
5. [ThreePrimeUntargeted — depth and genes](#3-threeprimeuntargeted--depth-and-genes)
6. [Imaging — inventory + area](#4-imaging--inventory--area)
7. [Legacy Transcript / Protein (older SDOs only)](#5-legacy-transcript--protein-older-sdos-only)
8. [Light plotting](#6-light-plotting)
9. [Soft interpret checklist](#7-soft-interpret-checklist)
10. [Common pitfalls](#common-pitfalls)
11. [Expected shapes (examples)](#expected-shapes-examples)

## 0. Setup helpers

```python
import numpy as np
import pandas as pd
from elembio_spatialdata_analysis.load import Loader
from elembio_spatialdata_analysis.qc import per_tile_cell_count

STORE = "/path/to/run.zarr.zip"  # local path or s3:// URI

loader = Loader(STORE)
run = next(iter(loader._runs.values()))
tables, types = set(loader.available_tables), set(loader.available_target_types)
specialized = run.table_suffixes_for_target_type(("OPS", "SpecializedTargeted"))
# Both are DISS SpecializedTargeted — split OPS assay vs QC/gene panel before Confidence
# (is_targeted_diss_qc defined below)
ops_assay_tables, targeted_diss_tables = [], []
for name in specialized:
    if "qc_panel" in name.lower() or name.lower().endswith("-qc") or "control_panel" in name.lower():
        targeted_diss_tables.append(name)
    else:
        ops_assay_tables.append(name)
present = [m for m, ok in [
    ("ThreePrimeUntargeted", "ThreePrimeUntargeted" in tables),
    ("OPS", bool(ops_assay_tables)),
    ("targeted_DISS", bool(targeted_diss_tables)),
    ("imaging", bool(tables & {"imaging", "CellPaint"}) or "imaging" in types),
    ("Transcript", "Transcript" in tables),
    ("Protein", "Protein" in tables),
] if ok]
print("Present:", present)
print("OPS assay:", ops_assay_tables)
print("Targeted DISS/QC:", targeted_diss_tables)

_CYTO_KEYS = ("c2s_runstats_cyto_wells", "runstats_cyto_wells")
_RAWREAD_KEYS = ("c2s_runstats_rawread_wells", "runstats_rawread_wells")

def uns_frame(run, candidates):
    uns = run._metadata
    for key in candidates:
        val = uns.get(key)
        if isinstance(val, pd.DataFrame) and not val.empty:
            return val.reset_index() if val.index.name else val.copy()
    return None

def attach_well_label(df, run, well_col="WellLocation"):
    """Add WellLabel next to WellLocation. Preserve a native WellLabel; backfill blanks."""
    wells = run.wells
    label_col = next((c for c in ("well_label", "WellLabel") if c in wells.columns), None)
    if df is None or df.empty or well_col not in df.columns or label_col is None:
        return df
    out = df.copy()
    mapped = out[well_col].astype(str).map(wells[label_col].astype(str).to_dict())
    if "WellLabel" in out.columns:
        existing = out["WellLabel"]
        out["WellLabel"] = existing.where(
            existing.notna() & (existing.astype(str).str.strip() != ""),
            mapped,
        )
    else:
        out.insert(out.columns.get_loc(well_col) + 1, "WellLabel", mapped)
    return out

def is_targeted_diss_qc(name: str, adata=None) -> bool:
    """True for targeted DISS QC / housekeeping gene panels — NOT OPS assay.

    These are still DISS (``DISSBatch`` / SpecializedTargeted), and their
    metadata ``target_type`` is ``OPS``, the same as OPS identity libraries
    (e.g. a ``B01-QC_Panel`` table with GAPDH+ACTB). Do **not** apply Confidence health gates or hand them to
    ops-pipeline.
    """
    n = name.lower()
    if "qc_panel" in n or n.endswith("-qc") or "control_panel" in n:
        return True
    if adata is not None and "target_gene" in getattr(adata, "var", pd.DataFrame()).columns:
        genes = {str(g).upper() for g in adata.var["target_gene"].dropna()}
        hk = {"GAPDH", "ACTB", "ACTIN", "B2M", "HPRT1"}
        if genes and (genes <= hk or (adata.n_vars <= 20 and bool(genes & hk))):
            return True
    return False

# Refine ambiguous specialized names with var content (housekeeping-only panels)
_refined_assay, _refined_diss = [], list(targeted_diss_tables)
for name in ops_assay_tables:
    ad = loader.load_tables(table=name, lazy=True)
    if is_targeted_diss_qc(name, ad):
        _refined_diss.append(name)
    else:
        _refined_assay.append(name)
ops_assay_tables, targeted_diss_tables = _refined_assay, _refined_diss
present = [m for m in present if m not in ("OPS", "targeted_DISS")]
if ops_assay_tables:
    present.append("OPS")
if targeted_diss_tables:
    present.append("targeted_DISS")
print("Present (refined):", present)
print("OPS assay:", ops_assay_tables)
print("Targeted DISS/QC:", targeted_diss_tables)
```

---

## 1. Universal — segmentation / cells per well

```python
seg = uns_frame(run, _CYTO_KEYS)
if seg is not None:
    seg = attach_well_label(seg, run).sort_values("WellLocation").reset_index(drop=True)
print(seg)
# Columns typically include: WellLocation, CellCount, PercentNucleatedCells,
# PercentConfluency, MedianCellDiameter
# MeanAssignedBarcodingCountPerCell / AssignedCountsPerMM2 are barcoding-assignment metrics;
# do not use them for OPS health.
```

Soft checks (compare wells to peers):

- `CellCount` orders of magnitude below median → dead / empty well
- `PercentNucleatedCells` ≪ ~85 → nucleation concern (cyto is 0–100 scale)
- `PercentConfluency` ≪ 20 or ≫ 70 → seeding / over-confluency concern

<!-- canonical: spatialdata-loading-and-access memory guard. Copy verbatim; edit only in that skill's SKILL.md -->
**Memory guard — before any `lazy=False`, `.to_memory()` or `.compute()`:**

1. Load lazily first and run `check_before` from `scripts/check_load.py` in `spatialdata-loading-and-access`. Run the script unchanged in the Python session that holds the data (`exec` its contents); don't rewrite it. It compares each table's estimated in-memory size with the memory this session can use, lists large session objects you could `del`, and raises on a FAIL.
2. Materialize only the tables and runs the plan needs. With several runs, load `concatenated=False` so it checks each run and their total. A WARN says what to `del` or to load one run at a time. A FAIL means subset before materializing: propose a subset to the user (runs, wells, or a set of features), then slice the lazy table and call `.to_memory()` on the slice, for example `adata[adata.obs["Well"].isin(wells)].to_memory()`. `.copy()` keeps a lazy table lazy. Feature subsets by modality: [ACCESS.md — Subset](../spatialdata-loading-and-access/ACCESS.md#subset-by-well-or-feature).
3. Keep `.X` in the form it was stored (`scipy.sparse.issparse(adata.X)`). Sparse tables (3′, OPS) are mostly zeros: never densify a whole one (`.toarray()`, `.todense()`, `np.asarray(adata.X)`); densify only a subset of cells or features. Dense tables (imaging, legacy barcoding) are float32 with `NaN` where a batch has no data: do not convert them to sparse, which turns missing into zero.
<!-- /canonical -->

Tile dropouts (any cell-shaped table):

```python
# Pick any present table for tile counts
tile_src = next(
    (t for t in ("ThreePrimeUntargeted", "imaging", "CellPaint", "Transcript")
     if t in loader.available_tables),
    None,
)
if tile_src is None:
    ops = run.table_suffixes_for_target_type(("OPS", "SpecializedTargeted"))
    tile_src = ops[0] if ops else None

if tile_src is not None:
    check_before(loader.load_tables(table=tile_src, concatenated=False))  # memory guard: lazy, raises on FAIL
    adata = loader.load_tables(table=tile_src, lazy=False)
    tiles = per_tile_cell_count(adata)
    ct = pd.Series(tiles, name="Cells").rename_axis("tile_id")
    if len(ct):
        print(f"{len(ct)} tiles; median={ct.median():.0f} min={ct.min()} max={ct.max()}")
```

### Spatial artifact detection status

```python
# Load any cell-aligned table to check for spatial_artifacts
_art_src = tile_src or next(iter(loader.available_tables), None)
if _art_src is not None:
    _art_adata = adata if tile_src else loader.load_tables(table=_art_src, lazy=True)
    if "spatial_artifacts" in _art_adata.obsm:
        art = _art_adata.obsm["spatial_artifacts"]
        params = _art_adata.uns.get("spatial_artifact_params", {})
        enabled = params.get("detection_enabled", None)
        reason = params.get("detection_disabled_reason", "")
        rollups = [c for c in art.columns if c.startswith("any")]
        flagged = {c: int(art[c].sum()) for c in rollups}
        print(f"Spatial artifacts: detection_enabled={enabled}"
              + (f" ({reason})" if reason else ""))
        print(f"  Flagged cells: {flagged}")
        # If any rollup has flagged cells, mention this in the QC summary
        # and recommend filter_spatial_artifacts=True in cell-quality-control.
    else:
        print("No obsm['spatial_artifacts'] (tiles-mode, older store, or --skip-spatial-artifacts)")
```

If artifact detection is enabled and cells are flagged, include the counts in the QC summary and recommend opting in with `filter_spatial_artifacts=True` during [cell-quality-control](../cell-quality-control/SKILL.md).

---

## 2a. OPS assay (DISS identity library) — Confidence by panel and well

Skip unless `ops_assay_tables` is non-empty (`"OPS" in present`).

```python
rows = []
for name in ops_assay_tables:
    ops = loader.load_tables(table=name, lazy=False, include_target_assignment=True)
    conf = ops.obs["Confidence"]
    total = ops.obs["TotalCount"]
    undetermined = (
        ops.obs["MostAbundantTarget"].astype(str).eq("Undetermined")
        if "MostAbundantTarget" in ops.obs
        else total.eq(0)
    )
    overall = {
        "panel": name,
        "kind": "ops_assay",
        "n_cells": ops.n_obs,
        "frac_conf_ge_0.5": float((conf >= 0.5).mean()),
        "frac_total_eq_0": float((total == 0).mean()),
        "frac_undetermined": float(undetermined.mean()),
        "median_total_count": float(total.median()),
        "n_vars": ops.n_vars,
    }
    rows.append(overall)
    by_well = (
        ops.obs.assign(conf_ge_0_5=conf >= 0.5, total_eq_0=total == 0)
        .groupby("Well", observed=True)
        .agg(
            n=("conf_ge_0_5", "size"),
            frac_conf_ge_0_5=("conf_ge_0_5", "mean"),
            frac_total_eq_0=("total_eq_0", "mean"),
            median_total=("TotalCount", "median"),
        )
    )
    print(f"\n=== OPS assay {name} ===")
    print(pd.Series(overall))
    print(by_well)

print(pd.DataFrame(rows))
```

**Interpret (assay only):** high `% Confidence>=0.5` and low undetermined is healthy. Near-zero on an assay panel → investigate. Hand off to [ops-pipeline](../ops-pipeline/SKILL.md).

## 2b. Targeted DISS QC / gene panels — depth only (no Confidence gate)

Skip unless `targeted_diss_tables` is non-empty. Example: a `B01-QC_Panel` table (GAPDH+ACTB) — still DISS SpecializedTargeted; metadata may say OPS, but it is a housekeeping gene panel, not an OPS identity library.

```python
for name in targeted_diss_tables:
    panel = loader.load_tables(table=name, lazy=False, include_target_assignment=True)
    assert is_targeted_diss_qc(name, panel)
    depth = panel.obs["TotalCount"] if "TotalCount" in panel.obs else np.asarray(panel.X.sum(axis=1)).ravel()
    print(f"\n=== Targeted DISS/QC {name} ===")
    print("targets:", panel.var_names.tolist())
    if "target_gene" in panel.var.columns:
        print(panel.var["target_gene"].value_counts())
    print(panel.obs.assign(depth=depth).groupby("Well", observed=True)["depth"].agg(["median", "mean", "count"]))
    # Do NOT compute or headline Confidence>=0.5 — low values are expected when probes compete
```

**Do not** hand these tables to ops-pipeline.

---

## 3. ThreePrimeUntargeted — depth and genes

Skip unless `"ThreePrimeUntargeted" in present`.

```python
if "ThreePrimeUntargeted" in present:
    rna = loader.load_tables(table="ThreePrimeUntargeted", lazy=False)
    well_col = next(c for c in ("Well", "well_name") if c in rna.obs.columns)

    tc = rna.obsm.get("total_counts")
    if isinstance(tc, pd.DataFrame):
        depth = tc.sum(axis=1).to_numpy()
    elif tc is not None:
        depth = np.asarray(tc).ravel()
    else:
        depth = np.asarray(rna.X.sum(axis=1)).ravel()

    genes = np.asarray((rna.X > 0).sum(axis=1)).ravel()
    summary = (
        pd.DataFrame({"well": rna.obs[well_col].to_numpy(), "depth": depth, "genes": genes})
        .groupby("well", observed=True)
        .agg(n_cells=("depth", "size"), median_depth=("depth", "median"), median_genes=("genes", "median"))
    )
    print(summary)

    # Optional richer uns (present on some builds)
    for key in ("f2p_3p_runstats_samples", "f2p_3p_runstats_sample_cell_metrics"):
        frame = uns_frame(run, (key,))
        if frame is not None:
            print(key, frame.head())
```

Raw reads (when populated):

```python
rr = uns_frame(run, _RAWREAD_KEYS)
if rr is not None:
    rr = attach_well_label(rr, run)
    # The raw-read column name depends on the store version; resolve it.
    rc = next((c for c in ("NumPoloniesBeforeTrimming", "PFReadCount", "ReadCount") if c in rr.columns), None)
    has = rc is not None and pd.to_numeric(rr[rc], errors="coerce").fillna(0).gt(0).any()
    print(f"raw {rc or 'read count'} populated:", has)
    if has:
        print(rr[["WellLocation", "WellLabel", rc]])
```

---

## 4. Imaging — inventory + area

Skip unless `"imaging" in present`.

```python
if "imaging" in present:
    # Alias: works whether available_tables lists "imaging" or "CellPaint"
    img = loader.load_tables(table="imaging", lazy=False)

    print(img.var["target_type"].value_counts(dropna=False))
    print(img.var["batch_name"].value_counts())
    # Morphology-only: CellPaint (+ NaN cellprofiler rows); IF present: ProteinIF / P* batches

    cg = img.obsm.get("cell_geometry")
    area = cg["Area"] if isinstance(cg, pd.DataFrame) and "Area" in cg.columns else img.obs.get("Area")
    wcol = next(c for c in ("Well", "well_name") if c in img.obs.columns)
    area_df = pd.DataFrame({"well": img.obs[wcol].to_numpy(), "Area": np.asarray(area)})
    print(area_df.groupby("well", observed=True)["Area"].agg(["median", "mean", "std"]))
```

Do **not** run intensity filtering here — that belongs in [cell-quality-control](../cell-quality-control/SKILL.md) (`passes_imaging_intensity`).

---

## 5. Legacy Transcript / Protein (older SDOs only)

Skip unless `"Transcript"` / `"Protein"` appear in `present`.

```python
for tbl in ("Transcript", "Protein"):
    if tbl not in present:
        continue
    ad = loader.load_tables(table=tbl, lazy=False)
    wcol = next(c for c in ("Well", "well_name") if c in ad.obs.columns)
    tc = ad.obsm.get("total_counts")
    depth = tc.sum(axis=1) if isinstance(tc, pd.DataFrame) else (
        np.asarray(tc).ravel() if tc is not None else np.asarray(ad.X.sum(axis=1)).ravel()
    )
    print(tbl, ad.obs.assign(depth=depth).groupby(wcol, observed=True)["depth"].median())
```

---

## 6. Light plotting

```python
import matplotlib.pyplot as plt

if seg is not None:
    metrics = [
        ("CellCount", "Cells per well"),
        ("PercentNucleatedCells", "Nucleated %"),
        ("PercentConfluency", "Confluency %"),
        ("MedianCellDiameter", "Median diameter (px)"),
    ]
    present = [(c, t) for c, t in metrics if c in seg.columns]

    labels = seg["WellLabel"].astype(str) if "WellLabel" in seg.columns else None
    if labels is not None:
        conds = sorted(labels.unique())
        import matplotlib as mpl
        cmap = mpl.colormaps["tab10" if len(conds) <= 10 else "tab20"].resampled(max(len(conds), 2))
        cond_color = {c: cmap(i) for i, c in enumerate(conds)}
        colors = [cond_color[v] for v in labels]
    else:
        colors = None

    fig, axes = plt.subplots(1, len(present), figsize=(4 * len(present), 3.2))
    for ax, (col, title) in zip(np.atleast_1d(axes), present):
        ax.bar(seg["WellLocation"].astype(str), pd.to_numeric(seg[col], errors="coerce"), color=colors)
        ax.set_title(title)
        ax.tick_params(axis="x", rotation=90)

    if labels is not None:
        from matplotlib.patches import Patch

        fig.legend(
            handles=[Patch(facecolor=cond_color[c], label=c) for c in conds],
            loc="lower center", ncol=min(len(conds), 6), fontsize=7,
        )
        fig.tight_layout(rect=[0, 0.12, 1, 1])
    else:
        fig.tight_layout()
    fig.savefig("run_qc_seg_per_well.png", dpi=120)
```

---

## 7. Soft interpret checklist

1. Universal cyto table present? Any dead wells / nucleation / confluency outliers?
2. For each **assay** OPS panel: is `% Confidence>=0.5` sane across wells? Any well-specific collapse?
3. For **QC_Panel** OPS: low single-assign expected — don't flag as assay failure.
4. If 3′ present: depth / genes/cell roughly consistent across wells?
5. If imaging present: morphology-only or +IF? Area distributions sane?
6. Then → [cell-quality-control](../cell-quality-control/SKILL.md) for actual filters.

## Common pitfalls

1. Assuming Transcript / 3′ / OPS always exist — discover first.
2. Treating listed `CellPaint` as “no imaging.”
3. Applying assay OPS thresholds to QC/control panels.
4. Using cyto barcoding-assignment columns (`MeanAssignedBarcodingCountPerCell`, `AssignedCountsPerMM2`) as OPS health metrics.
5. Reading raw-read counts without checking they are populated. If the column is all zeros, use count-modality depth. Resolve the column (`NumPoloniesBeforeTrimming`, `PFReadCount`, or `ReadCount`) rather than assuming one name.
6. Grouping lazy AnnData — always `lazy=False` for these summaries.
7. Filtering cells in this skill — don't; hand off.
8. Expecting `obs["well_label"]` — cells use `WellLabel` (sample-sheet name). `Well` is the plate position. `run.wells` stores the same strings as `well_label`. Cyto frames have `WellLocation` and no label until you backfill.

## Expected shapes (examples)

| Pattern | Expect |
|---|---|
| Full multimodal | 3′ + OPS panels + CellPaint-listed imaging; rich uns |
| Morphology-only | OPS assay + QC_Panel + `imaging` |
| 3′ + OPS + imaging | 3′ + OPS + CellPaint |
| 3′ + imaging | 3′ + CellPaint; no OPS |
| Legacy targeted | Transcript + Protein + CellPaint only |
