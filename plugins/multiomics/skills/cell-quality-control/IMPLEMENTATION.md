# Cell Quality Control — Implementation

Step-by-step workflows for QC filtering with `elembio_spatialdata_analysis.pp`.
Prefer [run-quality-overview](../run-quality-overview/SKILL.md) first for run-level context.

## Table of Contents

1. [Discover → Prepare QC Metrics](#1-discover--prepare-qc-metrics)
2. [Run `ebsd.pp.filter_cells`](#2-run-ebsdppfilter_cells)
3. [Examine Filter Diagnostics](#3-examine-filter-diagnostics)
4. [Multi-omic — `ebsd.pp.filter_modalities`](#4-multi-omic--ebsdppfilter_modalities)
5. [3′ Untargeted (cells → genes → normalize)](#5-3-untargeted-cells--genes--normalize)
6. [Imaging](#6-imaging)
7. [OPS assay only (DISS identity library)](#7-ops-assay-only-diss-identity-library)
8. [Filter by Morphology (Cell Area)](#8-filter-by-morphology-cell-area)
9. [Filter by Control Counts (NSB, GAPDH) — Legacy Transcript only](#9-filter-by-control-counts-nsb-gapdh--legacy-transcript-only)
10. [Assignment Rate (Exploratory) — Transcript / barcoding context only](#10-assignment-rate-exploratory--transcript--barcoding-context-only)
11. [Visualize Filtering Results](#11-visualize-filtering-results)
12. [Generate Filter Report](#12-generate-filter-report)
13. [Build a Structured QC Summary](#13-build-a-structured-qc-summary)
14. [Common Pitfalls](#common-pitfalls)
15. [Expected Outputs Validation](#expected-outputs-validation)
16. [Manual Fallback](#manual-fallback)

---

## 1. Discover → Prepare QC Metrics

<!-- canonical: spatialdata-loading-and-access memory guard. Copy verbatim; edit only in that skill's SKILL.md -->
**Memory guard — before any `lazy=False`, `.to_memory()` or `.compute()`:**

1. Load lazily first and run `check_before` from `scripts/check_load.py` in `spatialdata-loading-and-access`. Run the script unchanged in the Python session that holds the data (`exec` its contents); don't rewrite it. It compares each table's estimated in-memory size with the memory this session can use, lists large session objects you could `del`, and raises on a FAIL.
2. Materialize only the tables and runs the plan needs. With several runs, load `concatenated=False` so it checks each run and their total. A WARN says what to `del` or to load one run at a time. A FAIL means subset before materializing: propose a subset to the user (runs, wells, or a set of features), then slice the lazy table and call `.to_memory()` on the slice, for example `adata[adata.obs["Well"].isin(wells)].to_memory()`. `.copy()` keeps a lazy table lazy. Feature subsets by modality: [ACCESS.md — Subset](../spatialdata-loading-and-access/ACCESS.md#subset-by-well-or-feature).
3. Keep `.X` in the form it was stored (`scipy.sparse.issparse(adata.X)`). Sparse tables (3′, OPS) are mostly zeros: never densify a whole one (`.toarray()`, `.todense()`, `np.asarray(adata.X)`); densify only a subset of cells or features. Dense tables (imaging, legacy barcoding) are float32 with `NaN` where a batch has no data: do not convert them to sparse, which turns missing into zero.
<!-- /canonical -->

```python
from pathlib import Path
from elembio_spatialdata_analysis.load import Loader
import elembio_spatialdata_analysis as ebsd

STORE_PATH = Path("/path/to/run.zarr")
loader = Loader(STORE_PATH)
present = list(loader.available_tables)
print("present:", present)

# Pick a table that exists — never assume Transcript
# Prefer: ThreePrimeUntargeted → OPS suffixes → imaging/CellPaint → legacy
table = next(
    (t for t in ("ThreePrimeUntargeted", "imaging", "CellPaint", "Transcript", "Protein")
     if t in present),
    None,
)
run = next(iter(loader._runs.values()))
specialized = run.table_suffixes_for_target_type(("OPS", "SpecializedTargeted"))
ops_suffixes = [t for t in specialized if "qc_panel" not in t.lower()]  # OPS assay (DISS)
targeted_diss = [t for t in specialized if t not in ops_suffixes]  # QC/gene panel (also DISS)
if table is None:
    table = ops_suffixes[0] if ops_suffixes else (targeted_diss[0] if targeted_diss else None)
    assert table is not None, "No filterable table in this store"

check_before(loader.load_tables(table=table, concatenated=False))  # memory guard: lazy, raises on FAIL
adata = loader.load_tables(table=table, lazy=False)

# target_type for filter APIs: Imaging for imaging/CellPaint; OPS for OPS assay;
# otherwise use the table suffix as-is (QC/gene panels are DISS but not OPS Confidence QC)
if table in ("imaging", "CellPaint"):
    target_type = "Imaging"
elif table in ops_suffixes:
    target_type = "OPS"
else:
    target_type = table


ebsd.pp.compute_qc_metrics(adata, target_type)
print(adata.obsm["qc_metrics"].columns.tolist())
print(adata.obsm["qc_metrics"].describe())
```

Columns in `obsm["qc_metrics"]` are **modality-specific** and data-dependent (imaging is often sparse;
`nuclear_frac` needs `layers["nuclear"]`). **Discover at runtime** and guard every access.
See [REFERENCE.md](./REFERENCE.md#qc_metrics-columns).

---

## 2. Run `ebsd.pp.filter_cells`

Default: `inplace=True`, `diagnostics=False`. Prefer this form.

> **Guard against double-filtering.** After `filter_cells`, the AnnData has `obsm["filter_flags"]`.
> If that key is already present, the object has been filtered — do not call `filter_cells` again.
> Re-filtering already-filtered cells produces a trivially 100% pass rate.
>
> ```python
> if "filter_flags" not in adata.obsm:
>     ebsd.pp.filter_cells(adata, target_type=target_type)
> ```

```python
ebsd.pp.filter_cells(adata, target_type=target_type)
# Opt-in: drop cells in the modality rollup of obsm["spatial_artifacts"]
# (any-imaging / any-ThreePrimeUntargeted / any-OPS). Default is off.
# ebsd.pp.filter_cells(adata, target_type=target_type, filter_spatial_artifacts=True)
ebsd.pp.generate_filter_report(adata)  # reads from uns["filter_flags_full"] — safe after inplace filter
```

`inplace=False` + `diagnostics=True` when you need before/after copies:

```python
filtered_adata, diag = ebsd.pp.filter_cells(
    adata,
    target_type=target_type,  # "ThreePrimeUntargeted" | "Imaging" | "OPS" | "Transcript" | "Protein"
    inplace=False,
    diagnostics=True,
)
# Report from `diag` (or filtered_adata), not the pre-filter adata
ebsd.pp.generate_filter_report(diag)
```

After filtering counts modalities: **normalize only Transcript / Protein / ThreePrimeUntargeted**.
Do **not** normalize Imaging or OPS.

---

## 3. Examine Filter Diagnostics

```python
filter_flags = diag.flags  # or adata.obsm["filter_flags"] after inplace filter
print(filter_flags.head())

for col in filter_flags.columns:
    n_pass = filter_flags[col].sum()
    pct_pass = n_pass / len(filter_flags) * 100
    print(f"{col}: {n_pass} ({pct_pass:.1f}%)")

fails_by_filter = (~filter_flags).sum()
print("\nCells failing each filter:")
print(fails_by_filter.sort_values(ascending=False))
```

Imaging: look for `passes_imaging_intensity` (legacy alias `passes_imaging_counts` may still appear).

---

## 4. Multi-omic — `ebsd.pp.filter_modalities`

Only pass modalities that **exist**. Area owner is Transcript or ThreePrimeUntargeted when present;
imaging/protein area may be disabled for non-owners (`disable_area_for_non_owners`).

```python
tables = {}
if "ThreePrimeUntargeted" in present:
    tables["ThreePrimeUntargeted"] = loader.load_tables(table="ThreePrimeUntargeted", lazy=False)
if "Transcript" in present:
    tables["Transcript"] = loader.load_tables(table="Transcript", lazy=False)
if "Protein" in present:
    tables["Protein"] = loader.load_tables(table="Protein", lazy=False)
if "imaging" in present or "CellPaint" in present:
    key = "imaging" if "imaging" in present else "CellPaint"
    tables["Imaging"] = loader.load_tables(table=key, lazy=False)

results, pass_summary = ebsd.pp.filter_modalities(tables)
# results[mod].adata — intersected to cells passing all modalities
```

After `filter_modalities`, run `filter_genes` on each **RNA** result (3′ / Transcript) before normalize.

---

## 5. 3′ Untargeted (cells → genes → normalize)

**Required sequence for RNA.** Defaults first; print shapes; escalate `min_cells_per_gene` only with evidence + user OK.

```python
if "ThreePrimeUntargeted" in loader.available_tables:
    untargeted = loader.load_tables(table="ThreePrimeUntargeted", lazy=False)
    ebsd.pp.compute_qc_metrics(untargeted, "ThreePrimeUntargeted")
    ebsd.pp.filter_cells(untargeted, "ThreePrimeUntargeted")  # MT @ 0.10 default
    # Override cell MT only with evidence: filter_cells(..., max_pct_mt=0.2)

    print("after cells:", untargeted.shape)
    ebsd.pp.filter_genes(untargeted, "ThreePrimeUntargeted", min_cells_per_gene=10)
    print("after genes:", untargeted.shape)
    # Tighten if many near-zero genes remain / UMAP will be noisy:
    # ebsd.pp.filter_genes(untargeted, "ThreePrimeUntargeted", min_cells_per_gene=100)
    # Very high floors (hundreds) need large n_obs — ask before applying

    ebsd.pp.normalize(untargeted, "ThreePrimeUntargeted")
```

Same pattern for legacy Transcript (`filter_genes` after `filter_cells`). **Do not** call `filter_genes` on Imaging, OPS, or Protein.

---

## 6. Imaging

```python
img_key = "imaging" if "imaging" in loader.available_tables else (
    "CellPaint" if "CellPaint" in loader.available_tables else None
)
if img_key:
    img = loader.load_tables(table=img_key, lazy=False)
    ebsd.pp.compute_qc_metrics(img, "Imaging")
    ebsd.pp.filter_cells(img, "Imaging")  # passes_imaging_intensity on CP* channels
    # Do NOT normalize imaging
```

---

## 7. OPS assay only (DISS identity library)

`SpecializedTargeted` / `target_type=OPS` also covers **targeted DISS QC/gene panels** (e.g. `B01-QC_Panel`). Those are still DISS, but **not** OPS identity libraries — skip them here; summarize depth in [run-quality-overview](../run-quality-overview/SKILL.md).

```python
run = next(iter(loader._runs.values()))
specialized = run.table_suffixes_for_target_type(("OPS", "SpecializedTargeted"))
ops_tables = [t for t in specialized if "qc_panel" not in t.lower()]
for suffix in ops_tables:
    ops = loader.load_tables(table=suffix, lazy=False, include_target_assignment=True)
    ebsd.pp.compute_qc_metrics(ops, "OPS")
    ebsd.pp.filter_cells(ops, "OPS")  # area only by default; count quantiles off
    # Identity Confidence filter → ops-pipeline (assay panels only)
```

---

## 8. Filter by Morphology (Cell Area)

`ebsd.pp.filter_cells` resolves area from `obsm["cell_geometry"]["Area"]` first, then `obs["Area"]`. For plots, use the same source.

```python
import matplotlib.pyplot as plt
import seaborn as sns

if "cell_geometry" in adata.obsm and "Area" in adata.obsm["cell_geometry"].columns:
    area = adata.obsm["cell_geometry"]["Area"]
elif "Area" in adata.obs.columns:
    area = adata.obs["Area"]
else:
    raise KeyError("No Area in obsm['cell_geometry'] or obs")

fig, axes = plt.subplots(1, 2, figsize=(12, 4))
axes[0].hist(area, bins=50, edgecolor="black", alpha=0.7)
axes[0].set_title("Cell area")

plot_df = adata.obs.copy()
plot_df["_area"] = area.values
sns.boxplot(data=plot_df, x="Well", y="_area", ax=axes[1])
axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=45)
plt.tight_layout()
plt.show()
```

---

## 9. Filter by Control Counts (NSB, GAPDH) — Legacy Transcript only

**Skip this entire section** unless `Transcript` is in `available_tables`. Do not mention NSB/GAPDH for 3′, imaging, OPS, or Protein-only runs.

- **NSB:** columns in `obsm["control_variables"]` (names like `NSB*`).
- **Unassigned / background:** `obsm["unassigned_counts"]` (not mixed into `control_variables`).
- **GAPDH / positive controls:** expression in `.X` for `var["control_type"] == "Positive Control 1"`; `max_quant_gapdh` is **off by default** — set explicitly to filter.

```python
import matplotlib.pyplot as plt
import numpy as np

def plot_control_distributions(adata):
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    if "control_variables" in adata.obsm:
        ctrl_df = adata.obsm["control_variables"]
        nsb_cols = [c for c in ctrl_df.columns if "NSB" in c.upper()]
        if nsb_cols:
            axes[0, 0].hist(ctrl_df[nsb_cols].sum(axis=1), bins=50)
            axes[0, 0].set_title("NSB (sum)")
    mask = adata.var["control_type"] == "Positive Control 1"
    if mask.any():
        X = adata.X
        g = X[:, mask.values]
        gapdh = np.asarray(g.sum(axis=1)).ravel()
        axes[0, 1].hist(gapdh, bins=50)
        axes[0, 1].set_title("GAPDH / positive control (X)")
    if "unassigned_counts" in adata.obsm:
        ua = adata.obsm["unassigned_counts"].sum(axis=1)
        axes[1, 0].hist(ua, bins=50)
        axes[1, 0].set_title("Unassigned (sum across batches)")
    plt.tight_layout()
    plt.show()

plot_control_distributions(adata)
```

---

## 10. Assignment Rate (Exploratory) — Transcript / barcoding context only

Skip unless assignment metrics exist for the modality (typically legacy Transcript). Not OPS confidence. `passes_assignment_rate` applies only if `min_assignment_rate` is not `None` (off by default). You can still plot `total_assignment_rate` from `obsm["qc_metrics"]` after `ebsd.pp.compute_qc_metrics` when that column is present.

```python
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def calculate_assignment_rate(adata):
    tc = adata.obsm["total_counts"]
    ua = adata.obsm["unassigned_counts"]
    assigned = tc.sum(axis=1).values
    unassigned = ua.sum(axis=1).values
    total = assigned + unassigned
    return np.divide(assigned, total, where=total > 0, out=np.zeros_like(assigned, dtype=float))

assignment_rate = calculate_assignment_rate(adata)
adata.obs["assignment_rate"] = assignment_rate

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

axes[0].hist(assignment_rate, bins=50, edgecolor="black")
axes[0].axvline(0.5, color="red", linestyle="--", label="50% threshold")
axes[0].set_xlabel("Assignment Rate")
axes[0].set_ylabel("Count")
axes[0].set_title("Assignment Rate Distribution")
axes[0].legend()

sns.boxplot(data=adata.obs, x="Well", y="assignment_rate", ax=axes[1])
axes[1].axhline(0.5, color="red", linestyle="--")
axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=45)
axes[1].set_title("Assignment Rate by Well")

plt.tight_layout()
plt.show()

print(f"Cells with assignment rate >= 0.5: {(assignment_rate >= 0.5).sum()} / {len(assignment_rate)}")
```

---

## 11. Visualize Filtering Results

> ⚠️ **Never index `obsm["qc_metrics"]` by column name without checking availability first.**
> The schema is modality-dependent — and even within a modality some columns are only
> written when the underlying data exists (see [REFERENCE.md — `qc_metrics` columns](./REFERENCE.md#qc_metrics-columns)).
> Use `if col in df.columns` (or `df.get(col)`) — this eliminates the entire class of `KeyError`
> in multi-modality plotting and aggregation code.

```python
import matplotlib.pyplot as plt

# Plots require obsm["qc_metrics"] — call ebsd.pp.compute_qc_metrics first on both objects
# Each call returns a matplotlib Figure (per-well subplots internally)
fig1 = ebsd.vis.plot_qc_per_well(adata, "total_counts", title="Total counts (before)")
fig2 = ebsd.vis.plot_qc_per_well(filtered_adata, "total_counts", title="Total counts (after)")
plt.show()
```

Modality-agnostic plotting — guard every column, never assume the schema is uniform:

```python
# adata_by_modality: only keys that exist, e.g. {"ThreePrimeUntargeted": ..., "Imaging": ...}
# Each has had ebsd.pp.compute_qc_metrics(adata, name) called on it already.
CANDIDATE_METRICS = ["total_counts", "n_genes_by_counts", "total_assignment_rate", "nuclear_frac", "pct_counts_mt"]

for name, adata in adata_by_modality.items():
    qm = adata.obsm["qc_metrics"]
    present_cols = [col for col in CANDIDATE_METRICS if col in qm.columns]  # guard — schema varies
    if not present_cols:
        continue
    fig, axes = plt.subplots(1, len(present_cols), figsize=(4 * len(present_cols), 3), squeeze=False)
    for ax, col in zip(axes[0], present_cols):
        ax.hist(qm[col].dropna(), bins=50)
        ax.set_title(f"{name}: {col}")
    fig.tight_layout()
    plt.show()

# Overlay raw vs passing cells on the *unfiltered* AnnData using pre-subset flags
fig3 = ebsd.vis.plot_filtered_histograms_per_well(adata, diag.flags, metric="total_counts")
plt.show()
```


---

## 12. Generate Filter Report

```python
import pandas as pd

# Pass the AnnData directly — works after inplace or non-inplace filter_cells
report = ebsd.pp.generate_filter_report(adata)
print(report)

# Also accepts a FilterDiagnostics object or a flags DataFrame if you have them
# report = ebsd.pp.generate_filter_report(diagnostics)

# Custom summary
def summarize_filtering(adata_original, adata_filtered, filter_flags):
    summary = {
        "original_cells": adata_original.n_obs,
        "filtered_cells": adata_filtered.n_obs,
        "removed_cells": adata_original.n_obs - adata_filtered.n_obs,
        "removal_pct": (1 - adata_filtered.n_obs / adata_original.n_obs) * 100,
    }
    for col in filter_flags.columns:
        if col != "passes_all_filters":
            summary[f"fails_{col}"] = (~filter_flags[col]).sum()
    return pd.Series(summary)

summary = summarize_filtering(adata, filtered_adata, diagnostics.flags)
print(summary)
```

---

## 13. Build a Structured QC Summary

Build a canonical, structured QC summary with the helper below, then surface it
however your environment presents results. If the environment provides a
QC-rendering tool, hand it this payload (do not hand-write the JSON); otherwise
print or display it directly. Build the payload **once** with all modalities.
`modalities` is an object keyed by modality name (never a list).

```python
import json
import pandas as pd


def build_qc_render_payload(reports: dict[str, pd.DataFrame], run: str | None = None) -> dict:
    """Reshape per-modality `generate_filter_report` tables into a
    structured QC summary payload.

    reports: {modality_name: ebsd.pp.generate_filter_report(adata_or_diagnostics)}
    The TOTAL row gives survivors (cells_passed) and removed (cells_failed);
    each per-test row maps to a passes_<test> filter.
    """
    modalities: dict[str, dict] = {}
    for name, report in reports.items():
        total = report.loc[report["filter"] == "TOTAL"].iloc[0]  # include_total=True (default)
        n_after = int(total["cells_passed"])
        n_removed = int(total["cells_failed"])
        n_before = n_after + n_removed
        filters = {
            row["filter"]: {
                "cells_passed": int(row["cells_passed"]),
                "cells_failed": int(row["cells_failed"]),
                "pct_failed": round(float(row["pct_failed_total"]) * 100, 2),
            }
            for _, row in report.loc[report["filter"] != "TOTAL"].iterrows()
        }
        modalities[name] = {
            "n_cells_before": n_before,
            "n_cells_after": n_after,
            "n_cells_removed": n_removed,
            "pct_removed": round(n_removed / n_before * 100, 2) if n_before else 0.0,
            "filters": filters,
        }
    payload: dict = {"modalities": modalities}
    if run:
        payload["run"] = run
    return payload
```

Single modality:

```python
report = ebsd.pp.generate_filter_report(adata)            # TOTAL row included by default
payload = build_qc_render_payload({"ThreePrimeUntargeted": report}, run="<run-name>")
print(json.dumps(payload))   # emit the structured summary (one modality)
```

Multi-omic (`filter_modalities`) — only present modalities:

```python
results, _ = ebsd.pp.filter_modalities(tables)
reports = {name: ebsd.pp.generate_filter_report(res.diagnostics) for name, res in results.items()}
payload = build_qc_render_payload(reports, run="<run-name>")
print(json.dumps(payload))   # all present modalities in one payload
```

`generate_filter_report` already excludes `passes_all_filters` and hides disabled
(opt-in, threshold-`None`) filters, so the payload carries exactly the per-test
bars a summary chart should show.

---

## Common Pitfalls

1. **Too aggressive filtering** — Removing too many cells reduces statistical power. Start with default thresholds and adjust only if needed.

2. **Well-specific effects** — Quantile-based filtering is per-well by default. This handles well-to-well variation but may mask systematic issues.

3. **Narrating NSB / GAPDH / assignment-rate on non-Transcript modalities** — Those gates apply to legacy Transcript only. For 3′, imaging, or OPS, do not mention them.

4. **Missing control variables (Transcript)** — If `control_variables` is absent from `obsm`, NSB filtering is skipped silently. Check `.obsm` keys before assuming NSB ran.

5. **3′ data without MT annotation** — `ebsd.pp.compute_qc_metrics` handles MT detection internally, but if `var` lacks mitochondrial gene annotations, `pct_counts_mt` will be zero/absent.

6. **Unnecessary `inplace=False`** — Default is `inplace=True`; use it. Only pass `inplace=False` when you explicitly need to preserve the original for comparison or visualization alongside the filtered result.

7. **Assuming Transcript assignment or GAPDH filters are on by default** — They are **not**. Pass `min_assignment_rate` / `max_quant_gapdh` only for Transcript when needed.

8. **Assuming 3′ MT filtering is off** — For ThreePrimeUntargeted, `max_pct_mt=0.10` is **on by default**. Override only when you want a different threshold.

9. **Inventing absent modalities** — Filter only tables in `available_tables`. Do not run or narrate QC for modalities the store does not have.

10. **Normalizing imaging or OPS after filter** — Do not. Normalize only count modalities (3′, Transcript, Protein).

11. **Indexing `obsm["qc_metrics"]` by column name across modalities** — Schema varies: imaging is often sparse; `n_genes_by_counts` / MT / assignment columns may be absent. Guard with `if col in ...columns` — see [REFERENCE.md](./REFERENCE.md#qc_metrics-columns).

12. **OPS count quantiles** — Off by default. Do not treat OPS `X` like expression depth for filtering; use area (+ confidence elsewhere).

13. **Skipping `filter_genes` on 3′ / Transcript** — Required before normalize. Ultra-sparse genes often produce noisy PCA/UMAP even after cell QC. Gene filter ≠ HVG.

14. **Inventing custom thresholds before running defaults** — Start with package `filter_cells` defaults and `filter_genes(..., min_cells_per_gene=10)`. Inspect reports / shapes / distributions; escalate only with evidence and user approval (especially `min_cells_per_gene` ≫ 100).
15. **`sc.pp.calculate_qc_metrics` instead of `ebsd.pp.compute_qc_metrics`** — Scanpy QC is not the AVITI24 path. Use Element metrics + `filter_cells` so geometry, MT defaults, imaging intensity, and OPS area rules apply.

---

## Expected Outputs Validation

After filtering you should have:

1. **`filtered_adata`** — AnnData with low-quality cells removed
2. **`diagnostics.flags`** — DataFrame with per-cell filter flags
3. **QC visualizations** — Plots showing distributions before/after filtering

```python
# Check filtering was successful
assert filtered_adata.n_obs > 0, "All cells filtered out!"
assert filtered_adata.n_obs < adata.n_obs, "No cells were filtered"

# Check filter flags are stored
assert "filter_flags" in filtered_adata.obsm, "Filter flags not stored"

# Check remaining cells pass all filters
remaining_flags = filtered_adata.obsm["filter_flags"]
assert remaining_flags["passes_all_filters"].all(), "Some remaining cells fail filters"
```

---

## Manual Fallback

If `ebsd.pp.filter_cells` is unavailable (import error, version mismatch, etc.), implement the core filters directly with pandas/numpy. This replicates the two most important defaults — count outliers and area outliers — using the same per-well quantile logic:

```python
import numpy as np
import pandas as pd

# Total counts per cell
total = np.asarray(adata.X.sum(axis=1)).ravel()

keep = pd.Series(True, index=adata.obs_names)

# Per-well quantile filter on total counts (mirrors min_quant_count / max_quant_count)
for well, idx in adata.obs.groupby("Well").groups.items():
    pos = adata.obs_names.get_indexer(idx)
    t = total[pos]
    lo, hi = np.quantile(t, 0.025), np.quantile(t, 0.975)
    keep[idx] &= (t >= lo) & (t <= hi)

# Per-well quantile filter on cell area (mirrors min_quant_area / max_quant_area)
if "cell_geometry" in adata.obsm and "Area" in adata.obsm["cell_geometry"].columns:
    area = adata.obsm["cell_geometry"]["Area"]
elif "Area" in adata.obs.columns:
    area = adata.obs["Area"]
else:
    area = None

if area is not None:
    for well, idx in adata.obs.groupby("Well").groups.items():
        a = area.loc[idx]
        lo, hi = np.quantile(a, 0.01), np.quantile(a, 0.99)
        keep[idx] &= (a >= lo) & (a <= hi)

adata = adata[keep].copy()
print(f"Cells remaining: {adata.n_obs} / {keep.shape[0]}")
```

This does **not** replicate NSB, GAPDH, assignment rate, or MT filters — use `ebsd.pp.filter_cells` when available for the full suite.
