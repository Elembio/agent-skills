---
name: code-generation-guidelines
description: Import conventions, discovery-first pipeline order (Loader → run QC → ebsd cell QC → normalize only when appropriate), platform defaults vs generic Scanpy habits, and libraries for local/interactive AVITI24 analysis code. Use when generating elembio analysis code for a local machine or notebook. Never invent modalities absent from available_tables.
metadata:
  version: 0.5.2
  author: elembio
---

# Code Generation Guidelines

## Purpose

Guidelines for Python analysis code for AVITI24 data **for local or interactive (notebook) use** — imports, pipeline order, and performance. Follow these when producing code the user will read and run.

> **Scope.** Human-facing local/interactive workflow. Managed execution environments have their own conventions — follow those there instead.

## Platform defaults (follow before writing code)

AVITI24 is **DISS** (and imaging), not droplet-seq. Prefer Element APIs over generic Scanpy substitutes for load and QC.

| Do | Do not |
|---|---|
| `Loader(path)` + `load_tables` for tables / QC / normalize / DR | `spatialdata.read_zarr("…zarr.zip")` for table analysis |
| `ebsd.pp.compute_qc_metrics` → `ebsd.pp.filter_cells` | `sc.pp.calculate_qc_metrics` as the AVITI24 cell-QC path |
| `ebsd.pp.normalize(adata, target_type)` | `sc.pp.normalize_total` / `sc.pp.log1p` — Scanpy normalization ignores per-well sequencing depth and composition bias; use `ebsd.pp.normalize` for all count modalities |
| Say depth / `total_counts` / polony counts | Call DISS counts “UMIs” |
| Open the specialist skill for the step (load, run QC, cell QC, …) | Improvise from generic Scanpy tutorials after skimming the index |
| Mix Scanpy for HVG / PCA / neighbors / UMAP **after** Element QC + normalize | Mix Scanpy QC with `ebsd.qc` helpers inconsistently in the same workflow |

Scanpy remains appropriate **after** Element preprocessing (HVG, PCA, neighbors, UMAP, Leiden).

## Code Output Convention

Generate **inline notebook cells** (defs + calls). Do not:
- Write standalone `.py` scripts and ask the user to execute them
- Use `subprocess`, shell heredocs, or `python << 'PY'`
- Ask the user to run shell commands for analysis outputs

## Import Convention

```python
import elembio_spatialdata_analysis as ebsd   # all pp / tl / vis / qc via ebsd.*
from elembio_spatialdata_analysis.load import Loader  # exception — explicit import
```

Never `from elembio_spatialdata_analysis.pp import normalize` — call `ebsd.pp.normalize(...)`, `ebsd.pp.filter_cells(...)`, etc.

## Plotting convention

Pass colormaps as **strings** (`cmap="viridis"`) — `sc.pl.umap`, `plt.scatter`, `plt.imshow` all accept strings directly. If you need a colormap object (e.g. to call `cmap(value)` or build a colorbar): `matplotlib.colormaps["viridis"]`. For a discrete N-color palette: `matplotlib.colormaps["tab10"].resampled(N)`. Never use `matplotlib.cm.get_cmap` or `plt.get_cmap` — removed in matplotlib 3.9.

**If a library call fails for an API reason, diagnose it — do not reimplement the library's math.** Substitutes silently change results and drop the guards the real function has.

## Recommended Pipeline Order

**Discover first.** Print `available_tables` / `available_target_types`. Only load, QC, normalize, and recommend pipelines for modalities that exist. Do not mention absent modalities.

1. **Load** — `Loader` + `load_tables(..., lazy=False)` for table analysis. Use `spatialdata.read_zarr` **only** when images/labels/shapes are required (directory store; zip needs FUSE — see loading skill). Imaging: prefer `table="imaging"` (aliases listed `CellPaint`). OPS: discover suffixes via `table_suffixes_for_target_type(("OPS", "SpecializedTargeted"))` — never `table="OPS"`.
2. **Run QC** — [run-quality-overview](../run-quality-overview/SKILL.md): well/segmentation + per-**present**-modality summaries. No cell drops.
3. **Cell QC** — `ebsd.pp.compute_qc_metrics` → `filter_cells` / `filter_modalities` for **present** modalities ([cell-quality-control](../cell-quality-control/SKILL.md)). For **RNA** (3′ / Transcript): also **`filter_genes`** before normalize (`min_cells_per_gene=10` default; escalate only with evidence + approval). Imaging / OPS: cells only — no `filter_genes`. Defaults first; do not invent custom thresholds before inspecting reports.
4. **Normalize / scale** — `ebsd.pp.normalize` on **count** modalities only (3′, legacy Transcript/Protein). Imaging / IF → [imaging-pipeline](../imaging-pipeline/SKILL.md) **MAD-z** feature matrix (exclude Bkg/IgG/Cleave). **Do not** `pp.normalize` imaging or OPS. **Do not hard-code a per-well imaging reference** — check with [experiment-design](../experiment-design/SKILL.md) whether condition is constant within each well; if it is, scale against control wells or globally, because per-well centering deletes the condition.
5. **Downstream** — Counts: HVG (3′ only) → **subset to HVGs** → scale (`sc.pp.scale` densifies; scaling all genes on a large 3′ run OOMs) → [dimensionality-reduction](../dimensionality-reduction/SKILL.md) → [embedding-visualization](../embedding-visualization/SKILL.md). Imaging: MAD-z matrix → same DR + viz skills (no HVG). Treatment hits on imaging → well-level contrasts in imaging-pipeline, not UMAP alone.

## Quick Start

```python
import scanpy as sc
import elembio_spatialdata_analysis as ebsd
from elembio_spatialdata_analysis.load import Loader

loader = Loader("/path/to/run.zarr")
tables = set(loader.available_tables)
print("tables:", sorted(tables))
print("types:", sorted(loader.available_target_types))

# Example path when 3′ is present — skip this block if it is not
if "ThreePrimeUntargeted" in tables:
    adata = loader.load_tables(table="ThreePrimeUntargeted", lazy=False)
    ebsd.pp.compute_qc_metrics(adata, "ThreePrimeUntargeted")
    ebsd.pp.filter_cells(adata, target_type="ThreePrimeUntargeted")
    ebsd.pp.filter_genes(adata, "ThreePrimeUntargeted", min_cells_per_gene=10)  # pass explicitly
    ebsd.pp.normalize(adata, "ThreePrimeUntargeted")
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=False)
    # Scale the HVG subset only — sc.pp.scale densifies the whole matrix
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    sc.pp.pca(adata_hvg, n_comps=50)
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.5, flavor="igraph", n_iterations=2, directed=False)
```

## Detail Files

- [IMPLEMENTATION.md](./IMPLEMENTATION.md) — performance, structure, figures, errors
- [REFERENCE.md](./REFERENCE.md) — available libraries

## Related Skills

- [spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md)
- [run-quality-overview](../run-quality-overview/SKILL.md)
- [cell-quality-control](../cell-quality-control/SKILL.md)
- [normalization-strategies](../normalization-strategies/SKILL.md)
- Pipelines only if present: [3prime-untargeted-pipeline](../3prime-untargeted-pipeline/SKILL.md), [ops-pipeline](../ops-pipeline/SKILL.md), [imaging-pipeline](../imaging-pipeline/SKILL.md)

## Disclaimer

This content is provided "as-is" without warranty of any kind. AI-generated output may contain
inaccuracies. Review and verify all outputs before relying on them.
