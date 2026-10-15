---
name: imaging-pipeline
description: Imaging orchestrator — scope features by user intent (morphology vs protein), design-dependent MAD-z scaling, S/B and IgG diagnostics, well-level treatment contrasts. Hands off PCA/UMAP to dimensionality-reduction and coloring to embedding-visualization. Use when the user mentions imaging, morphology, CellPaint, immunofluorescence, IF protein, ProteinIF, cell painting, or imaging feature scaling. Not for legacy barcoding Protein. Never invent absent modalities.
metadata:
  version: 0.3.8
  author: elembio
---

# Imaging Pipeline

> **Confirm the SDO first.** Before running, make sure you are using the SDO the user selected in [locate-spatialdata-store](../locate-spatialdata-store/SKILL.md) — a run may have several (instrument-direct, `cells2stats`, `elembio-multiomics-nf`).

Orchestrator for the imaging table (`imaging` / `CellPaint` alias). Owns imaging-specific prep: inventory, feature scoping, S/B diagnostics, MAD-z scaling, and treatment contrasts. Hands off PCA/UMAP and plotting to specialist skills.

```text
Load → inventory → run QC → cell QC (Imaging)
  → scope features by intent (morphology / protein / both)
  → S/B + IgG diagnostics (obs only — never in the analysis matrix)
  → [experiment-design] is the biology constant within each well?
  → analysis matrix (annotate var, select view) → MAD-z against the design-implied reference
  → [dimensionality-reduction] PCA/neighbors/UMAP
  → [embedding-visualization] color plots
  → well-level treatment contrasts (this skill) when asking what differs
```

## Hard Rules

- Discover modalities first — never assume IF markers are present.
- Do **not** call `ebsd.pp.normalize` on imaging (`NotImplementedError`). Scaling is MAD-z.
- **Never scale per-well before checking confounding** — per-well centering removes every effect constant within a well, including the biology. See [Scaling](#scaling).
- Legacy barcoding `Protein` ≠ imaging IF.

## Scope by Intent

Everything lives in one imaging table. The user's question determines which features to scope. The `batch_name` split (CP* vs P*) is a **technical grouping**, not a semantic boundary — CP01 can contain antibodies, and P01 morphology panels contain structural/organelle markers that are morphology. Inventory the targets in the run and reason about them.

| User asks about | Default scope | Then offer |
|---|---|---|
| **Morphology** | CP01 targets + AreaShape geometry | Include morphology-panel P01 targets? Other proteins? |
| **Protein / IF** | All P01 targets + CellProfiler for those targets | Clarify which panel or subset if multiple are present |
| **Imaging** (broad) | Both | Present each scope separately and labeled. Clarify if ambiguous. |

**Always state which scope you are using and why.** After addressing the primary scope, offer to check the other — morphology changes alongside protein changes are biologically meaningful.

Which CellProfiler feature types matter most depends on what each target stains — punctate markers carry signal in Texture/Granularity, membrane markers in RadialDistribution, cytoskeletal markers in Texture. See [REFERENCE.md](./REFERENCE.md#feature-relevance-by-target-biology).

**P01 NaN gap:** `filter_cells("Imaging")` scopes QC to CP01 only. NaN in P01 columns survives QC. When P01 targets are in scope, check for and handle NaN before building the matrix. See [IMPLEMENTATION.md §4a](./IMPLEMENTATION.md#4a-p01-nanzero-check).

**Spatial artifacts:** After loading, check if `obsm["spatial_artifacts"]` exists and whether any cells are flagged in `any-imaging`. If so, report the count and ask the user whether to opt in with `filter_spatial_artifacts=True`. See [cell-quality-control](../cell-quality-control/SKILL.md).

## Controls and Diagnostics

These are for validation and channel-health checks. **Never include them in the analysis matrix.**

| Control | What it is | Relevant scope | Use |
|---|---|---|---|
| **Background (`Bkg_*`)** | Paired autofluorescence measurement per channel | Both (CP + Protein) | Compute S/B = `signal / max(bkg, ε)` into `obs`. Plot per-well medians to spot dead or saturated channels. |
| **IgG** | Antibody isotype control | Any antibody target | Validates IF staining — is the antibody signal above the control? Applies to all P01 targets and any antibody-based CP01 targets. |
| **Nuclear fraction** | `Nuclear / (Nuclear + Cytoplasmic)` | Both | Localization diagnostic. Collinear with Nuclear and Cytoplasmic columns already in the data. |
| **Cleave / FinalCleave** | Process markers | Neither | Exclude from everything. |

## Analysis Matrix

Annotate `img.var` with role flags, then select via AnnData views — no separate AnnData.

**Default (morphology scope):** base + CellProfiler + AreaShape, minus Bkg / IgG / Cleave / `Location_*`.

**Always exclude from the matrix:**

| Exclude | Why |
|---|---|
| `Bkg_*` / `is_background` | Background — S/B diagnostic only |
| IgG | Antibody control — validation only |
| Cleave / FinalCleave | Process markers |
| `Location_*` (absolute) | FOV pixel coordinates — positional noise. Keep `RadialDistribution_*` (relative within-cell localization). |

## Scaling

MAD-z (median / 1.4826·MAD). The transform is fixed; the **reference group** is a design decision.

| Is biology constant within each well? | Reference | Why |
|---|---|---|
| **No** — wells internally heterogeneous (pooled / OPS) | Per-well | Well is purely technical |
| **Yes** + shared control wells exist | Control wells (run-level) | Corrects drift, preserves condition differences |
| **Yes** + no shared control | Global (all cells) | Preserves biology; leaves per-well variation in — say so |

Per-well on a confounded plate is a **sensitivity check only** — label it as such.

**Guardrails:** floor MAD relative to the feature's global MAD (not absolute `1e-6`), fall back to std when MAD ≈ 0, clip to ±10, clean `inf`/`NaN` after scaling. See [IMPLEMENTATION.md §6](./IMPLEMENTATION.md#6-build-analysis-matrix-and-scale).

## Treatment Contrasts

When the question is "what differs by treatment?" — stay here (not DESeq2):

1. Confirm groupings with [experiment-design](../experiment-design/SKILL.md).
2. Scope features by intent (morphology / protein / both).
3. Aggregate to well medians per feature.
4. Rank features by effect size. UMAP is a soft check, not the primary evidence.

## Handoffs

| Next question | Skill |
|---|---|
| PCA / neighbors / UMAP | [dimensionality-reduction](../dimensionality-reduction/SKILL.md) — no HVG, no second `sc.pp.scale` |
| Color / plot the embedding | [embedding-visualization](../embedding-visualization/SKILL.md) — always Well + WellLabel/condition, then S/B |
| Is the well structure technical? | [batch-effect-assessment](../batch-effect-assessment/SKILL.md) |
| Count-style `pp.normalize`? | **No** — [normalization-strategies](../normalization-strategies/SKILL.md) is for counts only |

## Quick Start

<!-- canonical: spatialdata-loading-and-access memory guard. Copy verbatim; edit only in that skill's SKILL.md -->
**Memory guard — before any `lazy=False`, `.to_memory()` or `.compute()`:**

1. Load lazily first and run `check_before` from `scripts/check_load.py` in `spatialdata-loading-and-access`. Run the script unchanged in the Python session that holds the data (`exec` its contents); don't rewrite it. It compares each table's estimated in-memory size with the memory this session can use, lists large session objects you could `del`, and raises on a FAIL.
2. Materialize only the tables and runs the plan needs. With several runs, load `concatenated=False` so it checks each run and their total. A WARN says what to `del` or to load one run at a time. A FAIL means subset before materializing: propose a subset to the user (runs, wells, or a set of features), then slice the lazy table and call `.to_memory()` on the slice, for example `adata[adata.obs["Well"].isin(wells)].to_memory()`. `.copy()` keeps a lazy table lazy. Feature subsets by modality: [ACCESS.md — Subset](../spatialdata-loading-and-access/ACCESS.md#subset-by-well-or-feature).
3. Keep `.X` in the form it was stored (`scipy.sparse.issparse(adata.X)`). Sparse tables (3′, OPS) are mostly zeros: never densify a whole one (`.toarray()`, `.todense()`, `np.asarray(adata.X)`); densify only a subset of cells or features. Dense tables (imaging, legacy barcoding) are float32 with `NaN` where a batch has no data: do not convert them to sparse, which turns missing into zero.
<!-- /canonical -->

```python
import elembio_spatialdata_analysis as ebsd
from elembio_spatialdata_analysis.load import Loader
from elembio_spatialdata_analysis.tl.morphology import select_morphology_channels

loader = Loader("/path/to/run.zarr")
check_before(loader.load_tables(table="imaging", concatenated=False))  # memory guard: lazy, raises on FAIL
img = loader.load_tables(table="imaging", lazy=False)

# QC
ebsd.pp.compute_qc_metrics(img, "Imaging")
ebsd.pp.filter_cells(img, "Imaging")

# Annotate var roles and scope — IMPLEMENTATION.md §3-4
# S/B diagnostics — IMPLEMENTATION.md §5
# Build analysis matrix and scale — IMPLEMENTATION.md §6
# DR handoff — IMPLEMENTATION.md §7
```

## Detail Files

- [IMPLEMENTATION.md](./IMPLEMENTATION.md) — var annotation, scope selection, diagnostics, scaling, treatment contrasts
- [REFERENCE.md](./REFERENCE.md) — var schema, feature classes, scaling reference table

## Related Skills

- [spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md)
- [run-quality-overview](../run-quality-overview/SKILL.md) / [cell-quality-control](../cell-quality-control/SKILL.md)
- [experiment-design](../experiment-design/SKILL.md) — condition labels, well-confounding check
- [dimensionality-reduction](../dimensionality-reduction/SKILL.md) — PCA / UMAP on the scaled matrix
- [embedding-visualization](../embedding-visualization/SKILL.md) — color embeddings
- [batch-effect-assessment](../batch-effect-assessment/SKILL.md) — well-vs-condition structure
- [normalization-strategies](../normalization-strategies/SKILL.md) — counts only; imaging is MAD-z here
- [spatial-multimodal-visualization](../spatial-multimodal-visualization/SKILL.md)
- Joint with RNA/OPS: [3prime-untargeted-pipeline](../3prime-untargeted-pipeline/SKILL.md), [ops-pipeline](../ops-pipeline/SKILL.md)
- Legacy barcoding: [legacy-targeted-pipeline](../legacy-targeted-pipeline/SKILL.md)

## Disclaimer

This content is provided "as-is" without warranty of any kind. AI-generated output may contain
inaccuracies. Review and verify all outputs before relying on them.
