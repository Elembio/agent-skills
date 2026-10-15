---
name: normalization-strategies
description: Count-modality normalization for AVITI24 (3′, Transcript, Protein) via ebsd.pp.normalize — aggregated counts, median-of-ratios, log1p. Imaging / IF uses MAD-z in imaging-pipeline, scaled against a control-well, global, or per-well reference depending on the design (pp.normalize Imaging is unimplemented). Use when normalizing counts or asking how to scale imaging features.
metadata:
  version: 0.3.2
  author: elembio
---

# Normalization Strategies

## Purpose

**Count modalities** use `ebsd.pp.normalize` (agg-count → median-of-ratios → log1p for RNA; Protein is lighter). **Imaging / IF / morphology** do **not** use this API — see [imaging-pipeline](../imaging-pipeline/SKILL.md) for MAD-robust z (PyCytominer-style), whose **reference group depends on the experiment design**. Always finish the right transform before dimensionality reduction or DE.

## Prerequisites

- Completed: [Cell Quality Control](../cell-quality-control/SKILL.md)
- A filtered AnnData object

## Key Concepts

### What normalization does (conceptual overview)

Each step addresses a specific source of technical variation:

| Step | What it corrects | Applied to |
|---|---|---|
| **Aggregated-count scaling** | Differences in total counts per cell per batch/well (sequencing depth) — scales each cell to the median total counts for its well+batch | Transcript, ThreePrimeUntargeted, Protein (optional) |
| **Median-of-ratios** | Composition bias — a few highly expressed genes dominating normalization (DESeq2-style) | Transcript, ThreePrimeUntargeted |
| **log1p** | Variance stabilization — compresses dynamic range so high-count genes don't dominate PCA/clustering | Transcript, ThreePrimeUntargeted, Protein |

Raw counts are always snapshotted to `layers["raw"]` before any transformation. Use `layers["raw"]` for DE — never `X`.


### Pipeline Order

**ThreePrimeUntargeted / Transcript** (`ebsd.pp.normalize(adata, "ThreePrimeUntargeted")` or `"Transcript"`):
1. Aggregated-count scaling — scale cells to median total counts per batch/well
2. Median-of-ratios — DESeq2-like correction for composition bias
3. log1p — variance-stabilizing log(1+x) transformation

**Protein** (`ebsd.pp.normalize(adata, "Protein")`):
1. Aggregated-count scaling (optional, controlled by `apply_agg_norm`)
2. log1p — variance-stabilizing transformation

**Imaging** (`"Imaging"` / `"CellPaint"` alias) — `ebsd.pp.normalize` raises `NotImplementedError`. Do **not** call it. Scale imaging features with **MAD-robust z** on a matrix that excludes Bkg / IgG / Cleave — full recipe in [imaging-pipeline](../imaging-pipeline/SKILL.md). OPS: do not normalize.

**There is no per-well default for imaging.** The reference group that supplies median / MAD is a design decision, and picking wrong silently deletes the effect being measured:

| Is the biology constant within a well? | Reference |
|---|---|
| No (pooled / OPS / mixed lines per well) | Per-well |
| Yes, shared control wells exist | Control wells, run-level |
| Yes, no shared control | Global across the run |

Per-well robust z subtracts each well's median, so on a plate where line or treatment is applied per well it removes exactly the biology — replicate wells do not help. Decision procedure and numerical guardrails: [imaging-pipeline → Scaling](../imaging-pipeline/SKILL.md#scaling).

### `inplace` return value

`ebsd.pp.normalize` returns `None` when `inplace=True` (the default). **Do not capture the return value** — use `adata` directly after the call:

```python
ebsd.pp.normalize(adata, "ThreePrimeUntargeted")  # in place — do NOT assign the return value
```

Pass `inplace=False` to get a normalized copy without modifying the original.

### Layer Storage

| Layer | Description |
|---|---|
| `layers["raw"]` | Snapshot of raw counts before any normalization step (always written) |
| `layers["norm_agg_counts"]` | After aggregated-count scaling (opt-in via `agg_norm_layer_name`) |
| `layers["norm_med_of_ratios"]` | After median-of-ratios, before log1p (opt-in via `mor_norm_layer_name`) |
| `X` | Final normalized values |

### DE Note

Pseudobulk DE tools (PyDESeq2, DESeq2) require raw integer counts. Always use `adata.layers["raw"]` for DE — never pass `adata.X`.

## Quick Start

```python
import elembio_spatialdata_analysis as ebsd

# 3′ (preferred) or legacy Transcript: agg-count → median-of-ratios → log1p
ebsd.pp.normalize(adata, "ThreePrimeUntargeted")  # returns None — adata modified in place

# Protein: agg-count scaling → log1p
ebsd.pp.normalize(protein_adata, "Protein")

# To get a copy instead of modifying in place:
normalized = ebsd.pp.normalize(adata, "Transcript", inplace=False)

# Verify
assert "raw" in adata.layers
print(f"Raw range : {adata.layers['raw'].max():.0f}")
print(f".X range  : {adata.X.max():.2f}")

# Imaging: do NOT call ebsd.pp.normalize — use imaging-pipeline MAD-z
```

## Detail Files

- [IMPLEMENTATION.md](./IMPLEMENTATION.md) — count-modality `normalize` workflows; imaging points to imaging-pipeline
- [REFERENCE.md](./REFERENCE.md) — `normalize()` parameter table, expected outputs by modality

## Related Skills

- [cell-quality-control](../cell-quality-control/SKILL.md) — filter before normalizing
- [imaging-pipeline](../imaging-pipeline/SKILL.md) — imaging feature matrix + MAD-z and its design-dependent reference (this skill is counts-only)
- [experiment-design](../experiment-design/SKILL.md) — reports whether condition is confounded with Well
- [dimensionality-reduction](../dimensionality-reduction/SKILL.md) — PCA/UMAP after the right transform
- [pseudobulk-differential-expression](../pseudobulk-differential-expression/SKILL.md) — use `layers["raw"]` for count DE

## Disclaimer

This content is provided "as-is" without warranty of any kind. AI-generated output may contain
inaccuracies. Review and verify all outputs before relying on them.
