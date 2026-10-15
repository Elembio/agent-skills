---
name: integrated-analysis
description: Joint analysis of aligned multimodal AVITI24 data — 3′ + imaging, 3′ + OPS identity, or other present combinations. Concatenated PCA/UMAP, modality-weighted integration, cross-modality correlation, CCA, and MOFA+. Use for joint embedding or cross-modality stats after multimodal-data-preparation. Imaging inputs should already be MAD-z profiling matrices from imaging-pipeline.
metadata:
  version: 0.2.2
  author: elembio
---

# Integrated Analysis

## Purpose

Joint analysis for **present** modalities (typically 3′ ± imaging ± OPS; legacy Transcript/Protein on older SDOs). Assumes **aligned `AnnData`** from [multimodal-data-preparation](../multimodal-data-preparation/SKILL.md); imaging should be the MAD-z profiling matrix from [imaging-pipeline](../imaging-pipeline/SKILL.md).

## Prerequisites

- [cell-quality-control](../cell-quality-control/SKILL.md)
- [multimodal-data-preparation](../multimodal-data-preparation/SKILL.md)
- Imaging (if used): [imaging-pipeline](../imaging-pipeline/SKILL.md) MAD-z `prof` AnnData

## Key Concepts

### Integration Strategies

| Strategy | Description | When to Use |
|----------|-------------|-------------|
| **Early fusion** | Concatenate features, run joint analysis | Modalities measured on same cells |
| **Late fusion** | Analyze separately, integrate results | Different noise structures |
| **Intermediate fusion** | Project to shared space, then integrate | Factor analysis approaches |

### Scaling and dominance

Before concatenating heterogeneous modalities, scale each to unit variance (`StandardScaler` on that modality’s `X`). Count RNA and imaging MAD-z features have different scales; unscaled concatenation lets one modality dominate. Weights can be proportional to explained variance.

## Quick Start

```python
import scanpy as sc
import anndata as ad
import numpy as np
from sklearn.preprocessing import StandardScaler

# `prepared` from multimodal-data-preparation — keys that exist only
# e.g. {"ThreePrimeUntargeted": rna, "imaging": img_prof}
modalities = prepared

X_list, obs_ref = [], None
for name, adata in modalities.items():
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    X_list.append(StandardScaler().fit_transform(X))
    if obs_ref is None:
        obs_ref = adata.obs.copy()

combined = ad.AnnData(X=np.hstack(X_list), obs=obs_ref)
sc.pp.pca(combined, n_comps=min(50, combined.n_vars - 1))
sc.pp.neighbors(combined, n_neighbors=15, n_pcs=min(30, combined.obsm["X_pca"].shape[1]))
sc.tl.umap(combined)
# Plot → embedding-visualization
```

## Detail Files

- [IMPLEMENTATION.md](./IMPLEMENTATION.md) — Full implementations: concatenated PCA/UMAP, modality-weighted integration, cross-modality correlation, Procrustes alignment, CCA, factor analysis, modality-specific clusters, joint Leiden clustering, integration report.
- [REFERENCE.md](./REFERENCE.md) — Integration strategy comparison, expected outputs, common pitfalls.

## Related Skills

- [multimodal-data-preparation](../multimodal-data-preparation/SKILL.md) — align present modalities (3′ / OPS / imaging)
- [imaging-pipeline](../imaging-pipeline/SKILL.md) — prepare imaging MAD-z matrix before joint work
- [dimensionality-reduction](../dimensionality-reduction/SKILL.md) — compute joint or per-modality embeddings
- [spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md) — load SpatialData when images/masks are needed

## Disclaimer

This content is provided "as-is" without warranty of any kind. AI-generated output may contain
inaccuracies. Review and verify all outputs before relying on them.
