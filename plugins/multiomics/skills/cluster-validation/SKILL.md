---
name: cluster-validation
description: Validating cluster robustness, building multi-class classifiers to evaluate cluster reconstruction, assessing technical biases and batch contamination in clusters. Silhouette scores, cross-validation, and stability analysis. Use when the user asks whether clusters are good, robust, reliable, or biologically meaningful, or wants to validate, assess, or sanity-check clustering results. Applies to any AnnData (AVITI24, scRNA-seq-like, or AnnData built from CSV/parquet) with cluster labels in `obs` and a low-dimensional embedding in `obsm`.
metadata:
  version: 0.2.2
  author: elembio
---

# Cluster Validation

## Purpose

This skill covers validating cluster robustness, building classifiers to evaluate cluster reconstruction, assessing technical biases in clusters, and ensuring clusters represent meaningful biological variation. Use it whenever you need confidence that your clusters reflect genuine biology rather than algorithm artifacts or technical noise.

## Prerequisites

- Completed: [cell-quality-control](../cell-quality-control/SKILL.md) and [normalization-strategies](../normalization-strategies/SKILL.md) — clusters should reflect biology, not unfiltered noise
- Completed: [clustering-and-annotation](../clustering-and-annotation/SKILL.md)
- An AnnData object with cluster assignments in `adata.obs["leiden"]` (or any categorical column — pass the column name to the helpers below)

## Key Concepts

### Validation Approaches

| Approach | What It Tests | Primary Metric |
|---|---|---|
| **Logistic regression CV** | How well-defined cluster boundaries are | CV accuracy, F1 |
| **PC coefficients** | Which principal components drive each cluster | Coefficient heatmap |
| **Silhouette scoring** | Compactness and separation of clusters | Silhouette score |
| **Stability (bootstrap)** | Robustness to subsampling | ARI across subsamples |
| **Technical bias check** | Whether clusters are confounded with technical factors | NMI, chi-square |
| **Replicate comparison** | Same clusters reproduce across replicates | Proportion correlation |
| **Biological markers** | Whether clusters align with known cell-type markers | Recall of marker genes |

### When to Validate

Validate clusters:
- Before any downstream analysis (differential expression, trajectory)
- When NMI with a technical variable (batch, well, sequencing depth) is > 0.3
- When a cluster is unexpectedly small (< 50 cells) or very large
- When clustering resolution was changed post-hoc

## Quick Start

Logistic regression cross-validation is the fastest primary check:

```python
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score, StratifiedKFold

X = adata.obsm["X_pca"][:, :30]
y = adata.obs["leiden"].values

clf = LogisticRegression(max_iter=1000, solver="lbfgs", random_state=42)
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
scores = cross_val_score(clf, X, y, cv=cv, scoring="accuracy")

print(f"CV accuracy: {scores.mean():.3f} ± {scores.std():.3f}")
# > 0.9  → well-defined clusters
# 0.7–0.9 → acceptable, investigate low-recall clusters
# < 0.7  → fuzzy/overlapping — consider resolution change or batch correction
```

## Detail Files

- **[IMPLEMENTATION.md](./IMPLEMENTATION.md)** — Full step-by-step code: logistic regression CV, confusion matrix, PC coefficient analysis, technical bias check, replicate comparison, stability (bootstrap), silhouette scoring, borderline cell identification, summary report generator, and recommendations engine.
- **[REFERENCE.md](./REFERENCE.md)** — Metrics interpretation table: accuracy, silhouette, NMI, ARI, correlation — value ranges, what they mean, when each metric is appropriate, and threshold guidance.

## Related Skills

- [clustering-and-annotation](../clustering-and-annotation/SKILL.md) — Perform clustering
- [batch-effect-assessment](../batch-effect-assessment/SKILL.md) — Investigate batch effects before re-clustering
- [cell-quality-control](../cell-quality-control/SKILL.md) — Upstream QC

## Disclaimer

This content is provided "as-is" without warranty of any kind. AI-generated output may contain
inaccuracies. Review and verify all outputs before relying on them.
