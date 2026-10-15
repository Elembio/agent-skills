# Batch Effect Assessment — Reference

## Parameters and Thresholds

| Metric | Low concern | Moderate | High concern |
|---|---|---|---|
| Technical factor R² (variance explained) | < 5% | 5–15% | > 15% |
| Technical silhouette score | < 0.1 | 0.1–0.3 | > 0.3 |
| Within/between replicate correlation ratio | > 1.2 | 1.0–1.2 | < 1.0 |

**PCA components used for variance calculation:** 30 (adjust `n_pcs` for datasets with fewer cells).

**Harmony `max_iter_harmony`:** 20 is sufficient for most datasets; increase to 50 for large, heterogeneous experiments.

**ComBat:** Avoid when conditions are fully confounded with wells (one condition per well). The model will fail or produce meaningless results.

## When to Correct

| Scenario | Recommendation |
|---|---|
| Technical R² > biological R² | Correct — batch dominates |
| Cells cluster by well/run in UMAP | Correct |
| Control samples across wells fail to co-cluster | Correct |
| Batch effect < 5% variance explained | Skip correction |
| Correction removes biological signal | Skip correction |
| One condition per well (confounded) | Cannot correct — redesign experiment |

## Correction Method Comparison

| Method | What it corrects | Output location | Notes |
|---|---|---|---|
| **Harmony** | Embedding-level (PCA) | `obsm["X_pca_harmony"]` | Non-destructive; use for UMAP/clustering only |
| **ComBat** | Expression-level (`.X`) | `adata_combat.X` | Modifies values; use when DE is the goal |

## Expected Outputs

After running batch assessment you should have:

1. **Visualizations** — PCA and UMAP plots colored by `Well`, `tile_id`/`Tile`, `run_id`, total counts, and experimental condition.
2. **Variance metrics** — Mean R² per technical and biological factor across the top 30 PCs.
3. **Silhouette scores** — Per-well and per-run separation on UMAP coordinates.
4. **Replicate correlations** — Within- vs. between-condition pseudobulk Pearson correlations with a heatmap.
5. **A correction decision** — Whether to apply Harmony (embedding-level) or ComBat (expression-level).
6. **Corrected embeddings or expression** (if applicable) — `adata.obsm["X_pca_harmony"]` or `adata_combat.X`.
7. **Batch report** — Console summary with flagged metrics exceeding concern thresholds.

## Common Pitfalls

**Confounded design** — If each biological condition occupies exactly one well, you cannot distinguish batch from biology. No correction method can recover this; redesign the experiment.

**Over-correction** — Aggressive correction (especially ComBat) can remove genuine biological signal. Always compare corrected and uncorrected UMAP/clustering side by side before proceeding.

**Wrong correction variable** — Correcting for `Well` when the true batch is `run_id` (or vice versa) may be ineffective. Check which variable has the highest R² and silhouette score first.

**Post-hoc correction on confounded data** — Batch correction cannot rescue fundamentally confounded experiments. Check design balance before applying any correction.

**`sparse_output` deprecation** — `OneHotEncoder(sparse=False, ...)` is deprecated in recent scikit-learn; use `sparse_output=False` instead.
