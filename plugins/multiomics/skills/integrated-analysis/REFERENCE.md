# Integrated Analysis — Reference

## Integration Strategy Comparison

| Strategy | Description | When to Use |
|---|---|---|
| **Early fusion** | Concatenate features after per-modality scaling, run joint PCA/UMAP | Modalities measured on same cells; fast; interpretable |
| **Weighted fusion** | Scale and weight each modality by explained variance before concatenation | Modalities differ in feature count (e.g., 5k transcript + 50 protein) |
| **Late fusion** | Analyze modalities separately, integrate results (e.g., Procrustes, CCA) | Different noise structures; compare embeddings |
| **Intermediate fusion** | Project to shared latent space (MOFA+, WNN) | Factor analysis; shared vs modality-specific variation |

### Scaling and Dominance Warning

Always scale each modality to unit variance with `StandardScaler` **before** concatenating. A high-dim 3′ matrix can dominate a smaller imaging feature set (or vice versa) without per-modality scaling. Use weights proportional to per-modality explained variance for better balance. Imaging should already be MAD-z’d via imaging-pipeline before this step.

## Key Parameter Guidance

| Parameter | Recommended | Notes |
|---|---|---|
| `n_comps` for per-modality PCA | 30–50 | Reduce before CCA to avoid overfitting |
| `n_pcs` for joint neighbors | 30 | Use joint PCA space |
| CCA components | ≤ min(50, n_cells / 20) | CCA overfits when `n_features ≫ n_cells` |
| MOFA+ factors | 10–20 | Interpret factors with highest R² per modality |

## Expected Outputs

| Output | Description |
|---|---|
| Combined embedding | Joint PCA/UMAP across all modalities (`obsm["X_pca"]`, `obsm["X_umap"]`) |
| Cross-modality correlations | Feature-level transcript ↔ protein Spearman correlations |
| Factor loadings | Proportion of each factor explained by each modality (MOFA+) |
| Cluster comparisons | ARI scores and UMAP panels per modality |
| Integration report | `results/integration_report.md` summary |

## Common Pitfalls

1. **Scale mismatch** — Always scale each modality independently to unit variance before concatenation. Skipping this gives control of PCA axes to the modality with the largest numeric range.

2. **Modality dominance by feature count** — A 20,000-gene transcript panel will dwarf a 50-plex protein panel in concatenated PCA even after scaling. Use weighted concatenation with weights proportional to per-modality explained variance.

3. **CCA overfitting** — CCA maximises correlation and will find spurious structure when `n_features ≫ n_cells`. Reduce dimensionality with PCA first (keep ~50 PCs per modality) before fitting CCA, or use regularized CCA.

4. **Correlation ≠ causality** — High transcript-protein Spearman correlation does not imply a causal link; confounders (cell type, cell cycle, batch) often drive shared variance.

5. **Neighbors computed on wrong space** — `sc.tl.leiden` uses the neighbor graph in `combined.obsp["connectivities"]`. Make sure `sc.pp.neighbors` was called on `combined` after PCA, not on a modality-specific copy.
