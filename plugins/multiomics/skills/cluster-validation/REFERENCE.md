# Cluster Validation — Reference

## Metrics Interpretation Table

| Metric | Range | Interpretation | Action |
|---|---|---|---|
| **CV Accuracy** | > 0.90 | Well-defined, separable clusters | Proceed with downstream analysis |
| **CV Accuracy** | 0.70–0.90 | Acceptable; some cluster overlap | Inspect low-recall clusters in confusion matrix |
| **CV Accuracy** | < 0.70 | Poorly defined clusters | Lower resolution, check batch effects, revisit QC |
| **CV Accuracy** | > 0.98 | May indicate over-clustering | Review at lower resolution; look for biologically trivial splits |
| **Silhouette Score** | 0.50–1.00 | Strong cluster structure | High confidence in cluster assignments |
| **Silhouette Score** | 0.25–0.50 | Reasonable structure with overlap | Acceptable; note boundary clusters |
| **Silhouette Score** | 0.00–0.25 | Weak structure / significant overlap | Re-examine clustering parameters |
| **Silhouette Score** | < 0.00 | Cells closer to another cluster than their own | Likely mis-assigned; consider re-clustering |
| **NMI (technical)** | < 0.10 | Minimal technical confounding | Clusters likely represent biology |
| **NMI (technical)** | 0.10–0.30 | Moderate technical correlation | Acceptable; report for transparency |
| **NMI (technical)** | > 0.30 | Strong technical confounding | Apply batch correction before re-clustering |
| **ARI (stability)** | > 0.80 | Highly stable clusters | Robust to subsampling |
| **ARI (stability)** | 0.60–0.80 | Moderately stable | Clusters may shift slightly with different seeds |
| **ARI (stability)** | < 0.60 | Unstable clusters | Resolution too high or data too sparse |
| **ARI (replicates)** | > 0.80 | Consistent biology across replicates | Strong reproducibility |
| **ARI (replicates)** | 0.50–0.80 | Partial consistency | Investigate which clusters differ between replicates |
| **ARI (replicates)** | < 0.50 | Poor reproducibility | Biological or technical confounders present |
| **Proportion correlation** | > 0.90 | Excellent replicate consistency | |
| **Proportion correlation** | 0.70–0.90 | Good consistency | |
| **Proportion correlation** | < 0.70 | Poor consistency | Check for condition-specific sub-populations or batch effects |

---

## When Each Metric Is Appropriate

### CV Accuracy / F1 (Logistic Regression)

- **Use when:** You want a single summary number for how well-separated clusters are in PCA space.
- **Not appropriate when:** Clusters are intentionally defined on non-PCA embeddings (e.g., protein-only space).
- **Complementary metric:** Confusion matrix for identifying specific cluster pairs that are confused.

### Silhouette Score

- **Use when:** Assessing compactness and separation without a reference labelling.
- **Best computed on:** PCA coordinates (`X_pca[:, :30]`). UMAP coordinates are acceptable for visualisation but should not be used for thresholding decisions.
- **Not appropriate when:** Cluster shapes are non-convex (silhouette assumes convex clusters); Leiden clusters often fail this assumption for very complex topologies.

### NMI (Normalized Mutual Information)

- **Use when:** Comparing clustering against a categorical technical variable (batch, well, lane).
- **Symmetric:** NMI(A, B) = NMI(B, A), so direction does not matter.
- **Not appropriate when:** Comparing against a continuous variable (use correlation or linear regression instead).

### ARI (Adjusted Rand Index)

- **Use when:** Comparing two labellings of the same cells (stability across bootstrap runs, replicate agreement at cell level).
- **Adjusted for chance:** Random labellings produce ARI ≈ 0, making it directly interpretable.
- **Not appropriate when:** The two labellings cover different cell sets (use proportion correlation instead).

### Proportion Correlation

- **Use when:** Comparing cluster composition across samples or replicates that do not share individual cells.
- **Requires aligned cluster labels:** Cluster names must be consistent across samples (e.g., both clustered together or labels matched post-hoc).

---

## Thresholds at a Glance

| Decision point | Threshold | Recommended action |
|---|---|---|
| Proceed to DE analysis | CV accuracy > 0.80 | Safe to proceed |
| Flag for review | CV accuracy 0.70–0.80 | Note in methods; check confusion matrix |
| Re-cluster | CV accuracy < 0.70 | Adjust resolution or run batch correction |
| Flag technical bias | NMI > 0.30 | Apply Harmony / scVI before re-clustering |
| Merge small clusters | Cluster size < 50 | Review biology; merge if no marker support |
| Stable enough for trajectory | ARI (stability) > 0.75 | Proceed |
| Investigate borderline cells | Silhouette < 0.10 | Inspect on UMAP; consider soft-clustering |
