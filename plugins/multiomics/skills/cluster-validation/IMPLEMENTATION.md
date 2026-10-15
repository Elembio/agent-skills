# Cluster Validation — Implementation

## Table of Contents

1. [Logistic Regression Cross-Validation](#1-logistic-regression-cross-validation)
2. [Confusion Matrix and Per-Cluster Performance](#2-confusion-matrix-and-per-cluster-performance)
3. [PC Coefficient Analysis](#3-pc-coefficient-analysis)
4. [Technical Bias Check](#4-technical-bias-check)
5. [Replicate Comparison](#5-replicate-comparison)
6. [Stability Analysis (Bootstrap/Subsampling)](#6-stability-analysis-bootstrapsubsampling)
7. [Silhouette Scoring](#7-silhouette-scoring)
8. [Borderline Cell Identification](#8-borderline-cell-identification)
9. [Summary Report Generator](#9-summary-report-generator)
10. [Recommendations Engine](#10-recommendations-engine)
11. [Expected Outputs](#expected-outputs)
12. [Common Pitfalls](#common-pitfalls)

---

## 1. Logistic Regression Cross-Validation

```python
import scanpy as sc
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.metrics import classification_report, confusion_matrix
import matplotlib.pyplot as plt
import seaborn as sns

# Prepare data
X = adata.obsm["X_pca"][:, :30]  # Use PCA for features
y = adata.obs["leiden"].values

# Train logistic regression classifier
clf = LogisticRegression(
    max_iter=1000,
    solver="lbfgs",
    random_state=42,
)

# Cross-validation
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
scores = cross_val_score(clf, X, y, cv=cv, scoring="accuracy")

print(f"Cross-validation accuracy: {scores.mean():.3f} ± {scores.std():.3f}")
print(f"Per-fold scores: {scores}")

# > 0.9  → well-defined clusters
# 0.7–0.9 → acceptable; investigate low-recall clusters
# < 0.7  → fuzzy/overlapping clusters
```

---

## 2. Confusion Matrix and Per-Cluster Performance

```python
from sklearn.model_selection import train_test_split

# Train-test split
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, stratify=y, random_state=42
)

# Fit classifier
clf.fit(X_train, y_train)
y_pred = clf.predict(X_test)

# Per-cluster metrics
print("\nClassification Report:")
print(classification_report(y_test, y_pred))

# Confusion matrix
fig, ax = plt.subplots(figsize=(10, 8))
cm = confusion_matrix(y_test, y_pred)
sns.heatmap(cm, annot=True, fmt="d", cmap="Blues", ax=ax,
            xticklabels=clf.classes_, yticklabels=clf.classes_)
ax.set_xlabel("Predicted")
ax.set_ylabel("True")
ax.set_title("Cluster Classification Confusion Matrix")
plt.tight_layout()
plt.show()

# Clusters with low recall may be unstable or overlapping with neighbours
```

---

## 3. PC Coefficient Analysis

Identify which principal components and genes drive each cluster's separation.

```python
# Logistic regression coefficients indicate PC importance per cluster
coef_df = pd.DataFrame(
    clf.coef_,
    columns=[f"PC{i+1}" for i in range(clf.coef_.shape[1])],
    index=clf.classes_
)

# Plot top 20 PCs for all clusters
fig, ax = plt.subplots(figsize=(12, 6))
sns.heatmap(coef_df.iloc[:, :20], cmap="RdBu_r", center=0, ax=ax)
ax.set_title("Logistic Regression Coefficients (PCs)")
ax.set_xlabel("Principal Component")
ax.set_ylabel("Cluster")
plt.tight_layout()
plt.show()


def get_important_genes_for_cluster(adata, clf, cluster_id, n_genes=10):
    """Return genes most important for classifying a specific cluster."""
    cluster_idx = list(clf.classes_).index(cluster_id)
    pc_importance = clf.coef_[cluster_idx]

    # Get PC loadings stored in adata.varm["PCs"]
    pc_loadings = adata.varm["PCs"][:, :len(pc_importance)]

    # Weight loadings by classifier importance
    gene_importance = np.abs(pc_loadings @ pc_importance)

    top_idx = np.argsort(gene_importance)[-n_genes:][::-1]
    top_genes = adata.var_names[top_idx]
    top_scores = gene_importance[top_idx]

    return list(zip(top_genes, top_scores))


# Get important genes for each cluster
for cluster in clf.classes_[:3]:  # First 3 clusters as example
    genes = get_important_genes_for_cluster(adata, clf, cluster)
    print(f"\nCluster {cluster} — top classifier genes: {[g for g, _ in genes[:5]]}")
```

---

## 4. Technical Bias Check

Detect whether clusters are confounded with technical variables (well, batch, sequencing depth).

```python
from scipy.stats import chi2_contingency, entropy
from sklearn.metrics import normalized_mutual_info_score


def assess_technical_bias(adata, cluster_col, tech_col):
    """Assess if clusters are confounded with a technical variable."""
    crosstab = pd.crosstab(adata.obs[cluster_col], adata.obs[tech_col])

    # Chi-square test
    chi2, p_value, dof, expected = chi2_contingency(crosstab)

    # Normalized mutual information
    nmi = normalized_mutual_info_score(
        adata.obs[cluster_col].astype(str),
        adata.obs[tech_col].astype(str)
    )

    # Per-cluster entropy (higher = more evenly distributed across tech groups)
    cluster_entropy = {}
    for cluster in crosstab.index:
        props = crosstab.loc[cluster] / crosstab.loc[cluster].sum()
        cluster_entropy[cluster] = entropy(props)

    return {
        "chi2": chi2,
        "p_value": p_value,
        "nmi": nmi,
        "cluster_entropy": cluster_entropy,
        "crosstab": crosstab,
    }


# Assess bias with respect to Well
well_bias = assess_technical_bias(adata, "leiden", "Well")
print(f"Chi-square test p-value: {well_bias['p_value']:.2e}")
print(f"Normalized Mutual Information: {well_bias['nmi']:.3f}")
print("\nPer-cluster entropy (higher = less biased):")
for cluster, ent in well_bias["cluster_entropy"].items():
    print(f"  Cluster {cluster}: {ent:.3f}")

# Visualize
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Heatmap: normalized cluster composition by well
ax = axes[0]
crosstab_norm = well_bias["crosstab"].div(well_bias["crosstab"].sum(axis=1), axis=0)
sns.heatmap(crosstab_norm, cmap="YlOrRd", annot=True, fmt=".2f", ax=ax)
ax.set_title("Cluster Composition by Well (normalized)")

# Entropy bar plot per cluster
ax = axes[1]
ent_df = pd.Series(well_bias["cluster_entropy"])
ent_df.plot(kind="bar", ax=ax)
max_entropy = entropy(np.ones(len(adata.obs["Well"].unique())) /
                      len(adata.obs["Well"].unique()))
ax.axhline(max_entropy, color="red", linestyle="--", label="Max entropy (uniform)")
ax.set_xlabel("Cluster")
ax.set_ylabel("Entropy")
ax.set_title("Cluster Technical Entropy")
ax.legend()

plt.tight_layout()
plt.show()
```

---

## 5. Replicate Comparison

Check whether cluster proportions are consistent across biological replicates of the same condition.

```python
def compare_replicate_clustering(adata, condition_col, replicate_col="Well"):
    """Compare cluster proportions across replicates of the same condition."""
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

    conditions = adata.obs[condition_col].unique()
    results = []

    for cond in conditions:
        cond_mask = adata.obs[condition_col] == cond
        cond_adata = adata[cond_mask]

        replicates = cond_adata.obs[replicate_col].unique()
        if len(replicates) < 2:
            continue

        # Compare all pairs of replicates
        for i, rep1 in enumerate(replicates):
            for rep2 in replicates[i + 1:]:
                rep1_clusters = cond_adata[cond_adata.obs[replicate_col] == rep1].obs["leiden"]
                rep2_clusters = cond_adata[cond_adata.obs[replicate_col] == rep2].obs["leiden"]

                # Compare cluster proportions
                prop1 = rep1_clusters.value_counts(normalize=True)
                prop2 = rep2_clusters.value_counts(normalize=True)

                # Align indices
                all_clusters = set(prop1.index) | set(prop2.index)
                prop1 = prop1.reindex(all_clusters, fill_value=0)
                prop2 = prop2.reindex(all_clusters, fill_value=0)

                corr = np.corrcoef(prop1, prop2)[0, 1]

                results.append({
                    "condition": cond,
                    "rep1": rep1,
                    "rep2": rep2,
                    "proportion_correlation": corr,
                })

    return pd.DataFrame(results)


if "experimental_condition" in adata.obs.columns:
    rep_consistency = compare_replicate_clustering(adata, "experimental_condition")
    print("Replicate consistency (cluster proportion correlation):")
    print(rep_consistency)
    # > 0.8 → consistent clustering across replicates
```

---

## 6. Stability Analysis (Bootstrap/Subsampling)

Measure how reproducible clusters are when cells are subsampled.

```python
def assess_cluster_stability(adata, n_iterations=10, subsample_frac=0.8):
    """Assess cluster stability via repeated subsampling and re-clustering."""
    from sklearn.metrics import adjusted_rand_score

    original_clusters = adata.obs["leiden"].values
    stability_scores = []

    for i in range(n_iterations):
        n_sample = int(len(adata) * subsample_frac)
        idx = np.random.choice(len(adata), n_sample, replace=False)

        adata_sub = adata[idx].copy()
        sc.pp.neighbors(adata_sub, n_neighbors=15, n_pcs=30)
        sc.tl.leiden(adata_sub, resolution=0.5, random_state=i)

        ari = adjusted_rand_score(
            original_clusters[idx],
            adata_sub.obs["leiden"].values
        )
        stability_scores.append(ari)

    return {
        "mean_ari": np.mean(stability_scores),
        "std_ari": np.std(stability_scores),
        "scores": stability_scores,
    }


# Note: can be slow for large datasets; use n_iterations=3 for a quick check
stability = assess_cluster_stability(adata, n_iterations=5)
print(f"Cluster stability (ARI): {stability['mean_ari']:.3f} ± {stability['std_ari']:.3f}")
# > 0.8 → highly stable; 0.6–0.8 → moderately stable; < 0.6 → unstable
```

---

## 7. Silhouette Scoring

Measure how well each cell fits its assigned cluster vs. the nearest neighbouring cluster.

```python
from sklearn.metrics import silhouette_samples, silhouette_score

X_umap = adata.obsm["X_umap"]
labels = adata.obs["leiden"].astype("category").cat.codes

# Overall silhouette score
overall_silhouette = silhouette_score(X_umap, labels)
print(f"Overall silhouette score: {overall_silhouette:.3f}")

# Per-cell silhouette scores
sample_silhouettes = silhouette_samples(X_umap, labels)
adata.obs["silhouette"] = sample_silhouettes

# Per-cluster mean silhouette
cluster_silhouettes = adata.obs.groupby("leiden")["silhouette"].mean()
print("\nPer-cluster silhouette scores:")
print(cluster_silhouettes.sort_values())

# Visualize
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

sc.pl.umap(adata, color="silhouette", ax=axes[0], show=False,
           cmap="RdYlGn", vmin=-0.5, vmax=0.5,
           title="Silhouette Scores on UMAP")

cluster_silhouettes.plot(kind="bar", ax=axes[1])
axes[1].axhline(0, color="black", linestyle="--")
axes[1].set_xlabel("Cluster")
axes[1].set_ylabel("Mean Silhouette")
axes[1].set_title("Silhouette Score by Cluster")

plt.tight_layout()
plt.show()
```

> **Note:** Using UMAP coordinates for silhouette is an approximation. For more rigorous results use PCA: replace `X_umap` with `adata.obsm["X_pca"][:, :30]`.

---

## 8. Borderline Cell Identification

Find cells sitting at cluster boundaries — low silhouette scores indicate ambiguous cluster membership.

```python
threshold = 0.1
borderline_mask = adata.obs["silhouette"] < threshold
n_borderline = borderline_mask.sum()
print(f"Borderline cells (silhouette < {threshold}): "
      f"{n_borderline} ({n_borderline / len(adata):.1%})")

# Visualize borderline cells on UMAP
fig, ax = plt.subplots(figsize=(10, 8))
sc.pl.umap(adata[~borderline_mask], color="leiden", ax=ax, show=False,
           title="Cells with silhouette >= 0.1", alpha=0.8)
ax.scatter(
    adata[borderline_mask].obsm["X_umap"][:, 0],
    adata[borderline_mask].obsm["X_umap"][:, 1],
    c="red", s=5, alpha=0.3, label="Borderline",
)
ax.legend()
plt.show()

# Which clusters have the most borderline cells?
borderline_by_cluster = adata.obs[borderline_mask]["leiden"].value_counts()
total_by_cluster = adata.obs["leiden"].value_counts()
borderline_fraction = borderline_by_cluster / total_by_cluster
print("\nFraction borderline by cluster:")
print(borderline_fraction.sort_values(ascending=False))
```

---

## 9. Summary Report Generator

Aggregate all validation metrics into a single printable report.

```python
def generate_cluster_validation_report(adata, tech_vars=["Well"]):
    """Generate a comprehensive cluster validation report."""
    report = {
        "n_cells": adata.n_obs,
        "n_clusters": adata.obs["leiden"].nunique(),
    }

    # Cluster sizes
    sizes = adata.obs["leiden"].value_counts()
    report["min_cluster_size"] = sizes.min()
    report["max_cluster_size"] = sizes.max()
    report["median_cluster_size"] = sizes.median()

    # Classifier performance
    X = adata.obsm["X_pca"][:, :30]
    y = adata.obs["leiden"].values
    clf = LogisticRegression(max_iter=1000, random_state=42)
    cv_scores = cross_val_score(clf, X, y, cv=5, scoring="accuracy")
    report["classifier_accuracy"] = cv_scores.mean()
    report["classifier_std"] = cv_scores.std()

    # Silhouette (UMAP-based approximation)
    labels = adata.obs["leiden"].astype("category").cat.codes
    sil = silhouette_score(adata.obsm["X_umap"], labels)
    report["silhouette_score"] = sil

    # Technical bias
    for tech_var in tech_vars:
        if tech_var in adata.obs.columns:
            bias = assess_technical_bias(adata, "leiden", tech_var)
            report[f"nmi_with_{tech_var}"] = bias["nmi"]

    # Print summary
    print("=" * 60)
    print("CLUSTER VALIDATION REPORT")
    print("=" * 60)
    print(f"Dataset: {report['n_cells']} cells, {report['n_clusters']} clusters")
    print()
    print("Cluster Sizes:")
    print(f"  Min: {report['min_cluster_size']}, "
          f"Max: {report['max_cluster_size']}, "
          f"Median: {report['median_cluster_size']}")
    print()
    print("Cluster Quality:")
    print(f"  Classifier accuracy: {report['classifier_accuracy']:.3f} "
          f"± {report['classifier_std']:.3f}")
    print(f"  Silhouette score:    {report['silhouette_score']:.3f}")
    print()
    print("Technical Bias (NMI with technical variables):")
    for tech_var in tech_vars:
        key = f"nmi_with_{tech_var}"
        if key in report:
            marker = "WARN" if report[key] > 0.3 else "OK"
            print(f"  {tech_var}: {report[key]:.3f}  [{marker}]")
    print("=" * 60)

    return report


report = generate_cluster_validation_report(adata, tech_vars=["Well"])
```

---

## 10. Recommendations Engine

Translate numeric validation results into actionable guidance.

```python
def get_cluster_recommendations(report):
    """Generate recommendations based on a validation report dict."""
    recommendations = []

    if report["classifier_accuracy"] < 0.7:
        recommendations.append(
            "LOW CLASSIFIER ACCURACY: Clusters may be poorly defined. "
            "Consider lowering resolution or checking for batch effects."
        )

    if report["silhouette_score"] < 0.1:
        recommendations.append(
            "LOW SILHOUETTE SCORE: Clusters significantly overlap. "
            "Consider alternative clustering parameters or different n_pcs."
        )

    if report["min_cluster_size"] < 50:
        recommendations.append(
            f"SMALL CLUSTER: Minimum cluster has only {report['min_cluster_size']} cells. "
            "Consider merging small clusters or lowering resolution."
        )

    for key in report:
        if key.startswith("nmi_with_") and report[key] > 0.3:
            tech_var = key.replace("nmi_with_", "")
            recommendations.append(
                f"TECHNICAL BIAS: High correlation with {tech_var} "
                f"(NMI={report[key]:.2f}). "
                "Consider batch correction before clustering."
            )

    if not recommendations:
        recommendations.append("Clustering appears robust. No major issues detected.")

    print("\nRECOMMENDATIONS:")
    for i, rec in enumerate(recommendations, 1):
        print(f"{i}. {rec}")

    return recommendations


recommendations = get_cluster_recommendations(report)
```

---

## Expected Outputs

After running the full validation workflow you should have:

1. **CV accuracy** — single number summarising cluster separability (stored in `report["classifier_accuracy"]`)
2. **Confusion matrix plot** — reveals which clusters are confused with each other
3. **PC coefficient heatmap** — shows which PCs differentiate each cluster
4. **Technical bias NMI** — per-technical-variable contamination score (stored as `report["nmi_with_<var>"]`)
5. **Replicate proportion correlation** — consistency score across biological replicates
6. **Stability ARI** — robustness score from bootstrap subsampling
7. **Silhouette scores** — per-cell scores added to `adata.obs["silhouette"]`; per-cluster means printed
8. **Borderline cell mask** — UMAP overlay showing ambiguous cells
9. **Printed report** — all metrics in one place
10. **Recommendations list** — actionable next steps

---

## Common Pitfalls

1. **Over-interpreting very high accuracy**: CV accuracy > 0.98 often indicates over-clustering (too many clusters, each with slight expression differences). Re-examine at a lower resolution.

2. **Ignoring small clusters**: Clusters with fewer than 50 cells may be biologically important doublets or rare populations — investigate before merging.

3. **UMAP-based silhouette is an approximation**: UMAP distorts distances. Compute silhouette on `X_pca[:, :30]` for a more faithful result; use UMAP-based scores for visual intuition only.

4. **Confusing NMI with causation**: High NMI with a technical variable means correlation, not that the clusters are artifacts. They may still be biologically real — but batch correction is warranted before concluding either way.

5. **Stability analysis resolution mismatch**: The bootstrap re-clusters at the same resolution as the original. If you tuned resolution manually, ensure the `sc.tl.leiden` call in `assess_cluster_stability` uses the same value.
