# Integrated Analysis — Implementation

Implementation recipes for the [integrated-analysis skill](./SKILL.md).

## Table of Contents

1. [Concatenated PCA/UMAP (Early Fusion)](#1-concatenated-pcaumap-early-fusion)
2. [Modality-Weighted Integration](#2-modality-weighted-integration)
3. [Cross-Modality Correlation Analysis](#3-cross-modality-correlation-analysis)
4. [Procrustes Alignment — Compare Embeddings](#4-procrustes-alignment--compare-embeddings)
5. [Canonical Correlation Analysis (CCA)](#5-canonical-correlation-analysis-cca)
6. [Factor Analysis for Shared Variation](#6-factor-analysis-for-shared-variation)
7. [Identify Modality-Specific Signals](#7-identify-modality-specific-signals)
8. [Joint Leiden Clustering](#8-joint-leiden-clustering)
9. [Generate Integration Report](#9-generate-integration-report)
10. [Expected Outputs](#expected-outputs)
11. [Common Pitfalls](#common-pitfalls)

---

## 1. Concatenated PCA/UMAP (Early Fusion)

```python
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler

def create_concatenated_analysis(modalities, scale_features=True):
    """
    Concatenate modalities and run joint PCA/UMAP.
    """
    X_list = []
    var_info = []

    for name, adata in modalities.items():
        X = adata.X.copy()
        if hasattr(X, "toarray"):
            X = X.toarray()

        if scale_features:
            X = StandardScaler().fit_transform(X)

        X_list.append(X)

        for var_name in adata.var_names:
            var_info.append({"feature": var_name, "modality": name})

    X_combined = np.hstack(X_list)
    var_df = pd.DataFrame(var_info)
    var_df.index = var_df["modality"] + "_" + var_df["feature"]

    combined = ad.AnnData(
        X=X_combined,
        obs=list(modalities.values())[0].obs.copy(),
        var=var_df,
    )
    return combined


# Build modality dict from whatever was prepared (discover keys — do not hardcode)
# Example: ThreePrimeUntargeted + imaging (MAD-z prof)
# Legacy: Transcript / Protein; table alias CellPaint == imaging
modalities = {k: v for k, v in prepared.items()}

combined = create_concatenated_analysis(modalities, scale_features=True)
print(f"Combined shape: {combined.shape}")

sc.pp.pca(combined, n_comps=50)
sc.pp.neighbors(combined, n_neighbors=15, n_pcs=30)
sc.tl.umap(combined)

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
sc.pl.umap(combined, color="Well", ax=axes[0], show=False, title="Combined UMAP by Well")
if "experimental_condition" in combined.obs.columns:
    sc.pl.umap(combined, color="experimental_condition", ax=axes[1], show=False,
               title="Combined UMAP by Condition")
plt.tight_layout()
plt.show()
```

---

## 2. Modality-Weighted Integration

Use when one modality should contribute more signal than others. Weight proportional to explained variance is generally more principled than arbitrary values.

```python
def weighted_concatenation(modalities, weights=None):
    """
    Concatenate modalities with per-modality scalar weights.

    Parameters
    ----------
    modalities : dict
        Modality name -> AnnData
    weights : dict
        Modality name -> float weight (default: equal weights)
    """
    if weights is None:
        weights = {name: 1.0 for name in modalities}

    X_list = []
    for name, adata in modalities.items():
        X = adata.X.copy()
        if hasattr(X, "toarray"):
            X = X.toarray()
        X = StandardScaler().fit_transform(X)
        X = X * weights[name]
        X_list.append(X)

    return np.hstack(X_list)


# Example: up-weight imaging relative to RNA when both present
weights = {name: 1.0 for name in modalities}
for img_key in ("imaging", "CellPaint"):
    if img_key in weights:
        weights[img_key] = 2.0

X_weighted = weighted_concatenation(modalities, weights)
print(f"Weighted features shape: {X_weighted.shape}")
```

---

## 3. Cross-Modality Correlation Analysis

Compute pairwise Spearman (or Pearson) correlations between features in two modalities. Most useful for matched gene/protein panels where the same target appears in both.

```python
from scipy.stats import spearmanr, pearsonr


def cross_modality_correlation(adata1, adata2, name1="Modality1", name2="Modality2",
                                method="spearman"):
    """Compute all pairwise feature correlations across two modalities."""
    X1 = adata1.X.toarray() if hasattr(adata1.X, "toarray") else adata1.X
    X2 = adata2.X.toarray() if hasattr(adata2.X, "toarray") else adata2.X

    stat_fn = spearmanr if method == "spearman" else pearsonr
    correlations = []

    for i, feat1 in enumerate(adata1.var_names):
        for j, feat2 in enumerate(adata2.var_names):
            corr, pval = stat_fn(X1[:, i], X2[:, j])
            correlations.append({
                f"{name1}_feature": feat1,
                f"{name2}_feature": feat2,
                "correlation": corr,
                "p_value": pval,
            })

    corr_df = pd.DataFrame(correlations)
    corr_df["abs_corr"] = corr_df["correlation"].abs()
    return corr_df.sort_values("abs_corr", ascending=False)


# Correlate matching gene names between transcript and protein
if protein_adata is not None:
    common_genes = sorted(set(transcript_adata.var_names) & set(protein_adata.var_names))
    print(f"Common genes: {len(common_genes)}")

    if common_genes:
        correlations = []
        for gene in common_genes:
            t_expr = transcript_adata[:, gene].X
            p_expr = protein_adata[:, gene].X

            t_expr = t_expr.toarray().flatten() if hasattr(t_expr, "toarray") else t_expr.flatten()
            p_expr = p_expr.toarray().flatten() if hasattr(p_expr, "toarray") else p_expr.flatten()

            corr, pval = spearmanr(t_expr, p_expr)
            correlations.append({"gene": gene, "transcript_protein_correlation": corr, "p_value": pval})

        corr_df = pd.DataFrame(correlations).sort_values("transcript_protein_correlation", ascending=False)
        print("\nTop transcript-protein correlations:")
        print(corr_df.head(10))

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hist(corr_df["transcript_protein_correlation"], bins=30, edgecolor="black")
        ax.axvline(0, color="red", linestyle="--")
        ax.set_xlabel("Transcript-Protein Correlation (Spearman)")
        ax.set_ylabel("Count")
        ax.set_title("Distribution of Transcript-Protein Correlations")
        plt.show()
```

---

## 4. Procrustes Alignment — Compare Embeddings

Run PCA per modality, then use Procrustes disparity to quantify how similar the resulting embedding geometries are. Lower disparity means the modalities resolve similar cell structure.

```python
from scipy.spatial import procrustes


def compare_modality_embeddings(modalities, n_pcs=30):
    """Run PCA on each modality and compare via Procrustes analysis."""
    embeddings = {}

    for name, adata in modalities.items():
        adata_copy = adata.copy()
        sc.pp.scale(adata_copy, max_value=10)
        sc.pp.pca(adata_copy, n_comps=n_pcs)
        embeddings[name] = adata_copy.obsm["X_pca"]

    comparisons = []
    names = list(embeddings.keys())

    for i, name1 in enumerate(names):
        for name2 in names[i + 1:]:
            _, _, disparity = procrustes(
                embeddings[name1][:, :10],
                embeddings[name2][:, :10],
            )
            comparisons.append({
                "modality_1": name1,
                "modality_2": name2,
                "procrustes_disparity": disparity,
            })

    comparison_df = pd.DataFrame(comparisons)
    print("Embedding similarity (lower disparity = more similar):")
    print(comparison_df)

    return embeddings, comparison_df


embeddings, comparison = compare_modality_embeddings(modalities)
```

---

## 5. Canonical Correlation Analysis (CCA)

Find linear combinations of features across two modalities that are maximally correlated. Useful for identifying shared axes of variation without requiring matched features.

```python
from sklearn.cross_decomposition import CCA
from sklearn.preprocessing import StandardScaler


def canonical_correlation_analysis(adata1, adata2, n_components=10):
    """Find maximally correlated linear combinations across two modalities."""
    X1 = adata1.X.toarray() if hasattr(adata1.X, "toarray") else adata1.X.copy()
    X2 = adata2.X.toarray() if hasattr(adata2.X, "toarray") else adata2.X.copy()

    X1 = StandardScaler().fit_transform(X1)
    X2 = StandardScaler().fit_transform(X2)

    n_comp = min(n_components, X1.shape[1], X2.shape[1])
    cca = CCA(n_components=n_comp)
    X1_c, X2_c = cca.fit_transform(X1, X2)

    canonical_corrs = [
        np.corrcoef(X1_c[:, i], X2_c[:, i])[0, 1]
        for i in range(X1_c.shape[1])
    ]

    print("Canonical correlations:")
    for i, corr in enumerate(canonical_corrs):
        print(f"  Component {i + 1}: {corr:.3f}")

    return cca, X1_c, X2_c, canonical_corrs


if protein_adata is not None:
    cca, t_cca, p_cca, corrs = canonical_correlation_analysis(
        transcript_adata, protein_adata, n_components=10
    )

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for i, ax in enumerate(axes):
        ax.scatter(t_cca[:, i], p_cca[:, i], alpha=0.3, s=5)
        ax.set_xlabel(f"Transcript CCA{i + 1}")
        ax.set_ylabel(f"Protein CCA{i + 1}")
        ax.set_title(f"CCA Component {i + 1} (r={corrs[i]:.2f})")
    plt.tight_layout()
    plt.show()
```

> **Note:** CCA can overfit when `n_features >> n_cells`. Add regularization (`CCA(max_iter=1000)`) or reduce dimensionality with PCA first when dealing with full transcript panels.

---

## 6. Factor Analysis for Shared Variation

Fit a single factor model on the concatenated matrix to decompose variance into shared vs modality-specific factors.

```python
from sklearn.decomposition import FactorAnalysis


def multi_modality_factor_analysis(combined_adata, n_factors=20):
    """Find latent factors and attribute their loadings by modality."""
    X = combined_adata.X.toarray() if hasattr(combined_adata.X, "toarray") else combined_adata.X
    X = StandardScaler().fit_transform(X)

    fa = FactorAnalysis(n_components=n_factors, random_state=42)
    factors = fa.fit_transform(X)

    loadings = pd.DataFrame(
        fa.components_.T,
        index=combined_adata.var_names,
        columns=[f"Factor_{i + 1}" for i in range(n_factors)],
    )
    loadings["modality"] = combined_adata.var["modality"]

    # Variance attributed to each modality per factor
    records = []
    for factor in loadings.columns[:-1]:
        for modality in loadings["modality"].unique():
            mask = loadings["modality"] == modality
            records.append({
                "factor": factor,
                "modality": modality,
                "loading_variance": (loadings.loc[mask, factor] ** 2).sum(),
            })

    var_pivot = (
        pd.DataFrame(records)
        .pivot(index="factor", columns="modality", values="loading_variance")
    )
    var_pivot = var_pivot.div(var_pivot.sum(axis=1), axis=0)

    fig, ax = plt.subplots(figsize=(12, 6))
    var_pivot.iloc[:10].plot(kind="bar", stacked=True, ax=ax)
    ax.set_ylabel("Proportion of Loading Variance")
    ax.set_title("Factor Loading Variance by Modality")
    ax.legend(title="Modality", bbox_to_anchor=(1.05, 1))
    plt.tight_layout()
    plt.show()

    return fa, factors, loadings


if len(modalities) > 1:
    fa, factors, loadings = multi_modality_factor_analysis(combined, n_factors=20)
    combined.obsm["X_factors"] = factors
```

---

## 7. Identify Modality-Specific Signals

Cluster each modality independently and measure cluster agreement (Adjusted Rand Index). Low ARI between modalities means each captures distinct biology.

```python
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score


def identify_modality_specific_clusters(modalities, n_clusters=5):
    """Compare clusters found in each modality to reveal modality-specific signals."""
    cluster_labels = {}

    for name, adata in modalities.items():
        adata_copy = adata.copy()
        sc.pp.scale(adata_copy, max_value=10)
        sc.pp.pca(adata_copy, n_comps=30)

        kmeans = KMeans(n_clusters=n_clusters, random_state=42)
        cluster_labels[name] = kmeans.fit_predict(adata_copy.obsm["X_pca"][:, :20])

    print("Cluster agreement (Adjusted Rand Index):")
    names = list(cluster_labels.keys())
    for i, name1 in enumerate(names):
        for name2 in names[i + 1:]:
            ari = adjusted_rand_score(cluster_labels[name1], cluster_labels[name2])
            print(f"  {name1} vs {name2}: {ari:.3f}")

    return cluster_labels


if len(modalities) > 1:
    cluster_labels = identify_modality_specific_clusters(modalities, n_clusters=5)
```

---

## 8. Joint Leiden Clustering

Run Leiden on the joint UMAP neighborhood graph, then compare with per-modality Leiden labels transferred onto the same embedding.

```python
sc.tl.leiden(combined, resolution=0.5, key_added="leiden_combined",
             flavor="igraph", n_iterations=2, directed=False)

for name, adata in modalities.items():
    adata_copy = adata.copy()
    sc.pp.scale(adata_copy, max_value=10)
    sc.pp.pca(adata_copy, n_comps=30)
    sc.pp.neighbors(adata_copy, n_neighbors=15, n_pcs=20)
    sc.tl.leiden(adata_copy, resolution=0.5, key_added=f"leiden_{name}",
                 flavor="igraph", n_iterations=2, directed=False)
    combined.obs[f"leiden_{name}"] = adata_copy.obs[f"leiden_{name}"].values

# Visualize combined vs per-modality clusters on the same UMAP
n_panels = len(modalities) + 1
fig, axes = plt.subplots(1, n_panels, figsize=(5 * n_panels, 5))

sc.pl.umap(combined, color="leiden_combined", ax=axes[0], show=False,
           title="Combined Clustering")

for i, name in enumerate(modalities.keys()):
    sc.pl.umap(combined, color=f"leiden_{name}", ax=axes[i + 1], show=False,
               title=f"{name} Clustering")

plt.tight_layout()
plt.show()
```

---

## 9. Generate Integration Report

```python
import os


def generate_integration_report(modalities, combined):
    """Write a Markdown summary of the multimodal integration run."""
    lines = ["# Multimodal Integration Report\n\n"]

    lines.append("## Modalities\n\n")
    for name, adata in modalities.items():
        lines.append(f"- **{name}**: {adata.n_vars} features\n")
    lines.append(f"\n**Cells analyzed**: {combined.n_obs}\n")

    lines.append("\n## Combined Analysis\n\n")
    lines.append(f"- Total features: {combined.n_vars}\n")
    lines.append(f"- PCA components: {combined.obsm['X_pca'].shape[1]}\n")

    if "pca" in combined.uns:
        var_explained = combined.uns["pca"]["variance_ratio"].sum()
        lines.append(f"- Variance explained (top PCs): {var_explained:.1%}\n")

    if "leiden_combined" in combined.obs.columns:
        n_clusters = combined.obs["leiden_combined"].nunique()
        lines.append(f"- Combined clusters: {n_clusters}\n")

    report_text = "".join(lines)
    os.makedirs("results", exist_ok=True)

    with open("results/integration_report.md", "w") as f:
        f.write(report_text)

    print(report_text)
    return report_text


generate_integration_report(modalities, combined)
```

---

## Expected Outputs

| Output | Description |
|--------|-------------|
| **Combined embedding** | Joint PCA/UMAP across all modalities |
| **Cross-modality correlations** | Feature-level transcript ↔ protein relationships |
| **Factor loadings** | Proportion of each factor explained by each modality |
| **Cluster comparisons** | ARI scores and UMAP panels per modality |
| **Integration report** | `results/integration_report.md` summary |

---

## Common Pitfalls

1. **Scale mismatch** — Always scale each modality independently to unit variance before any concatenation step. Skipping this hands control of PCA axes to whichever modality has the largest numeric range.

2. **Modality dominance by feature count** — A 20,000-gene transcript panel will dwarf a 50-plex protein panel in concatenated PCA even after scaling. Use `weighted_concatenation` with weights set proportional to explained variance (e.g., top-PC variance ratio per modality) rather than arbitrary values.

3. **CCA overfitting** — CCA maximises correlation and will find spurious structure when `n_features ≫ n_cells`. Reduce dimensionality with PCA first (keep ~50 PCs per modality) before fitting CCA, or use regularized CCA.

4. **Correlation ≠ causality** — High transcript-protein Spearman correlation for a gene does not imply a causal or direct regulatory link; confounders (cell type, cell cycle, batch) often drive shared variance.

5. **Neighbors computed on wrong space** — `sc.tl.leiden` uses the neighbor graph in `combined.obsp["connectivities"]`. Make sure `sc.pp.neighbors` was called on `combined` after PCA, not on a modality-specific copy.
