# Cell Cycle Analysis — Implementation

## Table of Contents

1. [Scanpy cell cycle scoring](#1-scanpy-cell-cycle-scoring)
2. [Gene coverage check](#2-gene-coverage-check)
3. [UMAP visualization](#3-umap-visualization)
4. [Phase distribution plots](#4-phase-distribution-plots)
5. [Phase distribution by well / condition](#5-phase-distribution-by-well--condition)
6. [Regression to remove cell cycle signal](#6-regression-to-remove-cell-cycle-signal)
7. [Proliferation index scoring](#7-proliferation-index-scoring)
8. [Export results](#8-export-results)
10. [Expected Outputs](#expected-outputs)
11. [Common Pitfalls](#common-pitfalls)

---

## 1. Scanpy cell cycle scoring

Use the Tirosh et al. 2015 S/G2M gene lists (doi:10.1126/science.aad0501), filtered to genes present in your panel.

```python
import scanpy as sc
from elembio_spatialdata_analysis.load import Loader

adata = Loader("/path/to/run.zarr").load_tables(table="ThreePrimeUntargeted", lazy=False)
# Normalize before scoring — see normalization-strategies skill

# Tirosh et al. 2015 S/G2M gene lists (doi:10.1126/science.aad0501)
# Full lists in REFERENCE.md
S_GENES = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG",
    "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP",
    "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76",
    "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51",
    "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM",
    "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8",
]
G2M_GENES = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A",
    "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF",
    "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB",
    "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP",
    "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1",
    "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR",
    "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF",
    "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA",
]

s_genes   = [g for g in S_GENES   if g in adata.var_names]
g2m_genes = [g for g in G2M_GENES if g in adata.var_names]

sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

# adata.obs gains:
#   S_score   — continuous S-phase score
#   G2M_score — continuous G2/M-phase score
#   phase     — categorical: 'G1', 'S', or 'G2M'
print(adata.obs[["S_score", "G2M_score", "phase"]].head())
print(adata.obs["phase"].value_counts())
```

---

## 2. Gene coverage check

Targeted panels may only cover a subset of the standard lists. Check coverage before trusting scores.

```python
s_genes   = [g for g in S_GENES   if g in adata.var_names]
g2m_genes = [g for g in G2M_GENES if g in adata.var_names]

print(f"S genes found:   {len(s_genes)}/{len(S_GENES)}")
print(f"G2M genes found: {len(g2m_genes)}/{len(G2M_GENES)}")

# Check which genes are missing
s_missing   = [g for g in S_GENES   if g not in adata.var_names]
g2m_missing = [g for g in G2M_GENES if g not in adata.var_names]
print(f"Missing S genes:   {s_missing[:10]}")
print(f"Missing G2M genes: {g2m_missing[:10]}")
```

If fewer than ~5 genes per phase are found, consider using a single proliferation marker (e.g., MKI67) as a proxy instead (see section 8).

---

## 3. UMAP visualization

```python
import matplotlib.pyplot as plt
import scanpy as sc

if "X_umap" not in adata.obsm:
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata)

# Scanpy output
fig, axes = plt.subplots(1, 3, figsize=(15, 4))
sc.pl.umap(adata, color="phase",     ax=axes[0], show=False, title="Cell Cycle Phase")
sc.pl.umap(adata, color="S_score",   ax=axes[1], show=False, title="S Phase Score",
           cmap="RdBu_r", vcenter=0)
sc.pl.umap(adata, color="G2M_score", ax=axes[2], show=False, title="G2M Phase Score",
           cmap="RdBu_r", vcenter=0)
plt.tight_layout()
plt.savefig("cell_cycle_umap.png", dpi=150, bbox_inches="tight")
plt.close()
```

---

## 4. Phase distribution plots

```python
import matplotlib.pyplot as plt

phase_colors = {"G1": "#4CAF50", "S": "#2196F3", "G2M": "#FF9800"}

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

# Overall bar chart
ax = axes[0]
phase_counts = adata.obs["phase"].value_counts()
phase_counts.plot(kind="bar", ax=ax,
                  color=[phase_colors.get(p, "gray") for p in phase_counts.index])
ax.set_xlabel("Cell Cycle Phase")
ax.set_ylabel("Number of Cells")
ax.set_title("Cell Cycle Phase Distribution")
ax.tick_params(axis="x", rotation=0)

# Stacked bar by cluster
ax = axes[1]
phase_by_cluster = adata.obs.groupby(["leiden", "phase"]).size().unstack(fill_value=0)
phase_by_cluster_pct = phase_by_cluster.div(phase_by_cluster.sum(axis=1), axis=0)
phase_by_cluster_pct.plot(kind="bar", stacked=True, ax=ax,
                          color=[phase_colors.get(p, "gray") for p in phase_by_cluster_pct.columns])
ax.set_xlabel("Cluster")
ax.set_ylabel("Proportion")
ax.set_title("Cell Cycle by Cluster")
ax.legend(title="Phase", bbox_to_anchor=(1.05, 1))

plt.tight_layout()
plt.savefig("cell_cycle_distribution.png", dpi=150, bbox_inches="tight")
plt.close()
```

---

## 5. Phase distribution by well / condition

```python
from scipy.stats import chi2_contingency
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(10, 5))
phase_by_well = adata.obs.groupby(["Well", "phase"]).size().unstack(fill_value=0)
phase_by_well_pct = phase_by_well.div(phase_by_well.sum(axis=1), axis=0)
phase_by_well_pct.plot(kind="bar", stacked=True, ax=ax)
ax.set_xlabel("Well")
ax.set_ylabel("Proportion")
ax.set_title("Cell Cycle Phase by Well")
ax.legend(title="Phase", bbox_to_anchor=(1.05, 1))
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("cell_cycle_by_well.png", dpi=150, bbox_inches="tight")
plt.close()

chi2, p, dof, _ = chi2_contingency(phase_by_well)
print(f"Chi-square test for phase differences across wells: p = {p:.4e}")
```

---

## 6. Regression to remove cell cycle signal

> **Warning:** This permanently modifies `adata.X`. Only do this if you do **not** want to study cell cycle biology (see SKILL.md regression caveat).

```python
import scanpy as sc

# Preserve original counts
adata.layers["pre_cc_regression"] = adata.X.copy()

sc.pp.regress_out(adata, ["S_score", "G2M_score"])
sc.pp.scale(adata, max_value=10)

# Recompute embeddings
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
```

Compare before/after:

```python
import matplotlib.pyplot as plt

adata_before = adata.copy()
adata_before.X = adata_before.layers["pre_cc_regression"]
sc.pp.pca(adata_before)
sc.pp.neighbors(adata_before)
sc.tl.umap(adata_before)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
sc.pl.umap(adata_before, color="phase", ax=axes[0], show=False, title="Before CC Regression")
sc.pl.umap(adata,        color="phase", ax=axes[1], show=False, title="After CC Regression")
plt.tight_layout()
plt.savefig("cell_cycle_regression_comparison.png", dpi=150, bbox_inches="tight")
plt.close()
```

Alternatively, include phase as a covariate in differential expression:

```python
design = "~ experimental_condition + phase"
# See pseudobulk-differential-expression skill for the full DESeq2 workflow
```

---

## 7. Proliferation index scoring

```python
import scanpy as sc

# Combined S + G2M score as a continuous proliferation index
adata.obs["proliferation_score"] = adata.obs["S_score"] + adata.obs["G2M_score"]

threshold = adata.obs["proliferation_score"].quantile(0.75)
adata.obs["is_proliferating"] = adata.obs["proliferation_score"] > threshold

n = adata.obs["is_proliferating"].sum()
pct = adata.obs["is_proliferating"].mean() * 100
print(f"Proliferating cells: {n} ({pct:.1f}%)")

sc.pl.umap(adata, color=["proliferation_score", "is_proliferating"],
           save="_proliferation.png")
```

---

## 8. Export results

```python
# Save AnnData with all cell cycle columns
adata.write_h5ad("adata_with_cell_cycle.h5ad")

# Export phase table
adata.obs[["S_score", "G2M_score", "phase", "proliferation_score"]].to_csv(
    "cell_cycle_scores.csv"
)
```

---

## Expected Outputs

After `sc.tl.score_genes_cell_cycle`, `adata.obs` will contain:

| Column | Type | Values |
|---|---|---|
| `S_score` | float | continuous S-phase score |
| `G2M_score` | float | continuous G2/M score |
| `phase` | category | `G1`, `S`, `G2M` |
| `proliferation_score` | float | `S_score + G2M_score` (optional) |
| `is_proliferating` | bool | top-quartile flag (optional) |

Validation:

```python
assert "phase" in adata.obs.columns
assert set(adata.obs["phase"].unique()) <= {"G1", "S", "G2M"}
print(adata.obs["phase"].value_counts())
```

---

## Common Pitfalls

1. **Missing cell cycle genes** — Targeted panels may have very few S or G2M genes. Check coverage (section 2) before trusting scores. If < 5 genes per phase are found, results will be noisy.

2. **3′ untargeted vs. targeted panels** — Untargeted data generally has better gene coverage. For small targeted panels, verify cell cycle genes were included in the panel design.

3. **Over-regression** — Regressing out cell cycle removes real biology when cells of interest are proliferating (e.g., tumor, stem, immune activation). Prefer covariate inclusion in DE models.

4. **Phase imbalance** — Distribution > 95% G1 usually means few proliferating cells or mismatched gene list. Check coverage and normalization.

5. **Species differences** — `S_GENES` / `G2M_GENES` are human gene symbols. For mouse or other organisms, map to appropriate orthologs before scoring.
