---
name: cell-cycle-analysis
description: Cell cycle phase scoring and visualization for single-cell expression data. Uses Scanpy's score_genes_cell_cycle with Tirosh et al. 2015 S/G2M gene lists on any normalized AnnData (ThreePrimeUntargeted via Loader, scRNA-seq-like, or AnnData from CSV/parquet count matrix). Regressing out cell cycle effects for downstream analysis. Use when the user asks about cell cycle, proliferation, S phase, G2M, or cell state.
metadata:
  version: 0.3.3
  author: elembio
---

# Cell Cycle Analysis

## Purpose

Cell cycle phase scoring assigns G1, S, and G2M phases to individual cells based on gene expression signatures. Uses `sc.tl.score_genes_cell_cycle` with the Tirosh et al. 2015 S/G2M gene lists (doi:10.1126/science.aad0501, same as Seurat `cc.genes.updated.2019`). Produces `S_score`, `G2M_score`, and `phase` columns. Gene lists are defined inline in skill code — they are standard published reference data, not maintained in the analysis package.

## Prerequisites

- [cell-quality-control](../cell-quality-control/SKILL.md)
- [normalization-strategies](../normalization-strategies/SKILL.md)
- A normalized AnnData object with transcript data (AVITI24: `Loader.load_tables(table="ThreePrimeUntargeted", ...)` (legacy: `Transcript`); generic: any AnnData built from a CSV/parquet count matrix whose `var_names` are gene symbols)

## Key Concepts

### Approach

`sc.tl.score_genes_cell_cycle` with the standard Tirosh et al. 2015 S/G2M gene lists, filtered to genes present in the panel.

| Function | Columns written | Source |
|---|---|---|
| `sc.tl.score_genes_cell_cycle` | `S_score`, `G2M_score`, `phase` (G1/S/G2M) | Tirosh et al. 2015 (doi:10.1126/science.aad0501) |

### Regression Caveat

Regressing out cell cycle (`sc.pp.regress_out`) **removes** cell cycle signal from the expression matrix. Only do this if you want to study non-cell-cycle biology. If proliferating cells are biologically relevant (e.g., tumor cells), do not regress.

## Quick Start

```python
import scanpy as sc
from elembio_spatialdata_analysis.load import Loader

adata = Loader("/path/to/run.zarr").load_tables(table="ThreePrimeUntargeted", lazy=False)
# normalize first — see normalization-strategies skill

# Tirosh et al. 2015 S/G2M gene lists — full lists in REFERENCE.md
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

s_genes  = [g for g in S_GENES  if g in adata.var_names]
g2m_genes = [g for g in G2M_GENES if g in adata.var_names]

sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

# Writes to adata.obs: S_score, G2M_score, phase ('G1' | 'S' | 'G2M')
print(adata.obs[["S_score", "G2M_score", "phase"]].head())
print(adata.obs["phase"].value_counts())
```

## Detail Files

- [IMPLEMENTATION.md](./IMPLEMENTATION.md) — Full implementation: gene coverage check, UMAP visualization, regression, proliferation index, export, expected outputs, and common pitfalls.
- [REFERENCE.md](./REFERENCE.md) — Full S/G2M gene lists as copy-pasteable Python, with provenance citation.

## Related Skills

- [normalization-strategies](../normalization-strategies/SKILL.md) — Normalize before scoring
- [clustering-and-annotation](../clustering-and-annotation/SKILL.md) — Use phase in cluster annotation
- [pseudobulk-differential-expression](../pseudobulk-differential-expression/SKILL.md) — Include phase as covariate

## Disclaimer

This content is provided "as-is" without warranty of any kind. AI-generated output may contain
inaccuracies. Review and verify all outputs before relying on them.
