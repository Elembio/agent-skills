# Cell Cycle Analysis — Reference

## Table of Contents

1. [S-phase gene list](#s-phase-gene-list)
2. [G2M-phase gene list](#g2m-phase-gene-list)
3. [Provenance](#provenance)
4. [Usage notes](#usage-notes)

---

## S-phase gene list

43 genes from Tirosh et al. 2015 (doi:10.1126/science.aad0501). Same set as Seurat `cc.genes.updated.2019` S-phase partition.

```python
S_GENES = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG",
    "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP",
    "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76",
    "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51",
    "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM",
    "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8",
]
```

---

## G2M-phase gene list

54 genes from Tirosh et al. 2015 (doi:10.1126/science.aad0501). Same set as Seurat `cc.genes.updated.2019` G2M partition.

```python
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
```

---

## Provenance

These gene lists originate from:

> Tirosh, I. et al. "Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq." *Science* 352, 189–196 (2016). doi:[10.1126/science.aad0501](https://doi.org/10.1126/science.aad0501)

They are the same lists distributed by Seurat as `cc.genes.updated.2019` and recommended in [scanpy's cell cycle tutorial](https://scanpy.readthedocs.io/en/latest/how-to/cell-cycle.html). The lists are human gene symbols; for other species, map to orthologs before use.

These lists are maintained here as reference data for the skill layer.

---

## Usage notes

**Filter to panel genes before scoring:**

```python
s_genes   = [g for g in S_GENES   if g in adata.var_names]
g2m_genes = [g for g in G2M_GENES if g in adata.var_names]

print(f"S genes present:   {len(s_genes)}/{len(S_GENES)}")
print(f"G2M genes present: {len(g2m_genes)}/{len(G2M_GENES)}")
```

**Species:** These lists are human gene symbols. For mouse, map to mouse orthologs (e.g., via `biomart` or a curated conversion table) before using.

**Panel coverage threshold:** If fewer than ~5 genes per phase are found in the panel, cell cycle scoring will be unreliable. Consider skipping or using a proliferation marker (e.g., MKI67) as a proxy instead.
