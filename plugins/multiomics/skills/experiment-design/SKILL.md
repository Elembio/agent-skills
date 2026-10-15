---
name: experiment-design
description: Interpreting the experimental layout of an AVITI24 run — inferring conditions, replicates, and groupings from WellLabel (the sample-sheet name on each Well), and reporting whether condition is confounded with Well. Use when the user asks what conditions, cell lines, treatments, doses, timepoints, or replicates are in a run, what the wells represent, or before any analysis that requires grouping cells by condition (differential expression, pseudobulk, normalization, imaging feature scaling).
metadata:
  version: 0.3.2
  author: elembio
---

# Experiment Design

## Purpose

Interpret the experimental layout of an AVITI24 run by reading well-level labels and inferring conditions, replicates, and groupings. The output is a set of `obs` columns and recommended `groupby` / `reference_group` values that downstream skills (DE, pseudobulk, normalization) can use directly.

**This skill always verifies its interpretation with the user before proceeding.** Well labels often encode multi-factor designs in free-text strings — the agent must not silently assume a parsing is correct.

## Prerequisites

- [metadata-and-target-mapping](../metadata-and-target-mapping/SKILL.md) — `run.wells` when no cell table is loaded yet
- [spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md) — AnnData with `obs["Well"]` and `obs["WellLabel"]`

## Key Concepts

### `Well` vs `WellLabel`

| Column | Means | Where |
|---|---|---|
| `Well` | Plate position (`A1`) | `adata.obs` |
| `WellLabel` | Sample-sheet string (`Hela`, `HelaDMSO`) | `adata.obs` |
| `well_label` | Same strings as `WellLabel` | `run.wells` / `run.tiles` only |

On cells, color and groupby use **`WellLabel`** or **`condition`**. There is no `obs["well_label"]`.

### The inference → verify → commit loop

**Always follow this sequence:**

1. **Infer** — parse `WellLabel` strings, detect structure, propose condition/replicate columns
2. **Show** — print a summary table mapping each well to its inferred labels; ask the user to confirm or correct
3. **Commit** — only after explicit user confirmation, write columns to `adata.obs` and state the recommended `groupby` / `reference_group`

Never skip step 2. A wrong groupby in DE or pseudobulk produces silent garbage results.

### Source of truth: `WellLabel`, not `cell_type`

**When a modality table is loaded:** parse `adata.obs[["Well", "WellLabel"]].drop_duplicates()`. After confirmation, map derived columns onto cells with `adata.obs["Well"]`.

**Metadata-only fallback** (no AnnData yet): `run.wells["well_label"]` — same strings, different column name. Do not carry `well_label` onto `obs` or into `sc.pl.umap(color=...)`.

The `cell_type` column in `run.wells` may name the segmentation model selected for the run rather than the cell line in the well. Derive biological identity from `WellLabel`.

### What `WellLabel` typically encodes

`WellLabel` is a free-text field set at run configuration time. It can encode almost anything — the agent must not assume a specific design without asking. Common patterns:

| Pattern | Example | Possible meaning |
|---|---|---|
| Drug vs vehicle | `"Drug"`, `"DMSO"` | Treatment / control |
| Cell line | `"A375"`, `"HeLa"`, `"MCF7"` | Cell line comparison |
| Dose response | `"Drug_1nM"`, `"Drug_10nM"`, `"Drug_100nM"` | Dose series |
| Time course | `"Drug_2h"`, `"Drug_24h"` | Temporal variable |
| Multi-factor | `"Drug_10uM_24h"` | Condition + dose + timepoint |
| Patient / donor | `"Patient_001"`, `"Donor_A"` | Patient cohort or organoid |
| Co-culture | `"Mono_A375"`, `"CoC_A375_TIL"` | Mono vs co-culture |
| Control keywords | `"DMSO"`, `"Ctrl"`, `"Untreated"`, `"Vehicle"`, `"Mock"`, `"Neg"` | Vehicle / negative control |
| Replicate suffix | `"Drug_rep1"`, `"Drug_R2"`, `"A375_1"` | Replicate of same condition |

**Do not assume** — always show the parsed table to the user and ask. The same label pattern can mean different things in different experiments.

### Confounding between condition and Well

Every AVITI24 well is a separate physical unit — separately stained, separately imaged. So `Well` always carries some technical signal, and **whether `condition` is constant within a well determines what downstream analyses are allowed to do**.

Report this alongside the design, always:

| Quantity | Meaning |
|---|---|
| Conditions per well | `1` for every well ⇒ **condition is confounded with Well** |
| Wells per condition | `1` ⇒ no way to distinguish the condition from that one well's technique |
| Control wells | Which wells (if any) carry a shared reference — DMSO / vehicle / untreated / common line |

**State the consequences before the downstream skill runs, not after:**

- **Imaging scaling** — when condition is constant within wells, per-well MAD-z removes it. Tell [imaging-pipeline](../imaging-pipeline/SKILL.md) to scale against **control wells** if they exist, else **globally**; per-well becomes a sensitivity check only. Replicate wells do not change this.
- **Batch correction** — correcting on `Well` removes the condition too. Flag to [batch-effect-assessment](../batch-effect-assessment/SKILL.md).
- **Attribution** — with one well per condition, biology and well technique are aliased; no normalization separates them. Say "confounded, cannot attribute" rather than reporting a phenotype.

### What to produce

After confirmation, add to `adata.obs`:

| Column | Content |
|---|---|
| `condition` | Primary experimental grouping (e.g. `"Drug"`, `"DMSO"`) |
| `replicate` | Replicate index within condition (e.g. `"1"`, `"2"`) — `None` if not detected |
| `group` | `condition` + `replicate` combined (e.g. `"Drug_rep1"`) — useful as a pseudobulk key |

Then state clearly:
- Recommended `groupby` for DE / pseudobulk — `condition` or `group`, never `well_label`
- Recommended `reference_group` (the control condition)
- Number of replicates per condition
- **Whether `condition` is confounded with `Well`, the control wells (if any), and the resulting scaling / correction recommendation**

## Quick Start

```python
# Well = plate position (A1). WellLabel = sample-sheet string.
# Prefer the cell table when it is already loaded.
wells_df = (
    adata.obs[["Well", "WellLabel"]]
    .astype(str)
    .drop_duplicates()
    .set_index("Well")
    .sort_index()
)
print(wells_df.to_string())

# Color / groupby on cells: WellLabel or condition — never well_label
# See IMPLEMENTATION.md for inference, verification, and writing condition
```

## Detail Files

- [IMPLEMENTATION.md](./IMPLEMENTATION.md) — WellLabel parsing, multi-factor detection, replicate inference, verification prompt template, obs column writing, groupby recommendation

## Related Skills

- [metadata-and-target-mapping](../metadata-and-target-mapping/SKILL.md) — `run.wells` fallback when no table is loaded
- [pseudobulk-differential-expression](../pseudobulk-differential-expression/SKILL.md) — consumes groupby / reference_group
- [normalization-strategies](../normalization-strategies/SKILL.md) — may need condition column for grouped normalization
- [imaging-pipeline](../imaging-pipeline/SKILL.md) — consumes the confounding verdict + control wells to pick the MAD-z reference
- [batch-effect-assessment](../batch-effect-assessment/SKILL.md) — uses condition to separate biological from technical variation
- [embedding-visualization](../embedding-visualization/SKILL.md) — color `Well` (position) and `WellLabel` / `condition` (sample name)

## Disclaimer

This content is provided "as-is" without warranty of any kind. AI-generated output may contain
inaccuracies. Review and verify all outputs before relying on them.
