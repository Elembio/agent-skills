---
name: sandbox-explore
description: Explore what this mounted ElemBio run can answer — inspect available_tables, suggest only questions the present modalities support, then hand off to the matching multiomics skill. Use when the user asks "what else can I do with this run", "what can the sandbox do for me", "suggest analyses", "explore this run", "I don't know where to start", "give me sandbox examples", or "what should I ask next".
metadata:
  version: 0.1.0
  author: elembio
---

# Sandbox explore

Help the user go past the first listing. Inspect what this run actually contains, then offer questions it can answer. Demo proves the sandbox works. Explore turns that into a real starting set for this dataset.

This is the natural follow-on to `sandbox-demo`.

Do not write pipeline steps here. The Element Biosciences `multiomics` `index` skill owns order and specialists. This skill only picks questions and names the skill that answers them.

## When to use vs. other skills

- **sandbox-explore** (this skill) — modality-tailored questions for the mounted run.
- **sandbox-demo** — one run, `available_tables`, no analysis. Run this first if the user has not seen the sandbox work.
- **sandbox-onboard** — connect + companion plugin. Run this if tools or multiomics are missing.
- **sandbox-status** — health checks and audits.

If the user has not listed tables yet, route to **sandbox-demo** first. Explore works best after a win.

## Step 1: Inspect, do not interview

Zapier interviews because it cannot see the user's apps. We can see the run.

Require a live `Loader` (or create one sandbox, mount, and construct it — follow the server-published Quickstart and `MOUNTING.md`, and stay under the demo cost ceiling until the user picks a question). Then read `available_tables`.

Map table names to modalities. Do not invent a modality that is not present.

| Table signal | Treat as |
| --- | --- |
| `ThreePrimeUntargeted` | 3′ untargeted |
| `SpecializedTargeted` that is a true OPS identity library | OPS |
| Imaging / CellPaint / ProteinIF / morphology tables | Imaging |
| Legacy Transcript / Protein barcoding panels | Legacy targeted |
| Two or more of the above | Multimodal — also offer the multimodal questions |

If `available_tables` is empty or `Loader` failed, stop and diagnose via **sandbox-status**, not by guessing.

You may ask **one** follow-up only when the tables leave a real fork (which SDO, which well, which condition). Do not run a role survey.

## Step 2: Suggest questions this run supports

Pick 4–6 questions from the library below that match **present** modalities. Always include one cross-cutting question (run QC or experiment design) when those skills apply.

For each question, output three lines — not a table:

> **[What the question answers, in plain language]**
>
> Say to me: *"[exact prompt]"*. I'll [what comes back, one line].
>
> *Skill:* [multiomics skill name]

Do not list more than 6 at once. If you have more, give 4, then offer "want more?"

## Step 3: Hand off

When the user picks a question, open that skill and follow it. Start at `index` if you are unsure of order (run QC before cell QC, and so on). `index` encodes AVITI24 / DISS conventions that generic Scanpy defaults get wrong.

Do not paste pipeline steps from `index` into this chat. Name the skill and go.

## Question library

Do not read the user every entry. Pick the ones their tables support.

### Cross-cutting (any AVITI24 store)

> **Run health before you touch cells**
>
> Say to me: *"How does this run look — cells per well, tiles, any wells I should drop?"* I'll give well / tile health for the modalities that are present, with no cell filtering yet.
>
> *Skill:* `run-quality-overview`

> **What was actually plated**
>
> Say to me: *"What conditions, replicates, and cell lines are in this run?"* I'll infer design from `WellLabel` and flag when condition is confounded with well.
>
> *Skill:* `experiment-design`

> **Where things live in the object**
>
> Say to me: *"Where are counts, QC metrics, and spatial coordinates in this AnnData?"* I'll point at the canonical `obs` / `var` / `obsm` slots for the modalities we have.
>
> *Skill:* `anndata-structure-reference`

### 3′ untargeted (`ThreePrimeUntargeted`)

> **Cell QC on 3′ depth, not UMIs**
>
> Say to me: *"Filter low-quality cells on this 3′ transcriptome."* I'll run AVITI24 cell QC (`compute_qc_metrics` → `filter_cells`) and talk in depth / polony counts, not droplet-seq UMIs.
>
> *Skill:* `cell-quality-control`

> **Normalize counts, then reduce**
>
> Say to me: *"Normalize this 3′ table and show me a UMAP."* I'll normalize the count modality only, then PCA / UMAP.
>
> *Skill:* `normalization-strategies` then `dimensionality-reduction`

> **End-to-end 3′ pipeline**
>
> Say to me: *"Run the 3′ untargeted pipeline on this store."* I'll follow the specialist order for `ThreePrimeUntargeted`.
>
> *Skill:* `3prime-untargeted-pipeline`

> **Differential expression with replicates**
>
> Say to me: *"Which genes differ between [condition A] and [condition B]?"* I'll do pseudobulk DE if there are ≥2 biological replicates per condition.
>
> *Skill:* `pseudobulk-differential-expression`

### OPS (true Optical Pooled Screening identity library)

> **Confirm it is real OPS**
>
> Say to me: *"Is this SpecializedTargeted table a real OPS library, or a housekeeping panel?"* I'll distinguish identity libraries from QC gene panels that also label as OPS.
>
> *Skill:* `ops-pipeline`

> **OPS identity and area QC**
>
> Say to me: *"QC and assign identities for this OPS library."* I'll follow the OPS specialist (area-only cell QC, not 3′ MT filters).
>
> *Skill:* `ops-pipeline`

### Imaging / CellPaint / Protein IF

> **Morphology inventory**
>
> Say to me: *"What's in the imaging — channels, stains, any IF?"* I'll inventory morphology stains and protein IF channels before any scaling.
>
> *Skill:* `imaging-pipeline`

> **Imaging QC and MAD-z features**
>
> Say to me: *"Build the imaging feature matrix and tell me which scaling reference you used."* I'll run S/B, localization, and MAD-z with the design-dependent reference (control wells / global / per-well).
>
> *Skill:* `imaging-pipeline`

> **Well fluorescence and a cell close-up**
>
> Say to me: *"Show me well fluorescence and a representative cell."* I'll plot the plate and a `plot_cell` close-up.
>
> *Skill:* `spatial-multimodal-visualization`

### Legacy targeted (Transcript / Protein barcoding only)

> **Legacy panel path**
>
> Say to me: *"This looks like a barcoding panel — treat it as legacy targeted."* I'll use the legacy specialist, not the current 3′ / OPS / imaging-first pipelines.
>
> *Skill:* `legacy-targeted-pipeline`

### Multimodal (two or more present)

> **Align before you join**
>
> Say to me: *"Prepare the present modalities for joint analysis."* I'll apply the right transform per modality (count normalize vs imaging MAD-z) before concatenating.
>
> *Skill:* `multimodal-data-preparation`

> **Joint embedding**
>
> Say to me: *"Make a joint UMAP across the modalities we have."* I'll follow integrated analysis only after alignment.
>
> *Skill:* `integrated-analysis`

## Handling edge cases

- **Table name you do not recognize:** say you will not invent a modality. Open `spatialdata-loading-and-access` and discover. Never assume Scanpy defaults.
- **Several SDOs:** stop and use `locate-spatialdata-store` to confirm which store the user wants. Do not explore two stores at once.
- **No multiomics plugin:** route to `sandbox-onboard` **Companion missing**. Do not improvise pipelines.
- **User asks for a full analysis immediately:** pick the matching pipeline skill from `index` and go. Skip extra questions.

## Gotchas

- **Do not interview for role.** The tables are the source of truth. One follow-up is the cap.
- **Do not dump the whole library.** 4–6 questions, present modalities only.
- **Do not paste pipeline steps.** Name the skill. `index` owns order.
- **Do not offer 3′ UMI language, `sc.pp.calculate_qc_metrics`, or `read_zarr` on a zip for tables.** Those are the generic defaults `index` exists to prevent.
- **Do not call `load_tables` just to suggest questions.** `available_tables` is enough to choose. Load when the chosen skill needs it, and say the cost first.

## Tone

Concrete, never abstract. Avoid "the sandbox can help you analyze your data." Say instead "Try: 'How does this run look across wells?' and I will pull run QC." Show, do not pitch.
