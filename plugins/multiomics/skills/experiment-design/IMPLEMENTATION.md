# Experiment Design — Implementation

## Contents

- [Step 1: Load well labels](#step-1-load-well-labels)
- [Step 2: Infer design structure](#step-2-infer-design-structure)
- [Step 3: Show parsed table and verify with user](#step-3-show-parsed-table-and-verify-with-user)
- [Step 4: Commit — write obs columns](#step-4-commit--write-obs-columns)
- [Step 5: Recommend groupby and reference_group](#step-5-recommend-groupby-and-reference_group)
- [Step 6: Report condition–Well confounding](#step-6-report-conditionwell-confounding)
- [Multi-factor designs](#multi-factor-designs)
- [Edge cases](#edge-cases)
- [Common pitfalls](#common-pitfalls)

---

## Step 1: Load well labels

`Well` is the plate position (`A1`). `WellLabel` is the sample-sheet string. Prefer the cell table. `run.wells["well_label"]` is the same strings under a different name — use it only when no AnnData is loaded. Never copy `well_label` onto `obs` or into plot `color=`.

```python
from elembio_spatialdata_analysis.load import Loader

loader = Loader("/path/to/run.zarr")
run = next(iter(loader._runs.values()))

# Preferred — modality table already loaded
wells_df = (
    adata.obs[["Well", "WellLabel"]]
    .astype(str)
    .drop_duplicates()
    .set_index("Well")
    .sort_index()
)
print(wells_df.to_string())

# Metadata-only fallback (no AnnData yet):
# src = run.wells
# label_col = "well_label" if "well_label" in src.columns else "WellLabel"
# wells_df = src[[label_col]].rename(columns={label_col: "WellLabel"})
# wells_df.index = wells_df.index.astype(str)
# wells_df.index.name = "Well"
```

Example output:
```
     WellLabel
A1        DMSO
A2        DMSO
B1    Drug_1uM
B2    Drug_1uM
C1   Drug_10uM
C2   Drug_10uM
```

---

## Step 2: Infer design structure

Parse `WellLabel` to detect conditions, replicates, and factors. The logic below is heuristic — it must be verified with the user (Step 3) before use.

```python
import re
import pandas as pd

CONTROL_PATTERNS = re.compile(
    r"^(dmso|ctrl|control|untreated|vehicle|mock|neg|negative|parental|wt|wildtype|wild.type)$",
    re.IGNORECASE,
)
REPLICATE_PATTERNS = re.compile(
    r"[_\-]?(rep|r|replicate)[_\-]?(\d+)$", re.IGNORECASE
)
NUMERIC_SUFFIX = re.compile(r"^(.+)[_\-](\d+)$")


def infer_well_design(wells_df: pd.DataFrame) -> pd.DataFrame:
    """Parse WellLabel into condition, replicate, is_control columns.

    `wells_df` is indexed by Well and has a WellLabel column.
    Returns a copy with added columns. Always show to user before use.
    """
    df = wells_df.copy()
    conditions, replicates, is_controls = [], [], []

    for label in df["WellLabel"].fillna("unknown"):
        label = str(label).strip()

        # Check for explicit replicate suffix first
        rep_match = REPLICATE_PATTERNS.search(label)
        if rep_match:
            condition = label[: rep_match.start()].strip("_-")
            replicate = rep_match.group(2)
        else:
            # Check for numeric suffix as implicit replicate
            num_match = NUMERIC_SUFFIX.match(label)
            if num_match:
                condition = num_match.group(1)
                replicate = num_match.group(2)
            else:
                condition = label
                replicate = None

        is_control = bool(CONTROL_PATTERNS.match(condition))
        conditions.append(condition)
        replicates.append(replicate)
        is_controls.append(is_control)

    df["_inferred_condition"] = conditions
    df["_inferred_replicate"] = replicates
    df["_inferred_is_control"] = is_controls
    return df


inferred = infer_well_design(wells_df)
```

---

## Step 3: Show parsed table and verify with user

**Always do this before writing to `adata.obs`.** Print a clear summary and ask the user to confirm or correct.

```python
summary = inferred[["WellLabel", "_inferred_condition", "_inferred_replicate", "_inferred_is_control"]].copy()
summary.columns = ["WellLabel", "condition (inferred)", "replicate (inferred)", "is_control (inferred)"]
print("\n=== Inferred experiment design — please verify ===")
print(summary.to_string())
print()

# Summarise conditions found
unique_conditions = inferred["_inferred_condition"].unique().tolist()
control_conditions = inferred.loc[inferred["_inferred_is_control"], "_inferred_condition"].unique().tolist()
n_replicates = (
    inferred.groupby("_inferred_condition")["_inferred_replicate"]
    .nunique()
    .to_dict()
)

print(f"Conditions detected:  {unique_conditions}")
print(f"Control(s) detected:  {control_conditions or 'none detected — please specify'}")
print(f"Replicates per condition: {n_replicates}")
print()
print("Does this look correct? Please confirm:")
print("  1. Is the experiment type correct? (drug treatment / cell line comparison /")
print("     dose-response / time course / patient cohort / co-culture / other)")
print("  2. Are conditions and controls assigned correctly?")
print("  3. Are replicates identified correctly, or are these distinct conditions?")
print("  4. If labels are multi-factor (e.g. 'Drug_10uM_24h'), what does each part mean?")
print()
print("If anything is wrong, please provide the correct mapping or label convention.")
```

**Wait for user confirmation before proceeding to Step 4.** If the user corrects the mapping, apply their override:

```python
# Example: user-supplied override
override = {
    "A1": ("DMSO", "1"),
    "A2": ("DMSO", "2"),
    "B1": ("Drug_1uM", "1"),
    "B2": ("Drug_1uM", "2"),
}
# Apply override: dict of well_name -> (condition, replicate)
for well, (cond, rep) in override.items():
    if well in inferred.index:
        inferred.loc[well, "_inferred_condition"] = cond
        inferred.loc[well, "_inferred_replicate"] = rep
```

---

## Step 4: Commit — write obs columns

Only after user confirmation:

```python
adata = loader.load_tables(table="ThreePrimeUntargeted", lazy=False)  # or whichever table is present

well_to_condition = inferred["_inferred_condition"].to_dict()
well_to_replicate = inferred["_inferred_replicate"].to_dict()

adata.obs["condition"]  = adata.obs["Well"].map(well_to_condition)
adata.obs["replicate"]  = adata.obs["Well"].map(well_to_replicate)

# Composite group key — useful as pseudobulk sample identifier
adata.obs["group"] = adata.obs.apply(
    lambda r: f"{r['condition']}_{r['replicate']}" if pd.notna(r["replicate"])
    else r["condition"],
    axis=1,
)

print(f"obs columns added: condition, replicate, group")
print(adata.obs[["Well", "WellLabel", "condition", "replicate", "group"]].drop_duplicates().sort_values("Well").to_string())
# Plot / groupby on cells: WellLabel or condition — never well_label
```

---

## Step 5: Recommend groupby and reference_group

State recommendations explicitly so the user and downstream skills can reference them:

```python
control_conditions = inferred.loc[inferred["_inferred_is_control"], "_inferred_condition"].unique().tolist()
treatment_conditions = [c for c in inferred["_inferred_condition"].unique() if c not in control_conditions]

has_replicates = inferred["_inferred_replicate"].notna().any()

print("=== Recommended downstream settings ===")
if has_replicates:
    print(f"  groupby        = 'group'       # one pseudobulk sample per well")
    print(f"  condition_key  = 'condition'   # for DE contrast")
else:
    print(f"  groupby        = 'condition'   # no replicates detected — pseudobulk not recommended")

if control_conditions:
    print(f"  reference_group = {control_conditions[0]!r}")
else:
    print("  reference_group = ??? (no control detected — please specify)")

print(f"  treatment conditions: {treatment_conditions}")
print(f"  n_replicates per condition: { {c: (inferred['_inferred_condition'] == c).sum() for c in inferred['_inferred_condition'].unique()} }")
```

---

## Step 6: Report condition–Well confounding

Downstream skills need to know whether `condition` varies within a well. Report it with the design, every time — it decides the imaging scaling reference and whether `Well` batch correction is even legal.

`WellLabel` assigns exactly one label per well, so a label-derived `condition` is constant within a well **by construction**. The only escape is a well that physically contains a mixed population — a pooled OPS screen, mixed cell lines, multiplexed perturbations. That is not visible in `WellLabel`, so ask:

```python
print("Does any single well contain more than one cell line / perturbation")
print("(pooled or multiplexed)? If yes, per-well scaling IS the right correction.")
pooled_wells = []  # user-confirmed wells with internal heterogeneity
```

```python
def condition_well_confounding(inferred: pd.DataFrame, pooled_wells: list[str] | None = None) -> dict:
    """Well-level confounding summary. `inferred` is indexed by well name."""
    pooled = set(pooled_wells or [])
    conds = inferred["_inferred_condition"]
    wells_per_cond = conds.value_counts().to_dict()
    return {
        # Label-derived conditions never vary within a well unless the well is pooled
        "confounded": not pooled.issuperset(inferred.index),
        "pooled_wells": sorted(pooled),
        "wells_per_condition": wells_per_cond,
        "min_wells_per_condition": min(wells_per_cond.values()),
        "control_wells": inferred.index[inferred["_inferred_is_control"]].tolist(),
    }


conf = condition_well_confounding(inferred, pooled_wells)

print("=== condition vs Well ===")
print(f"  wells per condition : {conf['wells_per_condition']}")
print(f"  control wells       : {conf['control_wells'] or 'none detected'}")
print(f"  pooled wells        : {conf['pooled_wells'] or 'none'}")

if conf["confounded"]:
    print("  condition is CONSTANT within each well.")
    print("Consequences:")
    print("  - Imaging: do NOT scale per well — it subtracts the condition.")
    if conf["control_wells"]:
        print(f"    Scale against control wells {conf['control_wells']} (run-level reference).")
    else:
        print("    No shared control → scale globally across the run.")
    print("  - Batch correction on 'Well' would remove the condition. Do not apply blindly.")
    if conf["min_wells_per_condition"] < 2:
        print("  - WARNING: some conditions have a single well. Biology and well technique are")
        print("    aliased there — no normalization separates them. Report as confounded.")
else:
    print("  Wells are internally heterogeneous → per-well scaling is the correct correction.")
```

Pass the verdict on: [imaging-pipeline](../imaging-pipeline/SKILL.md) uses it to pick the MAD-z reference; [batch-effect-assessment](../batch-effect-assessment/SKILL.md) uses it to decide whether correction is possible.

---

## Multi-factor designs

When `WellLabel` encodes multiple factors (e.g. `"Drug_10uM_24h"`), split on `_` and ask the user to name each factor:

```python
def parse_multifactor(label: str, factor_names: list[str]) -> dict:
    """Split a label on '_' and map parts to named factors."""
    parts = str(label).strip().split("_")
    return {
        name: parts[i] if i < len(parts) else None
        for i, name in enumerate(factor_names)
    }

# Example: user says labels are "condition_dose_timepoint"
factor_names = ["condition", "dose", "timepoint"]   # confirmed with user

factor_df = pd.DataFrame([
    {"well_name": well, **parse_multifactor(row["WellLabel"], factor_names)}
    for well, row in wells_df.iterrows()
]).set_index("well_name")

print("=== Multi-factor design — please confirm ===")
print(factor_df.to_string())
print("Factor names assumed:", factor_names)
print("Correct? (yes / specify different names)")

# After confirmation, merge into obs:
for col in factor_names:
    adata.obs[col] = adata.obs["Well"].map(factor_df[col].to_dict())
```

---

## Edge cases

### All wells have the same label
```python
if inferred["_inferred_condition"].nunique() == 1:
    print("WARNING: All wells have the same label. This may be:")
    print("  - A single-condition run (no DE possible across conditions)")
    print("  - Well labels not set correctly at run time")
    print("Please confirm with the user what the intended design was.")
```

### Missing or blank well labels
```python
blank_wells = wells_df[wells_df["WellLabel"].fillna("").str.strip() == ""].index.tolist()
if blank_wells:
    print(f"WARNING: {len(blank_wells)} wells have no label: {blank_wells}")
    print("Please provide condition assignments for these wells before proceeding.")
```

### No replicates detected
```python
if not has_replicates:
    print("NOTE: No replicate structure detected.")
    print("Pseudobulk DE requires ≥2 replicates per condition.")
    print("With n=1 per condition, only cell-level DE (not pseudobulk) is possible.")
    print("Is this a pilot run, or were replicates not labeled?")
```

---

## Common pitfalls

1. **`cell_type` may not reflect the biological sample.** `run.wells["cell_type"]` (e.g. `"HeLa"`) may name the segmentation model selected for the run rather than the cell line in the well. Always derive biological identity from `obs["WellLabel"]` (or `run.wells["well_label"]` if no table is loaded).

2. **Coloring or grouping cells by `well_label`** — that column is only on `run.wells`. On AnnData use `obs["WellLabel"]` or `condition`. `Well` is the plate position (`A1`), not the sample-sheet name. Falling back to `Well` after a missing-key error plots 12 seats instead of the treatments.

3. **Silent wrong groupby** — the most dangerous failure mode. A `condition` column that maps all cells to the same group, or swaps treatment/control, produces results that look plausible but are wrong. Always print and confirm.

4. **Well label ≠ biological condition** — some runs use well labels for sample IDs (e.g. `"Patient_001"`) not treatment conditions. Ask the user what the label encodes before inferring.

5. **`Well` vs well index** — `adata.obs["Well"]` uses values like `"A1"`; `wells_df` is indexed by well name (same format). Verify they match before mapping.

6. **Replicate at well level, not cell level** — each well is a replicate, not each cell. Pseudobulk aggregates cells → wells → pseudobulk samples. Do not treat individual cells as replicates.

7. **Multi-run designs** — when loading multiple runs, `run_name` in `obs` is an additional grouping dimension. Confirm whether runs are biological or technical replicates before merging conditions across runs.

8. **Handing off `condition` without the confounding verdict** — downstream skills then pick a default that silently deletes the effect. Per-well imaging scaling and `Well` batch correction both subtract anything constant within a well. Always report Step 6 with the design, and note that **replicate wells do not remove the confounding** — only heterogeneity *within* a well does.
