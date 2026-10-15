# Imaging Pipeline — Reference

## `var` columns

| Column | Role |
|---|---|
| `feature_class` | `"base"` (derived intensities) or `"cellprofiler"` |
| `suffix` | `""` whole-cell; `_Nuclear`; `_Cytoplasmic`; `_FractionAboveThreshold`; `_MeanIntensityOfFractionAboveThreshold` |
| `is_background` | `True` for `Bkg_*` channels |
| `batch_name` | `CP*` = primary morphology group; `P*` = additional IF markers when present. Empty (`"nan"`) on `cellprofiler` features and backgrounds: a CellProfiler feature's stain is the end of its name (`Granularity_3_Actin` → `Actin`) |
| `target_type` | Often `"CellPaint"` for CP-group targets; `"ProteinIF"` for additional IF. Role ≠ chemistry (a histone in CP can be IF). |

## Scope: technical grouping vs biological meaning

The `batch_name` split (CP* vs P*) is a **technical grouping** — which targets were in the CP01 vs P01 staining batch. It is **not** a semantic boundary between "morphology" and "protein":

- CP01 can contain antibodies (e.g., a mitochondrial antibody may appear as "Mitochondria" in `var`). IgG validation applies to those targets.
- P01 morphology panels (Expanded, Neuro, Alternative) contain structural/organelle antibodies that describe cell organization — they are morphology.
- P01 functional panels (Immune Subtyping, Cell Signaling) contain markers for cell state or identity.
- Custom antibodies may appear that are not in any standard panel.

"Morphology" as a biological concept can span both batches. The agent must inventory what targets are present and reason about them, not assume CP = morphology and P = protein.

**Default scope for "morphology":** CP01 targets + AreaShape geometry. Offer to include morphology-panel P01 targets or other proteins.

**Default scope for "protein" / "IF":** all P01 targets. Clarify which panel or subset if multiple are present.

**IgG validation:** applies to any antibody target — all P01 targets plus any antibody-based CP01 targets. Not limited to the P01 batch.

## Product panels (ProteinIF)

Standard panels shipped with the product. Channel names in `var` may use abbreviated, alternate, or suffixed forms — do not rely on exact name matching. Custom antibodies may also appear. Use this catalog for reasoning context, not deterministic lookup.

| Panel | Targets |
|---|---|
| **Control** | Mouse Control IgG, Rabbit IgG isotype control |
| **Control (Mitochondria)** | TOMM20 (often labeled "Mitochondria" — may be in CP01 batch despite being an antibody) |
| **Expanded Morphology** | GM130 (Golgi), ERp57 (ER), α-Tubulin (cytoskeleton), LAMP2 (lysosomes) |
| **Neuro Morphology** | Nucleophosmin, GFAP, βIII-Tubulin, MAP2 |
| **Alternative Morphology** | β-Actin, Na+/K+-ATPase (membrane), Lamin A/C (nuclear envelope), Histone H3 |
| **Immune Subtyping** | CD8A, CD4, CD3E, CD19 |
| **Cell Signaling** | Phospho-ERK1/2, ERK1/2, Phospho-HSP27, Phospho-p38 MAPK |

**Handling unknown targets:** inspect `target_type`, `batch_name`, and name patterns. When unsure what a target stains or whether it is an antibody, ask the user.

## Feature relevance by target biology

Base intensity and CellProfiler features are extracted uniformly per target, but which features carry the most biological signal depends on what the target stains:

- **Punctate / foci markers** (DNA damage markers, phospho-proteins): Texture and Granularity capture spatial pattern — intensity alone misses the distribution of foci within the cell.
- **Membrane / organelle markers** (plasma membrane, lysosomes, ER, Golgi): RadialDistribution captures localization within the cell — where the signal concentrates relative to the nucleus.
- **Cytoskeletal markers** (tubulin, actin): Texture captures filament organization and structural integrity.

The specific targets in any run will vary. The agent should consider what each target stains when deciding which CellProfiler feature types are most relevant for the analysis question.

## Targets vs features

| Concept | Meaning |
|---|---|
| **Imaging target** | One channel/stain; each gets base + CellProfiler features (similar count per target) |
| **Primary morphology group** | `CP*` targets — QC intensity gate; may include IF reagents |
| **Additional IF markers** | `P*` / `ProteinIF` when present — not used for intensity QC |
| **Geometry** | Cell-level (`obsm["cell_geometry"]`), not per-target |

## P01 NaN/zero handling

`filter_cells("Imaging")` scopes NaN, empty, and intensity QC to CP01 morphology channels only (via `select_morphology_channels`). NaN or zero values in P01 targets survive QC. Before building an analysis matrix that includes P01 columns:

1. Check for NaN in the scoped P01 columns.
2. Report which targets have NaN and how many cells are affected.
3. Drop cells with NaN in scoped columns, or handle appropriately. NaN in a P01 column typically means the antibody was not measured for that cell — it is not zero signal.

See [IMPLEMENTATION.md §4a](./IMPLEMENTATION.md#4a-p01-nanzero-check) for the pattern.

## Helpers (`ebsd.tl.morphology`)

| Function | Use |
|---|---|
| `select_morphology_channels` | Channels for NaN, empty, and intensity QC gates (all scoped to these — CP01 only) |
| `get_morphology_feature_groups` | Optional inventory helper for CP-style name groups |

## Scaling reference (MAD-z: median / 1.4826·MAD)

The transform is fixed; the **reference group** is a design decision.

| Reference | When | Effect on biology |
|---|---|---|
| **Per-well** | Wells internally heterogeneous (pooled / OPS, mixed lines per well) | Removes everything constant within a well — correct only when that is technical |
| **Control wells** (run-level) | Biology constant within wells **and** shared controls exist | Corrects run drift, preserves condition differences (PyCytominer `samples=`) |
| **Global** (all cells in run) | Biology constant within wells, no shared control | Preserves biology; leaves per-well technical variation in |
| Per-well MAD + global median | Per-well spread differs but means must be kept | Corrects heteroscedasticity only |
| `ebsd.pp.normalize(..., "Imaging")` | **Do not** — `NotImplementedError` | — |

| Guardrail | Value |
|---|---|
| MAD floor | Relative to the feature's **global** MAD (~`1e-3`), never absolute `1e-6` |
| Degenerate MAD | Fall back to std, or drop the feature |
| Clip scaled values | ±10 (matches `sc.pp.scale(max_value=10)`) |
| `inf` / `NaN` cleanup | Before **and** after scaling |

## Analysis matrix exclusions

| Exclude | How to spot | Why |
|---|---|---|
| Background | `is_background` or `Bkg_*` (base **and** CP-on-Bkg) | S/B diagnostic only |
| IgG / antibody controls | Name contains `IgG` | Validation only |
| Cleave / process channels | `Cleave` / `FinalCleave` | Process markers |
| Absolute `Location_*` | Names starting with `Location_` | FOV pixel coordinates — positional noise |

Keep `RadialDistribution_*` (relative within-cell localization), Texture, Granularity, Intensity_*, and AreaShape.

## Derived metrics (obs-level diagnostics)

| Metric | Formula | Matrix? |
|---|---|---|
| S/B | `signal / max(bkg, ε)` | No — obs diagnostic only |
| Nuclear fraction | `Nuclear / (Nuclear + Cytoplasmic)` | No — obs diagnostic only |
| Nuclear area fraction | `NuclearArea / Area` from geometry | No — obs diagnostic only |
| Fraction above threshold | base suffix `_FractionAboveThreshold` | No — obs diagnostic only |

## Absolute `Location_*` (embedding noise)

| Example | What it is |
|---|---|
| `Location_CenterMassIntensity_X_<stain>` | Absolute X pixel of intensity-weighted centroid in the tile |
| `Location_MaxIntensity_Y_<stain>` | Absolute Y pixel of the brightest spot in the tile |

These vary with where the cell sits in the FOV and with well density / spatial layout. Do **not** confuse with `RadialDistribution_*` (relative within-cell localization), which belongs in the matrix.

## Not this skill

| Request | Skill |
|---|---|
| Legacy barcoding Protein / Transcript panels | [legacy-targeted-pipeline](../legacy-targeted-pipeline/SKILL.md) |
| Spatial image overlays | [spatial-multimodal-visualization](../spatial-multimodal-visualization/SKILL.md) |
| Run/cell QC details | [run-quality-overview](../run-quality-overview/SKILL.md), [cell-quality-control](../cell-quality-control/SKILL.md) |
