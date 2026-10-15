---
name: locate-spatialdata-store
description: >-
  Step zero for analyzing a run. Resolves a run name, a directory, or a given
  path to a validated SpatialData store (.zarr or .zarr.zip), and confirms
  which store to use when a run has several. Use before
  spatialdata-loading-and-access, even when the user gives a path.
metadata:
  version: 0.3.5
  author: elembio
---

# Locate SpatialData Store

Every skill in this stack operates on a SpatialData store (`.zarr` or `.zarr.zip`; SDO — spatialdata object). **Locate and validate a store before loading anything.**

> **Confirm before analysis.** A run often has more than one SDO: the one written directly off the instrument, plus one per downstream execution (e.g. `cells2stats` or `elembio-multiomics-nf`). When more than one exists, present the options **with their full paths** and **wait for the user to choose**. Do not load, mount at scale, or start any analysis until they confirm which SDO to use.

## Common mistakes

| Don't | Do instead |
|---|---|
| `elembio runs list -o json \| grep name` | `elembio runs list --max-items 0 --filter 'name:"..."' -o json` (quoted = exact match) |
| `elembio runs list --max-items 0` then sort, to find the latest run | `elembio runs list --max-items 5 -o json` (lists are newest first) |
| `find` or `rglob` on a FUSE mount | `ls` the known `Cytoprofiling/Instrument/SpatialData/` path directly |
| Mount a storage connection to search for a run | Search the runs API first; fall back to storage only if the run isn't registered |
| Offer to trigger / run `elembio-multiomics-nf` to produce missing 3′ | Tell the user to run it in ElemBio Cloud — the `elembio` CLI can't invoke a Catalyst pipeline |

## Decision tree

| User provides | Action |
|---|---|
| `.zarr` or `.zarr.zip` path | [Validate](#validate) |
| A directory (no `.zarr` suffix) | [Search a directory](#search-a-directory) |
| A run name — cloud CLI available / working in the cloud | [Discover SDOs and choose](#discover-sdos-and-choose) |
| A run name — working locally / no cloud CLI | [Discover a local SDO](#discover-a-local-sdo-no-cloud) |
| Nothing / ambiguous | Ask the user what to analyze |

**Cloud vs local for a run name:** if the `elembio` cloud CLI is installed and authed and the user is working in the cloud / sandbox, use the cloud flow; if the user is working from runs already on their own disk (or has no cloud CLI), use the local flow. When it's unclear, ask.

---

## Validate

Confirm the path exists and is a `.zarr.zip` (or `.zarr` directory):

```python
from pathlib import Path

p = Path(user_path)
assert p.exists(), f"Path not found: {p}"
assert p.name.endswith(".zarr.zip") or p.suffix == ".zarr", f"Not a .zarr store: {p}"
```

Then proceed to [spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md) to load with `Loader` and discover available tables.

---

## Search a directory

If the user gives a **local or downloaded** directory, search for `.zarr.zip` stores inside it:

```python
from pathlib import Path

stores = sorted(Path(user_path).rglob("*.zarr.zip")) + sorted(Path(user_path).rglob("*.zarr"))
```

> **Do not `rglob` or `find` on a FUSE mount.** FUSE directory walks issue one range-read per entry and will time out on a large run tree. On mounted runs, navigate directly to known paths:
> ```bash
> ls ./run/Cytoprofiling/Instrument/SpatialData/
> ```

**Typical locations within a run directory:**

| Path pattern | Source | Typical contents |
|---|---|---|
| `Cytoprofiling/Instrument/SpatialData/*.zarr.zip` | Instrument SDO | Imaging ± OPS; no 3′ |
| `Cytoprofiling/<analysis>/SpatialData/*.zarr.zip` | Analysis SDO | Imaging ± OPS; no 3′ |

If multiple stores exist, list them with their table inventories and let the user pick, or prefer the Instrument SDO as the default starting point.

After locating a store, [Validate](#validate).

---

## Discover SDOs and choose

When the user provides only a run name, that run may have **several** SDOs. Enumerate them, present every option **with its full path**, and let the user pick before doing anything else.

### 1. Instrument-direct SDO (from the run)

The run's `output.uri` plus the known relative path gives the instrument SDO URI directly — no mounting needed:

```bash
elembio runs list --filter 'name:"<run-name>"' --max-items 0 -o json | \
  python3 -c "import sys,json; [print((u+'Cytoprofiling/Instrument/SpatialData/') if (u:=(r.get('output') or {}).get('uri')) else '(no output.uri — use execution SDOs)') for r in json.loads(sys.stdin.read())]"
```

The instrument SDO is at `{output.uri}Cytoprofiling/Instrument/SpatialData/{run_name}.zarr.zip`. It is written **directly off the instrument** and is **missing any assigned 3′ data** (3′ polony points only — no `ThreePrimeUntargeted` counts table).

> **When `output.uri` is absent.** Some runs do not report `output` in `runs list` (the run's S3 root is then `storage.connection_name`), and the instrument SDO is not under a predictable run path. Use the execution SDOs below.

### 2. Execution SDOs (from downstream analyses)

List the run's executions and **generically** discover which produced a `.zarr.zip` — do not hard-code workflow names. [elembio-cloud-data-access](../elembio-cloud-data-access/SKILL.md) owns the CLI mechanics: see its [Listing a run's executions and finding their SDOs](../elembio-cloud-data-access/SKILL.md#listing-a-runs-executions-and-finding-their-sdos) section and the `executions_with_sdo()` helper in its [IMPLEMENTATION.md](../elembio-cloud-data-access/IMPLEMENTATION.md#executions).

Group the results by `workflow.name` for display (e.g. `cells2stats` and `elembio-multiomics-nf` executions), but treat the label as informational, not a filter. New SDO-producing workflows then appear automatically.

### 3. Present the options and confirm

Show the user **every** SDO with its full path, grouped by source, and wait for a choice. For example:

> I found several SDOs for **{run name}**:
> - **Instrument-direct** (missing any assigned 3′ data) — `s3://…/Cytoprofiling/Instrument/SpatialData/{run}.zarr.zip`
> - **cells2stats** executions:
>   - `{exec-id}` — {name}, {status}, {created} — `s3://…/{exec}/…/{run}.zarr.zip`
> - **elembio-multiomics-nf** executions:
>   - `{exec-id}` — {name}, {status}, {created} — `s3://…/{exec}/…/{run}.zarr.zip`
>
> Which would you like to analyze?

**Do not load, mount at scale, or start any analysis until the user picks one.** Once they choose, [validate](#validate) the path and hand off.

### Mount or download the chosen SDO

Use [elembio-cloud-data-access](../elembio-cloud-data-access/SKILL.md) to `mount` / `download` the chosen run or execution, then [validate](#validate) and load. The instrument SDO, once local, is at:

```
./<run-dir>/Cytoprofiling/Instrument/SpatialData/<run-name>.zarr.zip
```

---

## Discover a local SDO (no cloud)

Use this when the user is working from runs already on their own disk — no `elembio` cloud CLI (or they've said they're not using the cloud / sandbox). If the cloud CLI is present and authed and the user is working in the cloud, use [Discover SDOs and choose](#discover-sdos-and-choose) instead; when it's ambiguous, ask which.

### 1. Find where the run lives

There is **no default local runs location** — do not guess a path.

1. **Check what you already know.** If your memory or context records where this user keeps AVITI24 runs (a runs directory, or a specific SDO path for this run), start there.
2. **Otherwise, ask** the user for the run directory (or a direct `.zarr.zip` path).
3. **Remember it.** Once the user gives a runs directory, save it to memory so you can check there first next time. (Skip if your runtime has no memory.)

### 2. Search for the run's store

Search that location with the [directory search](#search-a-directory) (`find_stores(root)` in [IMPLEMENTATION.md](./IMPLEMENTATION.md#2-directory-search)), matching the run name. Prefer the instrument SDO path (`Cytoprofiling/Instrument/SpatialData/*.zarr.zip`). If several stores match, list them with their table inventories and let the user pick, then [validate](#validate).

### 3. Check for missing aligned 3′

Load the store, then apply the [aligned-3′-missing check](#if-aligned-3-is-missing--call-out-and-confirm) below **before any analysis** — a local instrument SDO commonly has 3′ in the panel but no assembled counts table.

---

## Instrument SDO vs 3′ execution

An instrument SDO (`Cytoprofiling/` path) always contains imaging data and possibly OPS. If 3′ was sequenced, the instrument SDO **will contain 3′ polony points** but **not** the aligned counts table (`ThreePrimeUntargeted`). The counts table is produced by a separate workflow and lives in an execution SDO.

### Detecting 3′ data

After loading the instrument SDO with `Loader` (see [spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md)):

```python
target_types = loader.available_target_types  # e.g. {'CellPaint', 'ThreePrimeUntargeted', 'OPS'}
tables = loader.available_tables
```

If `'ThreePrimeUntargeted' in target_types` but `'ThreePrimeUntargeted' not in tables`, the run had 3′ sequencing but the aligned counts table is in a separate SDO produced by a downstream workflow.

### If aligned 3′ is missing — call out and confirm

When 3′ is declared in the panel but the counts table is not in the store — i.e.

```python
has_3prime_config = "ThreePrimeUntargeted" in loader.available_target_types  # declared in panel
has_3prime_table  = "ThreePrimeUntargeted" in loader.available_tables         # counts assembled
if has_3prime_config and not has_3prime_table:
    ...  # STOP — do not start analysis; ask the user how to proceed
```

**stop and call it out to the user**, then let them choose — do not silently proceed or auto-hunt for another SDO:

> This SDO for **{run name}** has 3′ in the panel (`ThreePrimeUntargeted` is declared) but the **aligned 3′ counts table isn't assembled** in it. How would you like to proceed?
> 1. **Analyze without the aligned 3′ data** (imaging / OPS / whatever else is present), or
> 2. **Point me to the `elembio-multiomics-nf` SDO** that includes the aligned 3′ counts, and I'll use that instead.

If the user provides a path, **validate that it actually contains the 3′ table** before switching to it:

```python
from elembio_spatialdata_analysis.load import Loader
alt = Loader(user_provided_path)
if "ThreePrimeUntargeted" in alt.available_tables:
    ...  # use this SDO
else:
    ...  # warn: that store also lacks ThreePrimeUntargeted — ask again or proceed without
```

For **cloud** users, the aligned 3′ counts come from an `elembio-multiomics-nf` execution. First check the execution-SDO menu ([Discover SDOs and choose](#discover-sdos-and-choose)); **if such an SDO already exists, offer it.** If **none exists**, the 3′ counts have not been produced yet — tell the user to **run the `elembio-multiomics-nf` workflow in ElemBio Cloud (Catalyst)** to generate the 3′-assigned SDO, then come back and re-run this step once it completes. (This gate's "provide a path" option is the **local / no-cloud** equivalent.)

> **Never offer to trigger, launch, or run `elembio-multiomics-nf` yourself.** The `elembio` CLI here reads data (runs, executions, storage) but **cannot invoke a Catalyst pipeline** — only the user can start the workflow from ElemBio Cloud. State that they need to run it; do not imply you can kick it off.

---

## Next step

Once you have a validated `.zarr` path and know the available tables:

1. [spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md) — load tables with `Loader`
2. [index](../index/SKILL.md) — pick the right pipeline for the available modalities

---

## Detail files

- [IMPLEMENTATION.md](./IMPLEMENTATION.md) — step-by-step procedures, directory search logic, 3′ execution discovery
- [REFERENCE.md](./REFERENCE.md) — SDO path patterns, table names, CLI quick reference

## Disclaimer

This content is provided "as-is" without warranty of any kind. AI-generated output may contain
inaccuracies. Review and verify all outputs before relying on them.
