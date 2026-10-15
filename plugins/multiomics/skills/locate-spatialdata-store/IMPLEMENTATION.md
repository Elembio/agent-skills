# Locate SpatialData Store — Implementation

## Contents

1. [Local path validation](#1-local-path-validation)
2. [Directory search](#2-directory-search) — incl. [Local run discovery (no cloud CLI)](#local-run-discovery-no-cloud-cli)
3. [CLI-based run discovery](#3-cli-based-run-discovery)
4. [Verify and hand off](#4-verify-and-hand-off)

---

## 1. Local path validation

Determine what the user provided:

```python
from pathlib import Path

p = Path(user_input)

if not p.exists():
    raise FileNotFoundError(f"Path does not exist: {p}")

if p.suffix == ".zarr" or p.name.endswith(".zarr.zip"):
    sdo_path = p
elif p.is_dir():
    sdo_path = None  # proceed to directory search (§2)
else:
    raise ValueError(f"Not a .zarr store or directory: {p}")
```

---

## 2. Directory search

Search a directory for SpatialData stores:

```python
from pathlib import Path

def find_stores(root: Path) -> list[Path]:
    """Find .zarr / .zarr.zip stores under root."""
    stores = []
    for pattern in ("*.zarr.zip", "*.zarr"):
        stores.extend(root.rglob(pattern))
    seen = set()
    unique = []
    for s in stores:
        if s.name not in seen:
            seen.add(s.name)
            unique.append(s)
    return sorted(unique)
```

**Never `rglob` or `find` on a FUSE mount** — FUSE directory walks are one range-read per entry and will time out. Use `rglob` only on local/downloaded directories. On FUSE mounts, navigate directly to the known paths.

### Prioritization within a run directory

Prefer stores in this order:

1. `Cytoprofiling/Instrument/SpatialData/` — instrument SDO (primary)
2. `Cytoprofiling/<other>/SpatialData/` — analysis SDO
3. Other locations

If multiple stores exist, print each with its table inventory (§5) and ask the user which one to use.

### Local run discovery (no cloud CLI)

When the user names a run but is working locally (no `elembio` cloud CLI), there is **no default runs location** — resolve it, then reuse `find_stores`:

```python
from pathlib import Path

# 1. Where does this user keep runs? Prefer what you already know.
runs_root = remembered_runs_dir_from_memory_or_context()   # a path, or None
if runs_root is None:
    runs_root = ask_user_for_run_directory()               # or a direct .zarr.zip path
    remember_runs_dir(runs_root)                            # persist to memory for next time (if available)

# 2. Find the run's store(s) under that root, matching the run name.
stores = [s for s in find_stores(Path(runs_root)) if run_name in str(s)]
# Prefer the instrument SDO (Cytoprofiling/Instrument/SpatialData/); if several match,
# list them with table inventories and let the user pick.
```

`remembered_runs_dir_from_memory_or_context` / `ask_user_for_run_directory` / `remember_runs_dir` are **not library calls** — they describe agent behavior: check memory/context first, otherwise ask the user, then persist the answer to memory. Once a store is chosen, [validate](#4-verify-and-hand-off) and apply the aligned-3′-missing check (§4).

---

## 3. CLI-based run discovery

When the user provides a run name without a local path, use [elembio-cloud-data-access](../elembio-cloud-data-access/SKILL.md) to find the run by name and authenticate. That skill owns the CLI mechanics.

### Quick path — construct the S3 URI directly

The instrument SDO path is predictable. Combine the run's `output.uri` with the known relative path:

```python
import json, subprocess

result = subprocess.run(   # quoted = exact name; if several runs share it, runs[0] is the most recently updated
    ["elembio", "runs", "list", "--filter", f'name:"{run_name}"', "--max-items", "0", "-o", "json"],
    capture_output=True, text=True
)
runs = json.loads(result.stdout)
if runs:
    run_uri = (runs[0].get("output") or {}).get("uri")   # may be absent; see SKILL.md
    instrument_sdo = (f"{run_uri}Cytoprofiling/Instrument/SpatialData/{runs[0]['name']}.zarr.zip"
                      if run_uri else None)
```

This avoids mounting or downloading the entire run when you only need the SDO path. When `output.uri` is absent, skip the instrument option and rely on the execution SDOs (§ [Build the full SDO option list](#build-the-full-sdo-option-list)).

### Build the full SDO option list

A run may expose several SDOs. Combine the instrument-direct SDO with the SDO-bearing executions, then present them all **with full paths** and let the user choose. The execution enumeration is owned by [elembio-cloud-data-access](../elembio-cloud-data-access/SKILL.md) — reuse its `executions_with_sdo()` helper rather than re-implementing the probe:

```python
import json, subprocess
# from elembio-cloud-data-access/IMPLEMENTATION.md → Executions
# executions_with_sdo(run_id) -> [{id, name, workflow, status, sdo_path}, ...]

runs = json.loads(subprocess.run(
    ["elembio", "runs", "list", "--filter", f'name:"{run_name}"', "--max-items", "0", "-o", "json"],
    capture_output=True, text=True).stdout or "[]")
run = runs[0] if runs else None
# `output` may be absent (see SKILL.md), so resolve defensively.
run_uri = (run.get("output") or {}).get("uri") if run else None
instrument_sdo = (f"{run_uri}Cytoprofiling/Instrument/SpatialData/{run['name']}.zarr.zip"
                  if run_uri else None)

options = []
if instrument_sdo:
    options.append({"source": "instrument-direct (missing any assigned 3′ data)",
                    "path": instrument_sdo})   # omitted automatically when output.uri is absent
for e in (executions_with_sdo(run["id"]) if run else []):   # generic; groups naturally by e["workflow"]
    options.append({"source": f'{e["workflow"]} execution {e["id"]} '
                              f'({e["name"]}, {e["status"]})',
                    "path": e["sdo_path"]})

for o in options:
    print(f'- {o["source"]}\n  {o["path"]}')     # show the FULL path for every option
```

Then present these to the user and **wait for a selection** — do not load, mount at scale, or start analysis until they confirm which SDO to use.

### Full path — mount or download for local analysis

Use [elembio-cloud-data-access](../elembio-cloud-data-access/SKILL.md) to mount or download the run. After the run is local, search the resulting directory (§2). The instrument SDO is typically at:

```
./<run-dir>/Cytoprofiling/Instrument/SpatialData/<run-name>.zarr.zip
```

---

## 4. Verify and hand off

Once you have a `.zarr.zip` path, use [spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md) to open it with `Loader`. That skill owns the `Loader` API; the key check unique to this skill is 3′ detection:

```python
from elembio_spatialdata_analysis.load import Loader

loader = Loader(str(sdo_path))
target_types = loader.available_target_types
tables = loader.available_tables

# 3′ detection — unique to this skill's responsibility
if "ThreePrimeUntargeted" in target_types and "ThreePrimeUntargeted" not in tables:
    # Aligned 3′ is declared in the panel but not assembled in this store.
    # STOP and let the user choose (see SKILL.md → "If aligned 3′ is missing"):
    #   1) analyze without aligned 3′, or
    #   2) provide the elembio-multiomics-nf SDO that includes the counts.
    # If they give a path, validate BEFORE switching:
    #   alt = Loader(user_provided_path)
    #   assert "ThreePrimeUntargeted" in alt.available_tables  # else warn + re-ask
    ...
```

### Hand off

1. [spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md) — `Loader`, `load_tables`, table discovery
2. [index](../index/SKILL.md) — pick the pipeline for available modalities
