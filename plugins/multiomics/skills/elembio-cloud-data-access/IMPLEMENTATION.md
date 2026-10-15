# ElemBio Cloud Data Access — Implementation

## Contents

- [Preflight checklist](#preflight-checklist)
- [Finding a run by name](#finding-a-run-by-name)
- [Runs](#runs)
- [Executions](#executions)
- [Storage connections](#storage-connections)
- [Downloading data](#downloading-data)
- [Mounting data (browse without downloading)](#mounting-data-browse-without-downloading)
- [Connecting to downstream analysis skills](#connecting-to-downstream-analysis-skills)
- [Common pitfalls](#common-pitfalls)

---

## Preflight checklist

Run both checks before issuing any other `elembio` command. If either fails, fix it first — do not proceed.

```bash
# Check 1: CLI installed
which elembio && elembio version
# Expected: path to binary + version string

# Check 2: Authenticated
elembio whoami
# Expected: your email and organization name
```

**If `elembio whoami` fails:**

```bash
elembio login
# Paste your API key: elembio_xxxxxxxxxxxx
elembio whoami   # confirm before continuing
```

**If you have no API key:** stop and ask the user — *"Do you have an ElemBio Cloud API key? If not, please log in at elembio.io and generate one from account settings."* Do not continue without a confirmed authenticated session.

---

## Finding a run by name

Users often provide a run name (e.g. `<run-name>`) rather than a run ID (e.g. `run-abc123`). Always resolve the name to an ID before downloading.

```bash
# Exact name match server-side — default list is 100 rows, older runs will be missed by grep
elembio runs list --filter 'name:"<run-name>"' --max-items 0 -o wide
# Nothing found? Fall back to a case-insensitive substring match
elembio runs list --filter 'name:<partial-name>' --max-items 0 -o wide
```

`--max-items 0` lifts the default 100-row cap so no matches are missed; `-o json` gives no truncation notice. See [REFERENCE.md → Filter Syntax](./REFERENCE.md#filter-syntax) for quoting rules.

**If nothing is found after searching both runs and executions — stop immediately.** Do not retry variations or guess at IDs. Ask the user:

> "I could not find a run named 'X' in your account. Can you confirm the exact run name, or copy the run ID directly from [elembio.io](https://elembio.io)? It's also possible the run belongs to a different organization, or is still processing."

---

## Runs

A run represents a sequencing run from an AVITI instrument.

Lists are newest first — see [SKILL.md → Listing the most recent runs or executions](./SKILL.md#listing-the-most-recent-runs-or-executions).

```bash
elembio runs list                                                # newest 100
elembio runs list --max-items 10 -o wide                         # 10 most recently updated
elembio runs list --filter status:completed
elembio runs list --filter 'status:completed time_created>=7d'   # multiple filters combine with AND
elembio runs get run-abc123 -o json
elembio runs files run-abc123
elembio runs download run-abc123 --output-dir ./my-run --dry-run
```

See [Downloading data](#downloading-data) for the full set of download flags and [REFERENCE.md → Filter Syntax](./REFERENCE.md#filter-syntax) for filter fields.

---

## Executions

An execution is a secondary analysis derived from a run (e.g., a `cells2stats` cytoprofiling build, the `elembio-multiomics-nf` multiomics workflow, an alignment pipeline). Executions follow the same `list / get / files / download` verbs as runs, plus `logs`:

```bash
# List every execution for a run — input.run.id is always exact; quote input.run.name for exact
elembio executions list --filter 'input.run.id:cto_abc123' --max-items 0 -o wide
elembio executions list --filter 'input.run.name:"<run-name>"' --max-items 0 -o json
elembio executions list --max-items 10 -o wide                   # 10 most recently created, any run
elembio executions get exec-456 -o json
elembio executions files exec-456 --recursive
elembio executions download exec-456 --output-dir ./my-execution
elembio executions logs exec-456
```

> **The run→execution filter key is `input.run.id` / `input.run.name`** — not `run.id`, not `run_id`. Full valid keyword set: `flow.name`, `input.execution.id`, `input.execution.name`, `input.run.id`, `input.run.instrument_name`, `input.run.name`, `name`, `output.kind`, `platform`, `started_by`, `status`, `time_*`, `workflow.name`. An unknown keyword or bad value makes the CLI print the valid set — probe with a bogus value if unsure.

### Finding which executions produced an SDO (generic, workflow-agnostic)

A run may have several SDOs, one per secondary-analysis execution, in addition to the instrument SDO. **Discover them generically rather than by hard-coding workflow names**, so every SDO-producing workflow is picked up without code changes. For each execution, probe its outputs and keep any that contain a `.zarr.zip`:

```python
import json, subprocess

def executions_with_sdo(run_id: str) -> list[dict]:
    """Executions of a run that emitted a .zarr.zip: {id, name, workflow, status, sdo_path}."""
    execs = json.loads(subprocess.run(
        ["elembio", "executions", "list", "--filter", f"input.run.id:{run_id}",
         "--max-items", "0", "-o", "json"],
        capture_output=True, text=True).stdout or "[]")
    out = []
    for e in execs:
        files = subprocess.run(
            ["elembio", "executions", "files", e["id"], "--recursive", "-o", "json"],
            capture_output=True, text=True).stdout or "[]"
        for f in json.loads(files):
            path = str(f.get("path") or f.get("key") or f.get("uri") or "")
            if path.endswith(".zarr.zip"):
                out.append({
                    "id": e["id"], "name": e.get("name"),
                    "workflow": (e.get("workflow") or {}).get("name"),
                    "status": e.get("status"),
                    "sdo_path": f.get("uri") or f.get("path") or f.get("key"),
                })
    return out
```

`workflow.name` or `platform:catalyst` filters can *narrow* the probe (Catalyst runs the cytoprofiling workflows) to cut the number of `files` calls, but the source of truth is "the output contains a `.zarr.zip`." (`output.kind` only distinguishes `bam`/`cram`/`fastq`/`parquet`/`vcf` — it does not flag SDOs, so don't rely on it here.) Consume a chosen SDO via `executions mount` / `executions download` → `Loader` (see [Mounting](#mounting-data-browse-without-downloading) / [Downloading](#downloading-data)).

> Field names under `-o json` (`path` vs `key` vs `uri`) can vary by CLI version; the snippet checks the common ones. Confirm against `elembio executions files <id> -o json` for your version.

---

## Storage connections

A storage connection links an external bucket (e.g., S3) to your organization. Browse and download arbitrary files with the same verbs as runs:

```bash
elembio storage list
elembio storage get conn-789
elembio storage files conn-789 --prefix runs/2024/
elembio storage download conn-789 --prefix data/ --output-dir ./bucket-data
```

---

## Downloading data

All download commands (`runs download`, `executions download`, `storage download`) share a common set of flags. See [REFERENCE.md → Download Flags](./REFERENCE.md#download-flags) for the full table.

### Examples

Download only `.zarr` contents:

```bash
elembio runs download run-abc123 --include "*.zarr*" --output-dir ./my-run
```

Preview what would be downloaded:

```bash
elembio executions download exec-456 --dry-run
```

Resume an interrupted download (skip already-fetched files):

```bash
elembio runs download run-abc123 --skip-existing --output-dir ./my-run
```

Download with integrity verification:

```bash
elembio runs download run-abc123 --checksum --output-dir ./my-run
```

---

## Mounting data (browse without downloading)

The `mount` command presents remote files as a local read-only directory via FUSE — fetched on demand with two-tier caching.

```bash
elembio runs mount <run-id> <mountpoint>
elembio executions mount <execution-id> <mountpoint>
elembio storage mount <storage-connection-id> <mountpoint>
```

`mount` runs in the **foreground** by default, blocking until Ctrl+C or `unmount`. Pass `--detach` (or `-d`) to daemonize once FUSE is live — the command returns with the daemon's pid and log path, or non-zero on failure (auth, bad prefix, mountpoint conflict). `ELEMBIO_MOUNT_DETACH=1` flips the default environment-wide. `--log-file <path>` overrides the auto-generated log path (custom paths aren't discoverable by `elembio mount list`); `--detach-timeout` (default `15s`) bounds the parent's wait for child readiness.

### Inventory and cleanup

```bash
elembio mount list                       # live + stale mounts (supports -o wide / json / yaml)
elembio unmount <mountpoint>             # unmount a single path
elembio unmount <mountpoint> --force     # force-unmount, including stale mounts
```

### `--prefix` gotcha: trailing slash matters

For directory-style prefixes (the common case for `runs/<id>/`, `executions/<id>/`, or any S3 sub-prefix like `<prefix>/`), **always include the trailing slash**. The CLI passes the value through as a raw key prefix:

| Prefix value | Behavior |
|---|---|
| `<prefix>/` | Lists everything under the directory (correct for browsing) |
| `<prefix>` | Treated as a literal key prefix; matches only objects whose name *starts with* that string at the bucket root, often returning an empty mount |
| `<prefix>/<run-name>.zarr.zip` | Mounts a single object; the file appears at `<mountpoint>/<run-name>.zarr.zip` |

If `ls <mountpoint>/` is empty after a mount that otherwise looked successful, a missing trailing slash is the likely cause. `elembio unmount <mountpoint>` and re-mount with the slash.

### FUSE prerequisites

FUSE must be installed before using mount.

**Linux** — FUSE is a kernel module in all modern distributions:

| Distro | Install command |
|---|---|
| Debian / Ubuntu | `sudo apt install fuse3` |
| Fedora / RHEL | `sudo dnf install fuse3` |
| Arch | `sudo pacman -S fuse3` |
| Alpine | `sudo apk add fuse3` |

**macOS** — Requires macFUSE:

```bash
brew install --cask macfuse
```

First install on Apple Silicon requires allowing the kernel extension in **System Settings > Privacy & Security**, then a reboot.

**Windows** — Mount is not supported. Use `download` instead.

### Examples

Foreground mount — browse from a second shell:

```bash
elembio runs mount run-abc123 ./my-run
# in another shell: ls ./my-run/
```

Detached mount — the call returns once FUSE is live:

```bash
elembio runs mount run-abc123 ./my-run --detach
ls ./my-run/
elembio unmount ./my-run    # when finished
```

Scope to a subset and tune the cache:

```bash
elembio runs mount run-abc123 ./my-run --prefix "spatialdata/" --cache-size 1GB
```

See [REFERENCE.md → Mount Flags](./REFERENCE.md#mount-flags) for the full flag matrix.

### How caching works

The mount uses a two-tier cache: an in-memory L1 (256 MB LRU, instant repeated reads) and a persistent disk L2 (50 GB at `$TMPDIR/elembio-mount-cache`, survives across mount sessions so re-mounting the same data avoids re-downloading previously read blocks). Sequential access triggers readahead (up to 3 blocks prefetched); concurrent requests for the same block are deduplicated into a single HTTP request.

---

## Connecting to downstream analysis skills

A FUSE mount appears as a normal local directory, so the resulting local path — whether from `download` or `mount` — passes directly to the analysis package's `Loader`:

```python
from elembio_spatialdata_analysis.load import Loader

loader = Loader("./my-run/spatialdata/run.zarr")
print(loader.summary())
adata = loader.load_tables(table="ThreePrimeUntargeted", lazy=False)  # discover via available_tables
```

For mounted data the first access to each file region triggers an HTTP fetch with subsequent reads served from cache. Table and metadata loading (`Loader`'s first step) is fast because those payloads are small relative to images and point clouds.

After loading, proceed with the standard pipeline: [cell-quality-control](../cell-quality-control/SKILL.md) → [normalization-strategies](../normalization-strategies/SKILL.md) → downstream analysis.

---

## Common pitfalls

- **Authentication failures look like empty results.** If `elembio whoami` doesn't return an identity, list and download commands return nothing rather than erroring — the most common cause of a long failed session.
- **Mounted `.zarr` stores work with `Loader` but not with random-access-heavy tools.** For workloads that scan every file repeatedly (e.g., full image rendering), downloading first is faster than streaming through the mount.
