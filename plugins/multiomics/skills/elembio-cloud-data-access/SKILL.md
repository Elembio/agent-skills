---
name: elembio-cloud-data-access
description: Canonical CLI surface for retrieving data from ElemBio Cloud — listing runs, executions, and storage connections; downloading; mounting via FUSE; authentication. Use when the user needs to get data from ElemBio Cloud before running analysis skills.
metadata:
  version: 0.4.5
  author: elembio
---

# ElemBio Cloud Data Access

## Purpose

This skill describes the `elembio` CLI — what it does, how it authenticates, and how to use it to discover, download, and mount data from ElemBio Cloud. It is environment-agnostic: it documents CLI behavior, not where the CLI is run. Operating environments (sandboxes, IDE integrations) may pre-configure parts of this surface and provide their own usage guidance on top; consult that guidance when present.

---

## Preconditions

Before any `elembio` subcommand returns useful results, both must hold:

### CLI present

```bash
which elembio && elembio version
```

If absent, install per [Installation](#installation).

### Authenticated session

```bash
elembio whoami
```

Expected output: the caller's email and organization. If absent or "not authenticated," run `elembio login` and provide an API key, or set `ELEMBIO_API_KEY` (see [Authentication](#authentication)). Without a valid session, subcommands silently return empty results or fail with non-obvious errors.

In environments where the CLI is pre-installed or pre-authenticated, these checks are no-ops. Run them only when state is unknown.

---

## Installation

Recommended installation per platform. For the full matrix of architectures and direct-download URLs, consult the `elembio-cli` distribution.

### macOS and Linux — Homebrew

```bash
brew tap elembio/tap
brew install elembio
```

Upgrade:

```bash
brew update && brew upgrade elembio
```

### Windows — PowerShell

Mount is not currently supported on Windows; use `elembio runs download` instead. Install:

```powershell
$json = Invoke-RestMethod "https://elembio-cloud-downloads.s3.us-west-2.amazonaws.com/elembio-cli/releases/latest.json"
Invoke-WebRequest -Uri $json.builds."windows-amd64" -OutFile "$env:USERPROFILE\elembio.exe"
Move-Item "$env:USERPROFILE\elembio.exe" "C:\Windows\System32\elembio.exe"
elembio version
```

---

## Authentication

An API key is required. Generate one from account settings at [elembio.io](https://elembio.io).

```bash
elembio login
# Paste the API key when prompted: elembio_xxxxxxxxxxxx
elembio whoami    # verify
```

Alternatively, set the key via environment variable (scripts, CI, container bootstrap):

```bash
export ELEMBIO_API_KEY=elembio_xxxxxxxxxxxx
```

---

## Resource Model

| Resource               | Description                                             | Typical contents                                                |
| ---------------------- | ------------------------------------------------------- | --------------------------------------------------------------- |
| **Run**                | A sequencing run from an AVITI instrument               | Raw instrument outputs, FASTQ files, SpatialData `.zarr` stores |
| **Execution**          | A secondary analysis derived from a run                 | Processed outputs, SpatialData `.zarr` stores, analysis reports |
| **Storage connection** | An external bucket (e.g. S3) linked to the organization | Arbitrary files uploaded or synced to cloud storage             |

Runs and executions are the primary sources for SpatialData `.zarr` stores consumed by downstream analysis skills. Storage connection contents vary by source. **Any** downstream execution may emit an SDO (`cells2stats` and `elembio-multiomics-nf` are examples). Discover them generically by probing execution outputs for a `.zarr.zip` rather than by matching workflow names (see [Listing a run's executions and finding their SDOs](#listing-a-runs-executions-and-finding-their-sdos)).

> **This CLI reads, it does not launch.** `elembio` can list, get, download, and mount runs / executions / storage, but it **cannot start or trigger a workflow** (e.g. `elembio-multiomics-nf` and other Catalyst pipelines). When a needed execution doesn't exist yet, tell the user to run that workflow in ElemBio Cloud — never imply you can trigger it from here.

---

## Resolving a run by name

If the caller provides a run name (e.g. `<run-name>`) rather than a run ID (e.g. `run-abc123`):

```bash
# Exact name match — do not grep a paginated list (default is 100 rows, older runs will be missed)
elembio runs list --filter 'name:"<run-name>"' --max-items 0 -o wide
```

A double-quoted value is an **exact, case-sensitive** match; unquoted is a case-insensitive substring match, so `name:RUN-0001` also matches `RUN-00011`. `--max-items 0` removes the page-size cap so every match comes back.

If the exact match returns nothing, fall back once to an unquoted substring match for partial or differently-cased names:

```bash
# JSON output for scripting
elembio runs list --filter 'name:RUN-000' --max-items 0 -o json | jq '.[] | {id, name, status}'
```

To filter `list -o json` in Python, read stdin and loop over the array directly:

```bash
elembio storage list -o json | python3 -c '
import sys, json
for item in json.loads(sys.stdin.read()):   # bare array — parse and loop directly
    uri = item.get("uri", "")
    if "demo" in uri.lower():
        print(item["id"], item["name"], uri)
'
```

`list -o json` emits a **bare JSON array**, not a wrapped object — true for `runs`, `executions`, `storage`, and `instruments`. Iterate it directly (`jq '.[]'`, or `json.loads(out)` then loop); don't pipe `2>&1` into a JSON parser — a stray stderr line breaks `json.loads`.

If no match is found, stop and request clarification from the user. Do not guess at run IDs or attempt downloads against unconfirmed identifiers — silent failures and empty results consume significant time without surfacing the underlying mismatch.

---

## Listing the most recent runs or executions

Lists come back **newest first**, sorted server-side: runs by `time_updated`, executions by `time_created`. The order is fixed, and filters don't change it. Bound "recent" with a time window and ask for exactly the number you need:

```bash
elembio runs list --filter 'time_created>=7d' --max-items 10 -o wide                    # up to 10 runs created in the last 7 days
elembio runs list --filter 'type:multiomics time_created>=30d' --max-items 5 -o json    # up to 5 multiomics runs from the last 30 days
elembio executions list --filter 'time_created>=7d' --max-items 10 -o wide              # up to 10 executions from the last 7 days
elembio executions list --filter 'status:succeeded time_created>=30d' --max-items 1 -o json  # latest successful execution in the last 30 days
```

A relative date means "now minus that long": `>=7d` is *within* the last 7 days, and `<7d` is *older* than 7 days. The window makes "recent" explicit: an empty result means nothing in that window, so widen it (`30d`, `90d`) rather than dropping the filter.

**For "recent" or "latest", don't fetch the full list and sort it yourself** — `--max-items 0` without a narrowing filter pages through the organization's whole history for an answer the first page already holds. Keep `--max-items 0` for filters that must return **every** match (a name lookup, a run's executions): the default cap is 100 rows, and `-o json` prints no truncation notice, so a cut-off list looks complete.

Runs are ordered by last update, not creation, so an older run can rank above a newer one. For most recently *created*, *started*, or *completed*, narrow with a time filter first, then sort that small set:

```bash
elembio runs list --filter 'time_completed>=14d' --max-items 0 -o json \
  | jq 'sort_by(.time_completed) | reverse | .[:5] | .[] | {id, name, time_completed}'
```

---

## Listing a run's executions and finding their SDOs

A run can have several SDOs: the one written directly off the instrument, plus one per secondary-analysis **execution** (e.g. `cells2stats` or `elembio-multiomics-nf`). To offer the user a choice, enumerate the run's executions and discover which produced a `.zarr.zip`.

```bash
# Resolve the run first (see "Resolving a run by name"), then list its executions.
# Run→execution filter key is input.run.id / input.run.name — NOT run.id or run_id.
elembio executions list --filter 'input.run.id:<run-id>' --max-items 0 -o json
elembio executions list --filter 'input.run.name:"<run-name>"' --max-items 0 -o wide
```

Prefer `input.run.id` once the run is resolved — ID keywords always match exactly. `input.run.name` follows the quoting rule above; unquoted, it also returns executions of runs whose names contain yours. `workflow.name` labels which workflow produced each execution.

**Discover SDOs generically — do not hard-filter by workflow name.** For each execution, probe its outputs and keep any that contain a `.zarr.zip`:

```bash
elembio executions files <exec-id> --recursive        # look for *.zarr.zip in the output
```

Surface each SDO-bearing execution as an option, labeled with its `workflow.name` for readability. This picks up every SDO-producing workflow without hard-coded names. `workflow.name` or `platform:catalyst` filters can *narrow* the probe (e.g. `--filter 'input.run.id:<id> platform:catalyst'`), but the source of truth is "does the execution output contain a `.zarr.zip`."

Consume a chosen execution's SDO exactly like a run's — `executions mount` / `executions download` → `Loader` (see [Mounting](#mounting) / [Downloading](#downloading)).

See [IMPLEMENTATION.md → Executions](./IMPLEMENTATION.md#executions) for a scripted enumeration that builds the option list.

---

## Resolving a storage URI

A raw `s3://` URI is reachable through the CLI only if it falls under a registered storage connection. Before acting on one, confirm it:

```bash
elembio storage list -o json | python3 -c '
import sys, json
uri = "s3://bucket/key"   # the URI you were handed
hits = [i for i in json.loads(sys.stdin.read())
        if i.get("uri") and uri.startswith(i["uri"])]
print(hits or "NOT REGISTERED")
'
```

If the bucket is not in the list, it is **not accessible** through the CLI — stop and ask the user to register it as a storage connection (or supply a URI under one that already exists). Like an unconfirmed run ID, an unregistered bucket only burns turns on silent failures; do not attempt to reach it by other means.

---

## Downloading

```bash
# Preview without transferring
elembio runs download run-abc123 --dry-run

# Download to a local directory
elembio runs download run-abc123 --output-dir ./my-run
```

The same shape applies to `elembio executions download` and `elembio storage download`. Common flags: `--prefix <subpath>` to scope, `--include` / `--exclude` to filter, `--output-dir <path>` for destination. See REFERENCE.md for the full flag matrix.

### Vending temporary AWS credentials

`elembio storage download --mode credentials -o json` returns short-lived AWS STS credentials (access key, secret, session token) scoped to the storage connection. Useful for callers that need to drive boto3, the AWS CLI, or other AWS SDKs directly against the underlying bucket.

```bash
elembio storage download <connection-id> --mode credentials -o json
```

**Scope and restrictions.** Vended credentials reach only _that_ connection's bucket and prefix — they grant no access to any other bucket and must not be reused against one. Some operating environments disable credential vending entirely; consult environment-specific guidance (e.g. the `sandbox-environment` skill) before relying on this path.

---

## Mounting

Mounts expose remote data as a local read-only filesystem via FUSE. Files are fetched on demand via HTTP range requests with local caching.

```bash
elembio runs mount run-abc123 ./my-run
elembio executions mount exec-xyz789 ./my-exec
elembio storage mount conn-123 ./my-storage --prefix subdir/
```

By default, `elembio mount` runs in the **foreground** and blocks until SIGINT (Ctrl+C) or `unmount` — useful when the mount's lifetime is tied to an interactive shell. Pass `--detach` (or `-d`) to daemonize after the FUSE filesystem is live; the command then returns exit 0 with the daemon's pid and log path on success, or a non-zero exit code with the underlying error (authentication, bad prefix, mountpoint conflict, etc.) on failure. `ELEMBIO_MOUNT_DETACH=1` sets the detached mode as the default for the process environment, so callers (sandboxes, CI) can opt in without modifying every invocation.

For storage connections pointing at large or root buckets, scope with `--prefix <subpath>/` (note the trailing slash for directory-style prefixes — without it, the prefix is treated as a literal key and the mount usually lists empty).

To inventory or clean up mounts:

```bash
elembio mount list                # all live + stale mounts (table / JSON / YAML / wide)
elembio unmount <path>            # unmount a single path
elembio unmount <path> --force    # forcibly unmount, including stale mounts
```

Not currently supported on Windows. See IMPLEMENTATION.md for FUSE prerequisites and platform notes.

---

## Quick Start

```bash
elembio whoami                                                 # confirm session
elembio runs list --filter status:completed --max-items 10     # discover: 10 most recently updated
elembio runs mount run-abc123 ./my-run --detach                # mount in background (omit --detach to block in foreground)
ls ./my-run
```

```python
from elembio_spatialdata_analysis.load import Loader

loader = Loader("./my-run/spatialdata/run.zarr")
adata = loader.load_tables(table="ThreePrimeUntargeted", lazy=False)  # or imaging / OPS table if that is the goal
```

---

## Detail Files

- [IMPLEMENTATION.md](./IMPLEMENTATION.md) — Detailed workflows; FUSE prerequisites and mount internals; connecting downloaded/mounted data to downstream analysis skills
- [REFERENCE.md](./REFERENCE.md) — Command tables, download flags, mount flags, environment variables, output formats

## Related Skills

- [spatialdata-loading-and-access](../spatialdata-loading-and-access/SKILL.md) — Load `.zarr` stores into AnnData after downloading or mounting
- [raw-data-ingestion](../raw-data-ingestion/SKILL.md) — Load `.h5ad`, `.parquet`, or saved AnnData files
- [index](../index/SKILL.md) — Route to the right pipeline after data is available locally (3′ / OPS / imaging / legacy targeted)

## Disclaimer

This content is provided "as-is" without warranty of any kind. AI-generated output may contain
inaccuracies. Review and verify all outputs before relying on them.
