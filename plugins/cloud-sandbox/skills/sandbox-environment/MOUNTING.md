# Sandbox — Data Access & Mounting

How to reach the user's ElemBio Cloud data from inside the sandbox. Detail behind
[SKILL.md](SKILL.md).

**Boundaries — do not restate what these own:**

- The **canonical `elembio` CLI surface** (auth, the runs / executions / storage resource model,
  download/mount flag semantics, `elembio mount list`, `elembio unmount`) lives in the
  multiomics **`elembio-cloud-data-access`** skill. Consult it by name for a flag's existence or
  meaning — it installs as a separate plugin, so relative paths from here won't resolve.
- **Loading a store into AnnData** (`Loader`, `read_zarr`, modality tables, lazy vs
  materialized) lives in the multiomics **`spatialdata-loading-and-access`** skill.

This file carries only what the **sandbox runtime** forces or changes.

Contents: [What the sandbox already set up](#the-sandbox-already-did-the-setup) · [What you can reach](#what-you-can-reach)
· [Mounting](#mounting) · [Reading a `.zarr.zip`](#reading-a-zarrzip)
· [Getting local files in](#getting-local-files-in-request_upload)

## The sandbox already did the setup

- **The CLI is pre-installed and pre-authenticated** at bootstrap with the inbound API key,
  acting as the signed-in user. The `which elembio` / `elembio version` / `elembio whoami`
  preflights that `elembio-cloud-data-access` documents are **no-ops here — skip them**.
- **`ELEMBIO_MOUNT_DETACH=1` is set**, so mount commands daemonize **without** an explicit
  `--detach`: they return synchronously with exit 0 (plus a pid and log path) on success, or a
  non-zero exit with an actionable error on failure. No silent hangs.
- **No AWS credentials, by design.** Reach data only through `elembio` subcommands. Every
  bypass is a dead end: `boto3` / `s3fs(anon=True)` / `--no-sign-request`, IMDS
  (`169.254.169.254`), presigned-URL guessing, and `elembio storage download --mode credentials`
  to vend keys into the session (disabled here — do not use it). A subagent has no credentials
  either.

## What you can reach

`elembio-cli` reaches three kinds of ElemBio Cloud resource, all mounted the same way:

- **Runs / executions** — AVITI sequencing runs and the analyses derived from them; the usual
  source of SpatialData `.zarr` stores.
- **Storage connections** — S3 buckets registered in ElemBio Cloud, reached through the
  `elembio storage` verbs. Two flavors, accessed identically:
  - **Element-operated Catalyst storage** — managed buckets Element provisions for your runs and
    analysis outputs (named `ElemBio Catalyst …`).
  - **Your own linked bucket** — an external S3 bucket you have referenced into ElemBio Cloud;
    once registered it is reachable exactly like Catalyst storage.

A raw `s3://` URI works **only** if it falls under a registered connection; if it does not,
register it in ElemBio Cloud first — the sandbox has no credentials to reach an unregistered
bucket by other means. Objects that have aged into Deep Archive must be restored before they can
be read. Discovery, registration, and the archive-restore flow live in the multiomics
`elembio-cloud-data-access` skill.

## Mounting

**Mount and read in place — never copy a remote store into the writable budget.** Downloading,
`cp`-ing, or unzipping a multi-GB store into `/data/session` competes for ~900 MB and fails
mid-copy with `No space left on device`. Point your reader at the mount path directly.

Conventional (not pre-mounted) mount points, mounted on demand via `execute_command`:

```bash
elembio runs mount       <run-id>  /runs/<run-id>            --disk-cache-size 0
elembio executions mount <exec-id> /executions/<exec-id>     --disk-cache-size 0
elembio storage mount    <conn-id> /storage/<conn-id> --prefix <subpath>/ --disk-cache-size 0
```

Flags that the sandbox forces or that pay off here:

- **`--disk-cache-size 0` (mandatory).** The default 50 GB cache writes into the same ~900 MB
  volume and `ENOSPC`s the kernel. Zero keeps the mount as pure FUSE range-reads that never
  touch the writable budget.
- **`--prefix <subpath>/`** — scope a storage connection to just the path you need. Keep the
  **trailing slash**; without it the prefix is treated as a literal key and the mount usually
  lists empty.
- **`--block-size` / `--readahead`** — the in-AWS benchmark defaults (8 MB block, 8-block
  readahead) suit most reads; retune to the pattern: `--block-size 1MB --readahead 1` for random,
  chunk-sized reads (`read_zarr`, exploration); `--block-size 8MB` (or higher) for bulk
  sequential streaming.
- **`--cache-size` (in-memory read cache; default 1 GB per mount).** Because the disk cache is
  disabled above, this L1 cache is the *only* cache and it draws on the same RAM as your
  analysis. It is a ceiling, not a reservation — it only fills under heavy reads — but with
  several mounts open or a memory-tight workload, lower it to free RAM (keep it ≥ ~512 MB, or
  wide concurrent reads thrash). This is separate from the ~900 MB *disk* budget.

Mounts persist for the life of the session and are **lost on kernel death or a new sandbox** —
re-run your mount commands after any reset (see [RECOVERY.md](RECOVERY.md)). You rarely need to
unmount in a session-scoped sandbox; use `elembio unmount <path>` if you must
(`elembio-cloud-data-access` has the full lifecycle surface).

## Reading a `.zarr.zip`

- **Tables / points — no mount needed.** `Loader` opens a zip natively:
  `Loader("/storage/<conn-id>/run.zarr.zip")`. This is the common case; do not reach for
  ratarmount for it.
- **Rasters / a full `SpatialData` object — FUSE-mount the zip.** `read_zarr` cannot open a zip
  directly, so expose its contents read-only with **`ratarmount`, which is pre-installed in the
  sandbox** (no `pip install`, and none of the macOS FUSE setup the generic recipe mentions —
  this is Linux). The mount is a view over range-reads, not an extraction, so it does not
  consume the writable budget:

```bash
# via execute_command — mountpoint is an empty scratch dir under /tmp
mkdir -p /tmp/store
ratarmount /storage/<conn-id>/run.zarr.zip /tmp/store
```

```python
import spatialdata as sd
# point read_zarr at the nested .zarr path inside the mount
sdata = sd.read_zarr("/tmp/store/run.zarr")
```

Full recipe, fallbacks, and the `Loader(sdata)` handoff live in the multiomics
`spatialdata-loading-and-access` skill.

## Getting local files in (`request_upload`)

For a file **not** already in the Cloud, `request_upload(sandbox_id, path, size_bytes)` mints a
one-shot presigned PUT. `path` is the destination and must resolve under **`/data/session`**
(e.g. `/data/session/input.csv`); `size_bytes` must be the file's exact byte length. Upload with
the returned `url` / `method` / `headers` — replay the headers verbatim and PUT exactly
`size_bytes` bytes, or S3 rejects the signature.

**Where it lands:** once the PUT succeeds the object is at the `path` you gave, on the durable
S3-backed `/data/session` volume (the response echoes back `path` and its `s3_uri`). Read it
there directly with `execute_code` / `execute_command`, and it also shows up in `list_artifacts`.

The cap is **5 GiB**, but the local writable budget is only ~900 MB — so read a large upload in
place / stream it rather than materializing it all at once, and for very large sources prefer
mounting over uploading. Prefer mounting Cloud data over uploading it either way.
