# Sandbox — Working Reference

Detail behind [SKILL.md](SKILL.md). Read the section you need; you do not need all of it for a
given task. Everything here describes the `elembio-sandbox` MCP runtime
(`ELEMBIO_RUNTIME == "ebc-code-interpreter"` — the authoritative in-sandbox check).

Contents: [Limits at a glance](#limits-at-a-glance)
· [Pre-installed stack and installing packages](#pre-installed-stack-and-installing-packages)
· [Filesystem layout and the writable budget](#filesystem-layout-and-the-writable-budget)
· [Reusing kernel state across turns](#reusing-kernel-state-across-turns)
· [Executing code: budgets, detach, and polling](#executing-code-budgets-detach-and-polling)
· [Large stdout is spilled, not truncated](#large-stdout-is-spilled-not-truncated)
· [Generating analysis code](#generating-analysis-code-for-this-environment)

## Limits at a glance

| Limit | Value |
|---|---|
| Writable volume (`/`, `/tmp`, `/data/session`, shared) | ~900 MB |
| Output handed back off the sandbox | ≤ ~1 GB (downsample / summarize above) |
| `request_upload` inbound file | ≤ 5 GiB (mount larger sources) |
| `execute_code` / `execute_command` `timeout_seconds` | default 120 s · max 7200 s |
| Inline wait before a run detaches | ~45 s |
| `install_packages` build cap | 300 s (then detaches; poll `get_results`) |
| Inline stdout before it spills to a file | ~64 KB |
| Concurrent kernel calls | 1 (calls queue) |
| Idle hibernation | after ~1 h idle; auto-resumes on the next call |

## Pre-installed stack and installing packages

The image ships a scientific Python stack (`spatialdata`, `scanpy`, `anndata`, `squidpy`, plus
the usual numpy / pandas / scikit-learn / matplotlib, and `ratarmount`) alongside the `elembio`
analysis and CLI packages. **Don't assume an exact manifest — check before installing:**
`import <mod>`, `which <tool>`, or `execute_command` running `pip show <pkg>` / `pip list`.

- **`install_packages(sandbox_id, packages=[...])`** runs `pip install` in the persistent
  kernel. Use it only for a genuinely missing dependency, not to reinstall the stack. Pass
  **every package for the session in one call** — each call queues behind the kernel's in-flight
  run, and pip resolves them together. Version specifiers (`numpy==2.1.0`) are accepted.
- A failed install returns `success: false` with pip's stderr — re-check the name (a typo, or a
  package already present under a different import name).
- A long install detaches at the inline wait like any run (poll `get_results`); the install
  itself is capped at 300 s.
- Installs are **session state**: they persist across turns like variables and imports, and are
  **lost on a kernel reset or a new sandbox** — reinstall after `resume_state_lost` or an empty
  namespace (see [RECOVERY.md](RECOVERY.md)).

Prefer the preinstalled stack and the libraries the `multiomics` skills recommend over ad-hoc
additions.

## Filesystem layout and the writable budget

There is a **single writable volume, ~900 MB usable**, shared by `/`, `/tmp`, and
`/data/session`. This is the constraint most failures trace back to.

- **`$ELEMBIO_EXECUTION_OUTPUTS_DIR`** — a fresh per-call directory the kernel sets before each
  `execute_code`. Exactly the files you write there are returned as that call's `artifacts`
  (each with a `path`, `s3_uri`, `s3_status`, `size_bytes`). Write anything you want returned or
  kept beyond the call here, e.g.
  `os.path.join(os.environ["ELEMBIO_EXECUTION_OUTPUTS_DIR"], "plot.png")`.
- **`/data/session`** — writable and durable: it survives kernel restart, idle-reap/resume, and
  a fresh-provision reset. Use it for **checkpoints** (a stable path you control), not for
  per-call outputs.
- **`/tmp`, `/app`, other worker-local paths** — scratch only, not guaranteed to survive.
- **`/runs/<id>`, `/executions/<id>`, `/storage/<id>`** — conventional mount points (not
  pre-mounted; mount on demand — see below).

Rules that follow from per-call attribution:

- **Revise by writing a new file, never by overwriting a prior path.** Each past chat message
  keeps rendering the file it referenced, so overwriting rewrites history. A revised plot is a
  new file in the current call's output dir.
- **Read a prior call's output by the absolute path that call returned** — it stays valid.
- **Only `fetch_artifact` a path you received in a call's `artifacts`.** A failed
  `execute_code` writes no artifact, so a path "from" a failed call points at nothing.

Confirming durability and handing files back:

- Each artifact carries **`s3_status`** — `present` (landed durably), `syncing` (still
  uploading), or `unknown` / absent (no verdict; do **not** read as "gone"). Only `present`
  guarantees the file survives the session.
- **`fetch_artifact <path>`** mints a short-lived HTTPS URL for the user, and — when the
  bytes fit in the response — the image itself as a content block your client renders inline;
  describe what it shows rather than re-posting the link. When the bytes are omitted,
  `image_content_omitted` says why (`read_failed` is transient — retry; `too_large` and
  `unsupported_type` mean a different file is needed to see it), and `url` still delivers the
  original file untouched.

Budget hygiene:

- Check `df -h /data/session` before workloads expected to produce >100 MB.
- **Outputs larger than ~1 GB cannot be handed back** off the sandbox (there is no download path
  for arbitrary large files; `request_upload` is inbound only, ≤5 GiB). Downsample, chunk, or
  stream summary statistics back instead of materializing the full artifact.
- On `ENOSPC`, free space by deleting intermediates from `/tmp` and `/data/session`.
- Cache env vars (`MPLCONFIGDIR`, `NUMBA_CACHE_DIR`, `XDG_CACHE_HOME`, `TMPDIR`) draw on the same
  budget; for long sessions redirect them under `/data/session/.cache/`.
- Keep high-`fsync` workloads (e.g. SQLite committing in a loop) **off** `/data/session` — each
  commit blocks on a FUSE-to-S3 upload.

Getting data in — mount conventions, mandatory flags, the `ratarmount` recipe, and
`request_upload` — is its own topic: see **[MOUNTING.md](MOUNTING.md)**. The budget rule: mount
and read in place; never copy a remote store into `/data/session`.

## Reusing kernel state across turns

The Python kernel is long-lived within a session: variables, imports, loaded data, and
`Loader`/AnnData objects from earlier turns stay in the namespace.

- **Reuse what is already live** before reloading — re-creating a `Loader` re-reads the store
  over FUSE and re-parses metadata.
- **Check what survived with `execute_code`** (`dir()`, `'loader' in globals()`), **not**
  `get_status` — `get_status` reports liveness/resources, never namespace contents.
- **Reload only when** the kernel was reset (state lost — see [RECOVERY.md](RECOVERY.md)), you
  need a pristine copy after an in-place mutation, or a genuinely different store is required.

**Inherited state is inherited memory.** A sandbox reused from earlier work — especially one
resumed from hibernation — comes back holding whatever that work left in the namespace, which
can be many GB of the limit before you load anything. Reuse is still the right default — that
live namespace is the thing of value — but **check `get_status.resource_usage` against the size
of the store you are about to open**, and say what you found if headroom is thin. Two ways out,
in order of preference:

- **Free the specific objects** you no longer need (`del`, then `gc.collect()`) — keeps the
  mounts and installed packages.
- **`create_sandbox` fresh** when the old namespace is both large and irrelevant. Prefer this
  outright when a stale object could silently contaminate a result (a leftover `adata` or
  `loader` bound to a different run); a subprocess `execute_command` cannot be contaminated that
  way, but `execute_code` can.

### Idle hibernation and resume

An idle sandbox is **hibernated** after roughly an hour: its memory image (variables, imports,
loaded data) is snapshotted and its compute is released. The next tool call **auto-resumes** it
— nothing special to call; the call just takes longer while the image is restored (it may return
a provisioning response with `phase: "restoring"` — see [RECOVERY.md](RECOVERY.md)).
`/data/session` persists throughout, and a detached run survives the idle window whether or not
you poll. State is lost only if the restore **fails**, which surfaces as `resume_state_lost` —
the one case where you reload your checkpoint and re-run mounts. This is why checkpointing
expensive state (below) turns any reset into a one-line reload.

### Checkpointing

After a costly `load → QC → normalize` pipeline, write a checkpoint to a **stable** path on the
durable volume:

```python
adata.write_h5ad("/data/session/checkpoint/adata.h5ad")   # a fixed path you control
```

`/data/session` survives kernel restart, idle-reap/resume, and a fresh-provision reset, so any
state loss becomes a one-line reload instead of rerunning the pipeline. Do **not** checkpoint
into the per-call `$ELEMBIO_EXECUTION_OUTPUTS_DIR`.

**Reloading a large checkpoint is a cold S3 read — budget for it.** After a reset the file is no
longer in the mount cache, so it streams from S3 (~60–100 MB/s; a just-written file reads ~10×
faster from cache, which is why a reload can time out in production but not in a quick test). A
full `sc.read_h5ad(...)` of a multi-GB checkpoint materializes the whole matrix and can exceed
the default 120 s `execute_code` budget. Two ways through:

- Pass a larger **`timeout_seconds`** (up to 7200); the read detaches and you poll
  `get_results`.
- When you only need to slice or inspect it, open it lazily with
  **`adata = ad.read_h5ad(path, backed="r")`** — returns in under a second regardless of size;
  `.shape`, `.obs`, and `adata[mask].to_memory()` all work.

A reload that times out is a budget/throughput limit on a healthy file, **not** a corrupt
checkpoint — raise the timeout or use `backed="r"`. Only suspect the mount if
`get_status.kernel_status.session_storage_state` reads `degraded`.

## Executing code: budgets, detach, and polling

### Know the cost before you send it

Decide a call's expected duration *before* calling, because that decision is what you owe the
user out loud (see **Keeping the user in the loop** in [SKILL.md](SKILL.md)) and what sets
`timeout_seconds`. Rough costs over a `--disk-cache-size 0` mount:

| Operation | Expect |
|---|---|
| `elembio … mount` | ~1 s (detached; the FUSE mount is live on return) |
| `ls` / `find` over a mount | well under a second per listing — but each one is an S3 round trip |
| `Loader(<store>)` + `available_tables` | **seconds, even on a 45 GB `.zarr.zip`** — it opens the archive index, it does not read tables |
| **First `load_tables(..., lazy=True)` per store** | **tens of seconds to minutes — this is the expensive call.** Parses every table's metadata over FUSE, uncached. Measured at 98 s on a 45 GB `.zarr.zip`, which detached. |
| Reading from a table already loaded | seconds — the first `load_tables` is what paid for it |
| Reading `obs` / `obsm` columns, groupbys over them | seconds |
| A full pass over `X` | not measured here; a 6-of-95-chunk sample was ~3 s warm, so budget minutes cold and sample unless you need every cell |
| Cold checkpoint reload | throughput-bound, see [Checkpointing](#checkpointing) |

Two consequences worth internalizing:

- **Lazy is not free, and the first touch pays for the rest.** Constructing the `Loader` is
  cheap and tells you nothing about what follows — the bill arrives on the first
  `load_tables`, even with `lazy=True`. Budget *that* call generously and expect everything
  after it to be fast. Reversing the two is the easy mistake: a 1-second `Loader` reads as
  "this store is fast".
- **Reach for `obs` / `obsm` before `X`.** Most run-level QC (depth, missingness, per-tile
  counts, spatial coordinates) is precomputed in `obsm` and costs seconds. Going to `X` for a
  number that already exists in `obsm` turns a 3-second call into a full-store scan.

- **`execute_code` / `execute_command`** take an optional `timeout_seconds` (default **120**,
  max **7200**). This is the run's max wall-clock, not a transport limit.
- **One call runs at a time.** Batch independent steps into one cell.
- **Inline wait is ~45 s.** A run that outlasts it does not fail — it **detaches** and keeps
  running server-side. Branch on the result you got back:
  - A normal completed result (stdout + success/failure metadata) → treat as final.
  - **Detached metadata** `{status: "running", run_id, sandbox_id, poll_with: "get_results"}` →
    expected for a long run, not an error. The run survives past the idle window whether or not
    you poll. Then:
    - **Tell the user it detached before you poll.** One line naming the work and that it is
      running server-side. A detach is the single most disorienting thing the user cannot see:
      from their side a silent poll loop and a hung session look identical.
    - **Poll `get_results(sandbox_id, run_id)`** — returns `{status: "running"}` while in flight,
      then the same stdout/success/artifacts shape a short run returns. Space polls seconds
      apart; **retrieval is multi-turn** — submit in one turn, fetch in a later one if needed.
      **Do not narrate each poll** — one line at detach, one line when it lands. If the wait
      passes a couple of minutes, or you check `resource_usage` and something looks wrong
      (memory climbing toward the limit, a run far past its expected duration), say so then
      rather than at the end.
    - **Lead with the result when it lands.** The user has been waiting on this one — give them
      the headline number before you move on to the next call.
    - **A long wait is useful time, not dead time.** Read the skill you will need next, or
      check the docs for the step after this one, while the run is in flight — but say that is
      what you are doing, so the interleaved tool calls are legible.
    - **Recover a lost `run_id`** from `get_status.kernel_status.active_run_id` while the sandbox
      is still `busy` (it is absent when idle). `get_status.last_run` is the historical handle if
      you lost it entirely.
    - Use `get_status` alongside for resource progress (`resource_usage`, `busy`) — never for
      namespace contents.
- **`interrupt_execution(sandbox_id)`** stops the in-flight run and **keeps the kernel alive**
  (variables, imports, packages, mounts survive); `get_results` then returns the interrupted
  state (`success: false`). It is best-effort and Python-level: code stuck in a C extension
  (BLAS, native calls) or uninterruptible I/O yields only when it returns.

## Large stdout is spilled, not truncated

When a cell prints more than the inline limit (~64 KB), `execute_code` writes the **complete**
output to a spill file and returns a `stdout_spill` artifact naming the path.

- **Nothing is truncated.** Read the file instead of re-running — a follow-up `execute_code`
  can open the path and slice it (`open(p).read(200_000)`, `itertools.islice`), or
  `fetch_artifact` fetches it for the user.
- **Do not re-run the code to "get the output back"** — a re-run costs the same compute and
  spills again.
- **`size_bytes` is absent on the spill artifact** — absent means "not reported", not "empty";
  use `os.path.getsize(path)` if you need the length.

Printing a large result is usually the wrong shape anyway: write it to
`$ELEMBIO_EXECUTION_OUTPUTS_DIR` as a file and print a short summary instead.

## Generating analysis code for this environment

- **Deliver ordered cells through `execute_code`** (definitions then calls) and run them
  yourself. Do not write a standalone `.py` for a human to run, and do not launch Python via
  `subprocess` / heredocs — use `execute_code` (Python) and `execute_command` (shell) directly.
- **Use processes for parallelism.** `multiprocessing.Pool`, joblib's `loky` backend
  (`n_jobs`), and Numba `parallel=True` all use the session's vCPUs. `/dev/shm` is deliberately
  tiny (64 MiB) and is **not** scratch — put scratch in `/tmp`.
- **matplotlib is headless** (`MPLBACKEND=Agg`). Save figures into `$ELEMBIO_EXECUTION_OUTPUTS_DIR`
  and surface them with `fetch_artifact`.
- **Make runs reproducible** — set seeds (`np.random.seed`, `sc.settings.seed`) and print key
  library versions and parameters.

Library-specific import/call conventions belong to the `multiomics` skills, not here.
