---
name: sandbox-environment
description: "Use when working with Element Biosciences / AVITI / AVITI24 data — listing or resolving runs, executions, or cloud storage; downloading or mounting data; or running multiomics, single-cell, spatial, imaging, OPS, QC, or differential-expression analysis. Routes between the local `elembio` CLI (quick listing, metadata, small downloads) and the ElemBio Cloud sandbox MCP `elembio-sandbox` (for compute, or when no local CLI is available): create one sandbox and reuse it, mount cloud data in place with elembio-cli, and drive the analysis with the Element Biosciences `multiomics` skills (QC, normalization, and modality-specific pipelines) on the preinstalled stack (spatialdata / scanpy / anndata / squidpy)."
metadata:
  version: 0.7.2
  author: elembio
---

# ElemBio Cloud Sandbox — Environment & Routing

This skill gets you into the right environment with the user's data reachable. The analysis
itself is **driven by the Element Biosciences `multiomics` skills** — this skill hands off to
them once data is loaded and does not restate their pipelines.

Three companion files hold the detail that does not belong in this entry point — read them when
you reach the situation they cover, not before:

- **[MOUNTING.md](MOUNTING.md)** — reaching the user's Cloud data: sandbox mount conventions and
  mandatory flags, the `ratarmount` recipe for `.zarr.zip`, `request_upload`, and the
  no-credentials dead-ends.
- **[REFERENCE.md](REFERENCE.md)** — how to work efficiently in the sandbox: the preinstalled
  stack and installing packages, filesystem and the writable budget, where outputs go and how to
  confirm they are durable, reusing kernel state across turns, checkpointing, long-running
  detached runs, and what to say to the user while one is in flight.
- **[RECOVERY.md](RECOVERY.md)** — what to do when a call times out, errors, or returns an
  unexpected shape. The single most important rule lives there: a timeout is **not** a dead
  sandbox, and needlessly recreating one throws away loaded data.

## Choosing an environment

The sandbox (`elembio-sandbox` MCP) is a remote persistent-kernel Python environment with the
standard analysis stack (`spatialdata`, `scanpy`, `anndata`, `squidpy`, …) and `elembio-cli`
**preinstalled**. It acts as the calling user, so `elembio-cli` works inside it with no extra
credentials. A locally installed `elembio` CLI, when present, handles lightweight data access
with no spin-up. Check for it first — `which elembio && elembio whoami` — then route:

| Task | Use |
| --- | --- |
| List / resolve runs, executions, or storage; read metadata; small download | **Local `elembio` CLI** if present — fastest, no spin-up |
| Compute — QC, normalize, cluster, DE, imaging, or any multiomics / spatial analysis | **Sandbox** — the stack is preinstalled and kernel state persists across calls |
| No local CLI available | **Sandbox** — it ships `elembio-cli`; run `elembio …` via `execute_command` |
| Read cloud data in place (no copy) | **Either** — `elembio … mount` works from the local CLI or inside the sandbox |

**When both are viable, ask the user** whether to work locally or in the sandbox — for running
`elembio-cli` and for compute alike — rather than guessing. Skip the question only when the
choice is forced: no local CLI (use the sandbox), or a preference the user already stated.

## Quickstart (sandbox path)

The end-to-end shape for "analyze my run in the sandbox" — each step's detail is in the
companion files:

1. **`create_sandbox`** → keep the returned `sandbox_id` and reuse it on every call. If it comes
   back `INSTANCE_PROVISIONING`, poll `get_status` until `INSTANCE_READY`.
2. **Mount the data** with `execute_command` (see [MOUNTING.md](MOUNTING.md)), e.g.
   `elembio runs mount <run-id> /runs/<run-id> --disk-cache-size 0`.
3. **Load it** with `execute_code`:
   `from elembio_spatialdata_analysis.load import Loader; loader = Loader("<the .zarr under /runs/<run-id>>")`.
4. **Drive the analysis** with the `multiomics` skills (start at their `index`). Reuse live
   kernel state across turns; checkpoint expensive state to `/data/session` (see
   [REFERENCE.md](REFERENCE.md)).
5. **Return results** — write figures/tables to `$ELEMBIO_EXECUTION_OUTPUTS_DIR`, then
   `fetch_artifact` the path and describe what it shows to the user.

## Sandbox lifecycle

- **`create_sandbox` once** per body of work, and whenever you hold no `sandbox_id` — it is the
  only way in, as there is no tool that lists or adopts an existing sandbox. Then reuse the
  returned id on every subsequent call: the kernel keeps variables, imports, and installed
  packages between calls, so never create a second sandbox for the same task.
- **A `status: "INSTANCE_PROVISIONING"` response is success, not an error.** The `sandbox_id`
  is valid and reserved for you; compute is starting (usually seconds; a cold host boot can
  reach ~2–3 minutes). Wait `retry_after_seconds`, then poll `get_status` with that id until
  `status` is `"INSTANCE_READY"`. Do **not** call `create_sandbox` again — a second call starts
  a second sandbox. The response carries its own stopping rule (`max_wait_seconds`, `on_timeout`); honor
  it. See [RECOVERY.md](RECOVERY.md) for the full provisioning / restore handling.
- **`execute_code`** runs Python in the persistent kernel; **`execute_command`** runs bash
  (use it for all `elembio …` invocations). **One call runs at a time** — batch independent
  steps into a single cell rather than issuing many small calls back-to-back.
- To reuse a sandbox whose id you still hold, probe `get_status` first (a liveness/resource
  probe, **not** a namespace check). A hibernated sandbox auto-resumes on the next call; one
  that reports dead or is not found means `create_sandbox`.
- Only `destroy_sandbox` when the user explicitly asks to end the session.
- **Subagents get their own sandbox.** Never hand your `sandbox_id` to a subagent (it has no
  credentials for your session and a separate `/data/session`); pass mount metadata in the
  brief so it re-mounts in its own sandbox, and have it return results as text.

## Keeping the user in the loop

The sandbox is remote and invisible: the user sees your messages, not your tool calls. **They
cannot tell a 98-second metadata load from a hung session.** So narrate — a one-line
what-and-why *before* a call, a one-line outcome *after*. Skip it only for genuinely trivial or
internal steps (introspecting a shape, retrying a one-line syntax fix).

The three cases where silence actually costs the user something — narrate these even when the
step feels internal:

| Case | Say this |
|---|---|
| **Anything that will take a while** — a first `load_tables` over FUSE, a full pass over `X`, a big install (*not* a mount; that is ~1 s) | What you are about to do and roughly why it is slow, *before* you call it. Then the outcome. |
| **A detached run** (`status: "running"`) | That it detached and what it is working on — *before* the first `get_results`. Never poll in silence; see [REFERENCE.md](REFERENCE.md). |
| **A discovery result that changes the plan** | Say it when you learn it, not in the final summary. "This run has no execution, so the store ships inside the run itself" is a course correction the user should see happen. |

Two habits that make this cheap rather than chatty:

- **Predict the cost out loud when you raise `timeout_seconds`.** Passing `timeout_seconds=1800`
  is you declaring the call may run 30 minutes — tell the user that, not just the tool.
- **Say what you sampled.** When you subset to keep a call affordable — N of M chunks, a cell
  subsample, one well of many — state it in the same message as the number it produced, not only
  in a closing caveat. A median from 5% of the cells is a different claim than a median.

## Getting the user's Cloud data in

The sandbox has **no AWS credentials** — reach data only through `elembio-cli` (pre-installed
and pre-authenticated as the user). **Mount and read Cloud data in place; never copy a remote
store into `/data/session`** (the ~900 MB writable budget will `ENOSPC` mid-copy).

The mount commands and their mandatory flags, the `ratarmount` recipe for `.zarr.zip`,
`request_upload` for local files, and the credential dead-ends are in **[MOUNTING.md](MOUNTING.md)**.
The canonical `elembio` CLI surface lives in the multiomics `elembio-cloud-data-access` skill.

## Running the analysis

Analysis in the sandbox is **driven by the Element Biosciences `multiomics` skills** — do not
improvise with generic single-cell / Scanpy defaults. Once data is mounted, load the `.zarr`
store with `Loader`, then **start with the `multiomics` `index` skill**, which routes to the
right specialist and platform order (it encodes AVITI24 / DISS conventions that generic
defaults get wrong). Refer to it by name — it installs as a separate plugin, so relative file
paths from here will not resolve. If the multiomics skills aren't available, install the
Element Biosciences multiomics-skills plugin (or ask the user to) rather than substituting
ad-hoc analysis. This skill deliberately does not duplicate those steps.

## Outputs

- Write durable outputs to **`$ELEMBIO_EXECUTION_OUTPUTS_DIR`** (a fresh per-call directory);
  exactly those files come back as `artifacts`. Use `/tmp` for scratch.
- Confirm a file landed with its artifact **`s3_status: present`**; surface a file — or an
  inline image — to the user with **`fetch_artifact`**.
- The handling rules (revise-by-new-file, the full `s3_status` values, and when `fetch_artifact`
  omits an image's bytes) are in **[REFERENCE.md](REFERENCE.md)**.
