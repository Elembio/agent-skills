---
name: sandbox-environment
description: Use when working with Element Biosciences / AVITI / AVITI24 data — listing or resolving runs, executions, or cloud storage; downloading or mounting data; or running multiomics, single-cell, spatial, imaging, OPS, QC, or differential-expression analysis. Routes between the local `elembio` CLI (quick listing, metadata, small downloads) and the ElemBio Cloud sandbox MCP `elembio-sandbox` (for compute, or when no local CLI is available): create one sandbox and reuse it, mount cloud data in place with elembio-cli, and drive the analysis with the Element Biosciences `multiomics` skills (QC, normalization, and modality-specific pipelines) on the preinstalled stack (spatialdata / scanpy / anndata / squidpy).
metadata:
  version: 0.1.0
  author: elembio
---

# ElemBio Cloud Sandbox — Environment & Routing

This skill gets you into the right environment with the user's data reachable. The analysis
itself is **driven by the Element Biosciences `multiomics` skills** — this skill hands off to
them once data is loaded and does not restate their pipelines.

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

**When it is unclear which to use — ask the user.** If the request does not clearly favor one
side and both are viable (e.g. a local `elembio` CLI is available *and* remote compute is an
option), ask whether they want to work in the **local environment** or the **sandbox** —
both for running `elembio-cli` and for performing compute — rather than guessing. Only skip the
question when the choice is forced: no local CLI (use the sandbox) or an explicit user
preference already stated.

## Sandbox lifecycle

- **`create_sandbox` once** per body of work, then reuse the returned `sandbox_id` on every
  subsequent call. The kernel keeps variables, imports, and installed packages between calls —
  never create a second sandbox for the same task. A `status=provisioning` response is
  success; poll `get_status` with the returned `sandbox_id` until it reports READY.
- **`execute_code`** runs Python in the persistent kernel; **`execute_command`** runs bash
  (use it for all `elembio …` invocations). One call runs at a time — batch independent steps
  into a single cell rather than issuing many small calls.
- To reuse a sandbox from earlier context, probe `get_status` first; if it reports dead,
  `create_sandbox` again. A hibernated sandbox auto-resumes on the next call.
- Only `destroy_sandbox` when the user explicitly asks to end the session.

## Getting the user's Cloud data in

The sandbox has **no AWS credentials** — `boto3` / AWS SDK calls fail with
`NoCredentialsError`, and `elembio storage download --mode credentials` will not work here.
Reach data through `elembio-cli` instead, which runs as the user. Prefer **mounting in place**
(no copy, no size cap); sandbox mounts are detached by default:

```bash
elembio runs mount <run-id> /runs/<run-id> --disk-cache-size 0
elembio executions mount <exec-id> /executions/<exec-id> --disk-cache-size 0
elembio storage mount <conn-id> /storage/<conn-id> --prefix <subpath>/ --disk-cache-size 0
```

For files that are **not** already in the Cloud, `request_upload` mints a one-shot URL into
the session mount — but prefer mounting Cloud data over uploading it.

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

- Write durable outputs to `$ELEMBIO_EXECUTION_OUTPUTS_DIR`; they are returned as artifacts.
- `/data/session` is writable scratch; confirm durability via an artifact's `s3_status`
  (only `present` guarantees it landed) rather than assuming a file saved.
- `download_artifact <path>` mints a short-lived HTTPS URL to hand a file to the user; a small
  PNG / JPEG / GIF / WebP comes back inline.

## Gotchas

- No AWS credentials in the sandbox (see above) — always go through `elembio-cli`.
- One in-flight kernel call at a time; long runs detach and return a `run_id` to poll with
  `get_results`.
- `list_sandboxes` is not available on this endpoint — use `get_status` with a known
  `sandbox_id`, or `create_sandbox` if you have none.
