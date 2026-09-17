---
name: sandbox-environment
description: "Use when working with Element Biosciences / AVITI / AVITI24 data — listing or resolving runs, executions, or cloud storage; downloading or mounting data; or running multiomics, single-cell, spatial, imaging, OPS, QC, or differential-expression analysis. Routes between the local `elembio` CLI (quick listing, metadata, small downloads) and the ElemBio Cloud sandbox MCP `elembio-sandbox` (for compute, or when no local CLI is available): create one sandbox and reuse it, mount cloud data in place with elembio-cli, and drive the analysis with the Element Biosciences `multiomics` skills (QC, normalization, and modality-specific pipelines) on the preinstalled stack (spatialdata / scanpy / anndata / squidpy)."
metadata:
  version: 0.8.0
  author: elembio
---

# ElemBio Cloud Sandbox — Environment & Routing

This skill picks the environment and keeps you out of the three irreversible mistakes. Runtime detail lives on the MCP server. Read it when you reach that situation; do not copy it here.

| URI | Read when |
| --- | --- |
| `elembio-skill://sandbox/sandbox-environment/SKILL.md` | Quickstart, lifecycle, narration, outputs |
| `elembio-skill://sandbox/sandbox-environment/MOUNTING.md` | Mount flags, `.zarr.zip`, `request_upload` |
| `elembio-skill://sandbox/sandbox-environment/REFERENCE.md` | Stack, writable budget, detach / poll, cost |
| `elembio-skill://sandbox/sandbox-environment/RECOVERY.md` | Timeouts, provisioning, kernel and mount faults |

If the host does not support `resources/*`, call `read_skill` with the same URI.

## Choosing an environment

The sandbox (`elembio-sandbox` MCP) is a remote persistent-kernel Python environment with the analysis stack and `elembio-cli` preinstalled. It acts as the signed-in user. A local `elembio` CLI, when present, handles lightweight data access with no spin-up. Check for it first — `which elembio && elembio whoami` — then route:

| Task | Use |
| --- | --- |
| List / resolve runs, executions, or storage; read metadata; small download | **Local `elembio` CLI** if present — fastest, no spin-up |
| Compute — QC, normalize, cluster, DE, imaging, or any multiomics / spatial analysis | **Sandbox** — the stack is preinstalled and kernel state persists across calls |
| No local CLI available | **Sandbox** — it ships `elembio-cli`; run `elembio …` via `execute_command` |
| Read cloud data in place (no copy) | **Either** — `elembio … mount` works from the local CLI or inside the sandbox |

**When both are viable, ask the user** whether to work locally or in the sandbox. Skip the question only when the choice is forced: no local CLI, or a preference the user already stated.

The canonical `elembio` CLI surface is the multiomics `elembio-cloud-data-access` skill. Analysis is driven by the multiomics `index` skill. Do not restate either here.

## Three rules that cost real work if you break them

1. **Create one sandbox per body of work.** There is no list-or-adopt tool. Reuse the returned `sandbox_id`. A second `create_sandbox` starts a second sandbox and discards the loaded namespace.
2. **`status: "INSTANCE_PROVISIONING"` is success, not an error.** The id is valid. Poll `get_status`. Do not create again.
3. **A timeout is not a dead kernel.** Open `RECOVERY.md` and apply it before you act. Only `destroy_sandbox` when the user explicitly asks to end the session.

Subagents get their own sandbox. Never hand your `sandbox_id` to a subagent.

## Onboarding skills (plugin)

- **sandbox-onboard** — install, authenticate, companion plugin, empty project.
- **sandbox-demo** — one run, mount it, `available_tables` only.
- **sandbox-explore** — questions this run's modalities can answer.
- **sandbox-status** — health, audit, diagnose (diagnose routes into `RECOVERY.md`).
