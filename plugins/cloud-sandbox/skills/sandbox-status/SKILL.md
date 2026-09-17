---
name: sandbox-status
description: Check the health of the ElemBio Cloud sandbox setup. Three modes — health check (dashboard), audit (waste, duplicate sandboxes, CLI overlap), diagnose (router into server-published recovery). Use when asking "is my sandbox working?", "check my sandbox", "audit my setup", "what's broken?", "sandbox status", or "elembio status".
metadata:
  version: 0.1.0
  author: elembio
---

# Sandbox status

Three modes for monitoring and maintaining a sandbox setup. Determine the mode from context, or ask if unclear.

Runtime signals live in the server-published files. This skill reports and routes. It does **not** restate the recovery table.

- Health and resources: `elembio-skill://sandbox/sandbox-environment/SKILL.md` and `REFERENCE.md`
- Recovery: `elembio-skill://sandbox/sandbox-environment/RECOVERY.md`

If the host does not support `resources/*`, call `read_skill` with the same URI.

## Mode 1: Health check

**Trigger:** "check my sandbox", "sandbox status", "is everything working?", or any general status inquiry.

A quick dashboard of the current state.

### Steps

1. Inspect available `elembio-sandbox` tools. If none are available, report disconnected and suggest **sandbox-onboard**. Stop.
2. If you hold a `sandbox_id`, call `get_status`. If you do not, say so — there is no tool that lists other people's sandboxes. Do not call `create_sandbox` just to fill the dashboard.
3. Check whether the multiomics `index` skill is installed.
4. Check for a local `elembio` CLI (`which elembio && elembio whoami`).
5. If the sandbox is live, note `resource_usage`, `busy`, live mounts if reported, and `kernel_status.session_storage_state`.

**Format as a dashboard:**

```
ElemBio sandbox status
======================
MCP:              connected | not connected
Sandbox id:       <id or "none this session">
Instance:         INSTANCE_READY | INSTANCE_PROVISIONING | none
Kernel:           alive | not yet | unknown
Busy:             yes | no
Session storage:  mounted | degraded | pending | unknown
Multiomics:       present | missing
Local elembio CLI: yes | no

Memory / disk:    <from resource_usage, or "n/a">
```

6. End with "Everything looks good." or "Found [N] issues. Want me to diagnose them?"

Do not dump `get_status` JSON. Translate it.

## Mode 2: Audit

**Trigger:** "audit my setup", "clean up", "did I create extra sandboxes?", "what should I remove?"

Find waste. Status **detects and reports**. It does not own the local-CLI-versus-sandbox routing table — that stays in `sandbox-environment`.

### Look for

1. **More than one sandbox created this session.** Each `create_sandbox` is a new host and an empty namespace. If the transcript shows a second create without a recovery reason, flag it.
2. **Unretrieved detached runs.** `status: "running"` with no later `get_results`. Recover `run_id` from `get_status` as `RECOVERY.md` / `REFERENCE.md` describe.
3. **Inherited namespace bloat.** A resumed sandbox holding a large leftover `adata` / `loader` before the current task. Report `resource_usage` against the next store. Do not delete objects until the user agrees.
4. **`/data/session` budget pressure.** Check free space. The writable-budget number and the `ENOSPC` guidance live in `REFERENCE.md` — read them there; do not copy the limits into this report besides what `df` just showed.
5. **Orphaned checkpoints.** Files under `/data/session` that no current variable points at. List paths. Do not delete them unless the user asks.
6. **Stale or dropped mounts.** `Transport endpoint is not connected`, or a mount path that is empty after a reset. Re-mount is the fix; see `MOUNTING.md`.
7. **Local CLI + sandbox overlap.** Both a working local `elembio` and a live sandbox, used for the same listing/download. Report it:

   "You have `elembio` locally and a sandbox. For listing and small downloads the local CLI is faster. The sandbox is the compute path. I will keep using [whichever the routing table says / you already chose] unless you want to switch."

   Do not call both for the same listing. Do not restate the routing table.

### Report shape

```
Audit results
=============
Extra sandboxes:     <n>  <one-line evidence>
Detached runs:       <n>  <run_id if known>
Namespace bloat:     yes | no
Session disk:        <used / avail from df>
Orphan checkpoints:  <paths or none>
Mount issues:        <none or path>
CLI overlap:         yes | no

Recommended next step: <one sentence>
```

Ask before you destroy a sandbox, delete checkpoints, or unmount.

## Mode 3: Diagnose

**Trigger:** "what's broken?", "the sandbox isn't working", "debug my sandbox", or when a specific tool call has failed.

This mode is a **router**. Do **not** restate the signal table, the socket-state matrix, or the provisioning stopping rule. Those go stale. Open `RECOVERY.md` and apply it.

### Steps

1. Gather the failing call, the `sandbox_id` you hold, and any error text from this conversation.
2. Open `elembio-skill://sandbox/sandbox-environment/RECOVERY.md` (or `read_skill` with that URI).
3. Apply it in order: provisioning / restore first, then named signals, then ambiguous timeout, then filesystem / mount.
4. Report findings in plain language:

   "Here's what I found:
   - **Connection:** …
   - **Sandbox:** … (id, alive / provisioning / missing)
   - **Cause:** … (one sentence from RECOVERY.md, not the raw error)
   - **Action:** … (what I will do, or what you need to do)"

5. If the fix is user-side (reconnect MCP, install multiomics, reload the client), send them to **sandbox-onboard**. If the fix is in-session (poll, remount, reload a checkpoint), do it after you say so.

### What diagnose must not do

- Recreate on timeout without a `RECOVERY.md` verdict.
- Recreate on `INSTANCE_PROVISIONING`.
- Destroy a sandbox to "clear the error".
- Paste the recovery table into the user-facing reply.

If the problem is beyond the sandbox (ElemBio Cloud outage, empty project, no registered storage), say so and stop. Empty listing is **No reachable data** in `sandbox-onboard`, not a dead kernel.

## Gotchas

- **No `sandbox_id` is not "the MCP is down".** It means this chat has not created one yet. Health-check that as "none this session".
- **Do not create a sandbox to audit waste.** Creation is the waste you would be looking for.
- **`get_status` is not a namespace listing.** Use `execute_code` (`dir()`, `'loader' in globals()`) when you need to know what survived. That distinction is in `RECOVERY.md`.
- **Diagnose duplicates recovery if you copy the table.** Cite the file. Apply it. Move on.

## Tone

Plain language. "The kernel is still running; that call was just slow" beats a stack trace. Ask before anything destructive.
