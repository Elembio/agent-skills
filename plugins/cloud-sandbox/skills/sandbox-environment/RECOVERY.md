# Sandbox — Recovery

What to do when a call times out, errors, or returns an unexpected shape. Detail behind
[SKILL.md](SKILL.md).

Contents: [The one rule](#the-one-rule) · [Provisioning and restore](#provisioning-and-restore-before-the-sandbox-is-usable)
· [Error / signal → action](#error--signal--action) · [Ambiguous timeout / connection](#ambiguous-timeout--connection-error)
· [Filesystem / mount problems](#filesystem--mount-problems)

## The one rule

**A timeout or connection-shaped error does not by itself mean the kernel died.** `create_sandbox`
permanently discards in-memory state (loaded data, variables, installed packages); `/data/session`
files survive but everything else does not. Read the signal first, then act. Waiting out an
ambiguous error costs a minute; a needless recreate can throw away a long analysis. Nor does
`status: "hibernated"` — that sandbox is idle with its state intact, and the next call on the same
`sandbox_id` resumes it.

Two probes, and they answer different questions:

- **`get_status`** — liveness (`status`), resources (`resource_usage`), and mount/kernel verdicts.
  It does **not** describe the Python namespace.
- **`execute_code`** (`dir()`, `'loader' in globals()`) — the only way to check what survived in
  the namespace.

## Provisioning and restore (before the sandbox is usable)

A `status: "provisioning"` or `"resuming"` response (from `create_sandbox` **or** `get_status`)
is a success, not an error. The `sandbox_id` is valid and reserved.

- **Poll `get_status`** with that id, waiting `retry_after_seconds` between polls, until
  `status` is `"ready"`. Do **not** call `create_sandbox` — while `provisioning`, a second call
  starts a second sandbox (and can start a second host); while `resuming`, a new sandbox starts
  empty, without the state being restored.
- **Honor the stopping rule in the payload.** `max_wait_seconds` + `on_timeout` say how long to
  wait and what to do after; count from your *first* poll of that id, not `elapsed_seconds`
  (which restarts on each attempt). When `poll_budget_exhausted: true` (or you pass
  `max_wait_seconds`), stop, tell the user, and do `on_timeout` (`create_sandbox`) — the one
  case where a second create is correct.
- **`status: "resuming"`** means the session's memory image is being pulled back — variables and
  loaded data return with it, so waiting almost always beats starting over. `bytes_total` sizes
  the wait (tens of GiB = minutes). There is deliberately no ETA and no `poll_budget_exhausted`
  here; `max_wait_seconds` is the stop. Tell the user what is being restored if it passes a
  minute ("restoring 46 GiB of session state").
- **`at capacity`** → the fleet is full; tell the user rather than retrying in a loop.
- A session-scoped tool answering with this shape means the sandbox was hibernated and is being
  restored — same rules, then retry the same call on the same `sandbox_id`.

## Error / signal → action

Apply the message literally where it names a fix. "Invalid arguments" answers name the tool,
field, and correction — apply and re-send; they are **not** health signals, so do not
`get_status` or `create_sandbox` on one.

| Signal | What it means | Action | Recreate? |
|---|---|---|---|
| `sandbox_id is required` | dropped parameter, usually | re-send with the id you hold | only if the id is genuinely lost |
| `sandbox X not found` with a `reason` | the sandbox is gone, and `reason` says why: `destroy_requested`, `resume_window_expired`, `resume_attempts_exhausted`, `snapshot_incompatible`, `snapshot_missing`, `restore_failed`, `hibernate_failed`, `host_lost`, `kernel_out_of_memory`, `kernel_crashed`, `provision_failed` | tell the user why; `create_sandbox` only if more work is needed (after `destroy_requested`, only if the user wants it). A new sandbox has an empty kernel **and its own empty `/data/session`**: files the old one wrote are not there | yes |
| `sandbox X not found` with no `reason` | gone, cause not recorded | `create_sandbox`, retry | yes |
| `not authorized for session` | API-key / ownership mismatch | `create_sandbox` | yes |
| `registry lookup failed` | registry didn't answer; no verdict | keep the id, retry in ~10 s | **no** |
| `create_sandbox failed` | provisioning error | wait 15 s, retry once, then surface | (retry once) |
| result has `typed_busy_since_unix_nanos` | your call queued behind another run; never entered the kernel | wait briefly, retry; send fewer back-to-back calls | **no** |
| result has `reason: "kernel_out_of_memory"` / `"kernel_crashed"` | kernel died on this call (terminal) | surface the reason. A new sandbox has the same memory, so after an out-of-memory load less, work in chunks, `del` what you no longer need | user decides |
| result has `state_lost: {reason, detail}` | sandbox was reset (not resumed); this call ran against an **empty kernel**. `reason` says why the snapshot could not be restored (e.g. `snapshot_incompatible`, `snapshot_missing`, `restore_failed`) | stay on same id; re-run mounts, reinstall, reload `/data/session` checkpoint, re-run work; do not report a result that depended on earlier state | **no** |
| unexpected `NameError` / empty namespace (seen via `execute_code`) | state gone, not dead | same id; re-run mounts, reload checkpoint | **no** |
| `get_results` returns `payload_shed: true` | the task finished, but its stdout was dropped to bound memory after later tasks completed. The receipt has `cause`, `error_summary`, `output_paths`, `duration_seconds` and **no `stdout` or `success` field** | report the outcome from `cause` / `error_summary` — a missing `success` is unknown, **not** success; read written files from `output_paths`; re-run only if you need the printed output | **no** |
| result has `outputs_unavailable` | code ran but `/data/session` was unreachable, so nothing there was saved | treat outputs as lost; check `session_storage_state`; re-run once storage is healthy | **no** |
| `status` still `"busy"` after `interrupt_task` | interrupt is best-effort | poll `get_status`; if still busy, surface to user; do **not** `destroy_sandbox` yourself | **no** |

`state_lost` appears **at most once** — on the call whose own resume did the reset — so
act on it immediately rather than waiting for confirmation on the next call.

## Ambiguous timeout / connection error

Call `get_status` and read only the liveness fields:

- **`kernel_status.alive: true`** (even at high `resource_usage.memory_percent` or `status: "busy"`)
  → do **not** recreate. High `%` is a resident-RAM clue, not a death warrant. Recover in place:
  free unused objects, split the work into smaller cells, or raise `timeout_seconds`.
- **`status: "provisioning"` or `"resuming"` with no `kernel_status`** → not booted yet; an absent
  `kernel_status` is **not** `alive: false`. Wait `retry_after_seconds` and poll.
- **`status: "hibernated"`** → idle-hibernated, **not dead**: variables, imports, and loaded data
  are intact. The reply carries only `sandbox_id`, `status`, `region`, and (when known)
  `last_activity_at` — the missing `kernel_status` and `resource_usage` are not a fault. Call
  `execute_code` / `execute_command` on the same `sandbox_id`; it resumes automatically. Never
  `create_sandbox` for it.
- **`status: "kernel_lost"` or `"destroying"`** (or a dead-socket reply: `read response: EOF`,
  `the client session is not running`) → a verdict came back: the kernel cannot work. Tell the
  user the `reason`, then `create_sandbox` and re-run mounts.
- **`kernel_status.alive: false` while `status` is `"ready"` or `"busy"`** → the kernel is alive
  but not answering — usually stuck in its current task. Not a verdict yet: follow
  [Non-OK `socket_state`](#non-ok-socket_state) below (`interrupt_task`, poll), and only
  `create_sandbox` if it stays unresponsive after that.
- **`get_status` itself fails to connect** (`dial guest … after N attempt(s)`, `i/o timeout`,
  `no route to host`) → **no verdict**, only a broken path to a possibly-healthy sandbox. This is
  the one place recreate is the expensive guess. Poll again after ~10 s and keep polling for
  about a minute before concluding it is gone.
- **`status: "snapshotting"`** → the sandbox is refusing new work while it hibernates, moves
  host, or warms up right after a resume; it will come back. Wait and poll; do not replace it.

### Non-OK `socket_state`

`socket_state` names how the kernel socket is behaving. `execute_stalled`, `unanswered`,
`refused`, and `failed` read `responsive: false`; `connect_blocked` stays `responsive: true`
(a saturated backlog, not a wedge).

- **`execute_stalled`** → most often your last cell is *still running* and the sandbox gave up
  waiting first (a long native call into numpy/BLAS, numba, scanpy, or a mounted-store read).
  Read `process_state` on the same reply: **`R`** with `cpu_percent` near 100 is compute that
  will finish — wait; **`S`** with `cpu_percent` ~1–2 is blocked on a mounted store. Poll ~3
  minutes before giving up; do not recreate meanwhile. For the `S` case, try
  `interrupt_task` (non-destructive, reaches a stalled store read) before giving up.
  Re-run with a larger `timeout_seconds` if it was just slow.
- **`connect_blocked`** → usually a briefly saturated kernel that clears; poll again. If
  `stalled_seconds` keeps climbing, treat it as the `execute_stalled` case above (e.g.
  `sc.tl.leiden` on the default `leidenalg` backend holds the interpreter — `flavor="igraph"`
  honors the timeout). Do not use `cpu_percent` to disambiguate here.
- **`unanswered` / `refused` / `failed`** → terminal, no self-recovery. `create_sandbox`,
  re-run mounts.

## Filesystem / mount problems

- **An artifact reads as missing**, or a `/data/session` write fails with `EIO` / the dir reads
  empty → call `get_status` and read `kernel_status.session_storage_state`:
  - `degraded` → the mount is up but S3 is refusing its writes; nothing written to
    `/data/session` is durable right now. `session_storage_detail` carries the real AWS error
    (op, bucket, key, region); `session_storage_failed_files` counts what is stuck. This is
    **not** recoverable from inside the sandbox — report it rather than retrying the write. A
    `degraded` that flips back to `mounted` on the next poll was a transient S3 error.
  - `pending` → mount not attempted yet; wait and retry.
  - absent → no verdict (older worker / not probed); not a fault either way.
  - **Caveat:** `session_storage_state` reports on the *mount*, not on any specific write. To
    confirm a particular file is durable, check its artifact `s3_status == present`, `fsync` it,
    or read it back via `fetch_artifact`.
- **Idle and responsive but nothing finishes** → read `kernel_status.fuse_waiting`. `0` rules
  the filesystem out. Non-zero alone is not "stuck" (a healthy large read has requests in
  flight) — poll again: a count that does not drain while `cpu_percent` stays idle is a wedged
  mount, named by `fuse_waiting_mount`. A confirmed wedge does not self-recover; surface it.
- **`Transport endpoint is not connected`** under a mount path → the FUSE mount dropped; re-mount
  via `execute_command`.
- **`list_files` times out on a large FUSE directory** → read known paths directly, or re-mount
  with a narrower `--prefix`.
