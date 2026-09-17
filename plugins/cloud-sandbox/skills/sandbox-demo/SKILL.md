---
name: sandbox-demo
description: Walk a new user through the smallest possible first win in the ElemBio Cloud sandbox — one run, one mount, Loader plus available_tables, shown live. Use when the user asks "show me how the sandbox works", "give me a quick demo", "I want to try it", "what's the fastest way to see this work", "minimal setup", "hello world", "smallest example", or "onboard elembio demo".
metadata:
  version: 0.1.0
  author: elembio
---

# Sandbox demo

Walk a new user through the smallest first win: one run, mounted in place, tables listed. The whole flow should feel quick — a few minutes from "I'm curious" to "oh, that worked."

This is the natural next step after `sandbox-onboard` for users who want to see the sandbox work before they start analysis.

The end-to-end sandbox path lives on the server. Follow `elembio-skill://sandbox/sandbox-environment/SKILL.md` (Quickstart) and `MOUNTING.md` for how to create, poll, and mount. This skill adds run selection, the SDO trap, narration, and a **cost ceiling**. Do not copy those files into the chat.

## When to use vs. other skills

- **sandbox-demo** (this skill) — one run, one mount, `available_tables`. Fastest path to "oh, this works."
- **sandbox-onboard** — pitch + connect + companion plugin. Run this first if the server is not connected, or if multiomics is missing.
- **sandbox-explore** — questions this run can actually answer. Run this after the demo.
- **sandbox-status** — health checks and audits on an existing session.

If the user says "show me," lean here. If they have not connected yet, route to `sandbox-onboard` first.

## Cost ceiling (the demo rule)

This is the corollary to "never a write action for a first Zapier demo."

`Loader` plus `available_tables` returns in seconds, even on a large `.zarr.zip`. The first `load_tables` is tens of seconds to minutes. **Stop before `load_tables`.** Do not pass over `X`. Do not cluster, normalize, or plot.

If you want a count, take it from metadata the `Loader` already exposed (table names, shapes in `available_tables`). Do not materialize a table to prove the sandbox works.

Cost numbers and the rest of the runtime live in the server-published `REFERENCE.md`. Read them there. Do not restate the table here.

## Step 1: Confirm the server and the companion

Inspect available `elembio-sandbox` tools. If none exist, authenticate as `sandbox-onboard` does. Do not continue until sandbox tools are available.

Inspect for the multiomics `index` skill. If it is missing, route to `sandbox-onboard`'s **Companion missing** branch. The demo can list tables without it, but the handoff after the win cannot.

Set the tone:

> "Let's get one run into the sandbox so you can see this work. We will pick a run, mount it, and list the tables — should take a couple of minutes."

## Step 2: Create one sandbox

Follow the server-published Quickstart: `create_sandbox` once, keep the `sandbox_id`, poll `get_status` through `INSTANCE_PROVISIONING`. Do **not** create a second sandbox because provisioning is slow. That rule is in the lifecycle rule and in `RECOVERY.md`.

Narrate before you create ("Starting a sandbox — first boot can take a couple of minutes") and after it is `INSTANCE_READY`.

## Step 3: Pick a run

List runs with the local `elembio` CLI when it is present. Otherwise list from inside the sandbox via `execute_command`. Do not restate CLI flags; the multiomics **`elembio-cloud-data-access`** skill owns that surface.

Lead with a short list the user can react to — run name, id, date — rather than asking them to recall an id cold.

### Avoid several SDOs on a first demo

A run can ship more than one SpatialData object (instrument-direct, `cells2stats`, `elembio-multiomics-nf`). Picking among them is the exact friction this demo is trying to skip.

- Prefer a run that has a single obvious store.
- If every candidate has several SDOs, say so and pick the analysis output (`cells2stats` / `elembio-multiomics-nf`) only if the user does not object. Do not stop to teach SDO selection here — that belongs to the multiomics `locate-spatialdata-store` skill.

If listing is empty, stop and route to `sandbox-onboard`'s **No reachable data** branch.

## Step 4: Mount, then load metadata only

Mount the chosen run as `MOUNTING.md` requires. The mandatory `--disk-cache-size 0` and the conventional `/runs/<run-id>` path live there — follow that file, do not invent a copy recipe.

Then, via `execute_code`:

```python
from elembio_spatialdata_analysis.load import Loader
loader = Loader("<the .zarr under /runs/<run-id>>")
loader.available_tables
```

Narrate before the `Loader` call ("Opening the store index — this is metadata, not the matrices") and after, with the table names.

If `Loader` needs a nested path inside a `.zarr.zip`, follow `MOUNTING.md` rather than extracting the zip into `/data/session`.

## Step 5: Name the win and hand off

When `available_tables` returns, name the win:

> "There you go — the sandbox is working. I mounted [run] and the store has [table list]. Same pattern is how we start QC, clustering, or imaging later. We have not loaded any matrices yet; that is the next, slower step, and only when you ask."

Then offer one next move, not a buffet:

> "Want me to suggest questions this run can actually answer? Say so and I will run **sandbox-explore**. Or name a table and we start analysis through the multiomics skills."

## Progress checklist

- [ ] Server connection verified (or freshly authenticated)
- [ ] Companion `index` skill present, or user sent to onboard
- [ ] One sandbox created; provisioning polled, not recreated
- [ ] Run picked; several-SDO runs avoided when possible
- [ ] Mount followed `MOUNTING.md`
- [ ] `Loader` + `available_tables` only — no `load_tables`, no `X`
- [ ] Win named, explore offered

## Gotchas

- **Never call `load_tables` in this skill.** That is the first expensive touch. The demo is over at `available_tables`.
- **Do not create a second sandbox on `INSTANCE_PROVISIONING`.** The id is already yours. Poll `get_status`.
- **Do not pick a several-SDO run "to be thorough".** The selection prompt is the friction we are removing.
- **Do not copy the store into `/data/session`.** Mount and read in place. See `MOUNTING.md`.
- **Do not dump the Quickstart, the cost table, or recovery signals into the user-facing copy.** Cite the server files; keep this chat short.
- **If the mount or `Loader` fails, do not improvise.** Open `RECOVERY.md` / `MOUNTING.md` and apply them. Then come back.

## Tone

Friendly, low-pressure, action-oriented. This is someone's first impression. Say "sandbox" and "run", not "MCP server" and "FUSE" unless they bring those words up.

If something breaks, do not apologize at length. Say what to do next.
