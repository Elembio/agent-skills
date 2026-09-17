---
name: sandbox-onboard
description: Onboard a new user to the ElemBio Cloud sandbox MCP — introduce what it can do, walk through authentication, install the companion multiomics plugin, and route into the right flow based on the state of their setup. Use when getting started, troubleshooting connection issues, or when the user asks "what is the elembio sandbox", "onboard elembio", "how do I get started with the sandbox", "set me up", "what can I do with the sandbox", or "tell me about ElemBio Cloud sandbox".
metadata:
  version: 0.1.0
  author: elembio
---

# Sandbox onboard

Introduce the ElemBio Cloud sandbox, get the user authenticated, then guide them through the setup their client is actually in.

Runtime, mounts, and recovery live on the server. Do not restate them here. Open `elembio-skill://sandbox/sandbox-environment/SKILL.md` (or `read_skill` with that URI) when you need the how.

## Step 1: Introduction

Start by describing what the sandbox can do, then check the connection.

### Pitch

"The ElemBio Cloud sandbox is a remote Python environment with your AVITI / AVITI24 data in reach. It ships the analysis stack and `elembio-cli` already signed in as you. You mount a run in place, then ask questions in this chat. No local GPU, no copy of a multi-GB store onto your laptop."

### Check connection

Inspect whether any `elembio-sandbox` tools are available:

- **Tools are available.** Give a shorter pitch — "The sandbox MCP is connected. Let me check the rest of the setup." — then go to Step 2.
- **No sandbox tools at all.** Try `mcp_auth` on the ElemBio sandbox server. If that succeeds, re-check tools and go to Step 2.

  If `mcp_auth` fails or is unavailable, use the client-specific connect steps in **Branch: Not connected**. Detect the client from the environment. If unclear, give the generic instructions.

  Wait for the user to confirm ("done"), then re-check tools.

Sign-in is browser OAuth against the user's ElemBio account. There is no API key to paste.

## Step 2: Diagnose

Inspect the sandbox tools and the installed skills. The result picks the branch:

| Result | Branch |
| --- | --- |
| Sandbox tools work, multiomics skills are present, and at least one run or storage connection is reachable | **Healthy** |
| Sandbox tools work, but the multiomics skills are absent | **Companion missing** |
| Sandbox tools work and multiomics is present, but listing runs and storage returns nothing | **No reachable data** |
| Tools exist, but calls return `not authorized for session` or equivalent auth failure | **Auth broken** |
| No `elembio-sandbox` tools at all | **Not connected** |

Detect multiomics by looking for its `index` skill (also `spatialdata-loading-and-access`, `elembio-cloud-data-access`). Absence of those names is **Companion missing**, not a sandbox fault.

## Branch: Healthy

The server is connected, the companion plugin is installed, and data is in reach. Show a summary and offer next steps.

1. Report, in one block:
   - Signed-in identity (`elembio whoami` locally, or via `execute_command` in an existing sandbox — do not create a sandbox only to print whoami)
   - Multiomics present
   - Local `elembio` CLI present or not (`which elembio && elembio whoami`)
   - How many runs or storage connections you could see
2. Offer options:
   - "Try one live" → **sandbox-demo**
   - "What can I ask of a run" → **sandbox-explore**
   - "Run a health check" → **sandbox-status**
   - Or just start using the sandbox

## Branch: Auth broken

The server is configured, but authentication has expired or is invalid.

1. Tell the user the sandbox is configured but the connection is broken.
2. Give the same reconnect steps as **Not connected** for their client.
3. Wait for "done".
4. Retry a cheap sandbox tool.
5. If it works, re-diagnose. If it still fails, suggest removing and re-adding the MCP server in the client, then reconnect.

## Branch: Not connected

The plugin is installed, but the MCP server has not been authenticated. This is the common state on a fresh install — zero sandbox tools are visible.

1. Tell the user the plugin is installed but needs a connection.
2. Try `mcp_auth` on the ElemBio sandbox server. If that succeeds, skip to step 5.
3. If `mcp_auth` fails or is unavailable, give client-specific steps:

   - **Cursor:** Settings → Cursor Settings → Tools & MCP → **Connect** next to `elembio-sandbox`. Or press Cmd+Shift+P and search for "MCP".
   - **Claude Desktop:** Customize → Connectors → ElemBio Cloud sandbox → **Connect**. Complete the browser sign-in.
   - **Claude Code:** The client opens the browser sign-in on first use of a sandbox tool. `/mcp` shows connection status.
   - **Other clients:** Find `elembio-sandbox` in the client's MCP settings and connect it. The client redirects to ElemBio account sign-in.

4. Wait for "done".
5. Re-diagnose. Most often the next branch is **Companion missing** or **Healthy**.

## Branch: Companion missing

The sandbox is connected, but the analysis skills are not installed. Analysis in the sandbox is driven by the Element Biosciences `multiomics` plugin. Do not substitute ad-hoc Scanpy.

The public `elembio` marketplace ships `cloud-sandbox`. The analysis skills currently ship as `multiomics-preview` from `Elembio/agent-skills-preview`.

1. Tell the user they are connected, and that analysis guidance is a separate plugin.
2. Give the install for their client:

   - **Claude Code:**

     ```
     /plugin marketplace add Elembio/agent-skills-preview
     /plugin install multiomics-preview@elembio-preview
     ```

   - **Claude Desktop:** Customize → Plugins → `+` → Add from a repository → `Elembio/agent-skills-preview`. Then install **multiomics-preview**.
   - **Cursor:** Dashboard → Plugins → Add Marketplace → Import from Repo → `Elembio/agent-skills-preview`. Then install **multiomics-preview**.

3. Tell them to reload the client so the new skills appear (see **Reload instructions by client**).
4. Wait for "done", then re-check for the `index` skill.
5. If it is present, re-diagnose (usually **Healthy** or **No reachable data**). If it is still missing, they likely skipped the reload.

## Branch: No reachable data

The sandbox is connected and the companion plugin is installed, but listing runs and storage is empty.

Do not restate `elembio` CLI verbs. Hand off to the multiomics **`elembio-cloud-data-access`** skill for discovery, registration, and archive restore.

Tell the user:

- The sandbox can run, but it has nothing of theirs to mount yet.
- Runs and storage connections must exist in ElemBio Cloud first.
- A raw `s3://` URI works only after that bucket is a registered storage connection.

Then stop. Do not invent a demo against empty storage.

## Reload instructions by client

New skills and newly authenticated MCP tools do not appear until the client re-reads its config.

| Client | How to reload |
| --- | --- |
| Cursor | Cmd+Shift+P → "Reload Window" |
| Claude Desktop | Quit and reopen the app |
| Claude Code | Run `/mcp` to check status. Restart the session if tools or skills are still missing |
| Windsurf | Cmd+Shift+P → "Reload Window" |

## MCP config by client

Use this only when the user needs to inspect or recreate the server entry.

| Client | Config file location | Scope |
| --- | --- | --- |
| Cursor | `.cursor/mcp.json` (project) or `~/.cursor/mcp.json` (global) | Project / Global |
| Claude Desktop | `~/Library/Application Support/Claude/claude_desktop_config.json` (macOS) | Global |
| Claude Code | `.mcp.json` (project) or `~/.claude/mcp.json` (global) | Project / Global |
| Windsurf | `~/.codeium/windsurf/mcp_config.json` | Global |

## Gotchas

- **Do not create a sandbox only to diagnose onboard state.** Tool presence, `mcp_auth`, and `which elembio` are enough. Create on demo or first real task.
- **Companion missing is not "no actions on the server".** Sandbox tools are compute. Analysis recipes live in the other plugin. Installing `cloud-sandbox` a second time will not fix it.
- **Skills added mid-session are invisible until reload.** If the user says they installed multiomics and you still cannot see `index`, ask for a reload before you reinstall.
- **Do not suggest `sandbox-status` while disconnected.** Status needs tools.

## Tone

Casual and efficient. Do not explain MCP or protocol details. Get them to the right place. If something breaks, be direct: "That did not work. Let's try…"
