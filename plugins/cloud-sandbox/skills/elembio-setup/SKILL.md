---
name: elembio-setup
description: Use when a user wants to install, verify, or troubleshoot Element Biosciences' Claude Code tooling — the multiomics-preview skills plugin and/or the cloud-sandbox MCP connector. Detects what's already installed/configured instead of assuming a clean slate, gives the exact next command for whatever is missing, verifies the sandbox connection with a real tool call rather than trusting that tools merely appear, and can run an end-to-end proof-of-life check (storage connections, recent runs, sandbox create/destroy). Use when the user asks to "set up Element Biosciences tools," "install elembio plugins," "connect my sandbox," or reports the sandbox/MCP not working.
metadata:
  version: 0.1.0
  author: elembio
---

# Element Biosciences Tools — Setup & Verification

## Purpose

Two independent plugins make up the Element Biosciences Claude Code experience, distributed differently on purpose — don't collapse them into one flow or assume one implies the other:

|                                   | multiomics-preview                        | cloud-sandbox                                                      |
| --------------------------------- | ----------------------------------------- | ------------------------------------------------------------------ |
| Repo                              | private, `agent-skills-preview`           | public, `agent-skills`                                             |
| Distribution                      | Enterprise Teams marketplace (org-pushed) | self-service                                                       |
| `/plugin marketplace add` needed? | No — already known to every org member    | Yes                                                                |
| `/plugin install` needed?         | Yes (one command)                         | Yes (one command)                                                  |
| Secret                            | none                                      | `ELEMBIO_API_KEY` (`ebp_...`, from account settings at elembio.io) |

cloud-sandbox stays self-service deliberately: it carries a remote MCP connector, and org-pushed (team/managed) marketplace distribution would route it through the org's Connector policy, which requires admin-provisioned OAuth — incompatible with the per-user API-key model here. Don't suggest moving it to a team marketplace as a "fix" for install friction.

## Step 1 — detect current state before telling the user anything

Never assume a clean slate or that "not installed" is the problem — check what's actually true first:

- Look for any `mcp__plugin_cloud-sandbox_elembio-sandbox__*` tools already visible to you (loaded or in the deferred-tools listing). Their presence means cloud-sandbox is installed and its MCP server is registered — **not** proof the API key is valid; verify that for real in Step 3.
- Look for multiomics-preview skill names in your available-skills listing. Their presence means multiomics-preview is installed.
- State only what you can actually observe. If you can't tell from what's visible to you, say so rather than guessing either way.

## Step 2 — give the exact next command for whatever is missing

`/plugin` commands are interactive CLI state changes — you cannot run them yourself. Tell the user the exact line to paste and wait for them to confirm it ran, rather than narrating the whole doc at once.

**multiomics-preview missing:**

```
/plugin install multiomics-preview@elembio-preview
```

No marketplace-add step — it's an org-pushed team marketplace, already known to every member. If this reports no access, the fix is a team-membership assignment in the admin console, not a GitHub/repo permission — tell the user to ask IT for team assignment, not repo access.

**cloud-sandbox missing:**

```
/plugin marketplace add Elembio/agent-skills
/plugin install cloud-sandbox@elembio
```

After install, Claude Code should show an interactive prompt for `ELEMBIO_API_KEY`. The user needs a key (prefix `ebp_`) generated from account settings at elembio.io first. If no prompt appears, don't ask the user to paste the key into chat — tell them to run `/plugin configure` and pick cloud-sandbox instead; that's the path confirmed to store the key in the OS keychain, never in `settings.json`.

## Step 3 — verify the connection for real

Tools being visible only means the MCP server is registered, not that the key works. Confirm with an actual call — `list_sandboxes` or `get_status` is the lowest-risk choice (read-only, creates nothing). Interpret the result:

- Succeeds → connected, proceed to Step 4 if the user wants it.
- Auth-shaped error (401 / "invalid or expired API key" / "missing API key") → the key wasn't accepted. Tell the user to run `/plugin configure`, re-enter the key, and retry. Never ask the user to paste the key into chat, and never echo it back if they do anyway.
- No cloud-sandbox tools visible at all → not installed/enabled yet, or a freshly-installed MCP server hasn't registered — ask the user to restart Claude Code (new session) before troubleshooting further.

## Step 4 — end-to-end proof of life (only once Step 3 succeeds, and only if the user asks for it)

Narrate every tool call as you make it — what it does, not just what it returns.

1. `create_sandbox`
2. Discover reachable storage / data connections; summarize as a table (name, type or bucket, purpose if available).
3. List the user's most recent instrument runs, newest first; summarize as a table (run name, instrument, date, run type).
4. If either comes back empty or errors, say specifically which one failed and why — never fabricate or guess placeholder data.
5. `destroy_sandbox` when done — including when a step failed. Don't leave a sandbox running because the check didn't fully succeed.

## Security

- Never print, log, or ask the user to paste `ELEMBIO_API_KEY` into chat. The only supported input path is the native `/plugin configure` prompt (keychain-backed).
- Never write the key into a file, `settings.json`, or a tool-call argument you construct — it belongs only in the `x-api-key` header the MCP client attaches automatically.

## Troubleshooting quick reference

| Symptom                                                 | Likely cause                                                                                        | Fix                                                                            |
| ------------------------------------------------------- | --------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------ |
| `/plugin install multiomics-preview@...` says no access | Not assigned to the right Enterprise team                                                           | Ask IT for team assignment — not repo access                                   |
| cloud-sandbox installs but no key prompt appears        | Known Claude Code flakiness: `userConfig` prompts don't always fire on non-interactive enable paths | Run `/plugin configure`, pick cloud-sandbox                                    |
| Key entered but tools still 401                         | Stale/incorrect key, or the prompt didn't persist it                                                | Re-run `/plugin configure`; regenerate the key at elembio.io if it still fails |
| Tools don't appear after configuring                    | MCP server registration needs a fresh session                                                       | Restart Claude Code                                                            |
