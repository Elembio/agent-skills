# ElemBio agent skills — preview

Preview channel for Element Biosciences plugins and agent skills, for Claude Code and Cursor.

Content here is **not released yet**. It graduates to the public `elembio` marketplace via a
deliberate publish; this repo is where it is authored and tested first.

Marketplace name: **`elembio-preview`**

## Install

You need read access to this repo. Claude Code and Cursor authenticate with your existing git
credentials — if `gh auth status` works, so will this.

### Claude Code

```
/plugin marketplace add Elembio/agent-skills-preview
/plugin install cloud-sandbox-preview@elembio-preview
```

You will be prompted for your **ElemBio API key** (production, prefix `ebp_`). Generate one from
your account settings at [cloud.elembio.io](https://cloud.elembio.io). Claude Code marks the field
`sensitive`, so it is stored in your system keychain rather than in `settings.json`.

For background auto-updates on this private repo, make sure git has a credential helper:

```
gh auth setup-git
```

### Cursor

Cursor has no CLI or slash command for this — marketplaces are added through the dashboard.

Add the marketplace (once per team):

1. **Dashboard → Plugins**, then **Add Marketplace** under *Team Marketplaces*.
2. Choose **Import from Repo** and give it `Elembio/agent-skills-preview`.
3. Under **Marketplace Settings**, set access and save. Enable **Auto Refresh** if you want new
   commits picked up automatically — this requires the **Cursor GitHub App** to be installed on
   this repo, and Cursor re-indexes at most once every 10 minutes.

Install the plugin:

1. Open **Customize** in the sidebar and find `cloud-sandbox-preview`.
2. **Install**, choosing project or user scope.
3. Enter your **ElemBio API key** when prompted for `ELEMBIO_API_KEY` — same production key as
   above, from [cloud.elembio.io](https://cloud.elembio.io).

Cursor auto-discovers `mcp.json` at the plugin root, so no manifest wiring is needed.

**Cursor has no secure-storage flag.** Its `variables` block supports only standard JSON Schema
keywords — there is no `sensitive` or `secret` equivalent to Claude Code's. A key configured in
Cursor is therefore not keychain-backed. Treat it accordingly.

## What is here

| Plugin | Contains | Becomes on publish |
| --- | --- | --- |
| `cloud-sandbox-preview` | ElemBio Cloud sandbox MCP (production: `mcp.usw2.elembio.io/sandbox`) | `cloud-sandbox@elembio` |

Not here yet:

- **The skills pack** — name still being settled.
- **`cloud`** — the `elembio` CLI data-access skill. This is a *separate* plugin from the sandbox:
  the two share only the API-key requirement, one is remote compute while the other reaches your
  data from wherever you already work, and a customer can want either without the other.

## Conventions

### The MCP server key is not the plugin name, deliberately

The plugin is named `cloud-sandbox-preview`, but the MCP server it defines is `elembio-sandbox`:

    { "mcpServers": { "elembio-sandbox": { "url": "..." } } }

**Do not "fix" this.** Server keys become part of every tool name the agent sees
(`mcp__elembio-sandbox__execute_code`); plugin names do not appear in tool names at all. Holding the
key fixed is what keeps tool names stable as this plugin moves from `elembio-marketplace-internal`
to preview to the public marketplace. Renaming it silently renames every tool and breaks anything
that references them.

### Versioning

**Bump a plugin's `version` for any change you want people to receive.** Claude Code decides whether
to update from the version field, not from the git commit, so a fix pushed without a bump reaches
nobody. Cursor's update trigger has not been verified here — bump regardless.

There are two independent version numbers. They do not track each other:

| Field | Meaning |
| --- | --- |
| `metadata.version` in `marketplace.json` | version of the marketplace manifest itself |
| `version` in a plugin's `plugin.json` | the plugin's released version — the one users receive |

### Two harnesses, one tree

`.claude-plugin/` and `.cursor-plugin/` hold the same entries in each harness's schema, at both the
marketplace and the plugin level. Change one, change the other.

Known schema differences, so they are not mistaken for drift:

| | Claude Code | Cursor |
| --- | --- | --- |
| Config declaration | `userConfig` | `variables` (JSON Schema object) |
| Secure storage | `sensitive: true` | not supported |
| MCP config | `mcpServers: "./claude-mcp.json"` | `mcp.json` auto-discovered at plugin root |
| Value reference | `${user_config.KEY}` | `${KEY}` |

No `$schema` field is set on any manifest: no schema URL is published for either harness, and
neither `anthropics/life-sciences`, `anthropics/claude-plugins-community` nor
`hashicorp/agent-skills` sets one.
