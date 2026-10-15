# Element Biosciences agent skills & plugins

Public marketplace of Element Biosciences plugins and agent skills for **Claude** (Desktop
and Claude Code) and **Cursor**.

Marketplace name: **`elembio`**

## Plugins

| Plugin | What it does |
| --- | --- |
| `cloud-sandbox` | Connects your agent to the **ElemBio Cloud sandbox MCP** — remote compute plus your Element Biosciences cloud data for multiomics, QC, and spatial analysis. Sign in with your ElemBio account via browser (OAuth). Production endpoint: `mcp.usw2.elembio.io/sandbox`. |
| `multiomics` | **Analysis skills** for AVITI24 multiomic and spatial data: loading, run and cell QC, normalization, clustering, differential expression, and the 3′ transcriptome, OPS, and imaging pipelines. The code they describe uses the `elembio-spatialdata-analysis` Python package. See [plugins/multiomics](./plugins/multiomics/README.md). |

### What `cloud-sandbox` is for

Use it when you want to **explore or analyze your Element Biosciences / AVITI run data** and
need compute for it:

- List and resolve runs, executions, and storage connections; download or mount your data.
- Run multiomics, single-cell, spatial, imaging, OPS, QC, or differential-expression analysis
  in a remote Python kernel with the analysis stack (`spatialdata`, `scanpy`, `anndata`,
  `squidpy`, …) and `elembio-cli` **preinstalled**.

The sandbox acts as you (the signed-in user) and mounts your cloud data in place — no copying,
no local setup. Everything runs in the sandbox, from listing runs to the full analysis, so you
don't need the `elembio` CLI installed locally. The bundled `sandbox-environment` skill tells
the agent how to use the sandbox for these scenarios, and hands the analysis itself off to the Element
Biosciences `multiomics` skills. The sandbox provides those skills, so you don't need to install
them separately. Install the `multiomics` plugin only for analysis on your own compute.

## Install

### Claude Desktop (and claude.ai)

1. **Customize → Plugins → `+` → Add from a repository**, and give it `Elembio/agent-skills`.
2. Install **`cloud-sandbox`** to analyze in a sandbox, or **`multiomics`** to analyze on your
   own compute.
3. For `cloud-sandbox`, open the connector and **Connect** to complete the browser sign-in
   (OAuth).

The bundled skills then load automatically and are used when relevant.

### Claude Code

```
/plugin marketplace add Elembio/agent-skills
/plugin install cloud-sandbox@elembio
```

For analysis on your own compute, install `multiomics` instead:

```
/plugin install multiomics@elembio
```

For `cloud-sandbox`, authenticate on first use of a sandbox tool — Claude Code opens the
browser sign-in automatically (`/mcp` shows connection status).

### Cursor

Marketplaces are added from the dashboard:

1. **Dashboard → Plugins → Add Marketplace** (under *Team Marketplaces*), choose **Import from
   Repo**, and give it `Elembio/agent-skills`.
2. Open **Customize**, find `cloud-sandbox` (sandbox) or `multiomics` (your own compute), and
   **Install**.
3. For `cloud-sandbox`, connect the `elembio-sandbox` MCP server and complete the browser
   sign-in.

## Authentication

`cloud-sandbox` uses your **ElemBio account** via browser OAuth — there is no API key to paste.
Generate/manage your account access at [elembio.io](https://elembio.io). Nothing secret is
stored in this repository.

## Disclaimer

This content is provided "as-is" without warranty of any kind. AI-generated output may contain
inaccuracies. Review and verify all outputs before relying on them.

## License

Licensed under the Apache License, Version 2.0. See [LICENSE](./LICENSE) and [NOTICE](./NOTICE).
