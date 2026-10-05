# ElemBio Cloud Sandbox

Connects your agent to the **ElemBio Cloud sandbox MCP**: remote compute plus your Element
Biosciences cloud data for multiomics, QC, and spatial analysis.

> **Disclaimer:** This content is provided "as-is" without warranty of any kind. AI-generated
> output may contain inaccuracies. Review and verify all outputs before relying on them.

## What it does

Use it to explore or analyze your Element Biosciences / AVITI run data:

- List and resolve runs, executions, and storage connections; download or mount your data.
- Run multiomics, single-cell, spatial, imaging, OPS, QC, or differential-expression analysis
  in a remote Python kernel with the analysis stack (`spatialdata`, `scanpy`, `anndata`,
  `squidpy`, …) and `elembio-cli` preinstalled.

The bundled `sandbox-environment` skill tells the agent how to use the sandbox and hands the
analysis itself to the Element Biosciences `multiomics` skills.

## What it runs, sends, and fetches

- **Nothing runs on your machine.** The plugin contains a skill (instructions for the agent)
  and the configuration for one remote MCP server. It has no hooks, scripts, or local servers.
- **Remote MCP server:** `https://mcp.usw2.elembio.io/sandbox`. The agent sends it tool calls
  (code, commands, and file paths to run in the sandbox) and receives their results and
  output files.
- **Sign-in:** browser OAuth with your ElemBio account. There is no API key, and nothing
  secret is stored in the plugin.
- **Your data:** the sandbox acts as you and reads your ElemBio Cloud data in place. Only data
  your account can already access is reachable.

## Install

See the [repository README](https://github.com/Elembio/agent-skills#install) for Claude
Desktop, Claude Code, and Cursor instructions.

## License

Licensed under the Apache License, Version 2.0. See
[LICENSE](https://github.com/Elembio/agent-skills/blob/main/LICENSE) and
[NOTICE](https://github.com/Elembio/agent-skills/blob/main/NOTICE).
