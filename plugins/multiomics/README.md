# ElemBio Multiomics

Analysis skills for Element Biosciences AVITI24 multiomic and spatial data.

> **Disclaimer:** This content is provided "as-is" without warranty of any kind. AI-generated
> output may contain inaccuracies. Review and verify all outputs before relying on them.

## What it does

Skills that guide an agent through AVITI24 analysis: loading SpatialData stores and other AnnData
sources, run and cell QC, normalization, dimensionality reduction, clustering, differential
expression, imaging, and optical pooled screening (OPS). Start from the
[index](./skills/index/SKILL.md) skill, which routes to the right specialist.

The code the skills describe uses the `elembio-spatialdata-analysis` Python package. On your
own compute, install it in the environment your agent uses. In an ElemBio Cloud sandbox it is
preinstalled, and the sandbox provides these skills, so you don't need this plugin there.

## What it runs, sends, and fetches

- **The plugin runs nothing on install.** It contains skills (instructions for the agent) and
  has no hooks or MCP servers.
- **Code the agent writes:** following these skills, the agent writes and runs Python analysis
  code wherever you are working, such as your machine or a notebook.
- **Your data:** the `elembio-cloud-data-access` skill tells the agent how to use the `elembio`
  CLI to list, download, or mount your ElemBio Cloud data, signed in as you.

## Install

See the [repository README](https://github.com/Elembio/agent-skills#install) for Claude
Desktop, Claude Code, and Cursor instructions.

## Changelog

See [CHANGELOG.md](./CHANGELOG.md).

## License

Licensed under the Apache License, Version 2.0. See
[LICENSE](https://github.com/Elembio/agent-skills/blob/main/LICENSE) and
[NOTICE](https://github.com/Elembio/agent-skills/blob/main/NOTICE).
