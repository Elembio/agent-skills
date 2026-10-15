# ElemBio Cloud Data Access — Reference

## Table of Contents

1. [CLI Commands](#cli-commands)
2. [Download Flags](#download-flags)
3. [Mount Flags](#mount-flags)
4. [Environment Variables](#environment-variables)
5. [Global Flags](#global-flags)
6. [Output Formats](#output-formats)
7. [Filter Syntax](#filter-syntax)

## CLI Commands

| Command | Description |
|---|---|
| `elembio login` | Authenticate with an API key |
| `elembio logout` | Clear stored credentials |
| `elembio whoami` | Show current authenticated identity |
| `elembio version` | Print CLI version |
| `elembio update` | Update the CLI to the latest release |
| `elembio completion <shell>` | Generate shell completion (`bash`, `zsh`, `fish`, `powershell`) |
| `elembio mount list` | List live and stale FUSE mounts started by the CLI |
| `elembio unmount <mountpoint>` | Unmount a single path (`--force` to clear stale mounts) |

### Runs

| Command | Description |
|---|---|
| `elembio runs list` | List runs, newest first (default 100; `--max-items N` to change) |
| `elembio runs get <run-id>` | Get run details |
| `elembio runs files <run-id>` | List files in a run |
| `elembio runs download <run-id>` | Download run files |
| `elembio runs mount <run-id> <mountpoint>` | Mount run files as a local read-only directory |

### Executions

| Command | Description |
|---|---|
| `elembio executions list` | List executions, newest first (default 100; `--max-items N` to change) |
| `elembio executions get <exec-id>` | Get execution details |
| `elembio executions files <exec-id>` | List execution files |
| `elembio executions download <exec-id>` | Download execution files |
| `elembio executions logs <exec-id>` | View execution logs |
| `elembio executions mount <exec-id> <mountpoint>` | Mount execution files as a local read-only directory |

### Instruments

| Command | Description |
|---|---|
| `elembio instruments list` | List all instruments |
| `elembio instruments get <serial>` | Get instrument details |

### Storage

| Command | Description |
|---|---|
| `elembio storage list` | List storage connections |
| `elembio storage get <conn-id>` | Get connection details |
| `elembio storage files <conn-id>` | List files in a storage connection |
| `elembio storage download <conn-id>` | Download files from a storage connection |
| `elembio storage mount <conn-id> <mountpoint>` | Mount storage connection as a local read-only directory |

## Download Flags

Shared across `runs download`, `executions download`, and `storage download`.

| Flag | Default | Description |
|---|---|---|
| `--output-dir, -o` | `./<resource-id>/` | Local directory for downloads |
| `--prefix` | *(all files)* | Only download files under this key prefix |
| `--include` | `*` | Glob pattern for files to include |
| `--exclude` | *(none)* | Glob pattern for files to exclude |
| `--parallel, -p` | `4` | Number of parallel downloads |
| `--dry-run` | `false` | List files without downloading |
| `--overwrite` | `false` | Overwrite existing files |
| `--skip-existing` | `false` | Skip files that already exist locally |
| `--checksum` | `false` | Verify checksums after download |

## Mount Flags

Shared across `runs mount`, `executions mount`, and `storage mount`.

| Flag | Default | Description |
|---|---|---|
| `-d, --detach` | `false` | Daemonize after FUSE is live; default can be flipped via `ELEMBIO_MOUNT_DETACH=1` |
| `--log-file` | *(auto)* | Custom log path for detached mounts (not discoverable by `elembio mount list`) |
| `--detach-timeout` | `15s` | Parent's wait for the detached child to signal readiness |
| `--prefix` | | Only mount files under this key prefix (include trailing `/` for directory-style prefixes) |
| `--include` | | Glob pattern for files to include (repeatable) |
| `--exclude` | | Glob pattern for files to exclude (repeatable) |
| `--cache-size` | `256MB` | Max in-memory L1 cache size |
| `--block-size` | `4MB` | HTTP range request block size |
| `--cache-dir` | `$TMPDIR/elembio-mount-cache` | Persistent disk cache directory |
| `--disk-cache-size` | `50GB` | Max disk L2 cache size (`0` to disable) |
| `--allow-other` | `false` | Allow other users to access the mount |
| `--debug-fuse` | `false` | Enable FUSE debug logging |

## Environment Variables

| Variable | Description |
|---|---|
| `ELEMBIO_API_KEY` | API key (overrides stored config) |
| `ELEMBIO_API_URL` | API endpoint override |
| `ELEMBIO_CONFIG` | Config file path |
| `ELEMBIO_CONFIG_DIR` | Config directory path |
| `ELEMBIO_MOUNT_DETACH` | When `1`, makes `--detach` the default for mount commands |

## Global Flags

Available on all commands.

| Flag | Description |
|---|---|
| `--api-url` | Override API endpoint |
| `--api-key` | Override API key |
| `-o, --output` | Output format (`table`, `json`, `yaml`, `wide`) |
| `--no-headers` | Omit table headers |
| `-q, --quiet` | Suppress non-essential output |
| `-v, --verbose` | Verbose output |
| `--debug` | Debug output (includes HTTP traces) |
| `-h, --help` | Help for command |

## Output Formats

All list commands support multiple output formats via `-o`:

| Format | Description |
|---|---|
| `table` | Default columnar output |
| `wide` | Extended table with additional columns |
| `json` | JSON output |
| `yaml` | YAML output |

Use `--no-headers` with `table` or `wide` for script-friendly output.

## Filter Syntax

List commands accept `--filter` with `field:value` pairs. Multiple filters in a single string combine with AND logic.

**Quoting.** An unquoted text value is a case-insensitive substring match (`value*` starts with, `*value` ends with). A double-quoted value is an **exact, case-sensitive** match; wrap the whole filter in single quotes so the shell keeps the double quotes: `--filter 'name:"RUN-0001"'`. ID keywords (`id`, `storage.connection_id`, `input.run.id`, `input.execution.id`) and enum keywords (`status`, `type`, `platform`, …) always match exactly.

**Result order.** Newest first and not configurable: runs by `time_updated`, executions by `time_created`.

Filter **keywords differ between resources** — an unknown keyword or an out-of-range value makes the CLI print the valid set, so probe with a bogus value (e.g. `--filter type:ZZ`) when unsure. Values verified against `elembio-cli/0.2.2`.

**Runs** (`elembio runs list`) — keywords include `name`, `id`, `type`, `status`, `instrument.*`, `storage.connection_id`, `tags`, `time_created`, …

| Example | Meaning |
|---|---|
| `status:completed` | Runs with status "completed" (valid: `completed, created, discarded, failed, paused, pausing, planned, primed, priming, queued, running, sequenced, sequencing, stopped, stopping`) |
| `type:multiomics` | Cytoprofiling / multiomics runs (valid `type`: `multiomics`, `sequencing`, `wash`) |
| `name:"<run-name>"` | Run by exact name (case-sensitive) |
| `name:<partial>` | Runs whose name contains the value (case-insensitive) |
| `status:failed time_created>=7d` | Failed runs created in the last 7 days |

**Executions** (`elembio executions list`) — valid keywords: `flow.name`, `input.execution.id`, `input.execution.name`, `input.run.id`, `input.run.instrument_name`, `input.run.name`, `name`, `output.kind`, `platform`, `started_by`, `status`, `time_completed`, `time_created`, `time_started`, `time_updated`, `workflow.name`.

| Example | Meaning |
|---|---|
| `input.run.id:cto_abc123` | Executions for a specific run, by run id (exact) — **not** `run.id` or `run_id` |
| `input.run.name:"<run-name>"` | Executions for a specific run, by exact run name |
| `workflow.name:cells2stats` | Executions produced by a given workflow (substring) |
| `platform:catalyst` | Executions run on a platform (valid: `aws_omics`, `catalyst`, `dnanexus`, `seqera`) — Catalyst runs cytoprofiling workflows |
| `output.kind:parquet` | Executions by output kind (valid: `bam`, `cram`, `fastq`, `parquet`, `vcf`) |
| `status:failed time_created>=7d` | Failed executions created in the last 7 days (valid `status`: `canceled, cancelled, completed, failed, initialized, pending, running, submitted, succeeded`) |
