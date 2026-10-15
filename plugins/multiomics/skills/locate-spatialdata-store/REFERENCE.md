# Locate SpatialData Store — Reference

## SDO sources

SDOs are always `.zarr.zip` files.

| Source | Location | Contents |
|---|---|---|
| Instrument (run) | `{run}/Cytoprofiling/Instrument/SpatialData/<name>.zarr.zip` | Imaging ± OPS; 3′ points if sequenced (no counts table) |
| Execution (`cells2stats`, `elembio-multiomics-nf`, …) | Under the execution output — discover generically via `elembio executions files <exec-id>` (look for `*.zarr.zip`) | Depends on the workflow; the multiomics workflow adds assigned 3′ counts (`ThreePrimeUntargeted`) |

## Common table names

| Table | Modality | Notes |
|---|---|---|
| `imaging` | CellPaint / ProteinIF | Present in all instrument SDOs |
| `ThreePrimeUntargeted` | 3′ transcriptome | Only in SDOs where the 3′ workflow ran |
| `B01-*` (batch names) | OPS / targeted panels | Panel-specific; names vary by experiment |

## Local discovery (no cloud)

When the user works from runs on their own disk (no `elembio` cloud CLI):

- **There is no env var or default path** for a local runs location. Resolve it **memory-first → ask → remember**: use a runs directory you already know from memory/context; otherwise ask the user for the run directory; then save it to memory for next time (if the runtime has memory).
- Search that location with the [directory search](./IMPLEMENTATION.md#local-run-discovery-no-cloud-cli) mechanics (`find_stores`, run-name match), preferring the instrument SDO path.

The aligned-3′-missing decision (proceed without, or supply the `elembio-multiomics-nf` SDO) is owned by the [SKILL.md](./SKILL.md#if-aligned-3-is-missing--call-out-and-confirm) gate and detailed in [IMPLEMENTATION.md](./IMPLEMENTATION.md#4-verify-and-hand-off) — not repeated here.

## CLI reference

CLI mechanics (install, auth, list, mount, download) are owned by [elembio-cloud-data-access](../elembio-cloud-data-access/SKILL.md). This skill delegates to it when the data is in the cloud.
