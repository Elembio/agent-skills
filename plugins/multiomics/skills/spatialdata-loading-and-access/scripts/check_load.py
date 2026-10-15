"""Memory and sanity checks around materializing AVITI24 AnnData tables.

Run this, don't reimplement it. Load it into the Python session with
``exec(open(PATH).read())`` and call:

    check_before(tables)   before any lazy=False / .to_memory() / .compute()
    check_after(tables)    after loading, before handing off

``tables`` is one AnnData or a dict of them (``load_tables(..., concatenated=False)``).
Each call prints one verdict line per table and returns them as a list of dicts.
Any FAIL raises LoadCheckFailed after printing: change the plan, don't work around it.

check_before compares the estimated in-memory footprint of each table with the
memory the session can get:

    room = ceiling - (MemTotal - MemAvailable)

The ceiling is the most memory the session can grow to. On a VM with hot-pluggable
memory, the script reads it from the virtio-mem device, and uses a default elastic
ceiling when the device is present but unreadable; otherwise it uses MemTotal.
Memory is read from /proc/meminfo, or from psutil where that file does not exist.
Set ELEMBIO_MEMORY_CEILING_GIB or pass ceiling_bytes to override the ceiling.

The footprint counts X, layers, obsm, obsp and a fresh copy of obs / var. Sparse
arrays count their stored values (data + indices + indptr); for a lazy sparse array
the non-zero count comes from a few sampled chunks. A 10% margin covers what is not
counted (string objects, read buffers). With several tables, check_before also
checks their total: a default ``load_tables(..., lazy=False)`` concatenates every
run with an outer join on features, and holds the runs and the joined copy at once.

Before the verdict, check_before runs gc.collect() and lists large objects in the
session (``__main__`` or the namespace you pass) that you could release. It never
deletes them. When the load fits only if you release them, the verdict is WARN with
the names to ``del``.

OK: fits. WARN: fits only with a condition (release objects, near the ceiling, or
growth that can lag a fast read). FAIL: does not fit even at the ceiling after
releasing everything listed; subset by run, well or feature.
"""

import gc
import math
import os
import re
import sys

_GIB = 1024**3
# Elastic ceiling to assume when a virtio-mem device is present but its size cannot be read.
_DEFAULT_CEILING = 64 * _GIB
_CEILING_ENV = "ELEMBIO_MEMORY_CEILING_GIB"
# Uncounted overhead: object-dtype strings in obs / var, zarr decode buffers.
_MARGIN = 1.10
# Above this share of the room at the ceiling, WARN: the kernel takes memmap and
# reserve from each grow, and the guest can refuse a grow step.
_NEAR_CEILING = 0.85
# Smallest object worth listing as releasable.
_RELEASE_MIN = 256 * 1024**2
# Row chunks read to estimate the non-zero count of a lazy sparse array.
_SAMPLE_CHUNKS = 3
_VIRTIO_MEM = "/sys/bus/virtio/drivers/virtio_mem"


class LoadCheckFailed(RuntimeError):
    pass


def _meminfo():
    try:
        with open("/proc/meminfo") as f:
            kb = {k: int(v.split()[0]) * 1024 for k, v in (line.split(":", 1) for line in f)}
        return {"total": kb["MemTotal"], "available": kb["MemAvailable"], "swap_free": kb.get("SwapFree", 0)}
    except (OSError, KeyError, ValueError):
        pass
    try:
        import psutil

        vm, sw = psutil.virtual_memory(), psutil.swap_memory()
        return {"total": vm.total, "available": vm.available, "swap_free": sw.free}
    except Exception:
        return None


def _virtio_mem_room():
    """Return (bound virtio-mem device names, region bytes not yet plugged or None)."""
    try:
        devices = {d for d in os.listdir(_VIRTIO_MEM) if d.startswith("virtio")}
    except OSError:
        return None, None
    if not devices:
        return set(), None
    region = plugged = 0
    try:
        with open("/proc/iomem") as f:
            for line in f:
                m = re.match(r"\s*([0-9a-f]+)-([0-9a-f]+) : (.+)$", line)
                if not m:
                    continue
                size = int(m.group(2), 16) - int(m.group(1), 16) + 1
                if m.group(3) in devices:
                    region += size
                elif m.group(3) == "System RAM (virtio_mem)":
                    plugged += size
    except OSError:
        return devices, None
    # Without root, /proc/iomem shows every range as 0-0.
    return devices, (region - plugged if region > _GIB else None)


def _ceiling(mem, ceiling_bytes):
    if ceiling_bytes:
        return ceiling_bytes, "given"
    env = os.environ.get(_CEILING_ENV)
    if env:
        return float(env) * _GIB, _CEILING_ENV
    devices, unplugged = _virtio_mem_room()
    if devices is None:
        return mem["total"], "no elastic memory"
    if not devices:
        return mem["total"], "no elastic pool"
    if unplugged is not None:
        return mem["total"] + unplugged, "read from virtio-mem"
    return max(_DEFAULT_CEILING, mem["total"]), "default elastic ceiling"


def _is_lazy(x):
    # The Loader silently falls back to an eager read when lazy loading fails,
    # so test the array itself, never the lazy= argument.
    return type(x).__module__.startswith("dask")


def _is_sparse(x):
    probe = getattr(x, "_meta", x)  # a dask array's chunk type
    return type(probe).__module__.startswith("scipy.sparse")


def _gb(n):
    return f"{n / _GIB:.1f} GiB"


def _as_dict(tables):
    return tables if isinstance(tables, dict) else {"table": tables}


def _finish(rows):
    for r in rows:
        print(f"[{r['verdict']}] {r['table']}: {r['message']}")
    failed = [f"{r['table']}: {r['message']}" for r in rows if r["verdict"] == "FAIL"]
    if failed:
        raise LoadCheckFailed("; ".join(failed))
    return rows


def _sparse_nbytes(x):
    return sum(getattr(x, a).nbytes for a in ("data", "indices", "indptr", "row", "col", "offsets") if hasattr(x, a))


def _sampled_nnz(x):
    """Non-zero count of a lazy sparse array: exact if few chunks, else from sampled row chunks."""
    bounds = [0]
    for c in x.chunks[0]:
        bounds.append(bounds[-1] + c)
    picks = range(len(bounds) - 1)
    if len(picks) > _SAMPLE_CHUNKS:
        step = (len(picks) - 1) / (_SAMPLE_CHUNKS - 1)
        picks = sorted({round(i * step) for i in range(_SAMPLE_CHUNKS)})
    rows = nnz = 0
    for i in picks:
        nnz += x[bounds[i] : bounds[i + 1]].compute().nnz
        rows += bounds[i + 1] - bounds[i]
    exact = len(picks) == len(bounds) - 1
    return (nnz if exact else math.ceil(nnz / max(rows, 1) * x.shape[0])), exact


def _array_bytes(x):
    """Return (bytes once in memory, how it was estimated) for one array."""
    n_rows = x.shape[0]
    if not _is_lazy(x):
        if _is_sparse(x):
            return _sparse_nbytes(x), "sparse, stored values"
        if hasattr(x, "nbytes"):
            return x.nbytes, "dense"
    dense = math.prod(x.shape) * x.dtype.itemsize
    if not _is_sparse(x):
        return dense, "dense"
    try:
        nnz, exact = _sampled_nnz(x)
    except Exception:
        return dense, "sparse, dense bound (non-zeros unknown)"
    idx = 4 if nnz < 2**31 else 8
    return nnz * (x.dtype.itemsize + idx) + (n_rows + 1) * idx, (
        "sparse, non-zeros counted" if exact else f"sparse, non-zeros sampled from {_SAMPLE_CHUNKS} chunks"
    )


def _frame_bytes(df):
    try:
        if hasattr(df, "memory_usage"):
            return int(df.memory_usage(index=True, deep=False).sum())
        return int(df.ds.nbytes)  # a lazy anndata Dataset2D, sized without reading it
    except Exception:
        return 0


def _footprint(ad):
    """Estimate one table once materialized. Also return its dense feature-aligned arrays."""
    total, var_dense, x_how = 0, [], "no X"
    # Newer anndata also lists X in layers, under the key None.
    for key, x in [("X", ad.X), *((k, v) for k, v in ad.layers.items() if k is not None)]:
        if x is None:
            continue
        n, how = _array_bytes(x)
        total += n
        if key == "X":
            x_how = how
        if how == "dense":
            var_dense.append((n, x.dtype))
    for m in (ad.obsm, ad.obsp):
        for x in m.values():
            total += _array_bytes(x)[0] if hasattr(x, "dtype") and hasattr(x, "shape") else _frame_bytes(x)
    total += _frame_bytes(ad.obs) + _frame_bytes(ad.var)
    return {"bytes": total, "var_dense": var_dense, "how": x_how}


def _joined_bytes(fp, n_vars, n_union):
    """A run's footprint inside an outer join on features."""
    extra = 0
    if n_union > n_vars:
        # Sparse arrays pad with implicit zeros. Dense arrays pad with NaN, which
        # turns integer and bool arrays into float64.
        for n, dtype in fp["var_dense"]:
            promote = 8 / dtype.itemsize if dtype.kind in "iub" else 1
            extra += n * (n_union / n_vars * promote - 1)
    return fp["bytes"] + extra


def _object_bytes(x, seen):
    if id(x) in seen:
        return 0
    seen.add(id(x))
    mod = type(x).__module__
    if mod.startswith("numpy") and hasattr(x, "nbytes"):
        # A view frees nothing on its own; count the buffer it shares once.
        base = x
        while getattr(base, "base", None) is not None:
            base = base.base
        if base is not x:
            if id(base) in seen:
                return 0
            seen.add(id(base))
        return getattr(base, "nbytes", x.nbytes)
    if mod.startswith("scipy.sparse"):
        return _sparse_nbytes(x)
    if mod.startswith("pandas"):
        return _frame_bytes(x)
    if mod.startswith("dask") and hasattr(x, "__dask_graph__"):
        # Only persisted chunks hold memory; a plain dask graph holds tasks.
        return sum(_object_bytes(getattr(v, "value", v), seen) for v in dict(x.__dask_graph__()).values())
    if mod.startswith("anndata") and hasattr(x, "obs"):
        parts = [x.X, *x.layers.values(), *x.obsm.values(), *x.obsp.values(), x.obs, x.var]
        return sum(_object_bytes(p, seen) for p in parts if p is not None)
    return 0


def find_releasable(namespace=None, min_bytes=_RELEASE_MIN, skip=()):
    """List (name, bytes) of large objects in the session. Never deletes anything."""
    if namespace is None:
        namespace = vars(sys.modules["__main__"])
    seen = {id(s) for s in skip}
    found = []
    for name, obj in list(namespace.items()):
        if name.startswith("__"):
            continue
        items = [(name, obj)]
        if isinstance(obj, dict):
            items += [(f"{name}[{k!r}]", v) for k, v in obj.items()]
        elif isinstance(obj, (list, tuple)):
            items += [(f"{name}[{i}]", v) for i, v in enumerate(obj)]
        for label, value in items:
            try:
                n = _object_bytes(value, seen)
            except Exception:
                continue
            if n >= min_bytes:
                found.append((label, n))
    return sorted(found, key=lambda t: -t[1])


def _memory_verdict(name, n_obs, n_vars, need, how, mem, ceiling, source, release):
    size = f"{n_obs:,} cells x {n_vars:,} features, needs ~{_gb(need)} ({how})"
    if n_obs == 0:
        return {"table": name, "verdict": "FAIL", "message": f"{size}; no cells, check the store and any filters"}
    if mem is None:
        return {"table": name, "verdict": "WARN", "message": f"{size}; memory unknown, materialize one run at a time"}
    room = ceiling - (mem["total"] - mem["available"])
    freeable = sum(n for _, n in release)
    names = ", ".join(n for n, _ in release)
    where = f"{_gb(mem['available'])} free now, {_gb(room)} at the {_gb(ceiling)} ceiling ({source})"
    if need <= mem["available"]:
        verdict, msg = "OK", f"{size}; fits ({where})"
    elif need <= room:
        grow = f"memory must grow by ~{_gb(need - mem['available'])}; the VM adds it on demand"
        if need > _NEAR_CEILING * room:
            verdict, msg = "WARN", f"{size}; {where}; {grow}, but this is near the ceiling and a grow can be refused"
        elif need - mem["available"] > mem["swap_free"]:
            verdict, msg = "WARN", (
                f"{size}; {where}; {grow}, but growth can lag one fast read and swap cannot cover the gap; "
                "load one run at a time"
            )
        else:
            verdict, msg = "OK", f"{size}; {where}; {grow}"
        if verdict == "WARN" and release:
            msg += f"; to make room, del {names} ({_gb(freeable)})"
    elif need <= room + freeable:
        verdict, msg = "WARN", f"{size}; {where}; fits only after you del {names} ({_gb(freeable)}) and any other references"
    else:
        verdict, msg = "FAIL", f"{size}; {where}; does not fit even after releasing objects; propose a subset (runs, wells or features) to the user, then .to_memory() only that slice"
    return {"table": name, "verdict": verdict, "message": msg}


def check_before(tables, ceiling_bytes=None, namespace=None, memory=None):
    """Estimate what materializing each (lazily loaded) table, and all of them together, would cost.

    namespace: where to look for releasable objects (default ``__main__``).
    memory: a dict with total / available / swap_free bytes, for tests.
    """
    tables = _as_dict(tables)
    gc.collect()
    mem = memory or _meminfo()
    ceiling, source = (ceiling_bytes or 0, "given") if mem is None else _ceiling(mem, ceiling_bytes)
    release = find_releasable(namespace, skip=[a for a in tables.values() if _is_lazy(a.X)])
    if release:
        listed = ", ".join(f"{n} ({_gb(b)})" for n, b in release)
        print(f"[INFO] releasable objects in the session (not deleted): {listed}")
    rows, fps = [], {}
    for name, ad in tables.items():
        fps[name] = _footprint(ad)
        need = fps[name]["bytes"] * _MARGIN
        rows.append(_memory_verdict(name, ad.n_obs, ad.n_vars, need, fps[name]["how"], mem, ceiling, source, release))
    if len(tables) > 1:
        union = set()
        for ad in tables.values():
            union.update(ad.var_names)
        joined = sum(_joined_bytes(fps[k], ad.n_vars, len(union)) for k, ad in tables.items())
        runs = sum(fp["bytes"] for fp in fps.values())
        hows = sorted({fp["how"] for fp in fps.values()})
        total = _memory_verdict(
            f"all {len(tables)} combined (runs + outer-joined copy)",
            sum(ad.n_obs for ad in tables.values()),
            len(union),
            (runs + joined) * _MARGIN,
            "; ".join(hows),
            mem,
            ceiling,
            source,
            release,
        )
        if total["verdict"] == "FAIL":
            total["message"] += "; or keep concatenated=False and work one run at a time"
        rows.append(total)
    return _finish(rows)


def check_after(tables, expect="materialized", ceiling_bytes=None, memory=None):
    """Check loaded tables. expect is "materialized" or "lazy"."""
    mem = memory or _meminfo()
    room = None
    if mem is not None:
        ceiling, _ = _ceiling(mem, ceiling_bytes)
        room = ceiling - (mem["total"] - mem["available"])
    rows = []
    for name, ad in _as_dict(tables).items():
        lazy = _is_lazy(ad.X)
        if ad.n_obs == 0:
            verdict, msg = "FAIL", "no cells; stop and report, the store or filters are wrong"
        elif expect == "materialized" and lazy:
            verdict, msg = "FAIL", ".X is still lazy (dask); reload with lazy=False before pp / Scanpy"
        elif expect == "lazy" and not lazy:
            verdict, msg = "WARN", ".X is in memory although lazy was requested; the library fell back to an eager read"
        elif not lazy and not _is_sparse(ad.X) and room is not None and ad.X.nbytes > room:
            verdict, msg = "WARN", (
                f".X is dense and {_gb(ad.X.nbytes)}; one more copy does not fit in the {_gb(room)} left "
                "at the ceiling; subset before more work"
            )
        else:
            kind = "lazy" if lazy else ("sparse" if _is_sparse(ad.X) else "dense")
            verdict, msg = "OK", f"{ad.n_obs:,} cells x {ad.n_vars:,} features, {kind}"
        rows.append({"table": name, "verdict": verdict, "message": msg})
    return _finish(rows)
