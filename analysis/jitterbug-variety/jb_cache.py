"""Content-addressed memoisation + speculative parallel prefetch for the
expensive pure solves in this directory's research scripts.

WHY THIS EXISTS. `jb_z_quasistatic_array.py` spends 99.8% of a 6m14s run
inside ~30 `crank_run` continuation solves (measured 2026-08-22: z17 291.6s,
z15 49.0s, z13 22.0s, z14 10.2s; every geometry gate together 0.58s). Those
solves do not change when a print statement, a gate row, or a docstring
changes -- but they were recomputed anyway, on one core of sixteen, on every
single edit. That is an edit-loop defect, not a physics defect, and it is
fixed here without touching a single number.

TWO INDEPENDENT MECHANISMS, deliberately kept separable:

  1. MEMOISATION (`memoize`). Keyed on the SHA-256 of the transitive source
     closure of the wrapped function -- every function and every module
     constant it reads, followed recursively -- together with its bound
     arguments. Editing the stepper, the contact kernel, or any constant
     either reads invalidates every affected entry automatically. Editing a
     print statement invalidates nothing. The closure is computed from
     bytecode (`co_names`), NOT from a hand-maintained dependency list, so a
     newly added helper cannot be silently forgotten -- the failure mode a
     declared list would have.

  2. PREFETCH (`prefetch`). The solves are pure and mutually independent, so
     they can be computed in any order, in parallel, ahead of the serial pass
     that prints them. Which argument tuples a run will need is not knowable
     before the run, so it is LEARNED: every memoised call appends its bound
     arguments to a trace, and the next run replays that trace through a
     process pool before `main()` starts. The trace is keyed on ARGUMENTS
     ONLY, never on source, which is what makes it survive exactly the edit
     that needs it most: change the kernel and the cache is correctly cold,
     but the trace still knows all 30 argument tuples, so the cold run is
     parallel rather than serial.

     A trace miss is not an error. An argument tuple the trace has not seen
     is simply computed in the serial pass and recorded for next time, so
     the mechanism is self-healing and needs no maintenance when a grid is
     re-priced.

CORRECTNESS. Neither mechanism can change a computed value: a cache hit
returns what the same source computed from the same arguments, and a
prefetched value is computed by the same function in a fresh interpreter.
The guard against that claim being merely asserted is `--no-cache`, which
bypasses both and must produce byte-identical output; that equivalence is
what the caller's own gate should check when either mechanism is changed.

NOT A BUILD SYSTEM. This is exploration tooling. There is no dependency
graph to declare, no rule file, and no daemon; the cache is a directory of
pickles that can be deleted at any time with no consequence beyond one slow
run.
"""

from __future__ import annotations

import hashlib
import inspect
import os
import pickle
import sys
import tempfile
import types
from concurrent.futures import ProcessPoolExecutor

import numpy as np

#: Cache root. Sits beside the scripts rather than in a temp dir so it
#: survives reboots (a cold cache costs minutes, so losing it matters) and so
#: `rm -rf` on it is an obvious, local gesture.
CACHE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".jbcache")

#: Bumped only when the KEY DERIVATION below changes shape. Entries written by
#: an older derivation are then unreadable rather than mis-matched.
_KEY_VERSION = "v1"

#: Set by `--no-cache` / JB_NO_CACHE=1. When true, `memoize` degrades to a
#: transparent pass-through and `prefetch` is a no-op -- the reference path
#: any equivalence check runs against.
_DISABLED = bool(os.environ.get("JB_NO_CACHE"))

#: Worker count for `prefetch`. One per core less two, so a prefetch does not
#: make the machine unusable; capped because the solves are single-threaded
#: numpy and oversubscription only adds scheduler noise.
DEFAULT_WORKERS = max(1, min(16, (os.cpu_count() or 4) - 2))


def _resolve_module(module_name):
    """The live module object for `module_name`.

    A script run as `python jb_z_quasistatic_array.py` is `__main__`, not
    `jb_z_quasistatic_array`, so the importable name it declares to `memoize`
    is not in `sys.modules` under that name. Falling back to `__main__` when
    its file basename matches keeps the SAME name usable from both the script
    path and a prefetch worker's fresh import -- which is what lets a worker
    re-enter the function at all.
    """
    mod = sys.modules.get(module_name)
    if mod is not None:
        return mod
    main = sys.modules.get("__main__")
    main_file = getattr(main, "__file__", None)
    if main_file and os.path.splitext(os.path.basename(main_file))[0] == module_name:
        return main
    return __import__(module_name)


def disable():
    """Turn both mechanisms off for the rest of the process."""
    global _DISABLED
    _DISABLED = True


def enabled():
    return not _DISABLED


# ==========================================================================
# SOURCE FINGERPRINT: the transitive closure of what a function reads.
# ==========================================================================

def _const_digest(value):
    """A stable byte encoding for a module-level CONSTANT.

    `repr` is not enough for arrays (numpy elides large ones with `...`, which
    would make two different arrays fingerprint identically -- a silent stale
    hit), so arrays are hashed from their raw buffer plus dtype and shape.
    Anything without a stable repr is refused loudly rather than hashed as its
    memory address, which would make the cache miss on every run and look like
    the cache was broken.
    """
    if isinstance(value, np.ndarray):
        return b"nd|" + str(value.dtype).encode() + b"|" \
            + str(value.shape).encode() + b"|" + value.tobytes()
    if isinstance(value, (int, float, complex, bool, str, bytes, type(None))):
        return b"sc|" + repr(value).encode()
    if isinstance(value, (tuple, list, frozenset, set)):
        return b"sq|" + b",".join(_const_digest(v) for v in value)
    if isinstance(value, dict):
        return b"dc|" + b",".join(
            _const_digest(k) + b"=" + _const_digest(v)
            for k, v in sorted(value.items(), key=lambda kv: repr(kv[0])))
    r = repr(value)
    if "0x" in r and " object at " in r:
        raise TypeError(
            f"jb_cache: module global of type {type(value).__name__!r} has no "
            "stable repr, so it cannot be fingerprinted. Give it one, or keep "
            "it out of the memoised call graph.")
    return b"rp|" + r.encode()


def _walk_names(code, out):
    """Every global name a code object reads, including inside nested code
    objects (comprehensions, closures, lambdas) -- which is where a dependency
    would otherwise hide."""
    out.update(code.co_names)
    for c in code.co_consts:
        if isinstance(c, types.CodeType):
            _walk_names(c, out)


def _source_names(src):
    """Identifiers appearing in `src` as CODE.

    `co_names` alone is not sufficient: a constant used as a DEFAULT ARGUMENT
    (`def f(a, rtol=QP_NULL_RTOL)`) is evaluated once at definition time, so
    the name never appears in any caller's bytecode -- changing the constant's
    value would leave the source text identical and the cache stale. Found by
    asserting invalidation constant-by-constant; QP_NULL_RTOL was the one that
    read "no effect" when it should not have. Decorator expressions and
    annotations have the same shape.

    Tokenising rather than regexing the text is what keeps this an
    over-approximation only in the harmless direction: names inside comments
    and string literals are NOT picked up, so a docstring that merely mentions
    a constant does not bind this function's cache to that constant's value.
    """
    import io as _io
    import tokenize as _tok
    names = set()
    try:
        for t in _tok.generate_tokens(_io.StringIO(src).readline):
            if t.type == _tok.NAME:
                names.add(t.string)
    except (_tok.TokenError, IndentationError, SyntaxError):
        # A decorated or indented fragment can fail to tokenise standalone.
        # co_names still covers the common case, so degrade rather than fail.
        pass
    return names


def source_fingerprint(fn, module):
    """SHA-256 over the transitive source closure of `fn` within `module`.

    Walks global names read by `fn`, then by everything it reaches, to a fixed
    point. Functions contribute their source text; constants contribute
    `_const_digest`. Names that resolve to modules, classes from elsewhere, or
    nothing at all contribute nothing -- they are not part of this file's own
    editable surface.

    The returned digest changes if and only if something the function actually
    reads changed. A print statement elsewhere in the file does not move it.
    """
    seen = set()
    parts = []
    frontier = [fn.__name__]
    mod_globals = vars(module)
    while frontier:
        name = frontier.pop()
        if name in seen:
            continue
        seen.add(name)
        obj = mod_globals.get(name, None)
        if obj is None:
            continue
        if isinstance(obj, types.ModuleType):
            continue
        if isinstance(obj, types.FunctionType):
            # Unwrap BEFORE the ownership test: a memoised function's wrapper
            # is defined in THIS module, so testing the wrapper's `__module__`
            # rejects it and the closure silently comes back empty -- which
            # looks like a working cache that never invalidates. Found live,
            # by asserting that a constant the stepper reads moves the hash.
            target = getattr(obj, "__wrapped__", obj)
            if getattr(target, "__module__", None) != module.__name__:
                continue
            try:
                src = inspect.getsource(target)
            except (OSError, TypeError):
                src = repr(target)
            parts.append(b"fn|" + name.encode() + b"|" + src.encode())
            nxt = set()
            _walk_names(target.__code__, nxt)
            nxt |= _source_names(src)
            frontier.extend(nxt)
            continue
        if callable(obj) and not isinstance(obj, type):
            continue
        if isinstance(obj, type):
            continue
        try:
            parts.append(b"cn|" + name.encode() + b"|" + _const_digest(obj))
        except TypeError:
            raise
    h = hashlib.sha256()
    h.update(_KEY_VERSION.encode())
    for p in sorted(parts):
        h.update(p)
        h.update(b"\x00")
    return h.hexdigest()


# ==========================================================================
# KEY DERIVATION AND STORE
# ==========================================================================

def _bound_args(fn, args, kwargs):
    """Arguments normalised through the signature with defaults applied, so
    `f(x)` and `f(x, h0=H_STEP)` are the SAME cache entry rather than two."""
    sig = inspect.signature(getattr(fn, "__wrapped__", fn))
    b = sig.bind(*args, **kwargs)
    b.apply_defaults()
    return b.arguments


def _arg_digest(arguments):
    return hashlib.sha256(
        pickle.dumps(arguments, protocol=pickle.HIGHEST_PROTOCOL)).hexdigest()


def _entry_path(fn_key, src_hash, arg_hash):
    d = os.path.join(CACHE_DIR, fn_key, src_hash[:16])
    return d, os.path.join(d, arg_hash + ".pkl")


def _write_atomic(path, payload):
    d = os.path.dirname(path)
    os.makedirs(d, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=d, suffix=".tmp")
    try:
        with os.fdopen(fd, "wb") as fh:
            pickle.dump(payload, fh, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(tmp, path)          # atomic: a torn write is never read
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


# ==========================================================================
# TRACE: which argument tuples this script actually needs.
# ==========================================================================

def _trace_path(fn_key):
    return os.path.join(CACHE_DIR, fn_key, "trace.pkl")


def _trace_load(fn_key):
    try:
        with open(_trace_path(fn_key), "rb") as fh:
            return pickle.load(fh)
    except (OSError, EOFError, pickle.UnpicklingError):
        return {}


def _trace_record(fn_key, arg_hash, arguments):
    """Append one argument tuple to the trace. Source-INDEPENDENT by
    construction: the trace answers "what will be needed", never "what is
    still valid", which is why a kernel edit leaves it usable."""
    tr = _trace_load(fn_key)
    if arg_hash in tr:
        return
    tr[arg_hash] = arguments
    _write_atomic(_trace_path(fn_key), tr)


# ==========================================================================
# THE DECORATOR
# ==========================================================================

def memoize(module_name):
    """Memoise a pure function on (transitive source, bound arguments).

    `module_name` is the containing module's `__name__`, used both to resolve
    the source closure and to let a prefetch worker re-import and call the
    same function in a fresh interpreter.

    The wrapped function MUST be pure and its arguments and return value
    picklable. Nothing here checks that; it is the caller's claim, and for
    `crank_run` it was verified (arguments are a picklable `Topology` plus
    scalars; the return is a dict of picklables; the body reads only module
    constants).
    """
    def deco(fn):
        fn_key = f"{module_name}.{fn.__name__}"
        cached_src = []

        def _src_hash():
            if not cached_src:
                cached_src.append(
                    source_fingerprint(fn, _resolve_module(module_name)))
            return cached_src[0]

        def wrapper(*args, **kwargs):
            if _DISABLED:
                return fn(*args, **kwargs)
            arguments = _bound_args(fn, args, kwargs)
            ah = _arg_digest(arguments)
            _, path = _entry_path(fn_key, _src_hash(), ah)
            try:
                with open(path, "rb") as fh:
                    return pickle.load(fh)
            except (OSError, EOFError, pickle.UnpicklingError):
                pass
            result = fn(*args, **kwargs)
            _write_atomic(path, result)
            _trace_record(fn_key, ah, arguments)
            return result

        wrapper.__wrapped__ = fn
        wrapper.__name__ = fn.__name__
        wrapper.__doc__ = fn.__doc__
        wrapper._jb_fn_key = fn_key
        wrapper._jb_src_hash = _src_hash
        wrapper._jb_module = module_name
        return wrapper
    return deco


# ==========================================================================
# PREFETCH
# ==========================================================================

def _worker(module_name, fn_name, arguments):
    """Runs in a pool process: import the module, call the memoised function,
    let it populate the cache. Returns only a status string -- the result
    travels through the cache file, not back over the pipe, so a large return
    value is not paid for twice."""
    mod = _resolve_module(module_name)
    fn = getattr(mod, fn_name)
    with np.errstate(all="ignore"):
        fn(**arguments)
    return fn_name


def prefetch(fn, workers=None, quiet=False):
    """Compute every traced-but-uncached argument tuple for `fn` in parallel.

    Returns (n_hit, n_computed). A no-op when caching is disabled, when no
    trace exists yet (the very first run), or when nothing is missing.

    Failures inside a worker are REPORTED AND SWALLOWED on purpose: a prefetch
    is an optimisation, and the serial pass that follows will recompute the
    entry and surface any genuine error at its real call site, with its real
    traceback, rather than here in a pool process.
    """
    if _DISABLED:
        return (0, 0)
    fn_key = getattr(fn, "_jb_fn_key", None)
    if fn_key is None:
        raise TypeError("prefetch() needs a @memoize-wrapped function")
    trace = _trace_load(fn_key)
    if not trace:
        return (0, 0)
    src = fn._jb_src_hash()
    missing = []
    hits = 0
    for ah, arguments in trace.items():
        _, path = _entry_path(fn_key, src, ah)
        if os.path.exists(path):
            hits += 1
        else:
            missing.append(arguments)
    if not missing:
        return (hits, 0)
    n = workers or DEFAULT_WORKERS
    n = max(1, min(n, len(missing)))
    if not quiet:
        print(f"  [prefetch] {fn.__name__}: {hits} cached, "
              f"{len(missing)} to compute on {n} workers ...", flush=True)
    done = 0
    with ProcessPoolExecutor(max_workers=n) as ex:
        futs = [ex.submit(_worker, fn._jb_module, fn.__name__, a)
                for a in missing]
        for f in futs:
            try:
                f.result()
                done += 1
            except Exception as exc:                       # noqa: BLE001
                print(f"  [prefetch] one task failed ({type(exc).__name__}: "
                      f"{exc}); it will be recomputed serially.", flush=True)
    if not quiet:
        print(f"  [prefetch] {done}/{len(missing)} computed.", flush=True)
    return (hits, done)


def parse_argv(argv):
    """Consume this module's own flags from `argv`, returning the remainder.

    `--no-cache` disables both mechanisms (the reference path).
    `--clear-cache` deletes the store, trace included, and continues.
    """
    rest = []
    for a in argv:
        if a == "--no-cache":
            disable()
        elif a == "--clear-cache":
            import shutil
            shutil.rmtree(CACHE_DIR, ignore_errors=True)
        else:
            rest.append(a)
    return rest
