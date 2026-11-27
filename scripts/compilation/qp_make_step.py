#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
qp_make_step.py: the step counter of the QP make build.

Prints ninja-style progress for the make build system: one
`[k/N] <description>` line per build step, and a final summary
(the `done` command). State lives in one JSON file per scope --
config/.qp_make_steps.json for the root graph, <module>/
.qp_make_steps.json for each module -- holding {run_id, n, k, t0}.

Commands:

  next <state-file> <n> <description>
      Increment the counter of <state-file> and print
      `[k/N] <description>` (or `[k] <description>` when k > n:
      more steps ran than the static count predicted, which is
      normal for the root scope -- the ezfio sub-build compiles
      against it -- and for the module scope -- the pattern rules
      of irpf90.make fire once per object). A state file that is
      absent, corrupt, or stamped with a different QP_MAKE_RUN_ID
      is (re)initialized under the lock, so a fresh make run
      always starts at [1/N].

  done <root-state> [n-modules [rust-edge] [module-state ...]]
      Print the summary: `make: no work to do.` when this is a
      no-op run -- the root counter never advanced past its no-op
      baseline (n-modules + rust-edge: the .PHONY module edges
      always fire, even when nothing rebuilds, and so does the rust
      edge when the default goal is built) and no module counter
      advanced in this run -- else
      `make: <k> steps, <t>s elapsed`.

All read-modify-write cycles run under an exclusive fcntl lock on
the state file itself (macOS has no flock(1); python3 is a hard QP
dependency), and the single output line is printed -- and flushed --
while the lock is held, so concurrent steps never interleave.
"""

import fcntl
import json
import os
import sys
import time

STATE_KEYS = ("run_id", "n", "k", "t0")


def _run_id():
    return os.environ.get("QP_MAKE_RUN_ID", "")


def _load(lock_fd):
    """Read the state under the lock: None when absent or corrupt."""
    lock_fd.seek(0)
    content = lock_fd.read()
    if not content.strip():
        return None
    try:
        state = json.loads(content)
    except ValueError:
        return None
    if not isinstance(state, dict):
        return None
    return state


def _fresh():
    return {"run_id": _run_id(), "n": 0, "k": 0, "t0": time.time()}


def _lock(state_file):
    # "a+" creates the file when absent (the state file is also the
    # lock: flock it before reading or writing it).
    return open(state_file, "a+")


def cmd_next(args):
    if len(args) != 3:
        sys.stderr.write(
            "qp_make_step.py: usage: next <state-file> <n> <description>\n")
        return 1
    state_file, n_text, desc = args
    try:
        n = int(n_text)
    except ValueError:
        n = 0  # a failed count must never break the build
    with _lock(state_file) as lock_fd:
        fcntl.flock(lock_fd, fcntl.LOCK_EX)
        state = _load(lock_fd)
        if state is None or state.get("run_id") != _run_id():
            # First call in this scope, or a new make run: re-init.
            state = _fresh()
        state["n"] = n
        state["k"] = int(state.get("k", 0)) + 1
        k, total = state["k"], state["n"]
        lock_fd.seek(0)
        lock_fd.truncate()
        lock_fd.write(json.dumps(state))
        lock_fd.flush()
        if k <= total:
            sys.stdout.write("[{0}/{1}] {2}\n".format(k, total, desc))
        else:
            sys.stdout.write("[{0}] {1}\n".format(k, desc))
        sys.stdout.flush()
    return 0


def _active_k(state_file):
    """
    The k of `state_file` for the current run: 0 when the file is
    absent, corrupt, or stamped with another QP_MAKE_RUN_ID. Probed
    under the lock (a separate open, so it serializes with the
    writer of that scope); never creates the file.
    """
    try:
        with open(state_file, "r+") as lock_fd:
            fcntl.flock(lock_fd, fcntl.LOCK_EX)
            state = _load(lock_fd)
    except OSError:
        return 0
    if state is None or state.get("run_id") != _run_id():
        return 0
    try:
        return int(state.get("k", 0))
    except (TypeError, ValueError):
        return 0


def cmd_done(args):
    if len(args) < 2:
        sys.stderr.write("qp_make_step.py: usage: done <root-state> "
                         "[n-modules [rust-edge] [module-state ...]]\n")
        return 1
    state_file = args[0]
    n_modules = int(args[1]) if len(args) > 1 and args[1].isdigit() else 0
    rust_edge = int(args[2]) if len(args) > 2 and args[2].isdigit() else 0
    module_states = args[3:]
    with _lock(state_file) as lock_fd:
        fcntl.flock(lock_fd, fcntl.LOCK_EX)
        state = _load(lock_fd)
        k = (int(state.get("k", 0))
             if state is not None and state.get("run_id") == _run_id()
             else 0)
        # A no-op run: the root counter stopped at its no-op baseline
        # and no module did any work (see the module scopes above).
        baseline = n_modules + rust_edge
        no_work = (k <= baseline
                   and all(_active_k(ms) == 0 for ms in module_states))
        if no_work or k == 0:
            sys.stdout.write("make: no work to do.\n")
        else:
            elapsed = time.time() - float(state.get("t0", time.time()))
            sys.stdout.write("make: {0} steps, {1:.1f}s elapsed\n".format(
                k, elapsed))
        sys.stdout.flush()
    return 0


def main():
    if len(sys.argv) < 2:
        sys.stderr.write(__doc__)
        return 1
    if sys.argv[1] == "next":
        return cmd_next(sys.argv[2:])
    if sys.argv[1] == "done":
        return cmd_done(sys.argv[2:])
    sys.stderr.write(
        "qp_make_step.py: unknown command {0!r} (next|done)\n".format(
            sys.argv[1]))
    return 1


if __name__ == "__main__":
    sys.exit(main())
