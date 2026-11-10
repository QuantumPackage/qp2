#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Golden-corpus generator, replayer and differ for the QP2 Rust port.

What this tool does
-------------------
QP2's EZFIO directories and the output of the ``qp_*`` binaries are the
contract that the Rust port must reproduce.  This tool pins that contract:

``generate``
    Runs a *given* tool set once and writes the golden corpus:
    ``fixtures/`` (six EZFIO directories), ``transcripts/`` (the
    ``qp_edit -c`` transcript of each fixture, the files it created with
    their decompressed content, and the pre-fix coordinate-quantisation
    probe) and ``md5/basis-md5.jsonl`` (the ``ao_md5`` /
    ``ao_extra_md5`` of 27 bases).  It also writes ``manifest.json``,
    which is what the other subcommands replay.  The *normative* part of
    the check contract, ``contracts/``, is hand-written and is never
    touched by ``generate``.

``run``
    Replays every scenario of the manifest with a given tool set and
    writes a *snapshot*: stdout, stderr, exit status and the resulting
    EZFIO directory (arrays decompressed to text) for each
    ``(fixture, scenario)`` pair.

``diff``
    Compares two snapshots and exits 1 on any difference.

``verify``
    Replays the requested scenarios (``--only create,check`` by default)
    and gates them against the corpus.  ``create`` must reproduce
    ``fixtures/`` byte for byte, and ``check`` must reproduce its
    transcript and obey ``contracts/check-mutated.json``
    (``--contract post-fix``, the Rust target, or ``--contract pre-fix``,
    what the reference OCaml build does).  The contract's
    ``post-fix-changed`` list is a **one-way allowance** for the w12
    one-time ``nucl_coord`` normalisation of c2h2/b2: a toolset that
    applies the normalisation must produce exactly the stored expected
    payload in both scenarios (a difference outside the list, or a
    difference of the wrong bytes, is a failure); a toolset that does
    not apply it (the reference OCaml build, whose multiply pipeline is
    the corpus's own fixed point, as run by CI) leaves the file
    untouched and passes.  The negative side -- that the port does
    apply the normalisation -- is enforced by the port's own test suite
    (``qp_edit`` integration test
    ``c2_no_pre_existing_file_is_rewritten``, which also asserts the
    change is real), not by this tool.

Why the corpus is normalised
----------------------------
``ezfio/creation`` (timestamp), ``ezfio/user`` (login) and
``ezfio/library`` (absolute path of the EZFIO source tree) are written by
every tool run and carry no contract.  They are removed from the corpus
and from snapshots.  ``.version`` is *kept*: it is the EZFIO format
version and a Rust writer must still produce it.

gzip containers embed the compression time, so ``.gz`` files are
recompressed with ``mtime=0`` (level 9) before being written.  The
decompressed payload -- which is what the tools actually read and what
D6 of the plan makes the oracle -- is preserved byte for byte.  Snapshots
store the decompressed payload as plain text, so ``diff -ru`` on two
snapshot trees is also a valid comparison.

Machine-dependent lines of the tool output (dates, wall times, Git
stamps, task-server address, resource usage, memory budget) are dropped;
see ``NOISE``.  Every other byte is kept.

Units
-----
Coordinates in EZFIO are bohr; QP2's ``data/basis`` files store nothing
dimensional.  The script never converts units itself: it copies values
as text.

Cost
----
``generate`` runs ``qp_create_ezfio`` six times (plus one ``qp_run scf``
on HCN, a few seconds), ``qp_edit -c`` six times and
``qp_create_ezfio`` 27 times for the md5 corpus: about 20 s and ~30 MB
of scratch space on the reference laptop build.  ``run`` costs the same
minus the md5 part.  All work happens in a temporary directory under
``$TMPDIR`` and is copied into place only on success.

Usage
-----
    python3 golden_harness.py generate --qp-root ~/qp2 --corpus rust/tests/golden
    python3 golden_harness.py run      --qp-root ~/qp2 --corpus rust/tests/golden \
                                       --out $HOME/OpenCode/tmp/golden-a
    python3 golden_harness.py run      --qp-root ~/qp2 --corpus rust/tests/golden \
                                       --out $HOME/OpenCode/tmp/golden-b
    python3 golden_harness.py diff     $HOME/OpenCode/tmp/golden-a $HOME/OpenCode/tmp/golden-b
    python3 golden_harness.py verify   --qp-root ~/qp2 --corpus rust/tests/golden
    python3 golden_harness.py verify   --qp-root ~/qp2 --only check \\
                                       --contract pre-fix

``--qp-root`` may be any build tree: the reference OCaml build (binaries
in ``$QP_ROOT/ocaml``), the Rust build (``$QP_ROOT/rust/target/release``)
or the exported static package (``$QP_ROOT/bin``).  ``--qp-bin`` overrides
the search.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))

#: Files written by every EZFIO tool run that carry machine state, not data.
VOLATILE = ("ezfio/creation", "ezfio/user", "ezfio/library")

#: Directories dropped from snapshots.  ``work/`` is documented scratch
#: space (man qp_run) and may contain uncompressible Fortran files; the
#: committed fixtures keep it because it is small and part of a valid
#: directory.
SNAPSHOT_SKIP_DIRS = ("work",)

#: Output lines that depend on the machine, the clock or the build.
#: Dropped, not replaced: both sides of a diff lose the same lines.
NOISE = (
    r"^Date: ",
    r"^Wall time: ",
    r"^Git (Commit|Date  |SHA1  ): ",
    r"^Task server running : ",
    r"^\.\. >>>>> \[ ",
    r"^\* Total available memory \(GB\)",
    r"^\* Target maximum memory \(GB\)",
    r"^\s*Time to provide ",
)
NOISE_RE = tuple(re.compile(p) for p in NOISE)

#: Names of the two ``qp_edit -c`` contracts kept in
#: ``contracts/check-mutated.json``: the normative post-``32c9899d``
#: target and the pre-fix behaviour of the reference build (evidence).
CHECK_CONTRACTS = ("post-fix", "pre-fix")

#: The six fixtures.  ``create`` is the exact argv of the reference
#: ``qp_create_ezfio``; ``post`` is run in the same directory afterwards.
FIXTURES = (
    {
        "name": "h2",
        "notes": "H2, cc-pvdz, singlet; the smallest possible directory.",
        "inputs": ["h2_1.xyz"],
        "create": ["qp_create_ezfio", "-b", "cc-pvdz",
                   "-o", "ezfio", "h2_1.xyz"],
        "post": [],
    },
    {
        "name": "h3",
        "notes": "H3, cc-pvdz, multiplicity 2 (open shell, 3 electrons).",
        "inputs": ["h3_2.xyz"],
        "create": ["qp_create_ezfio", "-b", "cc-pvdz", "-m", "2",
                   "-o", "ezfio", "h3_2.xyz"],
        "post": [],
    },
    {
        "name": "c2h2",
        "notes": "C2H2, cc-pvdz_ecp_bfd with the bfd pseudopotential.",
        "inputs": ["c2h2.xyz"],
        "create": ["qp_create_ezfio", "-b", "cc-pvdz_ecp_bfd", "-p", "bfd",
                   "-o", "ezfio", "c2h2.xyz"],
        "post": [],
    },
    {
        "name": "b2",
        "notes": "B2 at 3.0 Angstrom, 6-31g, from a Gaussian Z-matrix "
                 "(the Z-matrix input path).",
        "inputs": ["b2_stretched.zmt"],
        "create": ["qp_create_ezfio", "-b", "6-31g",
                   "-o", "ezfio", "b2_stretched.zmt"],
        "post": [],
    },
    {
        "name": "h2g",
        "notes": "H2, cc-pvdz, one ghost atom X at the bond midpoint "
                 "(-g 5.0); exercises zero nuclear charge.",
        "inputs": ["h2_1.xyz"],
        "create": ["qp_create_ezfio", "-b", "cc-pvdz", "-g", "5.0",
                   "-o", "ezfio", "h2_1.xyz"],
        "post": [],
    },
    {
        "name": "hcn",
        "notes": "HCN, aug-cc-pvdz, converged SCF MOs.  Carries the "
                 "ao_md5 / mo_basis/ao_md5 pair of the D4 contract.",
        "inputs": ["hcn.xyz"],
        "create": ["qp_create_ezfio", "-b", "aug-cc-pvdz",
                   "-o", "ezfio", "hcn.xyz"],
        "post": [["qp_run", "scf", "ezfio"]],
    },
)

#: Bases of the md5 corpus.  For each file the *last* element listed in
#: the basis database is used; that is deterministic, exercises the tail
#: of the parser, and needs no per-basis table.  ``def2-svp``/``def2-tzvpp``
#: (RADON) and ``6-31g_star``/``6-31g_star_star`` (KRYPTON) end with
#: byte-identical blocks in the shipped data, so their md5s collide by
#: construction; see ``rust/tests/golden/README.md``.
BASES = (
    "sto-3g", "3-21g", "6-31g", "6-31g_star", "6-31g_star_star",
    "6-31+g_star", "6-311g", "6-311G_star", "6-311++g_3df_3pd", "m6-31g",
    "cc-pvdz", "cc-pvtz", "cc-pvqz", "aug-cc-pvdz", "aug-cc-pvtz",
    "d-aug-cc-pvdz", "cc-pcvdz", "cc-pwcvdz", "ano-rcc", "def2-svp",
    "def2-tzvpp", "lanl08", "lanl08+", "cc-pvdz_ecp_bfd",
    "cc-pvdz_ecp_ccecp", "aug-cc-pvtz_ecp_bfd", "aug-cc-pvtz_ecp_ccecp",
)

ELEMENT_RE = re.compile(r"^[A-Z][A-Z]+$")
PSEUDO_RE = re.compile(r"_ecp_([a-z0-9]+)$")


# --------------------------------------------------------------------------
# Tool discovery and process execution
# --------------------------------------------------------------------------

def find_tools(qp_root: str, qp_bin: str | None) -> str:
    """Return the directory holding the ``qp_*`` executables.

    Search order: ``--qp-bin``, then ``$QP_ROOT/ocaml`` (the reference
    OCaml build), then ``$QP_ROOT/rust/target/release`` (the Rust build),
    then ``$QP_ROOT/bin`` (the exported static package).  A directory
    qualifies when it has ``qp_create_ezfio`` or
    ``qp_create_ezfio.native``.
    """
    candidates = [qp_bin] if qp_bin else [
        os.path.join(qp_root, "ocaml"),
        os.path.join(qp_root, "rust", "target", "release"),
        os.path.join(qp_root, "bin"),
    ]
    for d in candidates:
        if not d:
            continue
        for name in ("qp_create_ezfio", "qp_create_ezfio.native"):
            if os.path.exists(os.path.join(d, name)):
                return os.path.abspath(d)
    raise SystemExit(
        "golden_harness: no qp_create_ezfio under "
        + ", ".join(repr(d) for d in candidates)
        + "; pass --qp-root or --qp-bin"
    )


def tool_env(qp_root: str, tools: str) -> dict:
    """Environment of every tool run.

    ``QP_ROOT`` is what the binaries use to find ``data/basis`` and
    ``data/pseudo``; ``OMP_NUM_THREADS=1`` keeps the SCF fixture
    reproducible (threaded BLAS reductions are not bit-reproducible).
    """
    env = dict(os.environ)
    env["QP_ROOT"] = os.path.abspath(qp_root)
    env["OMP_NUM_THREADS"] = "1"
    env["PATH"] = tools + os.pathsep + env.get("PATH", "/usr/bin:/bin")
    return env


def run_step(argv, cwd: str, env: dict) -> dict:
    """Run one command, capture stdout/stderr/exit, never raise on exit."""
    try:
        p = subprocess.run(argv, cwd=cwd, env=env,
                           capture_output=True, text=False, timeout=1800)
        return {
            "argv": list(argv),
            "stdout": normalise_output(p.stdout),
            "stderr": normalise_output(p.stderr),
            "exit": p.returncode,
        }
    except subprocess.TimeoutExpired:
        return {"argv": list(argv), "stdout": b"", "stderr": b"", "exit": -1}


def normalise_output(data: bytes) -> bytes:
    """Drop machine-dependent lines from a captured stream.

    Text is decoded as UTF-8 with ``surrogateescape`` so that a stray byte
    cannot crash the harness; the bytes are then re-encoded unchanged.
    Lines are split on ``\\n`` only, and a trailing newline is preserved.
    """
    text = data.decode("utf-8", "surrogateescape")
    if not text:
        return b""
    trailing = text.endswith("\n")
    lines = text.split("\n")
    if trailing:
        lines = lines[:-1]
    kept = [l for l in lines if not any(r.search(l) for r in NOISE_RE)]
    out = "\n".join(kept)
    if trailing:
        out += "\n"
    return out.encode("utf-8", "surrogateescape")


def git_info(qp_root: str) -> dict:
    """Commit and dirty paths of the tool tree, for provenance only."""
    def git(*args):
        p = subprocess.run(["git", "-C", qp_root, *args],
                           capture_output=True, text=True)
        return p.stdout if p.returncode == 0 else ""

    commit = git("rev-parse", "HEAD").strip()
    dirty = sorted(l for l in git("status", "--porcelain").splitlines() if l)
    return {"git_commit": commit or None, "git_dirty_paths": dirty}


# --------------------------------------------------------------------------
# EZFIO tree handling
# --------------------------------------------------------------------------

def is_gzip(path: str) -> bool:
    with open(path, "rb") as f:
        return f.read(2) == b"\x1f\x8b"


def recompress_deterministic(path: str) -> None:
    """Rewrite a ``.gz`` file with ``mtime=0`` and level 9.

    Content-preserving: the decompressed payload is written back
    unchanged, so every reader sees the same bytes.
    """
    with open(path, "rb") as f:
        payload = gzip.decompress(f.read())
    buf = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=buf,
                       mtime=0, compresslevel=9) as g:
        g.write(payload)
    with open(path, "wb") as f:
        f.write(buf.getvalue())


def strip_volatile(root: str) -> None:
    """Remove the files listed in :data:`VOLATILE` from *root*."""
    for rel in VOLATILE:
        p = os.path.join(root, *rel.split("/"))
        if os.path.isfile(p):
            os.remove(p)


def normalise_gz_tree(root: str) -> None:
    """Apply :func:`recompress_deterministic` to every ``.gz`` in *root*."""
    for dirpath, _, filenames in os.walk(root):
        for name in sorted(filenames):
            p = os.path.join(dirpath, name)
            if is_gzip(p):
                recompress_deterministic(p)


def read_tree(root: str, skip_dirs=SNAPSHOT_SKIP_DIRS) -> dict:
    """Map ``relative/path -> bytes`` for *root*, decompressed, sorted walk.

    Volatile files and *skip_dirs* are excluded.  The result is the
    comparison unit of this tool: identical dicts mean identical EZFIO
    state as far as any reader can tell.
    """
    out = {}
    root = os.path.abspath(root)
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames.sort()
        rel_dir = os.path.relpath(dirpath, root)
        if rel_dir == ".":
            rel_parts = []
        else:
            rel_parts = rel_dir.split(os.sep)
        if rel_parts and rel_parts[0] in skip_dirs:
            dirnames[:] = []
            continue
        for name in sorted(filenames):
            rel = "/".join(rel_parts + [name])
            if rel in VOLATILE:
                continue
            with open(os.path.join(dirpath, name), "rb") as f:
                data = f.read()
            if data[:2] == b"\x1f\x8b":
                data = gzip.decompress(data)
            out[rel] = data
    return out


def write_tree(tree: dict, dst: str, compress: bool = True) -> None:
    """Materialise a :func:`read_tree` map under *dst*.

    With *compress* (the default) a path ending in ``.gz`` is gzip
    compressed with ``mtime=0``, so the result is a valid EZFIO
    directory.  Snapshots pass ``compress=False`` and strip the suffix:
    they hold the decompressed payload, which is what D6 makes the oracle
    and what makes ``diff -ru`` on two snapshots readable.
    """
    for rel in sorted(tree):
        name = rel
        payload = tree[rel]
        if not compress and name.endswith(".gz"):
            name = name[:-3]
        p = os.path.join(dst, *name.split("/"))
        os.makedirs(os.path.dirname(p), exist_ok=True)
        if compress and rel.endswith(".gz"):
            buf = io.BytesIO()
            with gzip.GzipFile(filename="", mode="wb", fileobj=buf,
                               mtime=0, compresslevel=9) as g:
                g.write(payload)
            payload = buf.getvalue()
        with open(p, "wb") as f:
            f.write(payload)


def hash_tree(tree: dict) -> str:
    """Stable sha256 over a tree map (names and contents), for a quick id."""
    h = hashlib.sha256()
    for rel in sorted(tree):
        h.update(rel.encode())
        h.update(b"\0")
        h.update(tree[rel])
        h.update(b"\0")
    return h.hexdigest()


# --------------------------------------------------------------------------
# Corpus generation
# --------------------------------------------------------------------------

def create_fixture(fx: dict, corpus_inputs: str, work: str,
                   tools: str, env: dict) -> dict:
    """Run one fixture scenario and return its normalised tree."""
    for name in fx["inputs"]:
        shutil.copy2(os.path.join(corpus_inputs, name),
                     os.path.join(work, name))
    steps = [fx["create"]] + [list(s) for s in fx["post"]]
    for step in steps:
        r = run_step(step, work, env)
        if r["exit"] != 0:
            raise SystemExit(
                "golden_harness: fixture %s: %s exited %d\n%s"
                % (fx["name"], shlex.join(step), r["exit"],
                   r["stderr"].decode("utf-8", "replace"))
            )
    ezfio = os.path.join(work, "ezfio")
    if not os.path.isdir(ezfio):
        raise SystemExit("golden_harness: fixture %s produced no ezfio/"
                         % fx["name"])
    strip_volatile(ezfio)
    normalise_gz_tree(ezfio)
    return read_tree(ezfio)


def check_transcript(fixture_tree: dict, work: str,
                     tools: str, env: dict) -> dict:
    """Run ``qp_edit -c`` on a copy and report what it changed.

    Returns ``{"stdout", "stderr", "exit", "mutated", "added", "changed"}``
    where *mutated* lists ``"~ path"`` / ``"- path"`` / ``"+ path"``
    entries computed on decompressed content, *added* maps the ``"+ path"``
    entries (``.gz`` suffix kept) to their decompressed payload, and
    *changed* maps the ``"~ path"`` entries to their after payload.
    """
    dst = os.path.join(work, "ezfio")
    write_tree(fixture_tree, dst)
    r = run_step(["qp_edit", "-c", "ezfio"], work, env)
    after = read_tree(dst)
    before = fixture_tree
    mutated = []
    added = {}
    changed = {}
    for rel in sorted(set(before) | set(after)):
        if rel not in after:
            mutated.append("- " + rel)
        elif rel not in before:
            mutated.append("+ " + rel)
            added[rel] = after[rel]
        elif before[rel] != after[rel]:
            mutated.append("~ " + rel)
            changed[rel] = after[rel]
    shutil.rmtree(dst)
    return {"stdout": r["stdout"], "stderr": r["stderr"],
            "exit": r["exit"], "mutated": mutated, "added": added,
            "changed": changed}


def quantisation_probe(h2_tree: dict, work: str,
                       tools: str, env: dict) -> str:
    """Demonstrate the pre-fix ``qp_edit -c`` coordinate quantisation.

    Takes the H2 fixture, adds 5e-9 bohr to one coordinate -- half of the
    1.9e-8 bohr grid of the old ``%16.8f`` Angstrom rendering -- runs
    ``qp_edit -c`` and reports the value before and after.  This is the
    behaviour the post-``32c9899d`` contract forbids; it is recorded as
    evidence, not as a target.
    """
    dst = os.path.join(work, "ezfio")
    write_tree(h2_tree, dst)
    coord_rel = "nuclei/nucl_coord.gz"
    # decompressed text of the array; value index 4 is z of atom 2 for H2
    # (rank 2, dims (2,3), Fortran column-major: component-major).
    lines = h2_tree[coord_rel].decode().splitlines()
    dims = [int(x) for x in lines[1].split()]
    values = [float(x) for x in lines[2:2 + dims[0] * dims[1]]]
    idx = 4
    delta = 5e-9
    values[idx] += delta
    out = ["%3d" % 2,
           "".join("%20d " % d for d in dims)]
    out += ["%24.15e" % v for v in values]
    payload = ("\n".join(out) + "\n").encode()
    buf = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=buf,
                       mtime=0, compresslevel=9) as g:
        g.write(payload)
    with open(os.path.join(dst, *coord_rel.split("/")), "wb") as f:
        f.write(buf.getvalue())

    before = values[idx]
    r = run_step(["qp_edit", "-c", "ezfio"], work, env)
    after_tree = read_tree(dst)
    lines = after_tree[coord_rel].decode().splitlines()
    values_after = [float(x) for x in lines[2:2 + dims[0] * dims[1]]]
    after = values_after[idx]

    # Analytic prediction of the pre-fix round trip: bohr -> Angstrom
    # rendered with %16.8f -> parsed -> multiplied by 1/0.52917721067121.
    bohr_to_angstrom = 0.52917721067121
    angstrom_to_bohr = 1.0 / bohr_to_angstrom

    def round_trip(x: float) -> float:
        return float("%.8f" % (x * bohr_to_angstrom)) * angstrom_to_bohr

    predicted = round_trip(before)
    text = [
        "Pre-fix qp_edit -c coordinate quantisation (evidence, not target).",
        "",
        "Source: the H2 fixture.  One coordinate (nuclei/nucl_coord,",
        "dims (2,3), value index 4) is shifted by +5e-9 bohr, then",
        "`qp_edit -c` is run with the reference build (be4b71c1, i.e.",
        "before commit 32c9899d).",
        "",
        "command: qp_edit -c ezfio",
        "exit: %d" % r["exit"],
        "stdout: %r" % r["stdout"],
        "stderr: %r" % r["stderr"],
        "",
        "before:    %.17g bohr" % before,
        "predicted: %.17g bohr   (bohr -> Angstrom %.8f -> bohr)"
        % (predicted, before * bohr_to_angstrom),
        "after:     %.17g bohr" % after,
        "drift:     %.3g bohr" % (after - before),
        "",
        "The value written back is the prediction of the %16.8f Angstrom",
        "round trip, not the value that was put in: checking the input",
        "rewrote it.  32c9899d removes that write (and prints %.17g).",
    ]
    shutil.rmtree(dst)
    return "\n".join(text) + "\n"


def last_element(basis_path: str) -> str:
    """Last element name of a basis file (all upper case, first column)."""
    last = None
    with open(basis_path, "r", errors="replace") as f:
        for line in f:
            tok = line.strip()
            if ELEMENT_RE.match(tok):
                last = tok
    if last is None:
        raise SystemExit("golden_harness: no element in %s" % basis_path)
    return last


def md5_rows(qp_root: str, tools: str, env: dict, work: str) -> list:
    """``(basis, element, ao_md5, ao_extra_md5)`` for :data:`BASES`."""
    rows = []
    basis_dir = os.path.join(qp_root, "data", "basis")
    for basis in BASES:
        basis_path = os.path.join(basis_dir, basis)
        if not os.path.isfile(basis_path):
            raise SystemExit("golden_harness: missing basis %s" % basis_path)
        element = last_element(basis_path)
        m = PSEUDO_RE.search(basis)
        pseudo = m.group(1) if m else None
        d = os.path.join(work, "md5-" + basis.replace("/", "_"))
        os.makedirs(d, exist_ok=True)
        # Two atoms at 0.37 Angstrom on the z axis: enough to build the
        # basis and the md5, no integrals are computed by qp_create_ezfio.
        with open(os.path.join(d, "mol.xyz"), "w") as f:
            f.write("2\n%s2\n%s 0.0 0.0 -0.37\n%s 0.0 0.0 0.37\n"
                    % (element, element, element))
        argv = ["qp_create_ezfio", "-b", basis]
        if pseudo:
            argv += ["-p", pseudo]
        argv += ["-o", "ezfio", "mol.xyz"]
        r = run_step(argv, d, env)
        if r["exit"] != 0:
            raise SystemExit("golden_harness: md5 %s exited %d\n%s"
                             % (basis, r["exit"],
                                r["stderr"].decode("utf-8", "replace")))
        ez = os.path.join(d, "ezfio")

        def scalar(rel):
            with open(os.path.join(ez, *rel.split("/")), "r") as f:
                return f.read().strip()

        with open(basis_path, "rb") as f:
            basis_sha = hashlib.sha256(f.read()).hexdigest()
        rows.append({
            "basis": basis,
            "pseudo": pseudo,
            "element": element,
            "n_ao": int(scalar("ao_basis/ao_num")),
            "ao_md5": scalar("ao_basis/ao_md5"),
            "ao_extra_md5": scalar("ao_extra_basis/ao_extra_md5"),
            "basis_file_sha256": basis_sha,
        })
        shutil.rmtree(d)
    rows.sort(key=lambda r: r["basis"])
    return rows


def cmd_generate(args) -> int:
    """Write ``fixtures/``, ``transcripts/``, ``md5/`` and ``manifest.json``."""
    corpus = os.path.abspath(args.corpus)
    inputs = os.path.join(corpus, "inputs")
    if not os.path.isdir(inputs):
        raise SystemExit("golden_harness: no inputs/ in %s" % corpus)
    tools = find_tools(args.qp_root, args.qp_bin)
    env = tool_env(args.qp_root, tools)

    with tempfile.TemporaryDirectory(prefix="qp-golden-") as tmp:
        fixtures = {}
        for fx in FIXTURES:
            work = os.path.join(tmp, fx["name"])
            os.makedirs(work)
            fixtures[fx["name"]] = create_fixture(fx, inputs, work, tools, env)
            print("fixture  %-6s %4d files  %s"
                  % (fx["name"], len(fixtures[fx["name"]]),
                     hash_tree(fixtures[fx["name"]])[:16]))

        transcripts = {}
        for fx in FIXTURES:
            work = os.path.join(tmp, "check-" + fx["name"])
            os.makedirs(work)
            transcripts[fx["name"]] = check_transcript(
                fixtures[fx["name"]], work, tools, env)
            print("check    %-6s exit %d, %d file(s) rewritten"
                  % (fx["name"], transcripts[fx["name"]]["exit"],
                     len(transcripts[fx["name"]]["mutated"])))

        work = os.path.join(tmp, "probe")
        os.makedirs(work)
        probe = quantisation_probe(fixtures["h2"], work, tools, env)
        print("probe    quantisation: " + probe.splitlines()[-5].strip())

        rows = md5_rows(args.qp_root, tools, env, tmp)
        print("md5      %d bases" % len(rows))

        info = git_info(args.qp_root)

    # Everything succeeded: write into the corpus directory.
    fixtures_dir = os.path.join(corpus, "fixtures")
    transcripts_dir = os.path.join(corpus, "transcripts")
    md5_dir = os.path.join(corpus, "md5")
    for d in (fixtures_dir, transcripts_dir, md5_dir):
        if os.path.isdir(d):
            shutil.rmtree(d)
        os.makedirs(d)
    for name, tree in fixtures.items():
        write_tree(tree, os.path.join(fixtures_dir, name))
    for name, tr in transcripts.items():
        d = os.path.join(transcripts_dir, name)
        os.makedirs(d)
        with open(os.path.join(d, "stdout"), "wb") as f:
            f.write(tr["stdout"])
        with open(os.path.join(d, "stderr"), "wb") as f:
            f.write(tr["stderr"])
        with open(os.path.join(d, "exit"), "w") as f:
            f.write("%d\n" % tr["exit"])
        with open(os.path.join(d, "mutated.txt"), "w") as f:
            f.write("".join(m + "\n" for m in tr["mutated"]))
        if tr["added"]:
            # What the reference build created, decompressed and without
            # the .gz suffix (same convention as a snapshot tree).  The
            # normative list of files a *post-fix* build may create lives
            # in contracts/check-mutated.json; this tree is the content
            # evidence for it, including the two files that only the
            # section writer of a pre-fix build produced.
            write_tree(tr["added"], os.path.join(d, "added"), compress=False)
    with open(os.path.join(md5_dir, "basis-md5.jsonl"), "w") as f:
        for row in rows:
            f.write(json.dumps(row, sort_keys=True) + "\n")
    probe_dir = os.path.join(transcripts_dir, "pre-fix-quantisation")
    os.makedirs(probe_dir)
    with open(os.path.join(probe_dir, "probe.log"), "w") as f:
        f.write(probe)

    manifest = {
        "format": 1,
        "generator": "rust/tools/golden-harness/golden_harness.py",
        "reference": dict(info, ezfio_version=read_ezfio_version(
            os.path.join(corpus, "fixtures", "h2"))),
        "env": {"OMP_NUM_THREADS": "1"},
        "normalisation": {
            "excluded_paths": list(VOLATILE),
            "snapshot_skipped_dirs": list(SNAPSHOT_SKIP_DIRS),
            "gzip": "recompressed with mtime=0, level=9; the decompressed "
                    "payload is the oracle",
            "output_noise_patterns": list(NOISE),
        },
        "fixtures": [dict(f, scenarios=["create", "check"])
                     for f in FIXTURES],
        "contracts": {
            "check_mode": "contracts/check-mode.md",
            "check_mutated": "contracts/check-mutated.json",
        },
        "md5": {"file": "md5/basis-md5.jsonl", "bases": list(BASES),
                "count": len(rows)},
    }
    with open(os.path.join(corpus, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2, sort_keys=False)
        f.write("\n")
    print("wrote    %s" % corpus)
    return 0


def read_ezfio_version(fixture_dir: str) -> str | None:
    """Value of ``.version`` at the root of an EZFIO directory."""
    p = os.path.join(fixture_dir, ".version")
    if not os.path.isfile(p):
        return None
    with open(p) as f:
        return f.read().strip()


# --------------------------------------------------------------------------
# Scenario replay
# --------------------------------------------------------------------------

def load_manifest(corpus: str) -> dict:
    p = os.path.join(corpus, "manifest.json")
    if not os.path.isfile(p):
        raise SystemExit("golden_harness: no %s; run `generate` first" % p)
    with open(p) as f:
        return json.load(f)


def scenario_steps(fx: dict, scenario: str) -> list:
    """Command lines of a scenario, with ``ezfio`` already in place."""
    if scenario == "create":
        return [list(fx["create"])] + [list(s) for s in fx["post"]]
    if scenario == "check":
        return [["qp_edit", "-c", "ezfio"]]
    raise SystemExit("golden_harness: unknown scenario %r" % scenario)


def replay(corpus: str, fx: dict, scenario: str, tools: str,
           env: dict, work: str, snapshot: str) -> None:
    """Replay one ``(fixture, scenario)`` and write its snapshot."""
    if scenario == "check":
        src = os.path.join(corpus, "fixtures", fx["name"])
        if not os.path.isdir(src):
            raise SystemExit("golden_harness: no fixture %s" % src)
        shutil.copytree(src, os.path.join(work, "ezfio"))
    else:
        for name in fx["inputs"]:
            shutil.copy2(os.path.join(corpus, "inputs", name),
                         os.path.join(work, name))

    steps = scenario_steps(fx, scenario)
    results = [run_step(step, work, env) for step in steps]

    os.makedirs(snapshot)
    with open(os.path.join(snapshot, "cmd.txt"), "w") as f:
        for step in steps:
            f.write(shlex.join(step) + "\n")
    # A scenario may have several steps; concatenate their streams so that
    # a diff still points at the step that changed.  Exit is the first
    # non-zero status, or 0.
    stdout = b"".join(r["stdout"] for r in results)
    stderr = b"".join(r["stderr"] for r in results)
    exit_code = next((r["exit"] for r in results if r["exit"] != 0), 0)
    with open(os.path.join(snapshot, "stdout"), "wb") as f:
        f.write(stdout)
    with open(os.path.join(snapshot, "stderr"), "wb") as f:
        f.write(stderr)
    with open(os.path.join(snapshot, "exit"), "w") as f:
        f.write("%d\n" % exit_code)
    write_tree(read_tree(os.path.join(work, "ezfio")),
               os.path.join(snapshot, "tree"), compress=False)


def cmd_run(args) -> int:
    """Replay the manifest with one tool set into ``--out``."""
    corpus = os.path.abspath(args.corpus)
    out = os.path.abspath(args.out)
    manifest = load_manifest(corpus)
    tools = find_tools(args.qp_root, args.qp_bin)
    env = tool_env(args.qp_root, tools)
    wanted = set(args.fixtures.split(",")) if args.fixtures else None
    only = set(args.only.split(",")) if args.only else None

    if os.path.isdir(out):
        if not args.force:
            raise SystemExit("golden_harness: %s exists; pass --force" % out)
        shutil.rmtree(out)
    os.makedirs(out)
    with open(os.path.join(out, "meta.json"), "w") as f:
        json.dump({"qp_root": os.path.abspath(args.qp_root),
                   "qp_bin": tools,
                   "reference": manifest.get("reference", {}),
                   **git_info(args.qp_root)}, f, indent=2)
        f.write("\n")

    with tempfile.TemporaryDirectory(prefix="qp-harness-") as tmp:
        for fx in manifest["fixtures"]:
            if wanted and fx["name"] not in wanted:
                continue
            for scenario in fx["scenarios"]:
                if only and scenario not in only:
                    continue
                work = os.path.join(tmp, "%s-%s" % (fx["name"], scenario))
                os.makedirs(work)
                snapshot = os.path.join(out, fx["name"], scenario)
                replay(corpus, fx, scenario, tools, env, work, snapshot)
                print("%-6s %-7s exit %s"
                      % (fx["name"], scenario,
                         open(os.path.join(snapshot, "exit")).read().strip()))
    return 0


# --------------------------------------------------------------------------
# Diff
# --------------------------------------------------------------------------

def read_snapshot(root: str) -> dict:
    """Map ``fixture/scenario -> {"cmd", "stdout", "stderr", "exit", "tree"}``."""
    out = {}
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames.sort()
        if "tree" not in dirnames:
            continue
        rel = os.path.relpath(dirpath, root).replace(os.sep, "/")
        entry = {"tree": read_tree(os.path.join(dirpath, "tree"),
                                   skip_dirs=())}
        for name in ("cmd.txt", "stdout", "stderr", "exit"):
            p = os.path.join(dirpath, name)
            entry[name] = open(p, "rb").read() if os.path.isfile(p) else b""
        out[rel] = entry
    return out


def first_diff(a: bytes, b: bytes) -> str:
    """Human-readable location of the first difference between two blobs."""
    al = a.decode("utf-8", "replace").splitlines()
    bl = b.decode("utf-8", "replace").splitlines()
    for i, (x, y) in enumerate(zip(al, bl), 1):
        if x != y:
            return ("line %d:\n    A: %s\n    B: %s"
                    % (i, x[:200], y[:200]))
    if len(al) != len(bl):
        return "A has %d lines, B has %d lines" % (len(al), len(bl))
    return "bytes differ but text is identical (encoding?)"


def cmd_diff(args) -> int:
    """Compare two snapshots; exit 1 on any difference."""
    a = read_snapshot(os.path.abspath(args.snapshot_a))
    b = read_snapshot(os.path.abspath(args.snapshot_b))
    ignore = [re.compile(p) for p in args.ignore]
    keys = sorted(set(a) | set(b))
    failures = 0
    for key in keys:
        if any(r.search(key) for r in ignore):
            continue
        if key not in a:
            print("%s: only in B" % key)
            failures += 1
            continue
        if key not in b:
            print("%s: only in A" % key)
            failures += 1
            continue
        ea, eb = a[key], b[key]
        for field in ("cmd.txt", "exit", "stdout", "stderr"):
            if ea[field] != eb[field]:
                print("%s/%s differs: %s"
                      % (key, field, first_diff(ea[field], eb[field])))
                failures += 1
        ta, tb = ea["tree"], eb["tree"]
        for rel in sorted(set(ta) | set(tb)):
            full = "%s/tree/%s" % (key, rel)
            if any(r.search(full) for r in ignore):
                continue
            if rel not in tb:
                print("%s: only in A" % full)
                failures += 1
            elif rel not in ta:
                print("%s: only in B" % full)
                failures += 1
            elif ta[rel] != tb[rel]:
                print("%s differs: %s" % (full, first_diff(ta[rel], tb[rel])))
                failures += 1
    if failures:
        print("\n%d difference(s)" % failures)
        return 1
    print("identical (%d scenario(s))" % len(keys))
    return 0


def read_bytes(path: str) -> bytes:
    """Bytes of a file, or ``b""`` when it does not exist."""
    if not os.path.isfile(path):
        return b""
    with open(path, "rb") as f:
        return f.read()


def load_contract(corpus: str) -> dict:
    """The normative ``qp_edit -c`` contract (``contracts/check-mutated.json``).

    It is hand-written, not generated: it is the *target* the Rust port
    must meet, transcribed from the source changes of commit 32c9899d
    (``15e725aa`` on ``origin/qp-edit-fix``), not observed behaviour.
    Besides the created-file lists (``post-fix`` / ``pre-fix``), the
    ``post-fix-changed`` list names the existing files the w12 one-time
    ``nucl_coord`` normalisation of c2h2/b2 may change (see the file's
    ``description``); the expected decompressed payloads live in
    ``transcripts/<name>/changed/``.  The key is a format contract (its
    absence is a load error), but the pin itself is a **one-way
    allowance**: a run that touches one of the listed files must match
    the stored payload; a run that leaves it untouched passes (see
    ``verify``).
    """
    p = os.path.join(corpus, "contracts", "check-mutated.json")
    if not os.path.isfile(p):
        raise SystemExit("golden_harness: no %s; the corpus is incomplete"
                         % p)
    with open(p) as f:
        contract = json.load(f)
    for mode in CHECK_CONTRACTS:
        if mode not in contract:
            raise SystemExit("golden_harness: %s has no %r contract"
                             % (p, mode))
    if "post-fix-changed" not in contract:
        raise SystemExit("golden_harness: %s has no 'post-fix-changed' "
                         "list (the pinned normalisation)" % p)
    return contract


def expected_payload(tdir: str, rel: str, subdir: str):
    """The stored expected payload for *rel* under ``<tdir>/<subdir>/``,
    or ``None`` when nothing is stored there.

    The storage convention (both ``added/`` and ``changed/``): the
    decompressed payload, ``.gz`` suffix stripped.
    """
    rel_path = rel[:-3] if rel.endswith(".gz") else rel
    p = os.path.join(tdir, subdir, *rel_path.split("/"))
    return p if os.path.isfile(p) else None


def verify_create(corpus: str, fx: dict, work: str,
                  tools: str, env: dict, contract: dict) -> int:
    """Replay a ``create`` scenario and compare with ``fixtures/<name>``.

    The create output must equal the committed fixture, with the one-way
    allowance of ``contract["post-fix-changed"][<name>]`` (the w12
    one-time ``nucl_coord`` normalisation, c2h2/b2): when the
    normalisation is applied (the port's create tail runs the same
    divide-by-K pipeline as ``-c``), the difference must be in exactly
    those files and their produced payload must equal the stored
    expected file (``transcripts/<name>/changed/``); a toolset that does
    not apply it (the reference build in CI, whose multiply pipeline is
    the corpus's fixed point) produces the fixture bytes and passes.
    Any other difference is a failure.  The negative side -- that the
    port's create does apply the normalisation -- is enforced by the
    port's own test suite, not by this tool.
    """
    name = fx["name"]
    produced = create_fixture(fx, os.path.join(corpus, "inputs"),
                              work, tools, env)
    expected = read_tree(os.path.join(corpus, "fixtures", name))
    tdir = os.path.join(corpus, "transcripts", name)
    pinned = {}
    for rel in contract["post-fix-changed"].get(name, []):
        p = expected_payload(tdir, rel, "changed")
        pinned[rel] = read_bytes(p) if p else None
    failures = 0
    pinned_verified = 0
    for rel in sorted(set(expected) | set(produced)):
        if rel not in produced:
            failures += 1
            print("%s: missing from replay" % rel)
        elif rel not in expected:
            failures += 1
            print("%s: not in committed fixture" % rel)
        elif expected[rel] != produced[rel]:
            if pinned.get(rel) is not None and produced[rel] == pinned[rel]:
                pinned_verified += 1
                continue
            failures += 1
            print("%s: %s" % (rel, first_diff(expected[rel], produced[rel])))
    if not failures:
        line = "%-6s create  ok  %4d files  %s" \
               % (name, len(produced), hash_tree(produced)[:16])
        if pinned_verified:
            line += ", %d changed as pinned" % pinned_verified
        print(line)
        return 0
    print("%-6s create  FAILED" % name)
    return 1


def verify_check(corpus: str, fx: dict, work: str, tools: str,
                 env: dict, contract: dict, mode: str) -> int:
    """Replay ``qp_edit -c`` and compare it with the transcript and contract.

    Four things are checked, all on decompressed content:

    * stdout, stderr and the exit status against ``transcripts/<name>/``
      (the reference build's transcript);
    * the set of files created, against ``contract[mode][<name>]``;
    * the content of every created file for which
      ``transcripts/<name>/added/`` holds the reference payload;
    * the existing files that were touched: in ``post-fix`` mode a
      touched file is a failure unless it is in
      ``contract["post-fix-changed"][<name>]`` (the w12 one-time
      ``nucl_coord`` normalisation) with its payload equal to the stored
      expected file (``transcripts/<name>/changed/``, when one is
      stored).  The pin is a **one-way allowance**: a toolset that does
      not apply the normalisation (the reference build in CI, whose
      multiply pipeline is the corpus's own fixed point) touches nothing
      and passes -- the pin is not a requirement to change the file.
      The negative side (that the port does change it) is enforced by
      the port's own test suite
      (``c2_no_pre_existing_file_is_rewritten``, which also asserts the
      change is real).  In ``pre-fix`` mode, and for deletions, any
      touch is a failure.
    """
    name = fx["name"]
    before = read_tree(os.path.join(corpus, "fixtures", name))
    tr = check_transcript(before, work, tools, env)
    tdir = os.path.join(corpus, "transcripts", name)
    failures = 0

    for field in ("stdout", "stderr"):
        want = read_bytes(os.path.join(tdir, field))
        if tr[field] != want:
            failures += 1
            print("%s: %s differs: %s"
                  % (name, field, first_diff(want, tr[field])))
    want_exit = int(read_bytes(os.path.join(tdir, "exit")).decode().strip())
    if tr["exit"] != want_exit:
        failures += 1
        print("%s: exit %d, transcript says %d"
              % (name, tr["exit"], want_exit))

    expected = set(contract[mode][name])
    created = [m[2:] for m in tr["mutated"] if m.startswith("+ ")]
    touched = [m for m in tr["mutated"] if not m.startswith("+ ")]
    for rel in sorted(expected - set(created)):
        failures += 1
        print("%s: not created: %s" % (name, rel))
    for rel in sorted(set(created) - expected):
        failures += 1
        print("%s: created, but not in the %s contract: %s"
              % (name, mode, rel))

    # The pinned one-time normalisation (w12) is a one-way allowance: in
    # post-fix mode a "~" touch is allowed exactly for the files of
    # contract["post-fix-changed"][<name>], with the stored expected
    # payload (when one is stored).  An unexpected touch, a deletion or
    # an unexpected payload is a failure; a toolset that does not apply
    # the normalisation (the reference build in CI, whose multiply
    # pipeline is the corpus's own fixed point) touches nothing and
    # passes.  The negative side -- that the port DOES apply it -- is
    # enforced by the port's own test suite
    # (c2_no_pre_existing_file_is_rewritten), not by this tool.
    pinned = (contract["post-fix-changed"].get(name, [])
              if mode == "post-fix" else [])
    pinned_verified = 0
    for entry in touched:
        if not entry.startswith("~ ") or entry[2:] not in pinned:
            failures += 1
            print("%s: check mode touched an existing file: %s"
                  % (name, entry))
            continue
        rel = entry[2:]
        p = expected_payload(tdir, rel, "changed")
        if p is not None and tr["changed"].get(rel) != read_bytes(p):
            failures += 1
            print("%s: %s differs from the pinned payload: %s"
                  % (name, rel,
                     first_diff(read_bytes(p), tr["changed"].get(rel, b""))))
        else:
            pinned_verified += 1

    # Content of the created files, where the reference payload is known.
    compared = 0
    for rel in sorted(expected):
        p = expected_payload(tdir, rel, "added")
        if p is None:
            continue
        want = read_bytes(p)
        if tr["added"].get(rel) != want:
            failures += 1
            print("%s: %s differs: %s"
                  % (name, rel, first_diff(want, tr["added"].get(rel, b""))))
        compared += 1

    if not failures:
        line = "%-6s check   ok  %d created, %d compared byte for byte" \
               % (name, len(expected), compared)
        if pinned_verified:
            line += ", %d changed as pinned" % pinned_verified
        print(line)
    else:
        print("%-6s check   FAILED (%s contract)" % (name, mode))
    return 1 if failures else 0


def cmd_verify(args) -> int:
    """Replay the requested scenarios and gate them against the corpus."""
    corpus = os.path.abspath(args.corpus)
    manifest = load_manifest(corpus)
    tools = find_tools(args.qp_root, args.qp_bin)
    env = tool_env(args.qp_root, tools)
    scenarios = set(args.only.split(","))
    unknown = scenarios - {"create", "check"}
    if unknown:
        raise SystemExit("golden_harness: unknown scenario(s) %s"
                         % ", ".join(sorted(unknown)))
    # The contract gates both scenarios: ``check`` against the created
    # and the pinned changed files, ``create`` against the same pinned
    # changed files (create's tail runs the same normalisation).
    contract = load_contract(corpus) if scenarios else None
    failures = 0
    with tempfile.TemporaryDirectory(prefix="qp-verify-") as tmp:
        for fx in manifest["fixtures"]:
            for scenario in ("create", "check"):
                if scenario not in scenarios:
                    continue
                work = os.path.join(tmp, "%s-%s" % (fx["name"], scenario))
                os.makedirs(work)
                if scenario == "create":
                    failures += verify_create(corpus, fx, work, tools, env,
                                              contract)
                else:
                    failures += verify_check(corpus, fx, work, tools, env,
                                             contract, args.contract)
    return 1 if failures else 0


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def main(argv) -> int:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)

    def common(sp):
        sp.add_argument("--qp-root", default=os.environ.get("QP_ROOT"),
                        help="QP build tree holding data/basis and the "
                             "qp_* binaries (default: $QP_ROOT)")
        sp.add_argument("--qp-bin", default=None,
                        help="directory holding the qp_* binaries "
                             "(default: $QP_ROOT/ocaml, then $QP_ROOT/bin)")
        sp.add_argument("--corpus", default=os.path.join(
                            HERE, "..", "..", "tests", "golden"),
                        help="golden-corpus directory (default: "
                             "rust/tests/golden next to this tool)")

    g = sub.add_parser("generate", help="write the golden corpus")
    common(g)
    g.set_defaults(func=cmd_generate)

    r = sub.add_parser("run", help="replay scenarios into a snapshot")
    common(r)
    r.add_argument("--out", required=True, help="snapshot directory")
    r.add_argument("--force", action="store_true",
                   help="replace --out if it exists")
    r.add_argument("--fixtures", default=None,
                   help="comma-separated subset of fixture names")
    r.add_argument("--only", default=None,
                   help="comma-separated subset of scenarios (create,check)")
    r.set_defaults(func=cmd_run)

    d = sub.add_parser("diff", help="compare two snapshots")
    d.add_argument("snapshot_a")
    d.add_argument("snapshot_b")
    d.add_argument("--ignore", action="append", default=[],
                   help="regular expression on fixture/scenario/path to skip "
                        "(repeatable)")
    d.set_defaults(func=cmd_diff)

    v = sub.add_parser("verify", help="replay and gate against the corpus")
    common(v)
    v.add_argument("--only", default="create,check",
                   help="comma-separated subset of scenarios (default: "
                        "create,check)")
    v.add_argument("--contract", default="post-fix", choices=CHECK_CONTRACTS,
                   help="which qp_edit -c contract to gate against: "
                        "post-fix (the Rust target, default) or pre-fix "
                        "(what the reference OCaml build does)")
    v.set_defaults(func=cmd_verify)

    args = p.parse_args(argv)
    if args.func is not cmd_diff and not getattr(args, "qp_root", None):
        raise SystemExit("golden_harness: --qp-root (or $QP_ROOT) is required")
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
