#!/usr/bin/env python3
"""
qp_rag: build and interrogate the QP2 code-model RAG.

`build` wipes $QP_ROOT/share/rag (gitignored) and rebuilds the whole
corpus:

  1. picks the irpf90 binary ($RAG_IRPF90 or external/irpf90/bin/irpf90)
     and requires the `-J` capability (irpf90 >= 2.3.3);
  2. reads config/build.ninja for the (module, include_dir) build edges
     (22) and the exact $IRPF90_FLAGS of the build (-J output depends on
     the flags, so the verbatim build set is used);
  3. runs `irpf90 <include_dir> <IRPF90_FLAGS> -J` in each
     $QP_ROOT/src/<module> with the build environment (QP_ROOT, QP_EZFIO
     and the PYTHONPATH of etc/paths.rc, needed by the embedded
     `! shell` blocks) and stores the dump in
     share/rag/per-module/<module>.json;
  4. merges the per-module views into one tree-wide model (list fields
     are sorted unions, the view-dependent touched/freed/read/written
     flags are OR-merged, the identity fields must agree across views,
     file paths are tree-relative) and writes share/rag/qp2.json;
  5. disambiguates cross-module name collisions: the colliding records
     only are renamed name -> name__<module> (module = the dir of the
     record's own file); the rename table goes to the report;
  6. indexes the model in the SQLite FTS5/BM25 database
     share/rag/qp2_rag.db (one chunk per record + one project overview,
     doc_id = qp2:<kind>:<name>) and writes share/rag/merge-report.json.

`query`, `graph` and `stats` only read share/rag: query is BM25 over the
index, graph is a BFS over the entity needs/needed_by edges, stats gives
counts and the external API surface. If the corpus is missing: ninja rag.

The `ninja rag` target at the repo root runs `build`; ninja rebuilds the
corpus automatically when any module source (or this tool) changes. It
is not part of the default build.

Usage:
    python3 scripts/rag/qp_rag.py build
    python3 scripts/rag/qp_rag.py query "PT2 selection" -k 5
    python3 scripts/rag/qp_rag.py query "bitmask" --kind routine
    python3 scripts/rag/qp_rag.py graph <entity> --depth 3
    python3 scripts/rag/qp_rag.py stats

Stdlib only, Python 3.9+.
"""

import argparse
import difflib
import json
import os
import re
import shutil
import sqlite3
import subprocess
import sys
from collections import Counter
from pathlib import Path

QP_ROOT = Path(__file__).resolve().parents[2]
RAGDIR = QP_ROOT / "share" / "rag"
PER_MODULE = RAGDIR / "per-module"
MODEL = RAGDIR / "qp2.json"
DB = RAGDIR / "qp2_rag.db"
REPORT = RAGDIR / "merge-report.json"
PROJECT = "qp2"

# Merge semantics (verified over the 22-view QP2 corpus, 2026-10-31):
# identity fields must agree across views (a mismatch is a hard error);
# list fields are sorted unions; the view-dependent lifecycle flags are
# OR-merged (monotone w.r.t. visible-set inclusion, 1.2M run-pair check).
ENTITY_NONLIST = ["same_as", "type", "dims", "is_main", "doc"]
ENTITY_LISTS = ["needs", "needed_by", "calls"]
ENTITY_VIEWFLAGS = ["touched", "freed", "read", "written"]
ROUTINE_NONLIST = ["kind", "signature", "doc"]
ROUTINE_LISTS = ["needs", "to_provide", "touches", "calls",
                 "external_calls", "called_by"]

STOP = set("the and for with from that this are was were can may have has had not its"
           " into over under per via within between more less also than then using used"
           " use uses new some such other".split())


def die(msg, code=1):
    print("error: " + msg, file=sys.stderr)
    sys.exit(code)


# ~#~#~#~#~#~ #
# b u i l d   #
# ~#~#~#~#~#~ #

def parse_build_ninja(path):
    """Return (edges, flags) from config/build.ninja: one
    (module, include_dir) pair per build_irpf90.ninja edge and the
    verbatim $IRPF90_FLAGS value."""
    edge_re = re.compile(
        r"^build\s+\$QP_ROOT/src/(\w+)/irpf90_entities\s+\$QP_ROOT/src/\w+/tags\s+"
        r"\$QP_ROOT/src/\w+/IRPF90_temp/build\.ninja:\s+build_irpf90\.ninja\s")
    with open(path, encoding="utf-8") as f:
        lines = f.readlines()
    flags = None
    edges = []
    for i, line in enumerate(lines):
        m = re.match(r"^IRPF90_FLAGS\s*=\s*(\S.*)$", line)
        if m:
            flags = m.group(1).strip()
            continue
        m = edge_re.match(line)
        if not m:
            continue
        module = m.group(1)
        vars_ = {}
        for j in range(1, 6):
            idx = i + j
            if idx >= len(lines):
                die("edge %s: missing variable line %d" % (module, j))
            vm = re.match(r"^   +(\w+) = (.*)$", lines[idx])
            if not vm:
                die("edge %s: line %d is not an indented variable: %r"
                    % (module, idx + 1, lines[idx]))
            vars_[vm.group(1)] = vm.group(2).strip()
        if list(vars_.keys()) != ["module_abs", "module_rel", "SRC",
                                  "OBJ", "include_dir"]:
            die("edge %s: unexpected variable order %s"
                % (module, list(vars_.keys())))
        if vars_["module_rel"] != module or \
                vars_["module_abs"] != "$QP_ROOT/src/" + module:
            die("edge %s: inconsistent module_abs/module_rel" % module)
        include_dir = vars_["include_dir"].replace("$QP_ROOT", str(QP_ROOT))
        if not (QP_ROOT / "src" / module).is_dir():
            die("edge %s: module directory src/%s missing" % (module, module))
        edges.append({"module": module, "include_dir": include_dir})
    if not flags:
        die("config/build.ninja: no IRPF90_FLAGS variable")
    if not edges:
        die("config/build.ninja: no build_irpf90.ninja edges found")
    if len(edges) != 22:
        print("warning: expected 22 irpf90 edges, found %d" % len(edges),
              file=sys.stderr)
    return edges, flags.split()


def find_irpf90():
    """The irpf90 binary, capability-checked for -J (>= 2.3.3)."""
    irpf90 = os.environ.get("RAG_IRPF90") or str(
        QP_ROOT / "external" / "irpf90" / "bin" / "irpf90")
    if not (os.path.isfile(irpf90) and os.access(irpf90, os.X_OK)):
        die("irpf90 %s not found or not executable; set "
            "RAG_IRPF90=/path/to/bin/irpf90 or bump the external/irpf90 pin"
            % irpf90)
    p = subprocess.run([irpf90, "--help"], capture_output=True, text=True)
    help_text = (p.stdout or "") + (p.stderr or "")
    if "-J" not in help_text and "--json" not in help_text:
        die("irpf90 %s lacks -J (needs >= 2.3.3); set "
            "RAG_IRPF90=/path/to/bin/irpf90 or bump the external/irpf90 pin"
            % irpf90)
    return irpf90


def parse_env():
    """The parse environment: the caller's environment + QP_ROOT,
    QP_EZFIO and the PYTHONPATH of etc/paths.rc (the embedded `! shell`
    blocks run Python modules such as qp_path / generate_h_apply /
    perturbations during the parse)."""
    env = dict(os.environ)
    ezfio = str(QP_ROOT / "external" / "ezfio")
    env["QP_ROOT"] = str(QP_ROOT)
    env["QP_EZFIO"] = ezfio
    env["PYTHONPATH"] = ":".join([
        ezfio + "/Python",
        ezfio + "/Python",
        str(QP_ROOT / "scripts"),
        str(QP_ROOT / "scripts" / "ezfio_interface"),
        str(QP_ROOT / "scripts" / "utility"),
        str(QP_ROOT / "scripts" / "module"),
        str(QP_ROOT / "scripts" / "compilation"),
        str(QP_ROOT / "external" / "Python"),
        str(QP_ROOT / "external" / "Python" / "resultsFile"),
        os.environ.get("PYTHONPATH", ""),
    ])
    return env


def shell_leftovers():
    """The transient <file>_shell_<N> files irpf90 writes next to an
    embedded `! shell` block (removed by irpf90 itself on success)."""
    return sorted(p for p in (QP_ROOT / "src").rglob("*_shell_*")
                  if p.is_file())


def _cleanup_shell_leftovers(before):
    for p in set(shell_leftovers()) - before:
        try:
            p.unlink()
        except OSError:
            pass


def run_module(irpf90, edge, flags, env):
    """One `irpf90 <include_dir> <IRPF90_FLAGS> -J` in src/<module>;
    the dump goes to per-module/<module>.json."""
    m = edge["module"]
    out = PER_MODULE / (m + ".json")
    args = [irpf90] + edge["include_dir"].split() + flags + ["-J", str(out)]
    before = set(shell_leftovers())
    try:
        p = subprocess.run(args, cwd=str(QP_ROOT / "src" / m), env=env,
                           capture_output=True, text=True, timeout=600)
    except subprocess.TimeoutExpired:
        _cleanup_shell_leftovers(before)
        die("module %s: irpf90 -J timed out (600 s)" % m)
    if p.returncode != 0:
        with open(PER_MODULE / (m + ".stderr"), "w", encoding="utf-8") as fh:
            fh.write(p.stderr)
        _cleanup_shell_leftovers(before)
        tail = (p.stderr or "").strip().splitlines()[-20:]
        die("module %s: irpf90 -J failed (rc=%d)\n%s"
            % (m, p.returncode, "\n".join(tail)))
    with open(out, encoding="utf-8") as fh:
        d = json.load(fh)
    print("  %s: %d files, %d entities, %d routines"
          % (m, len(d.get("files") or []), len(d.get("entities") or []),
             len(d.get("routines") or [])))


def mangle(rel):
    """irpf90 fmodule formula applied to a path relative to src/
    (i.e. '<module>/<file>')."""
    return (rel.replace("/", "__").split(".irp.f")[0]
            .replace(".", "Dot") + "_mod")


def norm_file(run_module, entry):
    """Tree-relative path for a per-run files entry (entries are
    relative to src/<run_module>: own files are bare names, dep files
    carry their dir component)."""
    if "/" in entry:
        dep, f = entry.split("/", 1)
        return "src/%s/%s" % (dep, f)
    return "src/%s/%s" % (run_module, entry)


def module_of(rel):
    """Owning module dir of a tree-relative 'src/<module>/<file>' path."""
    parts = rel.split("/")
    if len(parts) >= 2 and parts[0] == "src":
        return parts[1]
    return None


def rename_collisions(store, coll, kind):
    """Cross-module name collisions (user rule, 2026-10-31): the
    colliding records only are renamed name -> name__<module>, module =
    the dir of the record's own file. Returns {(name, file, line): new}."""
    renames = {}
    for name, locs in sorted(coll.items()):
        for (file, line) in locs:
            mod = module_of(file)
            if mod is None:
                die("cannot disambiguate %s %r: no module dir in %r"
                    % (kind, name, file))
            renames[(name, file, line)] = "%s__%s" % (name, mod)
    existing = {k[0] for k in store}
    for new in set(renames.values()):
        if new in existing:
            die("disambiguation: renamed %s name %r already exists"
                % (kind, new))
    dup = sorted(n for n, c in Counter(renames.values()).items() if c > 1)
    if dup:
        die("disambiguation: duplicate renamed %s names %s" % (kind, dup))
    return renames


def merge_models(modules, docs, version):
    """Merge the per-module views into the tree-wide model; return
    (doc, report)."""
    warnings = []

    files = set()
    for m in modules:
        for entry in docs[m].get("files") or []:
            files.add(norm_file(m, entry))

    ents, ent_mismatch = {}, []
    n_viewflag_or = 0
    for m in modules:
        for e in docs[m].get("entities") or []:
            key = (e["name"], norm_file(m, e["file"]), e["line"])
            rec = ents.get(key)
            if rec is None:
                ents[key] = dict(e)
                ents[key]["__runs"] = [m]
                continue
            for f in ENTITY_NONLIST:
                if rec[f] != e[f]:
                    ent_mismatch.append(("entity", "%s:%s:%s" % key, f,
                                         repr(rec[f])[:80], repr(e[f])[:80],
                                         m))
            for f in ENTITY_VIEWFLAGS:
                if rec[f] != e[f]:
                    n_viewflag_or += 1
                rec[f] = rec[f] or e[f]
            if rec.get("source") != e.get("source"):
                ent_mismatch.append(("entity-source", "%s:%s:%s" % key,
                                     "source", len(rec.get("source") or ""),
                                     len(e.get("source") or ""), m))
            for f in ENTITY_LISTS:
                rec[f] = sorted(set(rec[f]) | set(e[f]))
            rec["__runs"].append(m)

    rts, rt_mismatch = {}, []
    for m in modules:
        for r in docs[m].get("routines") or []:
            key = (r["name"], norm_file(m, r["file"]), r["line"])
            rec = rts.get(key)
            if rec is None:
                rts[key] = dict(r)
                rts[key]["__runs"] = [m]
                continue
            for f in ROUTINE_NONLIST:
                if rec[f] != r[f]:
                    rt_mismatch.append(("routine", "%s:%s:%s" % key, f,
                                        repr(rec[f])[:80], repr(r[f])[:80],
                                        m))
            if rec.get("source") != r.get("source"):
                rt_mismatch.append(("routine-source", "%s:%s:%s" % key,
                                    "source", len(rec.get("source") or ""),
                                    len(r.get("source") or ""), m))
            for f in ROUTINE_LISTS:
                rec[f] = sorted(set(rec[f]) | set(r[f]))
            rec["__runs"].append(m)

    if ent_mismatch or rt_mismatch:
        for x in (ent_mismatch + rt_mismatch)[:10]:
            print("  identity mismatch: %s" % (x,), file=sys.stderr)
        die("identity fields disagree across module views (%d entity + "
            "%d routine mismatches); the tree changed mid-build, re-run"
            % (len(ent_mismatch), len(rt_mismatch)))

    # external_calls: the per-view definition ("a call not in this
    # view's routines") is made consistent over the merged routine set
    rt_names = {k[0] for k in rts}
    pruned = {}
    for key, rec in rts.items():
        before = set(rec["external_calls"])
        removed = sorted(before & rt_names)
        if removed:
            pruned[key[0]] = removed
        rec["external_calls"] = sorted(before - rt_names)

    # the module field is recomputed from the tree-relative file with
    # irpf90's own formula (the per-run value is mangled from the
    # run-relative filename, which legitimately differs)
    for key, rec in ents.items():
        rec["module"] = mangle(key[1][4:])
    for key, rec in rts.items():
        rec["module"] = mangle(key[1][4:])

    def collisions(store):
        byname = {}
        for (name, file, line) in store:
            byname.setdefault(name, []).append((file, line))
        return {n: sorted(v) for n, v in byname.items() if len(v) > 1}

    ent_coll = collisions(ents)
    rt_coll = collisions(rts)

    ent_renames = rename_collisions(ents, ent_coll, "entity")
    rt_renames = rename_collisions(rts, rt_coll, "routine")
    if ent_renames or rt_renames:
        warnings.append(
            "renamed records keep their pre-rename names in the edge lists"
            " (needs/needed_by/calls/...): a RAG-level id change, not a"
            " Fortran rename")

    out_ents = []
    for key in sorted(ents):
        rec = ents[key]
        rec["seen_in_runs"] = sorted(set(rec.pop("__runs")))
        for f in ENTITY_LISTS:
            rec[f] = sorted(set(rec[f]))
        if key in ent_renames:
            rec["name"] = ent_renames[key]
        out_ents.append(rec)

    out_rts = []
    for key in sorted(rts):
        rec = rts[key]
        rec["seen_in_runs"] = sorted(set(rec.pop("__runs")))
        for f in ROUTINE_LISTS:
            rec[f] = sorted(set(rec[f]))
        if key in rt_renames:
            rec["name"] = rt_renames[key]
        out_rts.append(rec)

    doc = {
        "format": "irpf90-model",
        "format_version": 1,
        "irpf90_version": version,
        "project": PROJECT,
        "files": sorted(files),
        "entities": out_ents,
        "routines": out_rts,
    }

    def renames_report(table):
        return [[k[0], k[1], k[2], v] for k, v in sorted(table.items())]

    report = {
        "out": str(MODEL),
        "irpf90_version": version,
        "modules": len(modules),
        "per_module": {
            m: {"files": len(docs[m].get("files") or []),
                "entities": len(docs[m].get("entities") or []),
                "routines": len(docs[m].get("routines") or [])}
            for m in modules},
        "files": len(doc["files"]),
        "entities": len(out_ents),
        "routines": len(out_rts),
        "entity_name_collisions": ent_coll,
        "routine_name_collisions": rt_coll,
        "entity_renames": renames_report(ent_renames),
        "routine_renames": renames_report(rt_renames),
        "external_calls_pruned_from_in_tree_names": pruned,
        "entity_viewflag_or_merges": n_viewflag_or,
        "warnings": warnings,
    }
    return doc, report


def cmd_build():
    # share/rag is a pure derived artifact: wipe and recreate
    if RAGDIR.exists():
        shutil.rmtree(RAGDIR)
    PER_MODULE.mkdir(parents=True)

    irpf90 = find_irpf90()
    edges, flags = parse_build_ninja(QP_ROOT / "config" / "build.ninja")
    env = parse_env()

    print("irpf90: %s" % irpf90)
    print("flags:  %s" % " ".join(flags))
    print("running %d module dumps..." % len(edges))
    for e in edges:
        run_module(irpf90, e, flags, env)

    modules = [e["module"] for e in edges]
    docs = {}
    for m in modules:
        with open(PER_MODULE / (m + ".json"), encoding="utf-8") as fh:
            docs[m] = json.load(fh)
    versions = {d.get("irpf90_version") for d in docs.values()}
    if len(versions) != 1:
        die("inconsistent irpf90 versions across the per-module dumps: %s"
            % sorted(versions))
    version = versions.pop()

    doc, report = merge_models(modules, docs, version)
    with open(MODEL, "w", encoding="utf-8") as fh:
        json.dump(doc, fh, ensure_ascii=False, indent=1)

    n_docs = index_model(doc)

    with open(REPORT, "w", encoding="utf-8") as fh:
        json.dump(report, fh, ensure_ascii=False, indent=1)

    print("qp2 RAG built: %d files, %d entities, %d routines"
          % (len(doc["files"]), len(doc["entities"]), len(doc["routines"])))
    print("  %d entity name collision(s) -> %d record(s) renamed;"
          % (len(report["entity_name_collisions"]),
             len(report["entity_renames"])))
    print("  %d routine name collision(s) -> %d record(s) renamed"
          % (len(report["routine_name_collisions"]),
             len(report["routine_renames"])))
    print("  indexed %d docs -> %s" % (n_docs, DB))
    print("  corpus: %s" % RAGDIR)
    print("interrogate with:")
    print('  python3 scripts/rag/qp_rag.py query "PT2 selection" -k 5')
    print('  python3 scripts/rag/qp_rag.py query "bitmask" --kind routine')
    print("  python3 scripts/rag/qp_rag.py graph <entity> --depth 3")
    print("  python3 scripts/rag/qp_rag.py stats")


# ~#~#~#~#~#~ #
# i n d e x   #
# ~#~#~#~#~#~ #

def meta_of(rec):
    return {k: v for k, v in rec.items() if k != "source"}


def _edges(lines, pairs):
    for label, vals in pairs:
        if vals:
            lines.append("%s: %s" % (label, " ".join(vals)))


def entity_body(e, by_name):
    if e.get("is_main", True):
        head = "Entity %s (%s" % (e["name"], e.get("type", "?"))
        dims = e.get("dims") or []
        if dims:
            head += " (%s)" % ", ".join(dims)
        head += ")"
    else:
        head = "Entity %s (alias of %s)" % (e["name"], e.get("same_as", "?"))
    out = [head, "%s:%s" % (e.get("file", "?"), e.get("line", 0))]
    if e.get("module"):
        out.append("module: " + e["module"])
    out += list(e.get("doc") or [])
    _edges(out, [("Needs", e.get("needs")), ("Needed by", e.get("needed_by")),
                 ("Calls", e.get("calls"))])
    life = [k for k in ("touched", "freed", "read", "written") if e.get(k)]
    if life:
        out.append("Lifecycle: " + ", ".join(life))
    if e.get("is_main", True):
        src = e.get("source") or ""
        if src:
            out.append(src.rstrip("\n"))
    else:
        rep = by_name.get(e.get("same_as") or "")
        ptr = "alias of " + str(e.get("same_as") or "?")
        if rep:
            ptr += " (defined at %s:%s)" % (rep.get("file", "?"),
                                            rep.get("line", 0))
        out.append(ptr)
    return "\n".join(out)


def routine_body(r):
    out = ["Routine %s (%s)" % (r["name"], r.get("kind", "?")),
           "%s:%s" % (r.get("file", "?"), r.get("line", 0))]
    if r.get("module"):
        out.append("module: " + r["module"])
    out += list(r.get("doc") or [])
    _edges(out, [("Needs", r.get("needs")), ("Calls", r.get("calls")),
                 ("Touches", r.get("touches")),
                 ("Called by", r.get("called_by"))])
    src = r.get("source") or ""
    if src:
        out.append(src.rstrip("\n"))
    return "\n".join(out)


def overview(doc):
    ents = doc.get("entities") or []
    routs = doc.get("routines") or []
    files = doc.get("files") or []
    progs = sorted(r["name"] for r in routs if r.get("kind") == "program")
    mods = {}
    for e in ents:
        mods.setdefault(e.get("module") or "?", []).append(e)
    mod_map = {}
    out = ["Project %s" % PROJECT,
           "%d entities, %d routines, %d files"
           % (len(ents), len(routs), len(files))]
    out += list(files)
    out.append("entry programs: " + (", ".join(progs) if progs else "(none)"))
    out.append("modules:")
    for m in sorted(mods):
        fs = sorted(set(e.get("file", "?") for e in mods[m]))
        mod_map[m] = {"files": fs, "n_entities": len(mods[m])}
        out.append("  %s (%s): %d entities" % (m, ", ".join(fs), len(mods[m])))
    meta = {"files": files, "n_entities": len(ents), "n_routines": len(routs),
            "entry_programs": progs, "modules": mod_map}
    return "\n".join(out), meta


def add_doc(con, kind, name, file_, line, meta, body):
    doc_id = "%s:%s:%s" % (PROJECT, kind, name)
    con.execute("INSERT INTO docs (doc_id, project, kind, name, file, line,"
                " meta) VALUES (?,?,?,?,?,?,?)",
                (doc_id, PROJECT, kind, name, file_ or "", line or 0,
                 json.dumps(meta, ensure_ascii=False)))
    con.execute("INSERT INTO chunks_fts (doc_id, body) VALUES (?,?)",
                (doc_id, body))


def index_model(doc):
    """Drop and rebuild qp2_rag.db from the merged model; return the
    number of indexed docs."""
    if os.path.exists(DB):
        os.remove(DB)
    con = sqlite3.connect(str(DB))
    con.execute("CREATE TABLE docs(doc_id TEXT PRIMARY KEY, project TEXT,"
                " kind TEXT, name TEXT, file TEXT, line INTEGER, meta TEXT)")
    con.execute("CREATE VIRTUAL TABLE chunks_fts"
                " USING fts5(doc_id UNINDEXED, body)")
    ents = doc.get("entities") or []
    routs = doc.get("routines") or []
    by_name = {e["name"]: e for e in ents}
    n = 0
    try:
        for e in ents:
            add_doc(con, "entity", e["name"], e.get("file"), e.get("line"),
                    meta_of(e), entity_body(e, by_name))
            n += 1
        for r in routs:
            add_doc(con, "routine", r["name"], r.get("file"), r.get("line"),
                    meta_of(r), routine_body(r))
            n += 1
        body, meta = overview(doc)
        add_doc(con, "overview", PROJECT, doc.get("project"), 0, meta, body)
        n += 1
    except sqlite3.IntegrityError as e:
        die("cannot index the model: %s (duplicate doc_id?)" % e)
    con.commit()
    con.close()
    return n


# ~#~#~#~#~#~ #
# q u e r y   #
# ~#~#~#~#~#~ #

def open_db():
    if not os.path.exists(DB):
        die("qp2_rag.db not found; run: ninja rag")
    return sqlite3.connect(str(DB))


def load_model():
    if not os.path.exists(MODEL):
        die("qp2.json not found; run: ninja rag")
    with open(MODEL, encoding="utf-8") as fh:
        return json.load(fh)


def fts_match(text):
    toks = [t for t in re.findall(r"[\w'’-]+", text.lower())
            if len(t) > 1 and t not in STOP]
    if not toks:
        die("no searchable terms in %r (need a word of 2+ characters)" % text)
    return " OR ".join('"%s"*' % t for t in dict.fromkeys(toks))


def cmd_query(text, k, kind, as_json):
    con = open_db()
    sql = ("SELECT d.doc_id, d.kind, d.name, d.file, d.line, d.meta, c.body,"
           " bm25(chunks_fts) AS score FROM chunks_fts c JOIN docs d"
           " ON d.doc_id = c.doc_id WHERE chunks_fts MATCH ?")
    args = [fts_match(text)]
    if kind:
        sql += " AND d.kind = ?"
        args.append(kind)
    sql += " ORDER BY score LIMIT ?"
    args.append(k)
    try:
        rows = con.execute(sql, args).fetchall()
    except sqlite3.OperationalError as e:
        die("FTS5 error: %s" % e)
    con.close()
    if as_json:
        out = [{"doc_id": r[0], "kind": r[1], "name": r[2], "file": r[3],
                "line": r[4], "score": r[7], "meta": json.loads(r[5]),
                "body": r[6]} for r in rows]
        json.dump(out, sys.stdout, indent=1, ensure_ascii=False)
        print()
        return
    if not rows:
        print("no results for %r" % text)
        return
    for i, r in enumerate(rows, 1):
        loc = "%s:%s" % (r[3], r[4]) if r[3] else r[2]
        print("%d. %s  (%s)" % (i, r[0], loc))
        print(r[6])
        print()


# ~#~#~#~#~#~ #
# g r a p h   #
# ~#~#~#~#~#~ #

def bfs(ents, start, edge, depth):
    """BFS levels {depth: sorted names}; each entity at its shortest
    depth."""
    seen, levels, frontier = {start}, {}, [start]
    for d in range(1, depth + 1):
        nxt = set()
        for n in frontier:
            for m in ents[n].get(edge) or []:
                if m in ents and m not in seen:
                    nxt.add(m)
        if not nxt:
            break
        levels[d] = sorted(nxt)
        seen |= nxt
        frontier = levels[d]
    return levels


def cmd_graph(entity, depth, direction):
    doc = load_model()
    ents = {e["name"]: e for e in doc.get("entities") or []}
    routs = {r["name"] for r in doc.get("routines") or []}
    if entity not in ents:
        close = difflib.get_close_matches(entity, list(ents), n=5, cutoff=0.3)
        msg = "unknown entity %r" % entity
        if entity in routs:
            msg += " (it is a routine; graph works on entities)"
        msg += "; closest: %s" % (", ".join(close) if close
                                  else ", ".join(sorted(ents)[:10]))
        die(msg)
    for d in (["needs", "needed_by"] if direction == "both" else [direction]):
        what = ("what %s depends on" % entity) if d == "needs" \
            else ("what breaks if %s changes" % entity)
        levels = bfs(ents, entity, d, depth)
        n = sum(len(v) for v in levels.values())
        print("%s: %s (depth <= %d)" % (entity, what, depth))
        if n == 0:
            print("  (none)")
        for d2 in sorted(levels):
            print("  depth %d: %s" % (d2, ", ".join(levels[d2])))
        print()


def cmd_stats():
    con = open_db()
    rows = con.execute(
        "SELECT kind, meta FROM docs WHERE kind IN ('entity','routine')"
        " ORDER BY kind, doc_id").fetchall()
    con.close()
    stats = {"entity": [0, 0], "routine": [0, 0], "ext": set()}
    for kind, meta in rows:
        m = json.loads(meta)
        stats[kind][0] += 1
        if m.get("doc"):
            stats[kind][1] += 1
        stats["ext"].update(m.get("external_calls") or [])
    print("%-14s %-14s %6s %6s %7s"
          % ("project", "kind", "total", "doc", "no-doc"))
    for kind in ("entity", "routine"):
        tot, docd = stats[kind]
        print("%-14s %-14s %6d %6d %7d"
              % (PROJECT, kind, tot, docd, tot - docd))
    ext = sorted(stats["ext"])
    print("%-14s %-14s %s"
          % (PROJECT, "external_calls", ", ".join(ext) if ext else "(none)"))


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("build",
                   help="wipe and rebuild the RAG corpus in share/rag")
    p = sub.add_parser("query",
                       help="BM25 full-text search over the indexed corpus")
    p.add_argument("text")
    p.add_argument("-k", type=int, default=5, help="max results (default 5)")
    p.add_argument("--kind", choices=["entity", "routine"])
    p.add_argument("--json", action="store_true",
                   help="machine-readable output")
    p = sub.add_parser("graph",
                       help="BFS over the entity needs/needed_by edges")
    p.add_argument("entity")
    p.add_argument("--depth", type=int, default=2,
                   help="BFS depth cap (default 2)")
    p.add_argument("--dir", dest="direction", default="both",
                   choices=["needs", "needed_by", "both"])
    sub.add_parser("stats",
                   help="counts, doc coverage, external API surface")
    a = ap.parse_args()
    if a.cmd == "build":
        cmd_build()
    elif a.cmd == "query":
        cmd_query(a.text, a.k, a.kind, a.json)
    elif a.cmd == "graph":
        cmd_graph(a.entity, a.depth, a.direction)
    else:
        cmd_stats()


if __name__ == "__main__":
    main()
