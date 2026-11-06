#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Check that the `qp_edit` config emission leaves the OCaml generator alone.

The script builds isolated `$QP_ROOT` trees from this checkout:

* **A** -- the current working tree (`ei_handler.py` + `qp_edit_config.py`);
* **B** -- the pristine revision given by ``--pristine-ref`` (default
  ``b0eae6e2``), extracted with ``git archive``;
* **C** -- the current working tree again, with the ``src/<module>``
  directories created in reverse order.

A and B receive the same ``src/*/EZFIO.cfg`` modules, created in the same
order.  Every module is run through ``ei_handler.py`` (all targets), then
``ocaml_global`` is run.  The script hashes every generated file under
``ocaml/``, ``data/``, ``external/ezfio/config/`` and every
``src/**/ezfio_interface.irp.f`` and requires that:

* a file present only in A lives under ``data/qp_edit_config/`` (the new
  build artifact);
* no file exists only in B;
* every common file is byte-identical.

C checks that ``data/qp_edit_config/`` is byte-identical to A: unlike the
OCaml ``Input_auto_generated.ml`` order, the Rust config must not depend on
``os.listdir`` order.

Absolute scratch paths embedded in generated headers are replaced by
``<QP_ROOT>`` before hashing, because A, B and C necessarily live in
different directories.  The replacement is applied to the snapshot only;
the config files themselves are additionally checked to contain no absolute
path at all.

Usage::

    TMPDIR=$HOME/OpenCode/tmp python3 scripts/ezfio_interface/test_dual_emission.py \
        --pristine-ref b0eae6e2 --qp-edit rust/target/debug/qp_edit

Exit status is 0 on success and 1 on the first failed check.  On failure the
scratch tree is kept and its path is printed.
"""

import argparse
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile
import tomllib
from pathlib import Path

#: Repository root, inferred from this file (`scripts/ezfio_interface/...`).
REPO_ROOT = Path(__file__).resolve().parents[2]

#: Generated subtrees compared between A and B.
SNAPSHOT_DIRS = ("ocaml", "data", "external/ezfio/config")

#: Prefix (inside `$QP_ROOT`) reserved for the new build artifact.
CONFIG_DIR = "data/qp_edit_config"

#: Handler sections and their order, mirrored from `qp_edit_config.py`.
SECTION_HEAD = ["nuclei_by_hand", "ao_basis", "ao_extra_basis", "electrons"]
SECTION_TAIL = ["mo_basis", "determinants_by_hand"]
HANDLER_MODULES = [
    "nuclei_by_hand",
    "ao_basis",
    "ao_extra_basis",
    "electrons",
    "mo_basis",
    "determinants_by_hand",
    "bitmasks",
]


class CheckError(Exception):
    """One assertion of the test failed."""


def require(condition, message):
    """Raise :class:`CheckError` with *message* when *condition* is false."""
    if not condition:
        raise CheckError(message)


# ------------------------------------------------------------------------- #
# Tree construction                                                         #
# ------------------------------------------------------------------------- #


def module_names():
    """Return the sorted names of the `src/<module>` directories with EZFIO.cfg."""
    source = REPO_ROOT / "src"
    names = sorted(
        path.name
        for path in source.iterdir()
        if path.is_dir() and (path / "EZFIO.cfg").is_file()
    )
    require(names, "no src/*/EZFIO.cfg found in {0}".format(source))
    return names


def copy_working_scripts(root):
    """Copy the working-tree generator scripts and the OCaml type source."""
    ignore = shutil.ignore_patterns("__pycache__", "*.pyc")
    for name in ("ezfio_interface", "utility"):
        shutil.copytree(
            REPO_ROOT / "scripts" / name,
            root / "scripts" / name,
            ignore=ignore,
        )
    (root / "ocaml").mkdir(parents=True, exist_ok=True)
    shutil.copy2(
        REPO_ROOT / "ocaml" / "qptypes_generator.ml",
        root / "ocaml" / "qptypes_generator.ml",
    )


def extract_pristine_scripts(root, rev):
    """Extract the generator scripts of *rev* into *root* with `git archive`."""
    paths = [
        "scripts/ezfio_interface",
        "scripts/utility",
        "ocaml/qptypes_generator.ml",
    ]
    archive = subprocess.run(
        ["git", "archive", rev] + paths,
        cwd=str(REPO_ROOT),
        capture_output=True,
        check=False,
    )
    require(
        archive.returncode == 0,
        "git archive {0} failed:\n{1}".format(
            rev, archive.stderr.decode("utf-8", "replace")
        ),
    )
    extracted = subprocess.run(
        ["tar", "-x", "-C", str(root)],
        input=archive.stdout,
        capture_output=True,
        check=False,
    )
    require(
        extracted.returncode == 0,
        "tar extraction failed:\n{0}".format(
            extracted.stderr.decode("utf-8", "replace")
        ),
    )
    require(
        (root / "scripts" / "ezfio_interface" / "ei_handler.py").is_file(),
        "pristine revision {0} does not contain scripts/ezfio_interface".format(rev),
    )


def make_layout(root):
    """Create the writable directories the generator expects."""
    for relative in ("src", "external/ezfio/config", "data/ezfio_defaults"):
        (root / relative).mkdir(parents=True, exist_ok=True)


def populate_modules(root, names, reverse=False):
    """Create `src/<module>/EZFIO.cfg` for every name, in the given order."""
    order = list(reversed(names)) if reverse else names
    for name in order:
        destination = root / "src" / name
        destination.mkdir(parents=True, exist_ok=True)
        shutil.copy2(
            REPO_ROOT / "src" / name / "EZFIO.cfg",
            destination / "EZFIO.cfg",
        )


def build_tree(root, names, rev=None, reverse=False):
    """Prepare one complete `$QP_ROOT` and populate the module directories."""
    root.mkdir(parents=True, exist_ok=True)
    if rev is None:
        copy_working_scripts(root)
    else:
        extract_pristine_scripts(root, rev)
    make_layout(root)
    populate_modules(root, names, reverse=reverse)


# ------------------------------------------------------------------------- #
# Running the generator                                                     #
# ------------------------------------------------------------------------- #


def generator_env(root):
    """Environment for one tree: `$QP_ROOT`, `$QP_EZFIO`, `$PYTHONPATH`."""
    env = os.environ.copy()
    env["QP_ROOT"] = str(root)
    env["QP_EZFIO"] = str(root / "external" / "ezfio")
    env["PYTHONPATH"] = os.pathsep.join(
        [
            str(root / "scripts" / "ezfio_interface"),
            str(root / "scripts" / "utility"),
        ]
    )
    return env


def run_python(root, args, label, cwd=None, verbose=False):
    """Run the tree's `ei_handler.py` and return its stdout."""
    script = root / "scripts" / "ezfio_interface" / "ei_handler.py"
    result = subprocess.run(
        [sys.executable, str(script)] + args,
        cwd=str(cwd if cwd is not None else root),
        env=generator_env(root),
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise CheckError(
            "generator failed for {0}\ncommand: {1}\nstdout:\n{2}\nstderr:\n{3}".format(
                label,
                " ".join([sys.executable, str(script)] + args),
                result.stdout,
                result.stderr,
            )
        )
    if verbose and result.stdout.strip():
        print("  {0}: {1}".format(label, result.stdout.strip()))
    return result.stdout


def run_generator(root, names, label, verbose=False):
    """Run every module, then `ocaml_global`, inside one tree."""
    for name in names:
        run_python(
            root,
            ["--path_module", str(root / "src" / name)],
            "{0} {1}".format(label, name),
            verbose=verbose,
        )
    run_python(
        root,
        ["ocaml_global"],
        "{0} ocaml_global".format(label),
        cwd=root / "ocaml",
        verbose=verbose,
    )


def check_template_consistency(root, verbose=False):
    """Run `qp_edit_config.check_template_consistency()` in the tree."""
    code = "import qp_edit_config; print(qp_edit_config.check_template_consistency())"
    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=str(root),
        env=generator_env(root),
        capture_output=True,
        text=True,
        check=False,
    )
    require(
        result.returncode == 0,
        "check_template_consistency failed:\nstdout:\n{0}\nstderr:\n{1}".format(
            result.stdout, result.stderr
        ),
    )
    if verbose and result.stdout.strip():
        print("  template consistency: {0}".format(result.stdout.strip()))
    return int(result.stdout.strip().splitlines()[-1])


# ------------------------------------------------------------------------- #
# Snapshots and comparison                                                  #
# ------------------------------------------------------------------------- #


def normalized_digest(path, root):
    """SHA-256 of *path* with the tree's absolute path replaced by <QP_ROOT>."""
    data = path.read_bytes()
    data = data.replace(str(root).encode("utf-8"), b"<QP_ROOT>")
    return hashlib.sha256(data).hexdigest()


def snapshot(root):
    """Map relative path -> digest for every generated file we compare."""
    digests = {}
    for relative in SNAPSHOT_DIRS:
        base = root / relative
        if not base.is_dir():
            continue
        for path in sorted(base.rglob("*")):
            if path.is_file():
                key = path.relative_to(root).as_posix()
                digests[key] = normalized_digest(path, root)
    source = root / "src"
    if source.is_dir():
        for path in sorted(source.rglob("ezfio_interface.irp.f")):
            key = path.relative_to(root).as_posix()
            digests[key] = normalized_digest(path, root)
    return digests


def config_snapshot(root):
    """Map config file name -> digest for `data/qp_edit_config/*.toml`."""
    base = root / CONFIG_DIR
    require(base.is_dir(), "{0} was not created".format(base))
    return {
        path.name: normalized_digest(path, root)
        for path in sorted(base.glob("*.toml"))
    }


def compare_emission(modified, pristine):
    """Assert that the only difference is the new config directory."""
    only_pristine = sorted(set(pristine) - set(modified))
    require(
        not only_pristine,
        "files generated by the pristine generator but not by the modified one:\n  "
        + "\n  ".join(only_pristine),
    )
    only_modified = sorted(set(modified) - set(pristine))
    foreign = [name for name in only_modified if not name.startswith(CONFIG_DIR + "/")]
    require(
        not foreign,
        "modified generator created files outside {0}/:\n  {1}".format(
            CONFIG_DIR, "\n  ".join(foreign)
        ),
    )
    common = sorted(set(modified) & set(pristine))
    changed = [name for name in common if modified[name] != pristine[name]]
    require(
        not changed,
        "files generated identically by both generators differ:\n  "
        + "\n  ".join(changed),
    )
    return common, only_modified


# ------------------------------------------------------------------------- #
# Config validation                                                         #
# ------------------------------------------------------------------------- #


def read_toml(path):
    """Parse one TOML file with the standard library."""
    with open(path, "rb") as stream:
        return tomllib.load(stream)


def validate_config(root):
    """Validate the generated config set; return a dict of measured counts."""
    config_dir = root / CONFIG_DIR
    registry_path = config_dir / "registry.toml"
    require(registry_path.is_file(), "{0} is missing".format(registry_path))
    registry = read_toml(registry_path)
    require(
        registry.get("schema_version") == 1,
        "registry.toml: schema_version is {0!r}".format(
            registry.get("schema_version")
        ),
    )

    entries = registry.get("module", [])
    require(entries, "registry.toml lists no module")
    names = [entry.get("name") for entry in entries]
    files = [entry.get("file") for entry in entries]
    require(
        len(set(names)) == len(names),
        "registry.toml has duplicate module names: {0}".format(names),
    )
    require(
        len(set(files)) == len(files),
        "registry.toml has duplicate file names: {0}".format(files),
    )

    # Order contract: hand-written sections first, then the generated modules
    # sorted by name, then the trailing handlers, then the non-sections.
    require(
        names[: len(SECTION_HEAD)] == SECTION_HEAD,
        "registry head is {0}, expected {1}".format(
            names[: len(SECTION_HEAD)], SECTION_HEAD
        ),
    )
    require(
        names[-len(SECTION_TAIL) - 1 : -1] == SECTION_TAIL,
        "registry tail is {0}, expected {1}".format(
            names[-len(SECTION_TAIL) - 1 :], SECTION_TAIL + ["<non-section>"]
        ),
    )
    require(
        names[-1] == "bitmasks",
        "registry must end with the non-section `bitmasks`, got {0!r}".format(
            names[-1]
        ),
    )
    middle = names[len(SECTION_HEAD) : -len(SECTION_TAIL) - 1]
    require(
        middle == sorted(middle),
        "generated modules are not sorted: {0}".format(middle),
    )
    for handler in HANDLER_MODULES:
        require(handler in names, "handler module {0!r} is missing".format(handler))

    registered = set(files)
    present = {path.name for path in config_dir.glob("*.toml")} - {"registry.toml"}
    require(
        registered == present,
        "registry/files mismatch: only in registry {0}, only on disk {1}".format(
            sorted(registered - present), sorted(present - registered)
        ),
    )

    stats = {
        "modules": len(entries),
        "sections": sum(1 for entry in entries if entry.get("section")),
        "fields": 0,
        "editable": 0,
        "formulas": 0,
        "files": len(present),
    }

    for entry in entries:
        name = entry["name"]
        expected_file = name + ".toml"
        require(
            entry["file"] == expected_file,
            "module {0!r}: file is {1!r}, expected {2!r}".format(
                name, entry["file"], expected_file
            ),
        )
        path = config_dir / entry["file"]
        module = read_toml(path)
        require(
            module.get("schema_version") == 1,
            "{0}: schema_version is {1!r}".format(
                path.name, module.get("schema_version")
            ),
        )
        require(
            module.get("name") == name,
            "{0}: name is {1!r}, registry says {2!r}".format(
                path.name, module.get("name"), name
            ),
        )
        require(
            bool(module.get("section")) == bool(entry.get("section")),
            "{0}: section disagrees with the registry".format(path.name),
        )
        require(
            module.get("handler") == entry.get("handler"),
            "{0}: handler {1!r} disagrees with the registry {2!r}".format(
                path.name, module.get("handler"), entry.get("handler")
            ),
        )
        title = module.get("title", "")
        require(
            bool(title) == bool(module.get("section")),
            "{0}: a title is required for sections and forbidden otherwise".format(
                path.name
            ),
        )

        source_cfg = module.get("source_cfg", "")
        source_path = Path(source_cfg)
        require(
            source_cfg
            and not source_path.is_absolute()
            and ".." not in source_path.parts,
            "{0}: source_cfg {1!r} is not a plain relative path".format(
                path.name, source_cfg
            ),
        )
        ezfio_dirs = module.get("ezfio_dir", [])
        require(
            len(ezfio_dirs) == 1,
            "{0}: expected one EZFIO directory, got {1!r}".format(
                path.name, ezfio_dirs
            ),
        )

        fields = module.get("field", [])
        providers = [field.get("provider") for field in fields]
        require(
            len(set(providers)) == len(providers),
            "{0}: duplicate providers {1}".format(path.name, providers),
        )
        handler = module.get("handler")
        for field in fields:
            provider = field.get("provider")
            require(provider, "{0}: a field has no provider".format(path.name))
            require(
                field.get("doc", "").strip(),
                "{0}/{1}: empty doc".format(path.name, provider),
            )
            require(
                field.get("type", "").strip(),
                "{0}/{1}: empty type".format(path.name, provider),
            )
            require(
                field.get("size", "").strip(),
                "{0}/{1}: empty size".format(path.name, provider),
            )
            interfaces = field.get("interfaces", [])
            require(
                interfaces,
                "{0}/{1}: empty interfaces".format(path.name, provider),
            )
            require(
                set(interfaces) <= {"ezfio", "provider", "ocaml"},
                "{0}/{1}: unknown interface in {2}".format(
                    path.name, provider, interfaces
                ),
            )
            editable = bool(field.get("editable"))
            if editable:
                require(
                    "ocaml" in interfaces,
                    "{0}/{1}: editable without the ocaml interface".format(
                        path.name, provider
                    ),
                )
                require(
                    handler is None,
                    "{0}/{1}: editable inside handler section {2!r}".format(
                        path.name, provider, handler
                    ),
                )
                stats["editable"] += 1
            has_default = "default" in field
            formula = bool(field.get("default_is_formula"))
            if formula:
                require(
                    has_default
                    and str(field["default"]).lstrip().startswith("="),
                    "{0}/{1}: default_is_formula set on a non-formula".format(
                        path.name, provider
                    ),
                )
                stats["formulas"] += 1
            elif has_default:
                require(
                    not str(field["default"]).lstrip().startswith("="),
                    "{0}/{1}: formula default without the flag".format(
                        path.name, provider
                    ),
                )
            stats["fields"] += 1

    # The build artifact must be relocatable: no absolute path at all.
    for path in sorted(config_dir.glob("*.toml")):
        text = path.read_text(encoding="utf-8")
        require(
            str(root) not in text,
            "{0}: contains the absolute build path {1!r}".format(
                path.name, str(root)
            ),
        )

    return stats


# ------------------------------------------------------------------------- #
# Main                                                                      #
# ------------------------------------------------------------------------- #


def parse_args(argv):
    """Parse the command line."""
    parser = argparse.ArgumentParser(
        description="Check that qp_edit config emission does not change the "
        "OCaml generator output.",
    )
    parser.add_argument(
        "--pristine-ref",
        default="b0eae6e2",
        help="git revision of the unmodified generator (default: b0eae6e2)",
    )
    parser.add_argument(
        "--qp-edit",
        type=Path,
        default=None,
        help="optional qp_edit binary; runs --check-config on the result",
    )
    parser.add_argument(
        "--scratch",
        type=Path,
        default=None,
        help="parent directory for the temporary trees (default: $TMPDIR)",
    )
    parser.add_argument(
        "--keep",
        action="store_true",
        help="keep the scratch trees after a successful run",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="print the generator output for every module",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Run the whole check; return a process exit status."""
    args = parse_args(argv)
    names = module_names()

    revision = subprocess.run(
        ["git", "rev-parse", "--verify", "--quiet", args.pristine_ref + "^{commit}"],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
        check=False,
    )
    require(
        revision.returncode == 0,
        "git revision {0!r} not found in {1}".format(
            args.pristine_ref, REPO_ROOT
        ),
    )

    parent = str(args.scratch) if args.scratch is not None else None
    scratch = Path(tempfile.mkdtemp(prefix="qp2-dual-emission-", dir=parent))
    keep = args.keep
    try:
        modified_root = scratch / "modified"
        pristine_root = scratch / "pristine"
        reverse_root = scratch / "reverse"

        build_tree(modified_root, names)
        build_tree(pristine_root, names, rev=args.pristine_ref)
        build_tree(reverse_root, names, reverse=True)

        run_generator(modified_root, names, "modified", verbose=args.verbose)
        run_generator(pristine_root, names, "pristine", verbose=args.verbose)
        run_generator(reverse_root, names, "reverse", verbose=args.verbose)

        titles = check_template_consistency(modified_root, verbose=args.verbose)

        modified = snapshot(modified_root)
        pristine = snapshot(pristine_root)
        common, added = compare_emission(modified, pristine)

        stats = validate_config(modified_root)
        require(
            config_snapshot(modified_root) == config_snapshot(reverse_root),
            "{0}/ depends on src/ directory creation order".format(CONFIG_DIR),
        )

        if args.qp_edit is not None:
            binary = Path(args.qp_edit)
            require(binary.is_file(), "qp_edit binary not found: {0}".format(binary))
            result = subprocess.run(
                [
                    str(binary.resolve()),
                    "--config",
                    str(modified_root / CONFIG_DIR),
                    "--check-config",
                ],
                capture_output=True,
                text=True,
                check=False,
            )
            require(
                result.returncode == 0,
                "qp_edit --check-config failed ({0}):\nstdout:\n{1}\nstderr:\n{2}".format(
                    result.returncode, result.stdout, result.stderr
                ),
            )
            require(
                "configuration OK" in result.stdout,
                "qp_edit --check-config did not report success:\n{0}".format(
                    result.stdout
                ),
            )
            if args.verbose:
                print(result.stdout.rstrip())

        print("dual emission: OK")
        print(
            "  generator runs: {0} modules + ocaml_global in 3 trees "
            "(modified, pristine {1}, reverse order)".format(
                len(names), args.pristine_ref
            )
        )
        print(
            "  compared: {0} identical generated files "
            "(ocaml/, data/, external/ezfio/config/, src/**/ezfio_interface.irp.f)".format(
                len(common)
            )
        )
        print(
            "  added: {0} files, all under {1}/".format(len(added), CONFIG_DIR)
        )
        print(
            "  config: {modules} modules ({sections} sections), "
            "{fields} fields ({editable} editable), "
            "{formulas} formula defaults".format(**stats)
        )
        print("  handler titles checked against qp_edit_template: {0}".format(titles))
        print("  qp_edit config is independent of src/ creation order")
        if args.qp_edit is not None:
            print("  qp_edit --check-config: OK")
        return 0
    except CheckError as error:
        keep = True
        print("dual emission: FAILED", file=sys.stderr)
        print(str(error), file=sys.stderr)
        print("scratch kept at {0}".format(scratch), file=sys.stderr)
        return 1
    finally:
        if not keep:
            shutil.rmtree(scratch, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
