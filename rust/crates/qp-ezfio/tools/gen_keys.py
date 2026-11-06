#!/usr/bin/env python3
"""Generate `src/keys.rs`, the static EZFIO attribute registry.

The source of truth is the set of `src/<module>/EZFIO.cfg` files of the QP2
tree. This script reproduces exactly what
`scripts/ezfio_interface/ei_handler.py` does when it generates the Fortran
`<module>.ezfio_interface_config` files:

* the group is `ezfio_dir` (default: the module directory name, lower-cased);
* the key is `ezfio_name` (default: the section name, lower-cased);
* the type is the fancy `type:` name resolved through the same table
  `ei_handler.get_type_dict()` builds from `ocaml/qptypes_generator.ml`;
* the dimensions are the comma-separated items of `size:` (default `1`, i.e.
  a scalar), with the outer parentheses and blanks removed. Ranges are kept as
  written (`0:pseudo.pseudo_lmax`) and evaluated at run time.

Usage, from `rust/crates/qp-ezfio`:

    python3 tools/gen_keys.py                 # write src/keys.rs
    python3 tools/gen_keys.py --check         # fail if src/keys.rs is stale
    python3 tools/gen_keys.py --cross-check ~/qp2/external/ezfio/config \
        --allow-missing qmcchem

`--cross-check` compares the generated table with reference-produced
`*.ezfio_interface_config` files (names and types only; the reference files
carry evaluated dimension expressions, not the raw ones). `--allow-missing`
names the groups that exist only in the reference tree; in QP2 that is
`qmcchem`, whose `src/qmcchem` is a symlink to an untracked plugin.
"""

import argparse
import configparser
import re
import subprocess
import sys
from pathlib import Path

# `<repo>/rust/crates/qp-ezfio/tools/gen_keys.py`
CRATE_DIR = Path(__file__).resolve().parents[1]
REPO_DIR = CRATE_DIR.parents[2]

# ocaml type (as written in qptypes_generator.ml / ei_handler.py) -> EZFIO type
OCAML_TO_SPEC = {
    "int": ("Integer", 0),
    "int64": ("Integer64", 0),
    "float": ("Real", 0),
    "bool": ("Logical", 0),
    "string": ("Character", 32),
    "sig": ("Character", 32),
}

# Raw Fortran type names accepted in `type:` (ei_handler.py's raw_types).
RAW_TYPES = {
    "integer": ("Integer", 0),
    "int": ("Integer", 0),
    "int64": ("Integer64", 0),
    "integer*8": ("Integer64", 0),
    "float": ("Real", 0),
    "double precision": ("Real", 0),
    "logical": ("Logical", 0),
    "bool": ("Logical", 0),
}

CHARACTER_RE = re.compile(r"^character\*\((\d+)\)$")
CHARACTER_NO_PARENS_RE = re.compile(r"^character\*(\d+)$")


def parse_qptypes(path):
    """Return {fancy_type: (rust_type, char_len)} from qptypes_generator.ml."""
    text = path.read_text()

    # `* Fancy_type : ocaml_type` lines of the `input_data` block.
    aliases = {}
    for line in text.splitlines():
        match = re.match(r"^\*\s+(\w+)\s*:\s*(\w+)\s*$", line.strip())
        if match:
            fancy, ocaml = match.groups()
            if ocaml not in OCAML_TO_SPEC:
                sys.exit(f"{path}: unknown OCaml type {ocaml!r} for {fancy!r}")
            aliases[fancy] = OCAML_TO_SPEC[ocaml]

    # `let untouched = "..."` block: modules whose OCaml type is opaque, mapped
    # to `character*(32)` by ei_handler.py.
    start = text.find('let untouched = "')
    end = text.find("let parse_input", start)
    if start < 0 or end < 0:
        sys.exit(f"{path}: cannot locate the `untouched` block")
    for line in text[start:end].splitlines():
        match = re.match(r"^module\s+(\w+)\s*:", line.strip())
        if match:
            aliases[match.group(1)] = ("Character", 32)

    # Hard-coded in ei_handler.py after the alias scan.
    aliases["MO_class"] = ("Character", 32)
    return aliases


def resolve_type(raw, aliases, where):
    """Map a `type:` value to (rust_type, char_len)."""
    name = raw.strip()
    if name in RAW_TYPES:
        return RAW_TYPES[name]
    match = CHARACTER_RE.match(name) or CHARACTER_NO_PARENS_RE.match(name)
    if match:
        return ("Character", int(match.group(1)))
    if name in aliases:
        return aliases[name]
    sys.exit(
        f"{where}: unknown EZFIO type {name!r}; add it to {Path(__file__).name} "
        f"and to scripts/ezfio_interface/ei_handler.py"
    )


def parse_size(raw, where):
    """Split a `size:` value into raw dimension expressions (None = scalar)."""
    if raw is None:
        return None
    size = raw.strip()
    if size == "1":
        return None
    if not (size.startswith("(") and size.endswith(")")):
        sys.exit(f"{where}: size {raw!r} is not parenthesised")
    inner = size[1:-1]
    if "(" in inner or ")" in inner:
        sys.exit(f"{where}: nested parentheses in size {raw!r}")
    dims = [item.strip() for item in inner.split(",")]
    if not dims or any(item == "" for item in dims):
        sys.exit(f"{where}: empty dimension in size {raw!r}")
    return dims


def read_configs(aliases):
    """Return the list of key specs sorted by (group, key)."""
    specs = []
    seen = {}
    for cfg_path in sorted(REPO_DIR.glob("src/*/EZFIO.cfg")):
        module = cfg_path.parent.name
        parser = configparser.RawConfigParser()
        parser.read_string(cfg_path.read_text())
        for section in parser.sections():
            where = f"{cfg_path.relative_to(REPO_DIR)}:[{section}]"
            group = parser.get(section, "ezfio_dir", fallback=module).strip().lower()
            key = parser.get(section, "ezfio_name", fallback=section).strip().lower()
            ty, char_len = resolve_type(parser.get(section, "type"), aliases, where)
            dims = parse_size(parser.get(section, "size", fallback=None), where)
            if (group, key) in seen:
                sys.exit(f"{where}: duplicate attribute {group}/{key} (also {seen[(group, key)]})")
            seen[(group, key)] = where
            specs.append(
                {"group": group, "key": key, "ty": ty, "char_len": char_len, "dims": dims}
            )
    specs.sort(key=lambda s: (s["group"], s["key"]))
    return specs


# A reference line is `  <key> <fortran type> [<dims>] [=<formula>]`, with the
# type starting at a fixed column but the dimension list, the default formula
# and the type itself possibly containing blanks (`double precision`) or
# parentheses (`character*(32)`). The type is therefore matched with an
# anchored alternation, longest first so that `integer*8` wins over `integer`.
REFERENCE_TYPE_RE = re.compile(
    r"^(character\*\s*\(?\s*\d+\s*\)?|double\s+precision|integer\s*\*\s*8"
    r"|integer|logical|real\s*\*\s*8|real)(.*)$"
)


def parse_reference_config(path):
    """Parse one reference `<module>.ezfio_interface_config`.

    Layout: the group name on line 1, then one line per attribute, with the key
    in the first column and the rest of the line holding the Fortran type, the
    parenthesised dimension list for arrays, and/or a `=formula` default.
    """
    lines = path.read_text().splitlines()
    if not lines:
        sys.exit(f"{path}: empty reference config")
    group = lines[0].strip()
    entries = {}
    for line in lines[1:]:
        if not line.strip():
            continue
        fields = line.split(None, 1)
        key = fields[0]
        rest = fields[1] if len(fields) > 1 else ""
        match = REFERENCE_TYPE_RE.match(rest.lstrip())
        if not match:
            sys.exit(f"{path}: cannot parse reference line {line!r}")
        fortran_type, tail = match.group(1), match.group(2).strip()
        if fortran_type.startswith("character"):
            ty, char_len = "Character", int(re.search(r"\d+", fortran_type).group())
        else:
            # Normalise the blanks around `*` so the RAW_TYPES lookup works.
            normalised = re.sub(r"\s*([*])\s*", r"\1", fortran_type)
            if normalised not in RAW_TYPES:
                sys.exit(f"{path}: unknown reference type {fortran_type!r}")
            ty, char_len = RAW_TYPES[normalised]
        # `tail` is either empty (scalar), `(dims)` (array), or `=formula`
        # (scalar with a computed default).
        if tail.startswith("("):
            is_array = True
        elif tail == "" or tail.startswith("="):
            is_array = False
        else:
            sys.exit(f"{path}: unexpected trailing text {tail!r} in {line!r}")
        entries[(group, key)] = (ty, char_len, is_array)
    return entries


def cross_check(specs, reference_dir, allowed_missing):
    """Compare the generated table with reference interface configs.

    `allowed_missing` is the set of groups that exist in the reference tree but
    are not tracked in this repository (currently `qmcchem`, a plugin symlink in
    the main tree). Their keys are reported but do not count as differences.
    """
    reference = {}
    for path in sorted(reference_dir.glob("*.ezfio_interface_config")):
        reference.update(parse_reference_config(path))

    ours = {
        (s["group"], s["key"]): (s["ty"], s["char_len"], s["dims"] is not None) for s in specs
    }

    known_groups = {name[0] for name in set(reference) - set(ours)}
    for group in sorted(allowed_missing - known_groups):
        sys.exit(f"--allow-missing {group}: no reference-only key with that group")

    problems = 0
    for name in sorted(set(ours) - set(reference)):
        print(f"only in the generated table: {name[0]}/{name[1]} = {ours[name]}")
        problems += 1
    for name in sorted(set(reference) - set(ours)):
        if name[0] in allowed_missing:
            print(f"known reference-only key: {name[0]}/{name[1]} = {reference[name]}")
        else:
            print(f"only in the reference configs: {name[0]}/{name[1]} = {reference[name]}")
            problems += 1
    for name in sorted(set(ours) & set(reference)):
        if ours[name][:2] != reference[name][:2] or ours[name][2] != reference[name][2]:
            print(f"mismatch {name[0]}/{name[1]}: generated {ours[name]}, reference {reference[name]}")
            problems += 1
    return problems


def render(specs, commit):
    """Render src/keys.rs."""
    lines = [
        "// This file is generated by `tools/gen_keys.py`; DO NOT EDIT.",
        f"// Source: the {len({s['group'] for s in specs})} tracked `src/*/EZFIO.cfg` files of commit {commit},",
        f"// {len(specs)} attributes, sorted by (group, key) for binary search.",
        "//",
        "// Re-generate: `python3 tools/gen_keys.py` from the crate directory.",
        "",
        "use crate::layout::ElementType::{Character, Integer, Integer64, Logical, Real};",
        "use crate::registry::KeySpec;",
        "",
        "/// Every EZFIO attribute declared by the QP2 modules, sorted by",
        "/// `(group, key)`.",
        "///",
        "/// Fields of each row are `(group, key, type, character length, raw",
        "/// dimension expressions)`; see [`KeySpec`].",
        "#[rustfmt::skip]",
        "pub static KEYS: &[KeySpec] = &[",
    ]
    for spec in specs:
        ty = spec["ty"]
        if spec["dims"] is None:
            dims = "None"
        else:
            quoted = ", ".join(f'"{item}"' for item in spec["dims"])
            dims = f"Some(&[{quoted}])"
        lines.append(
            f'    KeySpec::new("{spec["group"]}", "{spec["key"]}", {ty}, {spec["char_len"]}, {dims}),'
        )
    lines.append("];")
    lines.append("")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--check",
        action="store_true",
        help="do not write; exit 1 if src/keys.rs is out of date",
    )
    parser.add_argument(
        "--cross-check",
        metavar="DIR",
        help="compare with reference <module>.ezfio_interface_config files",
    )
    parser.add_argument(
        "--allow-missing",
        action="append",
        default=[],
        metavar="GROUP",
        help="reference-only group to accept (untracked module); repeatable",
    )
    args = parser.parse_args()

    aliases = parse_qptypes(REPO_DIR / "ocaml" / "qptypes_generator.ml")
    specs = read_configs(aliases)
    if args.cross_check:
        problems = cross_check(specs, Path(args.cross_check), set(args.allow_missing))
        if problems:
            print(f"{problems} difference(s) against {args.cross_check}", file=sys.stderr)
            return 1
        print(f"cross-check OK against {args.cross_check}")

    commit = subprocess.run(
        ["git", "-C", str(REPO_DIR), "rev-parse", "--short", "HEAD"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    rendered = render(specs, commit)

    out = CRATE_DIR / "src" / "keys.rs"
    if args.check:
        if not out.exists() or out.read_text() != rendered:
            print(f"{out} is out of date; run tools/gen_keys.py", file=sys.stderr)
            return 1
        print(f"{out} is up to date ({len(specs)} attributes)")
        return 0
    out.write_text(rendered)
    print(f"wrote {out} ({len(specs)} attributes from commit {commit})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
