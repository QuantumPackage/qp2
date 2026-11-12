#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Dual emission of the `qp_edit` configuration files (design decision D3).

The Rust `qp_edit` is a *fixed* binary: it loads a declarative description of
every editable EZFIO section from ``$QP_ROOT/data/qp_edit_config/``.  Those
files are **build artifacts** generated from the same ``EZFIO.cfg`` inputs that
already drive the Fortran generator, so editing an ``EZFIO.cfg``
regenerates the configuration without recompiling anything.

Layout produced here (see ``rust/DESIGN-qp_edit-config.md``)::

    data/qp_edit_config/registry.toml     # ordered list of modules
    data/qp_edit_config/<module>.toml     # one declarative table per module

``ei_handler.py`` used to emit the OCaml files byte-for-byte alongside these;
this module only *adds* files.  The OCaml emission and the dual-emission test
that checked the OCaml side did not move were removed in w11.

Config-vs-code split
--------------------

* Data (this module): field names, EZFIO keys, declared types, defaults, sizes,
  help texts, interface flags, section titles and order, and the name of the
  handler that owns a section.
* Code (Rust, ``crates/qp-input``): the semantics of each declared type
  (parsing, validation, printing) and the seven hand-written handlers ported
  from the tracked ``ocaml/Input_*.ml`` files.  The handler names emitted here
  must match ``crates/qp-input/src/handlers.rs``.

The handler table below is the only duplicated knowledge: the section titles
must equal the titles compiled into the Rust ``qp_edit`` binary (the
``TASKS``/registry).  The cross-check that used to parse the
``qp_edit_template`` file to assert that the two agree was removed in w11
with the template.
"""

import os

from qp_path import QP_ROOT

#: Version of the on-disk configuration format understood by ``qp-input``.
SCHEMA_VERSION = 1

#: Path of the generator, recorded in the file headers (never absolute: the
#: generated files must not depend on the build directory).
GENERATOR = "scripts/ezfio_interface/ei_handler.py"

#: Generated files are written here, relative to ``$QP_ROOT``.
CONFIG_DIR = "data/qp_edit_config"

#: Module directories (lowercase, as found under ``src/``) that are owned by a
#: hand-written handler instead of the generic engine.  Each entry is a list of
#: ``(module_id, title, handler, is_section)``:
#:
#: * ``module_id``  -- name of the ``<module_id>.toml`` file and of the section
#:   in the registry;
#: * ``title``      -- the RST header printed by ``qp_edit``; must equal the
#:   title compiled into the Rust ``qp_edit`` binary (``TASKS``/registry); the
#:   cross-check that used to parse ``qp_edit_template`` was removed in w11
#:   with the template; empty when the module is not a section;
#: * ``handler``    -- handler name compiled into the Rust binary;
#: * ``is_section`` -- whether ``qp_edit`` visits it as an editable section.
#:
#: One module directory can own several handlers: ``src/ao_basis/EZFIO.cfg``
#: declares both the ``ao_basis_*`` and the ``ao_extra_basis_*`` keys, so it
#: yields ``ao_basis.toml`` and ``ao_extra_basis.toml``.
HANDLER_MODULES = {
    "nuclei": [("nuclei_by_hand", "Molecule", "nuclei_by_hand", True)],
    "ao_basis": [
        ("ao_basis", "AO basis", "ao_basis", True),
        ("ao_extra_basis", "AO extra_basis", "ao_extra_basis", True),
    ],
    "electrons": [("electrons", "Electrons", "electrons", True)],
    "mo_basis": [("mo_basis", "MO basis", "mo_basis", True)],
    "determinants": [
        ("determinants_by_hand", "Determinants_by_hand", "determinants_by_hand", True)
    ],
}

#: Sections visited by ``qp_edit`` before the generated modules, in the order
#: of the ``tasks`` list of the removed ``qp_edit_template`` (now the ``TASKS``
#: list of the Rust binary).
SECTION_HEAD = ["nuclei_by_hand", "ao_basis", "ao_extra_basis", "electrons"]

#: Sections visited after the generated modules (same source).
SECTION_TAIL = ["mo_basis", "determinants_by_hand"]


# ------------------------------------------------------------------------- #
# Minimal TOML writer                                                       #
# ------------------------------------------------------------------------- #
# Only the constructs below are needed: bare keys, basic (single-line)
# strings, booleans, arrays of strings, tables and arrays of tables.  The
# `toml` Python package is not a build dependency, and the generated files
# must stay readable, so a hand-rolled writer is the smaller change.


def toml_string(value):
    """Render *value* as a single-line TOML basic string.

    Newlines and other control characters are escaped (TOML basic strings
    allow ``\\n``), so the result is always one line: the generated files stay
    diffable and greppable.
    """
    out = []
    for ch in str(value):
        if ch == '"':
            out.append('\\"')
        elif ch == "\\":
            out.append("\\\\")
        elif ch == "\n":
            out.append("\\n")
        elif ch == "\t":
            out.append("\\t")
        elif ch == "\r":
            out.append("\\r")
        elif ord(ch) < 0x20 or ord(ch) == 0x7F:
            out.append("\\u%04X" % ord(ch))
        else:
            out.append(ch)
    return '"' + "".join(out) + '"'


def toml_bool(value):
    """Render a Python bool as ``true`` or ``false``."""
    return "true" if value else "false"


def toml_string_array(values):
    """Render a list of strings as a single-line TOML array."""
    return "[" + ", ".join(toml_string(v) for v in values) + "]"


def _assign(items, key, rendered):
    """Append ``key = rendered`` to *items* unless *rendered* is ``None``."""
    if rendered is not None:
        items.append("{0} = {1}".format(key, rendered))


# ------------------------------------------------------------------------- #
# Field tables                                                              #
# ------------------------------------------------------------------------- #


def is_formula(default):
    """Tell whether a default is a Fortran expression rather than a literal.

    ``EZFIO.cfg`` uses ``default: = electrons.elec_alpha_num + ...`` for
    values that the Fortran side computes from other EZFIO entries.  Such a
    default cannot be parsed by the editor, so it is kept verbatim and marked
    as an expression: the Rust loader refuses to use it and the handler (or
    the Fortran side) must supply the value.
    """
    return default.startswith("=")


def field_text(provider, info):
    """Render one ``[[field]]`` table from a ``get_dict_config_file`` entry.

    *provider* is the section name of the ``EZFIO.cfg`` file; *info* is the
    dictionary built by :func:`ei_handler.get_dict_config_file`.  Fields are
    emitted in the order of the ``EZFIO.cfg`` file, which is the order of the
    generated OCaml record (the read/write snippets are sorted, but the record
    layout follows the file).
    """
    interfaces = sorted(info["interface"])
    items = []
    _assign(items, "provider", toml_string(provider))
    _assign(items, "ezfio_name", toml_string(info["ezfio_name"]))
    _assign(items, "type", toml_string(info["type_raw"]))
    _assign(items, "interfaces", toml_string_array(interfaces))
    _assign(items, "doc", toml_string(info["doc"]))
    _assign(items, "size", toml_string(info["size"]))
    if "default" in info:
        default = info["default"].ocaml
        _assign(items, "default", toml_string(default))
        if is_formula(default):
            _assign(items, "default_is_formula", toml_bool(True))
    # `editable` is filled in by `module_config_text`: it depends on whether
    # the module is owned by a handler.
    return items


def _module_fields(dict_ezfio_cfg, editable_from_interface):
    """Return the ``[[field]]`` blocks of a module."""
    blocks = []
    for provider, info in dict_ezfio_cfg.items():
        items = field_text(provider, info)
        editable = editable_from_interface and "ocaml" in info["interface"]
        _assign(items, "editable", toml_bool(editable))
        blocks.append("[[field]]\n" + "\n".join(items) + "\n")
    return blocks


# ------------------------------------------------------------------------- #
# Module and registry files                                                 #
# ------------------------------------------------------------------------- #


def module_config_text(
    module_id,
    title,
    is_section,
    handler,
    source_cfg,
    dict_ezfio_cfg,
    editable_from_interface,
):
    """Render one ``<module_id>.toml`` file.

    Arguments:

    * ``module_id`` -- file stem and registry key;
    * ``title`` -- RST header of the section (empty when not a section);
    * ``is_section`` -- whether ``qp_edit`` visits this module;
    * ``handler`` -- handler name, or ``None`` for the generic engine;
    * ``source_cfg`` -- path of the ``EZFIO.cfg`` relative to ``$QP_ROOT``;
    * ``dict_ezfio_cfg`` -- parsed configuration (see
      :func:`ei_handler.get_dict_config_file`);
    * ``editable_from_interface`` -- when true, a field is editable iff its
      ``interface`` list contains ``ocaml``; when false (handler modules) every
      field is non-editable because the handler owns the whole section.
    """
    ezfio_dirs = sorted({info["ezfio_dir"] for info in dict_ezfio_cfg.values()})

    lines = [
        "# Generated by {0} -- do not edit.".format(GENERATOR),
        "# Declarative description of the `qp_edit` section `{0}`.".format(module_id),
        "",
    ]
    items = []
    _assign(items, "schema_version", str(SCHEMA_VERSION))
    _assign(items, "name", toml_string(module_id))
    _assign(items, "title", toml_string(title))
    _assign(items, "section", toml_bool(is_section))
    _assign(items, "handler", toml_string(handler) if handler else None)
    _assign(items, "source_cfg", toml_string(source_cfg))
    # `create_ezfio_stuff` refuses several EZFIO directories in one config, so
    # a module has exactly one; the list keeps the schema honest should that
    # ever change.
    _assign(items, "ezfio_dir", toml_string_array(ezfio_dirs))
    lines.extend(items)
    lines.append("")

    for block in _module_fields(dict_ezfio_cfg, editable_from_interface):
        lines.append(block)

    return "\n".join(lines)


def registry_text(l_module_lower):
    """Render ``registry.toml`` from the list of generated modules.

    *l_module_lower* is the list returned by
    :func:`ei_handler.get_l_registry_modules`.

    Order: the handler sections that ``qp_edit`` visits first, then the
    generated modules **sorted by name**, then the handler sections visited
    last.  The OCaml generator emits
    the generated modules in ``os.listdir`` order, which is not reproducible
    across machines; sorting here makes the build artifact deterministic.
    Section order is not semantically significant -- sections are located by
    their header -- but it is visible in ``qp_edit -w`` output, see
    ``rust/DESIGN-qp_edit-config.md``.
    """
    generated = sorted(set(l_module_lower))
    by_id = {}
    for entries in HANDLER_MODULES.values():
        for module_id, title, handler, is_section in entries:
            by_id[module_id] = (title, handler, is_section)

    order = []
    order.extend((m, True) for m in SECTION_HEAD)
    order.extend((m, True) for m in generated)
    order.extend((m, True) for m in SECTION_TAIL)
    order.extend(
        (module_id, is_section)
        for module_id, (_, _, is_section) in sorted(by_id.items())
        if not is_section
    )

    lines = [
        "# Generated by {0} -- do not edit.".format(GENERATOR),
        "# Ordered list of the `qp_edit` configuration modules.",
        "",
    ]
    _assign(lines, "schema_version", str(SCHEMA_VERSION))
    lines.append("")

    seen = set()
    for module_id, is_section in order:
        if module_id in seen:
            # A handler module that is also a generated module (nuclei,
            # determinants) appears twice: once as the handler-owned
            # section coerced above, once from `generated`.  They have
            # different ids (`nuclei_by_hand` vs `nuclei`), so a duplicate
            # here is a bug in this function, not a legitimate case.
            raise ValueError("duplicate module in registry: {0}".format(module_id))
        seen.add(module_id)
        handler = by_id[module_id][1] if module_id in by_id else None
        items = [
            ("name", toml_string(module_id)),
            ("file", toml_string(module_id + ".toml")),
            ("section", toml_bool(is_section)),
        ]
        if handler:
            items.append(("handler", toml_string(handler)))
        lines.append("[[module]]")
        lines.extend("{0} = {1}".format(k, v) for k, v in items)
        lines.append("")

    return "\n".join(lines)


# ------------------------------------------------------------------------- #
# Emission (called from ei_handler.py)                                      #
# ------------------------------------------------------------------------- #


def _write(relative_path, text):
    """Write *text* to ``$QP_ROOT/<relative_path>``, creating directories."""
    path = os.path.join(QP_ROOT, relative_path)
    directory = os.path.dirname(path)
    # One `ei_handler.py` runs per module as a parallel ninja edge; all of
    # them create this directory, so creation must be race-free.
    os.makedirs(directory, exist_ok=True)
    with open(path, "w+") as f:
        f.write(text)
    return path


def save_module_config(module_id, text):
    """Write ``data/qp_edit_config/<module_id>.toml``; return its path."""
    return _write(os.path.join(CONFIG_DIR, module_id + ".toml"), text)


def save_registry(text):
    """Write ``data/qp_edit_config/registry.toml``; return its path."""
    return _write(os.path.join(CONFIG_DIR, "registry.toml"), text)


def emit_module(module_path, module_lower, dict_ezfio_cfg):
    """Emit every ``qp_edit`` configuration owned by one module directory.

    *module_path* is the path of the module's ``EZFIO.cfg`` (used to record
    the source file relative to ``$QP_ROOT``), *module_lower* its directory
    name, *dict_ezfio_cfg* the parsed configuration.

    Returns the list of written paths (empty when the module is neither a
    handler module nor an ``ocaml``-interface module).
    """
    source_cfg = os.path.relpath(module_path, QP_ROOT)
    written = []

    has_ocaml_field = any(
        "ocaml" in info["interface"] for info in dict_ezfio_cfg.values()
    )

    for module_id, title, handler, is_section in HANDLER_MODULES.get(
        module_lower, []
    ):
        text = module_config_text(
            module_id=module_id,
            title=title,
            is_section=is_section,
            handler=handler,
            source_cfg=source_cfg,
            dict_ezfio_cfg=dict_ezfio_cfg,
            editable_from_interface=False,
        )
        written.append(save_module_config(module_id, text))

    if has_ocaml_field:
        module_id = module_lower
        # One module id names exactly one file.  A directory that both owns a
        # handler section and exposes an `ocaml` field would make the generic
        # module shadow the handler file (written just above).  No directory
        # is in that situation today; fail loudly if one appears rather than
        # drop a handler from the build.
        handler_ids = {
            entry[0] for entries in HANDLER_MODULES.values() for entry in entries
        }
        if module_id in handler_ids:
            raise ValueError(
                "{0} is both a handler module and an ocaml-interface module; "
                "refusing to overwrite {0}.toml".format(module_id)
            )
        text = module_config_text(
            module_id=module_id,
            title=module_lower.capitalize(),
            is_section=True,
            handler=None,
            source_cfg=source_cfg,
            dict_ezfio_cfg=dict_ezfio_cfg,
            editable_from_interface=True,
        )
        written.append(save_module_config(module_id, text))

    return written
