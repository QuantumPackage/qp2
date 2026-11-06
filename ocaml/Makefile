# Build wrapper for the Rust tool layer.
#
# The OCaml sources in this directory are kept as the reference
# implementation; they are no longer compiled.
# This Makefile keeps the entry points the build already calls
#
#     make -C $QP_ROOT/ocaml executables .gitignore qp_edit.native qp_run.native
#
# and publishes the Rust binaries under their historical names, because
# `etc/paths.rc` puts this directory on PATH and `bin/qp_export_as_tgz`
# collects `*.native`.  The binaries still live in the Cargo target
# directory; the entries here are symlinks, so `cargo build` alone is enough
# to refresh them.
#
# It also regenerates `data/qp_edit_config/registry.toml`, the ordered module
# list the editor loads: that was a by-product of the old
# `Input_auto_generated.ml` rule, so it has to be requested explicitly here.

.NOPARALLEL:

# Check if QP_ROOT is defined
ifndef QP_ROOT
$(info  -------------------- Error --------------------)
$(info  QP_ROOT undefined. Source the quantum_package.rc script)
$(info  -----------------------------------------------)
$(error )
endif

CARGO ?= cargo

# Pinned with --target-dir so the symlinks cannot point at a stale tree when
# the caller exports CARGO_TARGET_DIR.
TARGET_DIR := $(QP_ROOT)/rust/target
RELEASE_DIR := $(TARGET_DIR)/release

# Published names: the seven Rust binaries, plus the `.native`/`.byte`
# spellings that `bin/qp_export_as_tgz` and the ninja rule still ask for.
BINARIES := element_create_db qp_create_ezfio qp_edit qp_print_basis qp_run \
            qp_set_mo_class qp_tunnel
LINKS := $(BINARIES) qp_edit.native qp_run.native element_create_db.byte

# The ordered list of modules read by `qp_edit` (`data/qp_edit_config/`
# is a build artifact, like the per-module tables).  This replaces the
# `Input_auto_generated.ml`/`qp_edit.ml` edge of the OCaml build: the
# same generator call still emits those two obsolete files, which the
# `.gitignore` rule below keeps out of `git status`.
REGISTRY := $(QP_ROOT)/data/qp_edit_config/registry.toml

.PHONY: default rust executables clean .gitignore

# The default goal is what `ninja` runs (the `rust_target` rule calls `make`
# without a goal).  It keeps the historical side effect of that call: the
# `data/executables` table, which `qp_run` uses to resolve PROGRAM names.
default: executables

# The build runs from `rust/` so that rustup picks up `rust/rust-toolchain.toml`
# (toolchain discovery starts at the current directory, not at the manifest).
rust:
	cd $(QP_ROOT)/rust && $(CARGO) build --release --target-dir $(TARGET_DIR)

# One rule for every published name.  The link target is the binary named by
# the target once a compatibility suffix is stripped, so `qp_edit.native` and
# `qp_edit` both point at `target/release/qp_edit`.
$(LINKS): rust
	ln -sfn $(RELEASE_DIR)/$(patsubst %.native,%,$(patsubst %.byte,%,$@)) $@

# The module list is recomputed from `src/*/EZFIO.cfg`, so the registry does
# not depend on which modules have already been built.  The module-local
# `ninja` default reaches this goal too, hence the dependency on `executables`
# below rather than on `rust`.
$(REGISTRY): $(QP_ROOT)/scripts/ezfio_interface/ei_handler.py \
             $(QP_ROOT)/scripts/ezfio_interface/qp_edit_config.py \
             $(wildcard $(QP_ROOT)/src/*/EZFIO.cfg)
	cd $(QP_ROOT) && ei_handler.py ocaml_global

# `data/executables` maps every Fortran binary to its path; the Rust tools are
# reached through PATH, not through this file.  Regenerated on every call
# because ninja asks for this phony target on every build.  It depends on
# every link so that a module-local `ninja` (whose default is `executables`)
# still publishes the whole tool set on PATH.
executables: $(LINKS) $(REGISTRY)
	$(QP_ROOT)/scripts/module/create_executables_list.sh
	$(QP_ROOT)/ocaml/element_create_db.byte

# The build asks for this goal by name.  The generated file ignores itself
# and the published symlinks (the `*.byte`/`*.native` globs also cover the
# stale artifacts of an old OCaml build), plus the two files that
# `ei_handler.py ocaml_global` still writes.
.gitignore:
	@for i in .gitignore *.byte *.native _build $(LINKS) element_create_db \
	           Input_*.ml qp_edit.ml Input_auto_generated.ml ; \
	do \
	   echo $$i ; \
	done | sort | uniq > .gitignore

clean:
	rm -rf _build $(LINKS)
