#!/bin/bash
#
# Build the Rust tool layer and the artifacts the Fortran build needs.
# Called by the ninja rules `executables` (root build.ninja) and `make_rust`
# (module-local `rust_target`); QP_ROOT must be set (source
# quantum_package.rc first).
#
# Artifacts:
#   rust/target/release/*                the six tool binaries (on PATH
#                                        through etc/paths.rc)
#   data/qp_edit_config/registry.toml    the ordered module list read by
#                                        qp_edit
#   data/executables                     the Fortran binary table read by
#                                        qp_run
#   data/list_element.txt                the element table (element_create_db)

set -e

if [[ -z ${QP_ROOT} ]] ; then
  echo "The QP_ROOT environment variable is not set."
  echo "Please reload the quantum_package.rc file."
  exit 1
fi

# The build runs from `rust/` so that rustup picks up `rust/rust-toolchain.toml`
# (toolchain discovery starts at the current directory, not at the manifest).
# The pinned `--target-dir` keeps the published path stable whatever
# CARGO_TARGET_DIR the caller exports.
cd ${QP_ROOT}/rust
cargo build --release --target-dir ${QP_ROOT}/rust/target
cd ${QP_ROOT}
scripts/ezfio_interface/ei_handler.py generate_registry
scripts/module/create_executables_list.sh
${QP_ROOT}/rust/target/release/element_create_db
