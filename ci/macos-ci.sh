 #!/usr/bin/env bash
 #
 # ci/macos-ci.sh — build QP2 on a macOS arm64 self-hosted runner (build only).
 #
 # Builds every C/Fortran target (`ninja all`). The OCaml tools target is NOT
 # built: the in-tree opampack ships OCaml 4.12 (2023) while the OCaml code
 # needs >= 5.0, so a fresh `./configure -i all` OCaml build fails on EVERY
 # platform (upstream gap, independent of macOS). `ocaml` is installed anyway
 # (via -i all) because configure's validation block requires `find_exe ocamlc`.
 #
 # How `./configure` behaves (matters for this script):
 #   - `-i` takes a SINGLE argument; `-i all` is the canonical package list
 #     (cmake zlib ninja zeromq f77zmq gmp ocaml docopt resultsFile bats
 #     trexio armadillo openorb simple_gpu). qmckl is NOT in `all`.
 #   - A failing package step is NON-fatal: the executed command list ends
 #     in `set +x`, which masks the failure; the final validation block
 #     (find_lib / find_exe / find_python_lib) is the real gate.
 #   - configure line 16 does `export CC=gcc` unconditionally, overriding any
 #     inherited CC.
 #
 # macOS-specific fixes applied here:
 #   1. qmckl: its autotools OpenMP check fails with Apple clang (the `gcc`
 #      that bare-CC resolves to on macOS). Run `./configure -i qmckl` with
 #      a PATH shim mapping gcc/g++ to Homebrew gcc-16/g++-16.
 #   2. the `cmake` package downloads a nonexistent Darwin tarball (upstream
 #      404, non-fatal): find_exe picks up Homebrew's cmake instead.
 #
 # simple_gpu: installed by `-i all` from the v1.0.3 release, which carries
 # the Apple Accelerate sdot fix (scemama/simple_gpu#16). Its own
 # `make check` (8/8, incl. test_sdot_cpu against the system Accelerate)
 # passes on this platform, so the step completes through `make install`.
 #
 # Prerequisites on the runner (Homebrew): gfortran, gcc-16/g++-16, cmake,
 # wget, autoconf toolchain.

set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"

export PATH=/opt/homebrew/bin:$PATH
export QP_ROOT="$ROOT"

echo "==> [1/4] qmckl with Homebrew gcc-16 (validation fails yet: deps not all installed)"
SHIM=$(mktemp -d)
ln -sf /opt/homebrew/bin/gcc-16 "$SHIM/gcc"
ln -sf /opt/homebrew/bin/g++-16 "$SHIM/g++"
PATH="$SHIM:$PATH" ./configure -i qmckl || echo "expected: remaining dependencies not installed yet"
rm -rf "$SHIM"

echo "==> [2/4] all in-tree dependencies (incl. simple_gpu v1.0.3; validation passes at the end)"
./configure -i all

echo "==> [3/4] final configuration"
./configure -c config/gfortran_macos.cfg

echo "==> [4/4] C/Fortran build (ninja all)"
bash -c 'source quantum_package.rc && ninja all'

echo "==> macOS build OK"
