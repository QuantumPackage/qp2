#!/usr/bin/env bash
#
# ci/macos-ci.sh — build QP2 on a GitHub-hosted macOS arm64 runner (build only).
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
#   1. Global gcc/g++ shim (below): on this platform a bare `gcc` resolves to
#      Apple's clang shim (/usr/bin/gcc), which breaks in-tree builds. zlib
#      1.2.11 in particular fires its `fdopen` hack because clang (not gcc)
#      predefines TARGET_OS_MAC. Homebrew installs only versioned gcc-NN/
#      g++-NN, and the glob below must exclude gcc-ar-NN/gcc-nm-NN/
#      gcc-ranlib-NN: with a plain `gcc-*` glob, sort -V picks gcc-ranlib-NN,
#      which "compiles" conftests as a silent no-op (autotools then report
#      "C compiler cannot create executables"). The shim must stay global:
#      it is needed by the zlib/zeromq/gmp/trexio/qmckl builds AND by the
#      gate's find_lib link tests, which invoke bare `gcc`.
#   2. qmckl is built AFTER `-i all`, not before: it links against trexio
#      (AC_CHECK_LIB trexio_open), which `-i all` installs.
#   3. The `cmake` package downloads a nonexistent Darwin tarball (upstream
#      404, non-fatal): find_exe picks up Homebrew's cmake instead.
#
# simple_gpu: installed by `-i all` from the v1.0.3 release, which carries
# the Apple Accelerate sdot fix (scemama/simple_gpu#16). Its own
# `make check` (8/8, incl. test_sdot_cpu against the system Accelerate)
# passes on this platform, so the step completes through `make install`.
#
# Prerequisites: gfortran, a versioned Homebrew gcc (gcc-NN/g++-NN), hdf5 and
# pkg-config, installed by the workflow; cmake, wget and the autoconf
# toolchain ship with the GitHub runner image.

set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"

export PATH=/opt/homebrew/bin:$PATH
export QP_ROOT="$ROOT"

# Global gcc/g++ shim (see header, item 1).
GCCC=$(ls -d /opt/homebrew/bin/gcc-[0-9]* | sort -V | tail -1)
GPPC=$(ls -d /opt/homebrew/bin/g++-[0-9]* | sort -V | tail -1)
SHIM=$(mktemp -d)
ln -sf "$GCCC" "$SHIM/gcc"
ln -sf "$GPPC" "$SHIM/g++"
export PATH="$SHIM:$PATH"

# configure masks per-package build failures (see header); only the trailing
# validation gate can make `./configure -i X` exit nonzero. That gate needs
# qmckl, which is built in [2/4] (after trexio), so the gate of [1/4] is
# expected to fail: tolerate it. The real gate is the one at the end of [2/4],
# where every dependency exists.

echo "==> [1/4] all in-tree dependencies (gate fails yet: qmckl not built)"
./configure -i all || :

echo "==> [2/4] qmckl (needs trexio from [1/4]); its trailing gate must pass"
./configure -i qmckl

echo "==> [3/4] final configuration"
./configure -c config/gfortran_macos.cfg

echo "==> [4/4] C/Fortran build (ninja all)"
bash -c 'source quantum_package.rc && ninja all'

echo "==> macOS build OK"
