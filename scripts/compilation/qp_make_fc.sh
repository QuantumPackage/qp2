#!/usr/bin/env bash
#
# qp_make_fc.sh: the FC (and CC) of the QP make build.
#
# The generated Makefiles set FC (and CC) to this script. It prints
# one ninja-style `[k/N] <description>` progress line through
# qp_make_step.py -- the full real command instead of the short
# description when V=1 (ninja's own default behavior) -- and then
# execs the real compiler with all its arguments untouched.
#
# Required environment:
#   QP_MAKE_REAL_FC or QP_MAKE_REAL_CC  the real compiler command;
#                                        it may carry flags (it is
#                                        word-split on exec)
#   QP_MAKE_COUNTER   the state file of the scope being counted
#   QP_MAKE_N         the total the counter prints
# Optional:
#   QP_MAKE_MODULE    description prefix ("<module>: <description>")
#   V=1               full real command lines instead of the short
#                     descriptions
#
# The generated makefiles wire both FC and CC to this same script,
# so the real compiler is picked by the language of the arguments:
# a C compile (a .c source) goes to QP_MAKE_REAL_CC, every Fortran
# compile and every link to QP_MAKE_REAL_FC (a link takes objects,
# and the QP builds link their programs with the Fortran compiler).
#
# It runs once per compile/link: dead simple on purpose.

set -u

real=${QP_MAKE_REAL_FC:-}
for a in "$@"; do
    case "$a" in
    *.c)
        real=${QP_MAKE_REAL_CC:-$real}
        ;;
    esac
done
if [ -z "$real" ]; then
    echo "qp_make_fc.sh: QP_MAKE_REAL_FC (or QP_MAKE_REAL_CC) is not set" >&2
    exit 1
fi

# Short description from the compiler arguments: -c with -o names
# the object (minus its .o), -c without -o names the source (the
# ezfio compiles), -o alone is a link, nothing is an unnamed compile.
desc=compiling
have_c=0
have_o=0
oarg=
carg=
prev=
for a in "$@"; do
    case "$prev" in
    -o)
        oarg=$a
        have_o=1
        ;;
    -c)
        carg=$a
        have_c=1
        ;;
    esac
    prev=$a
done
if [ "$have_c" -eq 1 ] && [ "$have_o" -eq 1 ]; then
    b=${oarg##*/}
    desc="compiling ${b%.o}"
elif [ "$have_c" -eq 1 ]; then
    desc="compiling ${carg##*/}"
elif [ "$have_o" -eq 1 ]; then
    desc="linking ${oarg##*/}"
fi

if [ -n "${QP_MAKE_MODULE:-}" ]; then
    desc="$QP_MAKE_MODULE: $desc"
fi
if [ "${V:-}" = "1" ]; then
    desc="$real $*"
fi

script_dir=$(cd "$(dirname "$0")" && pwd)
# A broken counter must never fail the build: the step line is
# cosmetic, the compile is not.
python3 "$script_dir/qp_make_step.py" next "$QP_MAKE_COUNTER" \
    "$QP_MAKE_N" "$desc" || :

# $real is deliberately unquoted: it may carry flags
# (e.g. "gfortran -ffree-line-length-none -I . -g -fPIC").
exec $real "$@"
