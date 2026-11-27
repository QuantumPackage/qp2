# quantum_package.rc is bash-only and /bin/sh is dash on some Linux
# distros (no `source` builtin): run all recipes under bash.
SHELL := $(shell command -v bash 2>/dev/null)
ifeq ($(SHELL),)
SHELL := /bin/bash
endif

QP_ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
# Total job budget of the build: the user's top-level -jN when given,
# else nproc; an env/command-line QP_BUILD_JOBS keeps overriding the
# nproc default (as before).
QP_BUILD_JOBS ?= $(shell nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 2)
# Capture the default before the lazy redefinition below (a recursive
# reference to QP_BUILD_JOBS inside its own definition would recurse).
QP_BUDGET_DEFAULT := $(QP_BUILD_JOBS)
# Extract the user's top-level -jN. $(MAKEFLAGS) must be read lazily,
# at recipe (execution) time: with GNU make 3.81 (macOS stock) it is
# still empty while the makefile is parsed. 3.81 records a parallel
# run as a bare `-j` (the number only lives in the jobserver's token
# pool), so when no explicit number is found the budget is probed from
# the jobserver: every free token is read out and handed back, plus
# one, because this recipe already holds its own slot when its line
# is expanded (verified on make 3.81).
QP_BUILD_JOBS = $(or $(strip $(shell python3 -c 'exec("import os,re,select,sys\nmf=sys.argv[1]\nj=re.search(r\"(?:^|\\s)-j(\\d+)\",mf)\njs=re.search(r\"(?:^|\\s)-j\\s+(\\d+)\",mf)\njb=re.search(r\"(?:^|\\s)--jobs=(\\d+)\",mf)\nfs=re.search(r\"--jobserver-fds=(\\d+),(\\d+)\",mf)\nn=int(j.group(1)) if j else (int(js.group(1)) if js else (int(jb.group(1)) if jb else None))\nif n is None and fs:\n    rf,wr=int(fs.group(1)),int(fs.group(2))\n    try:\n        toks=[]\n        while len(toks)<8192 and select.select([rf],[],[],0)[0]:\n            b=os.read(rf,1)\n            if not b: break\n            toks.append(b)\n        os.write(wr,b\"x\"*len(toks))\n        n=len(toks)+1\n    except OSError:\n        n=None\nprint(n if n else \"\")")' "$(MAKEFLAGS)")),$(QP_BUDGET_DEFAULT))
# JFLAGS must be expanded at recipe (execution) time, not at makefile parse
# time: with GNU make 3.81 (macOS stock) $(MAKEFLAGS) is still empty while
# the makefile is parsed, so a parse-time ifeq would never see a user-passed
# -j and would force -jN inside the sub-make, disabling jobserver mode.
# Empty when the user passed -j (the sub-make inherits the jobserver),
# otherwise -j<QP_BUILD_JOBS> so a bare `make` builds in parallel.
JFLAGS = $(if $(findstring j,$(MAKEFLAGS)),, -j$(QP_BUILD_JOBS))

ifneq ($(wildcard $(QP_ROOT)/config/Makefile),)
# ---- make build system ----
# As in the ninja branch below, the build runs inside a shell that has
# sourced quantum_package.rc: the per-module compiler commands carry only
# the cfg-file flags (same as the ninja compile rules) and rely on the
# QP environment for the include dirs (CPATH/C_INCLUDE_PATH ->
# $QP_ROOT/include, see etc/paths.rc), for PYTHONPATH (irpf90 ! shell
# blocks) and for the rest of the QP paths.
.PHONY: rust_target all clean tidy rag update
rust_target all clean tidy rag: update
	+cd $(QP_ROOT) && bash -c 'source quantum_package.rc && exec make $(JFLAGS) -f config/Makefile QP_BUILD_JOBS=$(QP_BUILD_JOBS) $@'
module_%: update
	+cd $(QP_ROOT) && bash -c 'source quantum_package.rc && exec make $(JFLAGS) -f config/Makefile QP_BUILD_JOBS=$(QP_BUILD_JOBS) $@'
update:
	cd $(QP_ROOT) && bash -c 'source quantum_package.rc && exec scripts/compilation/qp_create_ninja update'
else
# ---- ninja build system (or not configured) ----
.PHONY: rust_target all clean tidy rag
rust_target all clean tidy rag:
	cd $(QP_ROOT) && bash -c '[ -f config/build.ninja ] || { echo "QP is not configured yet. Please run ./configure -c config/<config_file>" ; exit 1 ; } ; source quantum_package.rc ; exec ninja $@'
module_%:
	cd $(QP_ROOT) && bash -c 'source quantum_package.rc ; exec ninja -f config/build.ninja $@'
endif
