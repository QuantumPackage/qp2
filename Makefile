QP_ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
QP_BUILD_JOBS ?= $(shell nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 2)
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
	cd $(QP_ROOT) && bash -c 'source quantum_package.rc && exec make $(JFLAGS) -f config/Makefile $@'
module_%: update
	cd $(QP_ROOT) && bash -c 'source quantum_package.rc && exec make $(JFLAGS) -f config/Makefile $@'
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
