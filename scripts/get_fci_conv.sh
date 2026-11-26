#!/bin/bash
# Usage: get_fci_conv.sh <EZFIO directory>
#
# Writes <EZFIO>.conv_fci with two columns, Ndet and E+PT2 (hartree,
# state 1) of the last CIPSI iteration, aligned by `column -s ' ' -t`
# as before.
#
# $1 is the EZFIO directory (or file name), no longer a log file:
# qp_run does not print the "# Ndet,E,E+PT2,E+RPT2,absPT2=" lines on
# stdout any more. The data now lives in the newest
# <EZFIO>/json/NNNNN.json result document, read through
# scripts/utility/qp_json.py.

set -euo pipefail

if [ $# -ne 1 ]; then
    echo "usage: $0 <EZFIO directory>" >&2
    exit 1
fi
ezfio="$1"

row="$(python3 - "$0" "$ezfio" <<'PY'
import os
import sys

# The wrapper passes this script's own path, so scripts/utility
# resolves the same way from any working directory.
__file__ = sys.argv[1]
ezfio = sys.argv[2]
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + "/utility")
import qp_json

try:
    doc = qp_json.load_last(ezfio)
except FileNotFoundError as exc:
    where = exc.args[0] if exc.args else "missing"
    sys.exit(f"get_fci_conv: no JSON result document: {where}")

cipsi = (doc.get("results") or {}).get("cipsi")
if cipsi is None:
    sys.exit("get_fci_conv: the newest JSON document has no /results/cipsi block")
if cipsi.get("n_det") is None:
    sys.exit("get_fci_conv: /results/cipsi/n_det is missing")
states = cipsi.get("states") or {}
if not states.get("energy") or not states.get("pt2"):
    sys.exit("get_fci_conv: /results/cipsi/states needs the energy and pt2 arrays")
if states["energy"][0] is None or states["pt2"][0] is None:
    sys.exit("get_fci_conv: state-1 energy or pt2 is null (not finite)")
print(cipsi["n_det"], states["energy"][0] + states["pt2"][0])
PY
)"

printf '%s\n' "$row" | column -s ' ' -t > "$ezfio.conv_fci"
