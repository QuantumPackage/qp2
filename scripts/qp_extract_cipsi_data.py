#!/usr/bin/env python3

import qp_json
import sys

if len(sys.argv) == 1:
  print(f"syntax: {sys.argv[0]} EZFIO_FILE")

d = qp_json.load_all(sys.argv[1])

print("# Energy  PT2   PT2_err  rPT2  rPT2_err  exFCI\n")
for f in sorted(d):
    try:
        j = d[f]['results']['cipsi']
    except KeyError:
        continue

    print(f"# {f}")
    st = j['states']
    n = len(st['energy'])
    out = f" {j['n_det']:8d}"
    ex = st.get('extrapolated_energy', [0.0]*n)
    for i in range(n):
        out += (f" {st['energy'][i]:16.8f}  {st['pt2'][i]:e}  {st['pt2_err'][i]:e}"
                f"  {st['rpt2'][i]:e}  {st['rpt2_err'][i]:e}  {ex[i]:16.8f}")
    print(out)
    print("\n")
