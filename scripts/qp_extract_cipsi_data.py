#!/usr/bin/env python3

import qp_json
import sys

if len(sys.argv) == 1:
  print(f"syntax: {sys.argv[0]} EZFIO_FILE")

d = qp_json.load_all(sys.argv[1])

k = [ x for x in d.keys() ]
k.sort()

print("# Energy  PT2   PT2_err  rPT2  rPT2_err  exFCI\n")
for f in k:
  try:
    j = d[f]["fci"]
  except:
     continue

  print(f"# {f}")
  for e in j:

       out = f" {e['n_det']:8d}"

       nstates = len(e["states"])
       for ee in e["states"]:
         try:
            exc_energy = ee['ex_energy'][0]
         except:
            exc_energy = 0.
         out += f" {ee['energy']:16.8f}  {ee['pt2']:e}  {ee['pt2_err']:e}  {ee['rpt2']:e}  {ee['rpt2_err']:e}  {exc_energy:16.8f}"
       print(out)

  print("\n")


