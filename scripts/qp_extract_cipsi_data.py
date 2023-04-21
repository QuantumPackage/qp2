#!/usr/bin/env python3

import re
import sys

# Read output file
with open(sys.argv[1], 'r') as file:
  output = file.read()


def extract_data(output):
    lines = output.split("\n")
    data = []

    n_det = None
    e = None
    pt2 = None
    err_pt2 = None
    rpt2 = None
    err_rpt2 = None
    e_ex = None


    reading = False
    for iline, line in enumerate(lines):
        if line.startswith("Summary at N_det"):
            reading = False

        if not reading and line.startswith(" N_det "):
            n_det = int(re.search(r"N_det\s+=\s+(\d+)", line).group(1))
            reading = True

        if reading:
            if line.startswith(" E "):
                e = float(re.search(r"E\s+=\s+(-?\d+\.\d+)", line).group(1))
            elif line.startswith(" PT2 "):
                pt2 = float(re.search(r"PT2\s+=\s+(-?\d+\.\d+E?.\d*)", line).group(1))
                err_pt2 = float(re.search(r"\+/-\s+(-?\d+\.\d+E?.\d*)", line).group(1))
            elif line.startswith(" rPT2 "):
                rpt2 = float(re.search(r"rPT2\s+=\s+(-?\d+\.\d+E?.\d*)", line).group(1))
                err_rpt2 = float(re.search(r"\+/-\s+(-?\d+\.\d+E?.\d*)", line).group(1))
            elif "minimum PT2 Extrapolated energy" in line:
                e_ex_line = lines[iline+2]
                e_ex = float(e_ex_line.split()[1])
                reading = False
                new_data = " {:8d}  {:16.8f}  {:e}  {:e}  {:e}  {:e}  {:16.8f}".format(n_det, e, pt2, err_pt2, rpt2, err_rpt2, e_ex)
                data.append(new_data)
                n_det = e = pt2 = err_pt2 = rpt2 = err_rpt2 = e_ex = None

    return data

data = extract_data(output)

for item in data:
    print(item)
