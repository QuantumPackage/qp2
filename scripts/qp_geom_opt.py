#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Usage:
       qp_geom_opt [-s state] [-r executable] [-f] [-t tolerance] <EZFIO_FILE>

Options:
      -s --state=<state>       Excited state to optimize
      -f --scf                 Perform an SCF after each geomety change
      -r --qp_run=executable   Excited state to optimize
      -t --tol=tolerance       Convergence criterion on the energy
"""


try:
    from docopt import docopt
    from module_handler import ModuleHandler, get_dict_child
    from module_handler import get_l_module_descendant
    from qp_path import QP_SRC, QP_PLUGINS, QP_DATA, QP_ROOT
except ImportError:
    print("Please check if you have sourced the ${QP_ROOT}/quantum_package.rc")
    print("(`source ${QP_ROOT}/quantum_package.rc`)")
    print(sys.exit(1))


import numpy as np
import subprocess
from scipy.optimize import minimize
from ezfio import ezfio

import sys


def set_unbuffered_output():
    """Ensure sys.stdout is unbuffered or line-buffered in a portable way."""
    if hasattr(sys.stdout, "reconfigure"):  # Python 3.7+
        sys.stdout.reconfigure(line_buffering=True)
    else:
        sys.stdout = open(sys.stdout.fileno(), mode='w', buffering=1)

set_unbuffered_output()




def get_energy(file, state, arguments):
    """Compute the energy of the given state by calling Quantum Package."""
    if not arguments["--qp_run"]:
        raise ValueError("--qp_run option missing")

    if arguments["--scf"]:
        executable = "scf"
    else:
        executable = "save_ortho_mos"

    result = subprocess.run(f"qp_run {executable} {file} > {file}.energy.out",
                            shell=True, capture_output=True, text=True, check=True
    )

    executable = arguments["--qp_run"]
    result = subprocess.run( f"qp_run {executable} {file} > {file}.energy.out",
                             shell=True)

    energy = None
    with open(f"{file}.energy.out", 'r') as f:
        for line in f:
            if "Energy of state" in line and f"{state}" in line:
                energy = float(line.split()[-1])  # Extracts the energy value

    return energy
    raise ValueError("Energy not found in Quantum Package output. Update script {sys.argv[0]}")

def set_coordinates(coord):
    """Update the nuclear coordinates in EZFIO."""
    ezfio.set_nuclei_nucl_coord(coord)


def get_coordinates():
    """Retrieve the current nuclear coordinates from EZFIO."""
    return np.array(ezfio.get_nuclei_nucl_coord())


memo_energy = {}
def energy_function(coord, file, state, arguments):
    """Wrapper for the energy calculation, ensuring coordinates are updated."""
    h = np.array_str(coord)
    if h in memo_energy:
        return memo_energy[h]

    set_coordinates(coord)
    energy = get_energy(file, state, arguments)
    memo_energy[h] = energy

    label = ezfio.get_nuclei_nucl_label()
    num_atoms = len(label)
    coord = coord.reshape(3, num_atoms).T  # Reshape into (num_atoms, 3)
    coord_angstrom = coord * 0.529177  # Convert atomic units to angstroms

    print(num_atoms)
    print(f"Energy: {energy:15.10f}")
    for i, (x, y, z) in enumerate(coord_angstrom):
        print(f"{label[i]:3s} {x:15.8f} {y:15.8f} {z:15.8f}")  # Replace 'X' with actual atomic symbols
    return energy


def optimize_geometry(file, state, arguments):
    """Perform geometry optimization using SciPy's minimize function."""

    x0 = get_coordinates().flatten()

    if arguments["--tol"]:
        tolerance = float(arguments["--tol"])
    else:
        tolerance = 1.e-3

    result = minimize(energy_function, x0, args=(file, state, arguments),
                      method='Powell',
                      tol=tolerance,
                      options={'xtol': tolerance, 'ftol': tolerance})

#    result = minimize(energy_function, x0, args=(file, state, arguments),
#                      method='BFGS',
#                      jac=None,
#                      tol=tolerance,
#                      options={'eps': 1.e-3})

    if result.success:
        print("Optimization successful!")
        print("Final energy:", result.fun)
        print("Optimized coordinates:", result.x)
    else:
        print("Optimization failed:", result.message)

    set_coordinates(result.x)  # Store the optimized geometry
    return result


def main(arguments):
    if arguments["--state"]:
        state=arguments["--state"]
    else:
        state=1
    ezfio_filename = arguments["<EZFIO_FILE>"]
    ezfio.set_file(ezfio_filename)

    optimize_geometry(ezfio_filename, state, arguments)


if __name__ == "__main__":
    ARG = docopt(__doc__)
    main(ARG)
