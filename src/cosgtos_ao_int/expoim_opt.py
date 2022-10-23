#!/usr/bin/env python

import sys
import os
import subprocess
from datetime import datetime
import time
import numpy as np
from modif_powell_imp import my_fmin_powell

QP_PATH=os.environ["QP_ROOT"]
sys.path.insert(0, QP_PATH+"external/ezfio/Python")
from ezfio import ezfio





#------------------------------------------------------------------------------
#
def get_expoim():

    expo_im = np.array(ezfio.get_cosgtos_ao_int_ao_expoim_cosgtos()).T
    #print(expo_im.shape)

    x = []
    for i in range(ao_num):
        for j in range(ao_prim_num[i]):
            x.append(expo_im[i,j])

    return x

# ---

def set_expoim(x):

    expo_im = np.zeros((ao_num, ao_prim_num_max))

    ii = 0
    for i in range(ao_num):
        for j in range(ao_prim_num[i]):
            expo_im[i,j] = x[ii]
            ii = ii + 1

    ezfio.set_cosgtos_ao_int_ao_expoim_cosgtos(expo_im.T)
#
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
#
def save_res(results, file_output):

    lines = results.splitlines()
    with open(file_output, "w") as f:
        for line in lines:
            f.write(f"{line}\n")

#
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
#
def get_scfenergy(results):

    scf_energy = 0.0

    lines = results.splitlines()
    for line in lines:
        if("SCF energy" in line):
            scf_energy = float(line.split()[-1])

    return scf_energy
#
#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
#
def run_scf():

    return subprocess.check_output( ['qp_run', 'scf', EZFIO_file]
                                  , encoding = "utf-8" )
#
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
#
def f_scf(x):

    global i_call
    i_call += 1

    #print(x)

    # set expo
    set_expoim(x)

    # run scf
    results = run_scf()
    #save_res(results, "scf_"+str(i_call))
    

    # get scf_energy
    scf_energy = get_scfenergy(results)
    print( scf_energy )
    sys.stdout.flush()

    return scf_energy
#
#------------------------------------------------------------------------------




if __name__ == '__main__':

    t0 = time.time()

    EZFIO_file = sys.argv[1]
    ezfio.set_file(EZFIO_file)
    print(" Today's date:", datetime.now() )
    print(" EZFIO file = {}".format(EZFIO_file))


    ao_num = ezfio.get_ao_basis_ao_num()
    print(f" ao_num = {ao_num}")

    ao_prim_num = ezfio.get_ao_basis_ao_prim_num()

    ao_prim_num_max = np.amax(ao_prim_num)
    print(f" ao_prim_num_max = {ao_prim_num_max}")

    ezfio.set_ao_basis_ao_prim_num_max(ao_prim_num_max)

    x = get_expoim()

    n_par = len(x)
    print(' nb of parameters = {}'.format(n_par))

    sys.stdout.flush()

    #x    = (np.random.rand(n_par) - 0.5) * 1.0
    x     = [ (+0.00) for _ in range(n_par)]

    x_min = [ (-10.0) for _ in range(n_par)]
    x_max = [ (+10.0) for _ in range(n_par)]

    i_call = 0
    memo_val = {'fmin': 100.}

    opt = my_fmin_powell( f_scf
                        , x, x_min, x_max
                        #, xtol        = 1e-1
                        #, ftol        = 1e-1
                        , maxfev      = 1e8
                        , full_output = 1
                        , verbose     = 1 )


    print(" x = " + str(opt))

    print(" end after {:.3f} minutes".format((time.time()-t0)/60.) )

    # !!!
# !!!
