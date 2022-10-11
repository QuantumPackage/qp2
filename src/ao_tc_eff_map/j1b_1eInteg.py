import sys, os
QP_PATH=os.environ["QP_EZFIO"]
sys.path.insert(0,QP_PATH+"/Python/")
from ezfio import ezfio
from datetime import datetime
import time
from math import exp, sqrt, pi
import numpy as np
import subprocess
from scipy.integrate import tplquad
import multiprocessing
from multiprocessing import Pool


# _____________________________________________________________________________
#
def read_ao():

    with open('ao_data') as f:
        lines = f.readlines()

    ao_prim_num = np.zeros((ao_num), dtype=int)
    ao_nucl     = np.zeros((ao_num), dtype=int)
    ao_power    = np.zeros((ao_num, 3))
    nucl_coord  = np.zeros((ao_num, 3))
    ao_expo     = np.zeros((ao_num, ao_num))
    ao_coef     = np.zeros((ao_num, ao_num))

    iline = 0
    for j in range(ao_num):
       
        line   = lines[iline]
        iline += 1
        ao_nucl[j] = int(line) - 1

        line   = lines[iline].split()
        iline += 1
        ao_power[j, 0] = float(line[0])
        ao_power[j, 1] = float(line[1])
        ao_power[j, 2] = float(line[2])

        line   = lines[iline].split()
        iline += 1
        nucl_coord[ao_nucl[j], 0] = float(line[0])
        nucl_coord[ao_nucl[j], 1] = float(line[1])
        nucl_coord[ao_nucl[j], 2] = float(line[2])

        line   = lines[iline]
        iline += 1
        ao_prim_num[j] = int(line)

        for l in range(ao_prim_num[j]):

            line   = lines[iline].split()
            iline += 1
            ao_expo[l, j] = float(line[0])
            ao_coef[l, j] = float(line[1])

    return( ao_prim_num
          , ao_nucl
          , ao_power
          , nucl_coord 
          , ao_expo 
          , ao_coef )
# _____________________________________________________________________________


# _____________________________________________________________________________
#
def Gao(X, i_ao):

    ii  = ao_nucl[i_ao]
    C   = np.array([nucl_coord[ii,0], nucl_coord[ii,1], nucl_coord[ii,2]])
    Y   = X - C
    dis = np.dot(Y,Y)

    ip  = np.array([ao_power[i_ao,0], ao_power[i_ao,1], ao_power[i_ao,2]])
    pol = np.prod(Y**ip)

    xi  = np.sum( ao_coef[:,i_ao] * np.exp(-dis*ao_expo[:,i_ao]) )

    return(xi*pol)
# _____________________________________________________________________________


# _____________________________________________________________________________
#
def grad_Gao(X, i_ao):

    ii = ao_nucl[i_ao]
    C  = np.array([nucl_coord[ii,0], nucl_coord[ii,1], nucl_coord[ii,2]])

    ix = ao_power[i_ao,0]
    iy = ao_power[i_ao,1]
    iz = ao_power[i_ao,2]

    Y   = X - C
    dis = np.dot(Y,Y)

    xm = np.sum(                ao_coef[:,i_ao]*np.exp(-dis*ao_expo[:,i_ao]))
    xp = np.sum(ao_expo[:,i_ao]*ao_coef[:,i_ao]*np.exp(-dis*ao_expo[:,i_ao]))

    ip  = np.array([ix+1, iy, iz])
    dx  = -2. * np.prod(Y**ip) * xp 
    if(ix > 0):
        ip  = np.array([ix-1, iy, iz])
        dx += ix * np.prod(Y**ip) * xm

    ip  = np.array([ix, iy+1, iz])
    dy  = -2. * np.prod(Y**ip) * xp 
    if(iy > 0):
        ip  = np.array([ix, iy-1, iz])
        dy += iy * np.prod(Y**ip) * xm

    ip  = np.array([ix, iy, iz+1])
    dz  = -2. * np.prod(Y**ip) * xp 
    if(iz > 0):
        ip  = np.array([ix, iy, iz-1])
        dz += iz * np.prod(Y**ip) * xm

    return(np.array([dx, dy, dz]))
# _____________________________________________________________________________


# _____________________________________________________________________________
#
#              3 x < XA |       exp[-gama r_C^2] | XB > 
#            - 2 x < XA | r_A^2 exp[-gama r_C^2] | XB >
#
def integ_lap(z, y, x, i_ao, j_ao):

    X = np.array([x, y, z])

    Gi = Gao(X, i_ao)
    Gj = Gao(X, j_ao)

    c = 0.
    for k in range(nucl_num):
        gama = j1b_gauss_pen[k]
        C    = nucl_coord[k,:]  
        Y    = X - C
        dis  = np.dot(Y, Y)
        arg  = exp(-gama*dis)
        arg  = exp(-gama*dis)
        c += ( 3. - 2. * dis * gama ) * arg * gama * Gi * Gj 

    return(c)
# _____________________________________________________________________________


# _____________________________________________________________________________
#
#
def integ_grad2(z, y, x, i_ao, j_ao):

    X = np.array([x, y, z])

    Gi = Gao(X, i_ao)
    Gj = Gao(X, j_ao)

    c = np.zeros((3))
    for k in range(nucl_num):
        gama = j1b_gauss_pen[k]
        C    = nucl_coord[k,:]  
        Y    = X - C
        c   += gama * exp(-gama*np.dot(Y, Y)) * Y

    return(-2*np.dot(c,c)*Gi*Gj)
# _____________________________________________________________________________


# _____________________________________________________________________________
#
#
def integ_nonh(z, y, x, i_ao, j_ao):

    X = np.array([x, y, z])

    Gi = Gao(X, i_ao)

    c = 0.
    for k in range(nucl_num):
        gama = j1b_gauss_pen[k]
        C    = nucl_coord[k,:]  
        Y    = X - C
        grad = grad_Gao(X, j_ao)
        c   += gama * exp(-gama*np.dot(Y,Y)) * np.dot(Y,grad)

    return(2*c*Gi)
# _____________________________________________________________________________


# _____________________________________________________________________________
#
def perform_integ( ind_ao ):

    i_ao = ind_ao[0]
    j_ao = ind_ao[1]

    a      = -15. #-np.Inf
    b      = +15. #+np.Inf
    epsrel = 1e-5

    res_lap, err_lap = tplquad( integ_lap
                              , a, b
                              , lambda x  : a, lambda x  : b
                              , lambda x,y: a, lambda x,y: b
                              , (i_ao, j_ao)
                              , epsrel=epsrel )

    res_grd, err_grd = tplquad( integ_grad2
                              , a, b
                              , lambda x  : a, lambda x  : b
                              , lambda x,y: a, lambda x,y: b
                              , (i_ao, j_ao)
                              , epsrel=epsrel )

    res_nnh, err_nnh = tplquad( integ_nonh
                              , a, b
                              , lambda x  : a, lambda x  : b
                              , lambda x,y: a, lambda x,y: b
                              , (i_ao, j_ao)
                              , epsrel=epsrel )
  
    return( [ res_lap, err_lap
            , res_grd, err_grd
            , res_nnh, err_nnh ])
# _____________________________________________________________________________


# _____________________________________________________________________________
#
def integ_eval():

    list_ind = []
    for i_ao in range(ao_num):
        for j_ao in range(ao_num):
            list_ind.append( [i_ao, j_ao] )

    nb_proc = multiprocessing.cpu_count()
    print(" --- Excexution with {} processors ---\n".format(nb_proc))

    p   = Pool(nb_proc)
    res = np.array( p.map( perform_integ, list_ind ) )

    ii = 0
    for i_ao in range(ao_num):
        for j_ao in range(ao_num):
            print(" {} {}     {:+e} {:+e}     {:+e} {:+e}".format( i_ao, j_ao
                 , res[ii][0], res[ii][1], res[ii][2], res[ii][3]) )
            ii += 1

    p.close()
# _____________________________________________________________________________



# _____________________________________________________________________________
#
if __name__=="__main__":
  
    t0 = time.time()

    EZFIO_file = sys.argv[1]
    ezfio.set_file(EZFIO_file)

    print(" Today's date:", datetime.now() )
    print(" EZFIO file = {}".format(EZFIO_file))

    nucl_num      = ezfio.get_nuclei_nucl_num()
    ao_num        = ezfio.get_ao_basis_ao_num()
    j1b_gauss_pen = ezfio.get_ao_tc_eff_map_j1b_gauss_pen() 

    ao_prim_num, ao_nucl, ao_power, nucl_coord, ao_expo, ao_coef = read_ao()

    #integ_eval()

    i_ao   = 0 
    j_ao   = 0

    a      = -5. 
    b      = +5. 
    epsrel = 1e-1
    res_grd, err_grd = tplquad( integ_nonh
                              , a, b
                              , lambda x  : a, lambda x  : b
                              , lambda x,y: a, lambda x,y: b
                              , (i_ao, j_ao)
                              , epsrel=epsrel )

    print(res_grd, err_grd)
 

    tf = time.time() - t0
    print(' end after {} min'.format(tf/60.))
# _____________________________________________________________________________



