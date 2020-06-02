#!/bin/env python

import gzip
import sys
from math import pi
inv_pi = 1.0/pi

def spec_dens(alpha,beta,z0,g_sign,e_shift):
    sze=len(alpha)
    sze_b=len(beta)
    if (sze != sze_b):
        print('Error: size mismatch',sze,sze_b)
        sys.exit(1)
    z=z0-g_sign*e_shift
    tmp=0.0+0.0j
    #for ai,bi in zip(reversed(a),reversed(b))
    for i in range(sze-1,0,-1):
        tmp=-(beta[i]**2)/(z+g_sign*alpha[i]+tmp)
    tmp=1.0/(z+g_sign*alpha[0]+tmp)
    return -1.0 * tmp.imag * inv_pi

def printspec(ezdir,wmin,wmax,nw,eps):
    gdir=ezdir+'/green/'
    with open(gdir+'n_green_vec') as infile:
        ngvec=int(infile.readline().strip())
    with open(ezdir+'/full_ci_zmq/energy') as infile:
        e0=float(infile.readline().strip())
    with open(gdir+'n_lanczos_complete') as infile:
        nlanc=int(infile.readline().strip())
    
    with gzip.open(gdir+'green_sign.gz') as infile:
        gsign0=infile.read().split()
    
    with gzip.open(gdir+'alpha_lanczos.gz') as infile:
        adata0=infile.read().split()
    with gzip.open(gdir+'beta_lanczos.gz') as infile:
        bdata0=infile.read().split()
    
    adim=int(adata0.pop(0))
    bdim=int(bdata0.pop(0))
    gsigndim=int(gsign0.pop(0))
    assert adim==2, 'dimension of alpha_lanczos should be 2'
    assert bdim==2, 'dimension of beta_lanczos should be 2'
    assert gsigndim==1, 'dimension of green_sign should be 1'
    
    ngvec_2=int(gsign0.pop(0))
    assert ngvec_2==ngvec, 'problem with size of green_sign.gz'
    
    ashape=tuple(map(int,adata0[:adim]))
    bshape=tuple(map(int,bdata0[:bdim]))
    assert ashape==(nlanc,ngvec), 'shape of alpha_lanczos should be (nlanc, ngvec)'
    assert bshape==(nlanc,ngvec), 'shape of  beta_lanczos should be (nlanc, ngvec)'
    
    amat=[]
    for xi in range(ngvec):
        amat.append(list(map(float,adata0[adim+xi*nlanc:adim+(xi+1)*nlanc])))
    
    bmat=[]
    b2mat=[]
    for xi in range(ngvec):
        #bmat.append(list(map(float,bdata0[bdim+xi*nlanc:bdim+(xi+1)*nlanc])))
        b_tmp=list(map(float,bdata0[bdim+xi*nlanc:bdim+(xi+1)*nlanc]))
        b2_tmp=[i*i for i in b_tmp]
        bmat.append(b_tmp)
        b2mat.append(b2_tmp)
    
    gsign=list(map(float,gsign0))
    dw=(wmax-wmin)/(nw-1)
    wlist = [wmin+iw*dw for iw in range(nw)]
    densmat=[]
    for ivec in range(ngvec):
        densmat.append([spec_dens(amat[ivec],bmat[ivec],iw+1.j*eps,gsign[ivec],e0) for iw in wlist])

    for i,dd in enumerate(zip(*densmat)):
        print(('{:15.6E}'+ngvec*'{:25.15E}').format(wlist[i],*dd))

if __name__ == '__main__':

    if len(sys.argv) != 6:
        print('bad args')
        print('USAGE: plot-spec-dens.py ezfio omega_min omega_max n_omega epsilon')
        sys.exit(1)
    ezfio=sys.argv[1]
    wmin=float(sys.argv[2])
    wmax=float(sys.argv[3])
    nw=int(sys.argv[4])
    eps=float(sys.argv[5]) 
    printspec(ezfio,wmin,wmax,nw,eps)


