#!/usr/bin/env python
from ezfio import ezfio
import h5py

import sys
import numpy as np
fname = sys.argv[1]
qph5name = sys.argv[2]

#qph5=h5py.File(qph5path,'r')

def convert_mol(filename,qph5path):
    ezfio.set_file(filename)
    ezfio.set_nuclei_is_complex(False)
    
    with h5py.File(qph5path,'r') as qph5:
        nucl_num = qph5['nuclei'].attrs['nucl_num']
        ao_num = qph5['ao_basis'].attrs['ao_num']
        mo_num = qph5['mo_basis'].attrs['mo_num']
        elec_alpha_num = qph5['electrons'].attrs['elec_alpha_num']
        elec_beta_num = qph5['electrons'].attrs['elec_beta_num']
    
    ezfio.set_nuclei_nucl_num(nucl_num)
    
    ezfio.set_ao_basis_ao_num(ao_num)
    ezfio.set_mo_basis_mo_num(mo_num)
    ezfio.electrons_elec_alpha_num = elec_alpha_num
    ezfio.electrons_elec_beta_num  = elec_beta_num
    
    
    
    ##ao_num = mo_num
    ##Important !
    #import math
    #nelec_per_kpt = num_elec // n_kpts
    #nelec_alpha_per_kpt = int(math.ceil(nelec_per_kpt / 2.))
    #nelec_beta_per_kpt = int(math.floor(nelec_per_kpt / 2.))
    #
    #ezfio.electrons_elec_alpha_num = int(nelec_alpha_per_kpt * n_kpts)
    #ezfio.electrons_elec_beta_num = int(nelec_beta_per_kpt * n_kpts)
    
    #ezfio.electrons_elec_alpha_num = int(math.ceil(num_elec / 2.))
    #ezfio.electrons_elec_beta_num = int(math.floor(num_elec / 2.))
    
    #ezfio.set_utils_num_kpts(n_kpts)
    #ezfio.set_integrals_bielec_df_num(n_aux)
    
    #(old)Important
    #ezfio.set_nuclei_nucl_num(nucl_num)
    #ezfio.set_nuclei_nucl_charge([0.]*nucl_num)
    #ezfio.set_nuclei_nucl_coord( [ [0.], [0.], [0.] ]*nucl_num )
    #ezfio.set_nuclei_nucl_label( ['He'] * nucl_num )
    
    
    with h5py.File(qph5path,'r') as qph5:
        nucl_charge=qph5['nuclei/nucl_charge'][()].tolist()
        nucl_coord=qph5['nuclei/nucl_coord'][()].T.tolist()
        nucl_label=qph5['nuclei/nucl_label'][()].tolist()
        nuclear_repulsion = qph5['nuclei'].attrs['nuclear_repulsion']
    
    ezfio.set_nuclei_nucl_charge(nucl_charge)
    ezfio.set_nuclei_nucl_coord(nucl_coord)
    ezfio.set_nuclei_nucl_label(nucl_label)
    
    ezfio.set_nuclei_io_nuclear_repulsion('Read')
    ezfio.set_nuclei_nuclear_repulsion(nuclear_repulsion)
    
    
    ##########################################
    #                                        #
    #                Basis                   #
    #                                        #
    ##########################################
    
    with h5py.File(qph5path,'r') as qph5:
        do_pseudo = qph5['pseudo'].attrs['do_pseudo']
        ezfio.set_pseudo_do_pseudo(do_pseudo)
        if (do_pseudo):
            ezfio.set_pseudo_pseudo_lmax(qph5['pseudo'].attrs['pseudo_lmax'])
            ezfio.set_pseudo_pseudo_klocmax(qph5['pseudo'].attrs['pseudo_klocmax'])
            ezfio.set_pseudo_pseudo_kmax(qph5['pseudo'].attrs['pseudo_kmax'])
            ezfio.set_pseudo_nucl_charge_remove(qph5['pseudo/nucl_charge_remove'][()].tolist())
            ezfio.set_pseudo_pseudo_n_k(qph5['pseudo/pseudo_n_k'][()].tolist())
            ezfio.set_pseudo_pseudo_n_kl(qph5['pseudo/pseudo_n_kl'][()].tolist())
            ezfio.set_pseudo_pseudo_v_k(qph5['pseudo/pseudo_v_k'][()].tolist())
            ezfio.set_pseudo_pseudo_v_kl(qph5['pseudo/pseudo_v_kl'][()].tolist())
            ezfio.set_pseudo_pseudo_dz_k(qph5['pseudo/pseudo_dz_k'][()].tolist())
            ezfio.set_pseudo_pseudo_dz_kl(qph5['pseudo/pseudo_dz_kl'][()].tolist())
    
    ##########################################
    #                                        #
    #                Basis                   #
    #                                        #
    ##########################################
    
    with h5py.File(qph5path,'r') as qph5:
        coeftmp = qph5['ao_basis/ao_coef'][()]
        expotmp = qph5['ao_basis/ao_expo'][()]
        ezfio.set_ao_basis_ao_basis(qph5['ao_basis'].attrs['ao_basis'])
        ezfio.set_ao_basis_ao_nucl(qph5['ao_basis/ao_nucl'][()].tolist())
        ezfio.set_ao_basis_ao_prim_num(qph5['ao_basis/ao_prim_num'][()].tolist())
        ezfio.set_ao_basis_ao_power(qph5['ao_basis/ao_power'][()].tolist())
        ezfio.set_ao_basis_ao_coef(qph5['ao_basis/ao_coef'][()].tolist())
        ezfio.set_ao_basis_ao_expo(qph5['ao_basis/ao_expo'][()].tolist())
    
    print(coeftmp)
    print(expotmp)
      
    ##########################################
    #                                        #
    #               MO Coef                  #
    #                                        #
    ##########################################
        
    
    with h5py.File(qph5path,'r') as qph5:
        mo_coef = qph5['mo_basis/mo_coef'][()].tolist()
    ezfio.set_mo_basis_mo_coef(mo_coef)
    #maybe fix qp so we don't need this?
    #ezfio.set_mo_basis_mo_coef([[i for i in range(mo_num)] * ao_num])
    
    return

convert_mol(fname,qph5name)
