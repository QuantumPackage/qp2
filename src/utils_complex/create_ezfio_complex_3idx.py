#!/usr/bin/env python
from ezfio import ezfio
import h5py

import sys
import numpy as np
filename = sys.argv[1]
qph5path = sys.argv[2]

ezfio.set_file(filename)
#qph5=h5py.File(qph5path,'r')


ezfio.set_nuclei_is_complex(True)

with h5py.File(qph5path,'r') as qph5:
    kpt_num = qph5['nuclei'].attrs['kpt_num']
    nucl_num = qph5['nuclei'].attrs['nucl_num']
    ao_num = qph5['ao_basis'].attrs['ao_num']
    mo_num = qph5['mo_basis'].attrs['mo_num']
    elec_alpha_num = qph5['electrons'].attrs['elec_alpha_num']
    elec_beta_num = qph5['electrons'].attrs['elec_beta_num']

ezfio.set_nuclei_kpt_num(kpt_num)
kpt_pair_num = (kpt_num*kpt_num + kpt_num)//2
ezfio.set_nuclei_kpt_pair_num(kpt_pair_num)

# don't multiply nuclei by kpt_num
# work in k-space, not in equivalent supercell
nucl_num_per_kpt = nucl_num 
ezfio.set_nuclei_nucl_num(nucl_num_per_kpt)

# these are totals (kpt_num * num_per_kpt)
# need to change if we want to truncate orbital space within pyscf
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
    ezfio.set_ao_basis_ao_basis(qph5['ao_basis'].attrs['ao_basis'])
    ezfio.set_ao_basis_ao_nucl(qph5['ao_basis/ao_nucl'][()].tolist())


#Just need one (can clean this up later)
ao_prim_num_max  = 5

d = [ [0] *ao_prim_num_max]*ao_num
ezfio.set_ao_basis_ao_prim_num([ao_prim_num_max]*ao_num)
ezfio.set_ao_basis_ao_power(d)
ezfio.set_ao_basis_ao_coef(d)
ezfio.set_ao_basis_ao_expo(d)



  
##########################################
#                                        #
#               MO Coef                  #
#                                        #
##########################################
    

with h5py.File(qph5path,'r') as qph5:
    mo_coef_reim = qph5['mo_basis/mo_coef_complex'][()].tolist()
ezfio.set_mo_basis_mo_coef_complex(mo_coef_reim)
#maybe fix qp so we don't need this?
#ezfio.set_mo_basis_mo_coef([[i for i in range(mo_num)] * ao_num])


##########################################
#                                        #
#            Integrals Mono              #
#                                        #
##########################################
    
with h5py.File(qph5path,'r') as qph5:
    if 'ao_one_e_ints' in qph5.keys():
        kin_ao_reim=qph5['ao_one_e_ints/ao_integrals_kinetic'][()].tolist()
        ovlp_ao_reim=qph5['ao_one_e_ints/ao_integrals_overlap'][()].tolist()
        ne_ao_reim=qph5['ao_one_e_ints/ao_integrals_n_e'][()].tolist()

        ezfio.set_ao_one_e_ints_ao_integrals_kinetic_complex(kin_ao_reim)
        ezfio.set_ao_one_e_ints_ao_integrals_overlap_complex(ovlp_ao_reim)
        ezfio.set_ao_one_e_ints_ao_integrals_n_e_complex(ne_ao_reim)
        
        ezfio.set_ao_one_e_ints_io_ao_integrals_kinetic('Read')
        ezfio.set_ao_one_e_ints_io_ao_integrals_overlap('Read')
        ezfio.set_ao_one_e_ints_io_ao_integrals_n_e('Read')
  
    
with h5py.File(qph5path,'r') as qph5:
    if 'mo_one_e_ints' in qph5.keys():
        kin_mo_reim=qph5['mo_one_e_ints/mo_integrals_kinetic'][()].tolist()
        #ovlp_mo_reim=qph5['mo_one_e_ints/mo_integrals_overlap'][()].tolist()
        ne_mo_reim=qph5['mo_one_e_ints/mo_integrals_n_e'][()].tolist()

        ezfio.set_mo_one_e_ints_mo_integrals_kinetic_complex(kin_mo_reim)
        #ezfio.set_mo_one_e_ints_mo_integrals_overlap_complex(ovlp_mo_reim)
        #ezfio.set_mo_one_e_ints_mo_integrals_n_e_complex(ne_mo_reim)
        ezfio.set_mo_one_e_ints_mo_integrals_e_n_complex(ne_mo_reim)
        
        ezfio.set_mo_one_e_ints_io_mo_integrals_kinetic('Read')
        #ezfio.set_mo_one_e_ints_io_mo_integrals_overlap('Read')
        #ezfio.set_mo_one_e_ints_io_mo_integrals_n_e('Read')
        ezfio.set_mo_one_e_ints_io_mo_integrals_e_n('Read')
  
##########################################
#                                        #
#               k-points                 #
#                                        #
##########################################

with h5py.File(qph5path,'r') as qph5:
    kconserv = qph5['nuclei/kconserv'][()].tolist()

ezfio.set_nuclei_kconserv(kconserv)
ezfio.set_nuclei_io_kconserv('Read')
  
##########################################
#                                        #
#             Integrals Bi               #
#                                        #
##########################################

# should this be in ao_basis? ao_two_e_ints?
with h5py.File(qph5path,'r') as qph5:
    if 'ao_two_e_ints' in qph5.keys():
        df_num = qph5['ao_two_e_ints'].attrs['df_num']
        ezfio.set_ao_two_e_ints_df_num(df_num)
        if 'df_ao_integrals' in qph5['ao_two_e_ints'].keys():
#            dfao_re0=qph5['ao_two_e_ints/df_ao_integrals_real'][()].transpose((3,2,1,0))
#            dfao_im0=qph5['ao_two_e_ints/df_ao_integrals_imag'][()].transpose((3,2,1,0))
#            dfao_cmplx0 = np.stack((dfao_re0,dfao_im0),axis=-1).tolist()
#            ezfio.set_ao_two_e_ints_df_ao_integrals_complex(dfao_cmplx0)
            dfao_reim=qph5['ao_two_e_ints/df_ao_integrals'][()].tolist()
            ezfio.set_ao_two_e_ints_df_ao_integrals_complex(dfao_reim)
            ezfio.set_ao_two_e_ints_io_df_ao_integrals('Read')

    if 'mo_two_e_ints' in qph5.keys():
        df_num = qph5['ao_two_e_ints'].attrs['df_num']
        ezfio.set_ao_two_e_ints_df_num(df_num)
#        dfmo_re0=qph5['mo_two_e_ints/df_mo_integrals_real'][()].transpose((3,2,1,0))
#        dfmo_im0=qph5['mo_two_e_ints/df_mo_integrals_imag'][()].transpose((3,2,1,0))
#        dfmo_cmplx0 = np.stack((dfmo_re0,dfmo_im0),axis=-1).tolist()
#        ezfio.set_mo_two_e_ints_df_mo_integrals_complex(dfmo_cmplx0)
        dfmo_reim=qph5['mo_two_e_ints/df_mo_integrals'][()].tolist()
        ezfio.set_mo_two_e_ints_df_mo_integrals_complex(dfmo_reim)
        ezfio.set_mo_two_e_ints_io_df_mo_integrals('Read')


#TODO: add check and only do this if ints exist
#dfmo_re=qph5['mo_two_e_ints/df_mo_integrals_real'][()].transpose((3,2,1,0)).tolist()
#dfmo_im=qph5['mo_two_e_ints/df_mo_integrals_imag'][()].transpose((3,2,1,0)).tolist()
#ezfio.set_mo_two_e_ints_df_mo_integrals_real(dfmo_re)
#ezfio.set_mo_two_e_ints_df_mo_integrals_imag(dfmo_im)
