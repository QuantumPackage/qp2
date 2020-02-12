#!/usr/bin/env python
from ezfio import ezfio
import h5py

import sys
import numpy as np
filename = sys.argv[1]
h5filename = sys.argv[2]
#num_elec, nucl_num, mo_num = map(int,sys.argv[2:5])

#nuclear_repulsion = float(sys.argv[5])
#ao_num = int(sys.argv[6])
#n_kpts = int(sys.argv[7])
#n_aux = int(sys.argv[8])
ezfio.set_file(filename)
qph5=h5py.File(h5filename,'r')


kpt_num = qph5['nuclei'].attrs['kpt_num']
ezfio.set_nuclei_kpt_num(kpt_num)
kpt_pair_num = (kpt_num*kpt_num + kpt_num)//2
ezfio.set_nuclei_kpt_pair_num(kpt_pair_num)

# should this be in ao_basis? ao_two_e_ints?
df_num = qph5['ao_two_e_ints'].attrs['df_num']
ezfio.set_ao_two_e_ints_df_num(df_num)

# these are totals (kpt_num * num_per_kpt)
# need to change if we want to truncate orbital space within pyscf
ezfio.electrons_elec_alpha_num = qph5['electrons'].attrs['elec_alpha_num']
ezfio.electrons_elec_beta_num = qph5['electrons'].attrs['elec_beta_num']
nucl_num = qph5['nuclei'].attrs['nucl_num']
nucl_num_per_kpt = nucl_num // kpt_num
ao_num = qph5['ao_basis'].attrs['ao_num']
mo_num = qph5['mo_basis'].attrs['mo_num']

ezfio.set_mo_basis_mo_num(mo_num)



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

ezfio.set_nuclei_nucl_num(nucl_num_per_kpt)

nucl_charge=qph5['nuclei/nucl_charge'][()].tolist()
ezfio.set_nuclei_nucl_charge(nucl_charge)

nucl_coord=qph5['nuclei/nucl_coord'][()].T.tolist()
ezfio.set_nuclei_nucl_coord(nucl_coord)

nucl_label=qph5['nuclei/nucl_label'][()].tolist()
ezfio.set_nuclei_nucl_label(nucl_label)


ezfio.set_nuclei_io_nuclear_repulsion('Read')
nuclear_repulsion = qph5['nuclei'].attrs['nuclear_repulsion']
ezfio.set_nuclei_nuclear_repulsion(nuclear_repulsion)

# Ao num
#ao_num = mo_num
ezfio.set_ao_basis_ao_basis("Dummy one. We read MO")
ezfio.set_ao_basis_ao_num(ao_num)
ezfio.set_ao_basis_ao_nucl([1]*ao_num) #Maybe put a realy incorrect stuff

#ezfio.set_ao_basis_ao_basis(qph5['ao_basis'].attrs['ao_basis'])
#ezfio.set_ao_basis_ao_nucl(qph5['ao_basis/ao_nucl'][()].tolist())


#Just need one (can clean this up later)
ao_prim_num_max  = 5

d = [ [0] *ao_prim_num_max]*ao_num
ezfio.set_ao_basis_ao_prim_num([ao_prim_num_max]*ao_num)
ezfio.set_ao_basis_ao_power(d)
ezfio.set_ao_basis_ao_coef(d)
ezfio.set_ao_basis_ao_expo(d)



ezfio.set_mo_basis_mo_num(mo_num)
#c_mo = [[1 if i==j else 0 for i in range(mo_num)] for j in range(ao_num)]
#ezfio.set_mo_basis_mo_coef([ [0]*mo_num] * ao_num)
##ezfio.set_mo_basis_mo_coef_real(c_mo)

mo_coef_re0 = qph5['mo_basis/mo_coef_real'][()].T
mo_coef_im0 = qph5['mo_basis/mo_coef_imag'][()].T
mo_coef_cmplx0 = np.stack((mo_coef_re0,mo_coef_im0),axis=-1).tolist()

#ezfio.set_mo_basis_mo_coef_real(qph5['mo_basis/mo_coef_real'][()].tolist())
#ezfio.set_mo_basis_mo_coef_imag(qph5['mo_basis/mo_coef_imag'][()].tolist())
ezfio.set_mo_basis_mo_coef_complex(mo_coef_cmplx0)

#maybe fix qp so we don't need this?
ezfio.set_mo_basis_mo_coef([[i for i in range(mo_num)] * ao_num])

ezfio.set_nuclei_is_complex(True)

# fortran-ordered re,im parts
kin_ao_re0=qph5['ao_one_e_ints/ao_integrals_kinetic_real'][()].T
kin_ao_im0=qph5['ao_one_e_ints/ao_integrals_kinetic_imag'][()].T
#test where to stack? (axis=0 or -1?)
kin_ao_cmplx0=np.stack((kin_ao_re0,kin_ao_im0),axis=-1).tolist()

ovlp_ao_re0=qph5['ao_one_e_ints/ao_integrals_overlap_real'][()].T
ovlp_ao_im0=qph5['ao_one_e_ints/ao_integrals_overlap_imag'][()].T
#test where to stack? (axis=0 or -1?)
ovlp_ao_cmplx0=np.stack((ovlp_ao_re0,ovlp_ao_im0),axis=-1).tolist()

ne_ao_re0=qph5['ao_one_e_ints/ao_integrals_n_e_real'][()].T
ne_ao_im0=qph5['ao_one_e_ints/ao_integrals_n_e_imag'][()].T
#test where to stack? (axis=0 or -1?)
ne_ao_cmplx0=np.stack((ne_ao_re0,ne_ao_im0),axis=-1).tolist()

kin_ao_re=kin_ao_re0.tolist()
kin_ao_im=kin_ao_im0.tolist()
ovlp_ao_re=ovlp_ao_re0.tolist()
ovlp_ao_im=ovlp_ao_im0.tolist()
ne_ao_re=ne_ao_re0.tolist()
ne_ao_im=ne_ao_im0.tolist()

#kin_ao_c = np.stack(kin_ao_re0,kin_ao_im0

#ezfio.set_ao_one_e_ints_ao_integrals_kinetic(kin_ao_re)
#ezfio.set_ao_one_e_ints_ao_integrals_kinetic_imag(kin_ao_im)
ezfio.set_ao_one_e_ints_ao_integrals_kinetic_complex(kin_ao_cmplx0)

#ezfio.set_ao_one_e_ints_ao_integrals_overlap(ovlp_ao_re)
#ezfio.set_ao_one_e_ints_ao_integrals_overlap_imag(ovlp_ao_im)
ezfio.set_ao_one_e_ints_ao_integrals_overlap_complex(ovlp_ao_cmplx0)

#ezfio.set_ao_one_e_ints_ao_integrals_n_e(ne_ao_re)
#ezfio.set_ao_one_e_ints_ao_integrals_n_e_imag(ne_ao_im)
ezfio.set_ao_one_e_ints_ao_integrals_n_e_complex(ne_ao_cmplx0)

dfao_re0=qph5['ao_two_e_ints/df_ao_integrals_real'][()].transpose((3,2,1,0))
dfao_im0=qph5['ao_two_e_ints/df_ao_integrals_imag'][()].transpose((3,2,1,0))
#ezfio.set_ao_two_e_ints_df_ao_integrals_real(dfao_re.tolist())
#ezfio.set_ao_two_e_ints_df_ao_integrals_imag(dfao_im.tolist())
dfao_cmplx0 = np.stack((dfao_re0,dfao_im0),axis=-1).tolist()
ezfio.set_ao_two_e_ints_df_ao_integrals_complex(dfao_cmplx0)


#TODO: add check and only do this if ints exist
#dfmo_re=qph5['mo_two_e_ints/df_mo_integrals_real'][()].transpose((3,2,1,0)).tolist()
#dfmo_im=qph5['mo_two_e_ints/df_mo_integrals_imag'][()].transpose((3,2,1,0)).tolist()
#ezfio.set_mo_two_e_ints_df_mo_integrals_real(dfmo_re)
#ezfio.set_mo_two_e_ints_df_mo_integrals_imag(dfmo_im)
