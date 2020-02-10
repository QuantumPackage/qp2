#!/usr/bin/env python
from ezfio import ezfio
import h5py

import sys
filename = sys.argv[1]
h5filename = sys.argv[2]
#num_elec, nucl_num, mo_num = map(int,sys.argv[2:5])

#nuclear_repulsion = float(sys.argv[5])
#ao_num = int(sys.argv[6])
#n_kpts = int(sys.argv[7])
#n_aux = int(sys.argv[8])
ezfio.set_file(filename)
qph5=h5py.File(h5filename.'r')

ezfio.electrons_elec_alpha_num = qph5['electrons'].attrs['elec_alpha_num']
ezfio.electrons_elec_beta_num = qph5['electrons'].attrs['elec_beta_num']

nucl_num = qph5['nuclei'].attrs['nucl_num']
kpt_num = qph5['nuclei'].attrs['kpt_num']
#df_num = qph5['???'].attrs['df_num']

mo_num = qph5['mo_basis'].attrs['mo_num']
ezfio.set_mo_basis_mo_num(mo_num)

#ao_num = mo_num
#Important !
import math
nelec_per_kpt = num_elec // n_kpts
nelec_alpha_per_kpt = int(math.ceil(nelec_per_kpt / 2.))
nelec_beta_per_kpt = int(math.floor(nelec_per_kpt / 2.))

ezfio.electrons_elec_alpha_num = int(nelec_alpha_per_kpt * n_kpts)
ezfio.electrons_elec_beta_num = int(nelec_beta_per_kpt * n_kpts)

#ezfio.electrons_elec_alpha_num = int(math.ceil(num_elec / 2.))
#ezfio.electrons_elec_beta_num = int(math.floor(num_elec / 2.))

#ezfio.set_utils_num_kpts(n_kpts)
#ezfio.set_integrals_bielec_df_num(n_aux)

#Important
ezfio.set_nuclei_nucl_num(nucl_num)
ezfio.set_nuclei_nucl_charge([0.]*nucl_num)
ezfio.set_nuclei_nucl_coord( [ [0.], [0.], [0.] ]*nucl_num )
ezfio.set_nuclei_nucl_label( ['He'] * nucl_num )

ezfio.set_nuclei_io_nuclear_repulsion('Read')
ezfio.set_nuclei_nuclear_repulsion(nuclear_repulsion)

# Ao num
#ao_num = mo_num
ezfio.set_ao_basis_ao_basis("Dummy one. We read MO")
ezfio.set_ao_basis_ao_num(ao_num)
ezfio.set_ao_basis_ao_nucl([1]*ao_num) #Maybe put a realy incorrect stuff

#Just need one
ao_prim_num_max  = 5

d = [ [0] *ao_prim_num_max]*ao_num
ezfio.set_ao_basis_ao_prim_num([ao_prim_num_max]*ao_num)
ezfio.set_ao_basis_ao_power(d)
ezfio.set_ao_basis_ao_coef(d)
ezfio.set_ao_basis_ao_expo(d)

#Dummy one
ao_md5 = '3b8b464dfc95f282129bde3efef3c502'
ezfio.set_ao_basis_ao_md5(ao_md5)
ezfio.set_mo_basis_ao_md5(ao_md5)


ezfio.set_mo_basis_mo_num(mo_num)
c_mo = [[1 if i==j else 0 for i in range(mo_num)] for j in range(ao_num)]
ezfio.set_mo_basis_mo_coef([ [0]*mo_num] * ao_num)
#ezfio.set_mo_basis_mo_coef_real(c_mo)

ezfio.set_nuclei_is_periodic(True)
