#!/bin/bash

ezfio=$1
h5file=$2
# Create the integral
echo 'Create Integral'

echo 'Create EZFIO'
#read nel nmo natom <<< $(cat param) 
#read e_nucl <<< $(cat e_nuc)
#read nao <<< $(cat num_ao)
#read nkpts <<< $(cat kpt_num)
#read ndf <<< $(cat num_df)
##./create_ezfio_complex_4idx.py $ezfio $nel $natom $nmo $e_nucl $nao $nkpts
./create_ezfio_complex_3idx.py $ezfio $h5file #$nel $natom $nmo $e_nucl $nao $nkpts $ndf
#Handle the orbital consitensy check
qp_edit -c $ezfio &> /dev/null
#cp $ezfio/{ao,mo}_basis/ao_md5 

qp_run import_kconserv $ezfio 
#qp_run import_ao_2e_complex $ezfio 
#qp_run dump_ao_2e_from_df $ezfio
#Read the integral
#echo 'Read Integral'


################################################
##  using AO mono, 4-idx from pyscf           ##
################################################
#qp_run import_integrals_ao_periodic $ezfio 


################################################
##  using AO mono, 3-idx, mo coef from pyscf  ##
################################################

#qp_run read_ao_mono_complex $ezfio 
#qp_run read_kconserv $ezfio
#qp_run read_ao_df_complex $ezfio
#qp_run read_mo_coef_complex $ezfio    #start from converged pyscf MOs
#
#qp_run save_mo_df_to_disk $ezfio
#qp_run save_mo_bielec_to_disk $ezfio

#qp_run mo_from_ao_orth $ezfio        #use canonical orthonormalized AOs as initial MO guess
#qp_run print_H_matrix_restart $ezfio > hmat.out


###############################################################
##  using AO mono, full 4-idx AO bielec, mo coef from pyscf  ##
###############################################################

#qp_run read_ao_mono_complex $ezfio 
#qp_run read_kconserv $ezfio
#qp_run read_ao_eri_chunk_complex $ezfio 
#qp_run read_mo_coef_complex $ezfio    #start from converged pyscf MOs
##qp_run mo_from_ao_orth $ezfio        #use canonical orthonormalized AOs as initial MO guess


######################################################
##  using MO mono, full 4-idx MO bielec from pyscf  ##
######################################################

#qp_run read_mo_mono_complex $ezfio 
#qp_run read_kconserv $ezfio
#qp_run read_mo_eri_chunk_complex $ezfio 

