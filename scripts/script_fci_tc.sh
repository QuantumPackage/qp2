 source ~/qp2/quantum_package.rc  
 alpha=1.8
 input=O
 basis=cc-pvdz
 mult=3
 output=${input}_${basis}_al_${alpha}
 qp create_ezfio -b ${basis} ${input}.xyz -m $mult
 qp run scf 
 qp set perturbation pt2_max 0.0001
 qp set_frozen_core 

########## FCI CALCULATION FOR REFERENCE 
 qp run fci | tee ${EZFIO_FILE}.fci.out
 qp run sort_wf 
 mv ${EZFIO_FILE}.wf_sorted ${EZFIO_FILE}_fci.wf_sorted
########### TC SCF CALCULATION
 qp reset -d 
 qp set ao_two_e_erf_ints mu_erf 0.87 
 qp set tc_keywords j1b_type 3
 qp set tc_keywords j1b_pen "[${alpha}]"
 qp set tc_keywords bi_ortho True 
 qp set tc_keywords test_cycle_tc True
 qp set tc_keywords write_tc_integ True
 qp set tc_keywords read_tc_integ False
 qp run tc_scf | tee ${EZFIO_FILE}.tc_scf.out
 qp set tc_keywords write_tc_integ False
 qp set tc_keywords read_tc_integ True
############ TC-FCI CALCULATION
 qp run fci_tc_bi_ortho | tee ${EZFIO_FILE}.fci_tc_bi_ortho.out
 grep "Ndet,E,E+PT2,E+RPT2,|PT2|=" ${EZFIO_FILE}.fci_tc_bi_ortho.out | cut -d "=" -f 2 > data_al_$alpha
 qp run sort_wf 
 mv ${EZFIO_FILE}.wf_sorted ${EZFIO_FILE}_tc_fci.wf_sorted
 
