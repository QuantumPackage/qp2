#!/bin/bash

source ~/qp2/quantum_package.rc

## Define the system/basis/charge/mult and genric keywords 
system=H2O
xyz=${system}.xyz
basis=6-31g
mult=1
charge=0
j2e_type="Boys_Handy"
thresh_tcscf=1e-10
io_tc_integ="Write"
nstates=4



##################### Function to create the EZFIO 
function create_ezfio (){
 qp create_ezfio -b $basis -m $mult -c $charge $xyz -o $ezfio
 qp run scf | tee ${EZFIO_FILE}.scf.out 
}

##################### Function to set parameters for BH9 jastrow 
function BH_9 (){
 j2e_type="Boys_Handy"   # type of correlation factor: Boys Handy type 
 env_type="None"         # Boys Handy J does not use our envelopes 
 j1e_type="None"         # Boys Handy J does not use our J1body
 tc_integ_type="numeric" # Boys Handy requires numerical integrals 
 jBH_size=9              # Number of parameters for the BH 

######## All parameters for the H2O and Boys Handy Jastrow 
 jBH_c=[[0.50000,-0.57070,0.49861,-0.78663,0.01990,0.13386,-0.60446,-1.67160,1.36590],[0.0,0.0,0.0,0.0,0.12063,-0.18527,0.12324,-0.11187,-0.06558],[0.0,0.0,0.0,0.0,0.12063,-0.18527,0.12324,-0.11187,-0.06558]]
 jBH_m=[[0,0,0,0,2,3,4,2,2],[0,0,0,0,2,3,4,2,2],[0,0,0,0,2,3,4,2,2]]
 jBH_n=[[0,0,0,0,0,0,0,2,0],[0,0,0,0,0,0,0,2,0],[0,0,0,0,0,0,0,2,0]]
 jBH_o=[[1,2,3,4,0,0,0,0,2],[1,2,3,4,0,0,0,0,2],[1,2,3,4,0,0,0,0,2]]
 jBH_ee=[1.0,1.0,1.0]
 jBH_en=[1.0,1.0,1.0]

 set_BH_J_keywords 
}


function set_BH_J_keywords (){
 qp set jastrow j2e_type $j2e_type # set the jastrow two-e type 
 qp set jastrow env_type $env_type 
 qp set jastrow j1e_type $j1e_type
 qp set jastrow jBH_size $jBH_size # set the number of parameters in Boys-Handy jastrow 
 qp set jastrow jBH_c "$jBH_c" # set the parameters which are lists for Boys-Handy 
 qp set jastrow jBH_m "$jBH_m" # 
 qp set jastrow jBH_n "$jBH_n" #
 qp set jastrow jBH_o "$jBH_o" # 
 qp set jastrow jBH_ee $jBH_ee # 
 qp set jastrow jBH_en $jBH_en #
 qp set tc_keywords tc_integ_type $tc_integ_type # set the analytical or numerical integrals 
 qp set tc_keywords thresh_tcscf $thresh_tcscf 
 qp set tc_keywords io_tc_integ $io_tc_integ # set the io 
 rm ${EZFIO_FILE}/tc_bi_ortho/psi_*
}

function run_ground_state (){
 qp set tc_keywords minimize_lr_angles True
 qp run tc_scf | tee ${EZFIO_FILE}.tc_scf.out
 qp set_frozen_core
 qp set determinants n_det_max 1e6
 qp run fci_tc_bi_ortho | tee ${EZFIO_FILE}.fci_tc_bi.out 
}

function run_excited_state (){
 qp set determinants n_states $nstates
 qp run cis | tee ${EZFIO_FILE}.cis.out 
 rm ${EZFIO_FILE}/tc_bi_ortho/psi_*
 qp run tc_bi_ortho | tee ${EZFIO_FILE}.tc_cis_nst_${nstates}.out
 qp set determinants read_wf True 
 qp run fci_tc_bi_ortho | tee ${EZFIO_FILE}.fci_tc_bi_nst_${nstates}.out 
 
}


## BH9 calculations 
ezfio=${system}_${charge}_${basis}_${j2e_type}
create_ezfio 
BH_9 
run_ground_state 
run_excited_state
