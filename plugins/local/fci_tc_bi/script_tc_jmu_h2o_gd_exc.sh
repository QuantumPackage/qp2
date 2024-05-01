#!/bin/bash
source ~/qp2/quantum_package.rc

## Define the system/basis/charge/mult and genric keywords 
system=H2O
xyz=${system}.xyz
basis=6-31g
mult=1
charge=0
j2e_type=Mu
thresh_tcscf=1e-10
io_tc_integ="Write"
nstates=4
nol_standard=False
tc_integ_type=numeric # can be changed for semi-analytic

if (( $nol_standard == "False" ))
then
 three_body_h_tc=True
else
 three_body_h_tc=False
fi


##################### Function to create the EZFIO 
function create_ezfio (){
 qp create_ezfio -b $basis -m $mult -c $charge $xyz -o $ezfio
 qp run scf | tee ${EZFIO_FILE}.scf.out 
}

function set_env_j_keywords (){
 
 qp set hamiltonian mu_erf 0.87
 qp set jastrow env_type Sum_Gauss 
 qp set jastrow env_coef "${coef}"
 qp set tc_keywords tc_integ_type $tc_integ_type
 qp set jastrow j1e_type $j1e_type
 qp set jastrow j2e_type $j2e_type
 qp set jastrow env_expo "${alpha}"
}

function run_ground_state (){
 qp set tc_keywords minimize_lr_angles True
 qp run tc_scf | tee ${EZFIO_FILE}.tc_scf.out
 qp set_frozen_core
 qp set determinants n_det_max 1e6
 qp set perturbation pt2_max 0.001 
 qp set tc_keywords nol_standard $nol_standard
 qp set tc_keywords three_body_h_tc $three_body_h_tc
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


# Define J(mu) with envelope and without j1e 
j2e_type=Mu
j1e_type=None
ezfio=${system}_${charge}_${basis}_${j2e_type}_${j1e_type}
create_ezfio 
alpha=[2.0,1000.,1000.] # parameters for H2O
coef=[1.,1.,1.] # parameters for H2O
set_env_j_keywords
run_ground_state 
run_excited_state 

# Define J(mu) with envelope and with a charge Harmonizer for J1e 
j2e_type=Mu
j1e_type=Charge_Harmonizer
ezfio=${system}_${charge}_${basis}_${j2e_type}_${j1e_type}
create_ezfio 
alpha=[2.5,1000.,1000.] # parameters for H2O
coef=[1.,1.,1.] # parameters for H2O
set_env_j_keywords
run_ground_state 
run_excited_state 
