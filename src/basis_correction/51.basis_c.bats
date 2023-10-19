#!/usr/bin/env bats

source $QP_ROOT/tests/bats/common.bats.sh
source $QP_ROOT/quantum_package.rc

function run() {
  thresh=$2
  test_exe fci || skip
  qp edit --check
  qp set perturbation do_pt2 False
  qp set determinants n_det_max 8000
  qp set determinants n_states  1
  qp set davidson_keywords threshold_davidson 1.e-10
  qp set davidson_keywords n_states_diag 8
  qp run fci 
  energy1="$(ezfio get fci energy | tr '[]' ' ' | cut -d ',' -f 1)"
  eq $energy1 $1 $thresh
}


function run_md() {
  thresh=$2
  qp set mu_of_r mu_of_r_potential cas_ful
  file_out=${EZFIO_FILE}.basis_corr.out
  qp run basis_correction | tee $file_out 
  energy1="$(grep 'ECMD SU-PBE-OT     , state    1  =' ${file_out} | cut -d '=' -f 2)"
  eq $energy1 $1 $thresh
}

function run_sd() {
  thresh=$2
  qp set mu_of_r mu_of_r_potential hf
  qp set_frozen_core
  file_out=${EZFIO_FILE}.basis_corr.out
  qp run basis_correction | tee $file_out 
  energy1="$(grep 'ECMD PBE-UEG       , state    1  =' ${file_out} | cut -d '=' -f 2)"
  eq $energy1 $1 $thresh
}

@test "O2 CAS" {  
  qp set_file o2_cas.gms.ezfio
  qp set_mo_class -c "[1-2]" -a "[3-10]" -d "[11-46]"
  run -149.72435425 3.e-4 10000
  qp set_mo_class -c "[1-2]" -a "[3-10]" -v "[11-46]"
  run_md -0.1160222327 1.e-6 
}


@test "LiF RHF" {  
  qp set_file lif.ezfio 
  run_sd -0.0649431665 1.e-6 
}

@test "F ROHF" {  
  qp set_file f.ezfio 
  run_sd -0.0355395041 1.e-6 
}

@test "Be RHF" {  
  qp set_file be.ezfio 
  run_sd -0.0325139011 1.e-6 
}

