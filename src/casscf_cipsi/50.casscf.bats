#!/usr/bin/env bats

source $QP_ROOT/tests/bats/common.bats.sh
source $QP_ROOT/quantum_package.rc


function run_stoch() {
  thresh=$2
  test_exe casscf || skip
  qp set perturbation do_pt2 True
  qp set determinants n_det_max $3
  qp set davidson_keywords threshold_davidson 1.e-10
  qp set davidson_keywords n_states_diag 4
  qp run casscf | tee casscf.out 
  energy1="$(ezfio get casscf energy_pt2 | tr '[]' ' ' | cut -d ',' -f 1)"
  eq $energy1 $1 $thresh
}

@test "F2" { # 18.0198s
  rm -rf f2_casscf
  qp_create_ezfio -b aug-cc-pvdz ../input/f2.zmt  -o f2_casscf
  qp set_file f2_casscf
  qp run scf 
  qp set_mo_class --core="[1-6,8-9]" --act="[7,10]" --virt="[11-46]"
  run_stoch  -198.773366970 1.e-4  100000
}

@test "N2" { # 18.0198s
  rm -rf n2_casscf
  qp_create_ezfio -b aug-cc-pvdz ../input/n2.xyz  -o n2_casscf
  qp set_file n2_casscf
  qp run scf 
  qp set_mo_class --core="[1-4]" --act="[5-10]" --virt="[11-46]"
  run_stoch  -109.0961643162 1.e-4  100000
}

@test "N2_stretched" {
  rm -rf n2_stretched_casscf
  qp_create_ezfio -b aug-cc-pvdz -m 7 ../input/n2_stretched.xyz  -o n2_stretched_casscf
  qp set_file n2_stretched_casscf 
  qp run scf  | tee scf.out 
  qp set_mo_class --core="[1-4]" --act="[5-10]" --virt="[11-46]"
  qp set electrons elec_alpha_num 7 
  qp set electrons elec_beta_num 7 
  run_stoch  -108.7860471300 1.e-4  100000
#   

}

