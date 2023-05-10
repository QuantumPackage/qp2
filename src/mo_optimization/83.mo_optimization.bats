#!/usr/bin/env bats

source $QP_ROOT/tests/bats/common.bats.sh
source $QP_ROOT/quantum_package.rc


function run() {
  thresh=2e-3
  test_exe scf || skip
  qp set_file $1
  qp edit --check
  qp reset -a
  qp run scf
  qp set_frozen_core
  qp set determinants n_states 2
  qp set determinants read_wf true
  qp set mo_two_e_ints io_mo_two_e_integrals None
  file="$(echo $1 | sed 's/.ezfio//g')"
  qp run cis
  qp run debug_gradient_list_opt > $file.debug_g.out
  err3="$(grep 'Max error:' $file.debug_g.out | awk '{print $3}')"
  qp run debug_hessian_list_opt > $file.debug_h1.out
  err1="$(grep 'Max error:' $file.debug_h1.out | awk '{print $3}')"
  qp run orb_opt > $file.opt1.out
  energy1="$(grep 'State average energy:' $file.opt1.out | tail -n 1 | awk '{print $4}')"
  qp set orbital_optimization optimization_method diag
  qp reset -d 
  qp run scf
  qp run cis
  qp run debug_hessian_list_opt > $file.debug_h2.out
  err2="$(grep 'Max error_H:' $file.debug_h2.out | awk '{print $3}')"
  qp run orb_opt > $file.opt2.out
  energy2="$(grep 'State average energy:' $file.opt2.out | tail -n 1 | awk '{print $4}')"
  qp set orbital_optimization optimization_method full
  qp reset -d
  qp run scf
  eq $energy1 $2 $thresh
  eq $energy2 $3 $thresh
  eq $err1 0.0 1e-12
  eq $err2 0.0 1e-12
  eq $err3 0.0 1e-12
}

@test "b2_stretched" {
run  b2_stretched.ezfio -48.9852901484277 -48.9852937541510
}

@test "h2o" {
run  h2o.ezfio -75.9025622449206 -75.8691844585879
}

@test "h2s" {
run  h2s.ezfio -398.576255809878 -398.574145943928
}

@test "hbo" {
run  hbo.ezfio -99.9234823022109 -99.9234763597840
}

@test "hco" {
run  hco.ezfio -113.204915552241 -113.204905207050
}
