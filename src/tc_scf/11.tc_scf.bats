#!/usr/bin/env bats

source $QP_ROOT/tests/bats/common.bats.sh
source $QP_ROOT/quantum_package.rc


function run_Ne() {
  rm -rf Ne_tc_scf
  echo Ne > Ne.xyz
  qp create_ezfio -b cc-pcvdz Ne.xyz -o Ne_tc_scf
  qp run scf 
  qp set tc_keywords bi_ortho True 
  qp set tc_keywords test_cycle_tc True
  qp set ao_two_e_erf_ints mu_erf 0.87 
  qp set tc_keywords j1b_pen [1.5]
  qp set tc_keywords j1b_type 3 
  qp run tc_scf | tee ${EZFIO_FILE}.tc_scf.out 
  eref=-128.552134
  energy="$(qp get tc_scf bitc_energy)"
  eq $energy $eref 1e-6
}


@test "Ne" {
 run_Ne 
}

