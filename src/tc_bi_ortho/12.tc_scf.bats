#!/usr/bin/env bats

source $QP_ROOT/tests/bats/common.bats.sh
source $QP_ROOT/quantum_package.rc


function run_Ne() {
  qp set_file Ne_tc_scf 
  qp run cisd 
  qp run tc_bi_ortho | tee Ne.ezfio.cisd_tc_bi_ortho.out  
  eref=-128.77020441279302
  energy="$(grep "eigval_right_tc_bi_orth =" Ne.ezfio.cisd_tc_bi_ortho.out)"
  eq $energy $eref 1e-6
}


@test "Ne" {
 run_Ne 
}

