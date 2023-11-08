#!/usr/bin/env bats

source $QP_ROOT/tests/bats/common.bats.sh
source $QP_ROOT/quantum_package.rc


function run_O() {
  qp set_file O_tc_scf 
  FILE=O_tc_scf/tc_bi_ortho/psi_l_coef_bi_ortho.gz
  if test -f "$FILE"; then
    rm O_tc_scf/tc_bi_ortho/psi*
  fi
  qp set determinants n_det_max 20000
  file=${EZFIO_FILE}.fci_tc_bi_ortho.out
  qp run fci_tc_bi_ortho | tee $file
  eref=-74.971188861115309
  energy="$(grep 'E(before) +rPT2   =' $file | tail -1 | cut -d '=' -f 2)"
  eq $energy $eref 1e-4
}


@test "O" {
 run_O 
}


