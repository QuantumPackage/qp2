#!/usr/bin/env bats

source $QP_ROOT/tests/bats/common.bats.sh
source $QP_ROOT/quantum_package.rc


function run() {
  thresh=1.e-7
  functional=$2
  qp set_file $1
  qp edit --check
  qp set determinants n_states 1
  qp set scf_utils thresh_scf 1.e-10
  qp set dft_keywords exchange_functional $functional
  qp set dft_keywords correlation_functional $functional
  qp set hamiltonian mu_erf 0.5 
  qp set becke_numerical_grid grid_type_sgn 1 
  qp_reset --mos $1 
  qp run rs_ks_scf 
  energy="$(ezfio get kohn_sham_rs energy)"
  eq $energy $3 $thresh
}

@test "H3COH" {
  run h3coh.ezfio sr_pbe -115.502382305049
}

@test "HCN" {
  run hcn.ezfio sr_pbe -93.26674673761752
}

@test "N2" {
  run n2.ezfio  sr_pbe -109.404692225719
}

@test "SiH2_3B1" {
  run sih2_3b1.ezfio sr_lda -289.439873304975
}


