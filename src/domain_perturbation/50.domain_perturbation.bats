#!/usr/bin/env bats

source $QP_ROOT/tests/bats/common.bats.sh
source $QP_ROOT/quantum_package.rc

# Deterministic evaluation, compared to reference values
# The tests need converged orbitals: the SCF is run first. The
# perturbative energies depend on the orbitals (not only on the HF energy),
# for example on the rotations of degenerate orbitals, hence the tolerances.
function scf() {
  test_exe scf || skip
  qp set_file $1
  qp reset --mos
  qp set scf_utils n_it_scf_max 50
  qp run scf
}

function run() {
  thresh=1.e-5
  test_exe domain_pert || skip
  qp edit --check
  qp set domain_perturbation domain_chebyshev false
  qp set domain_perturbation n_pert $2
  qp set domain_perturbation n_det_domain $3
  qp set domain_perturbation n_walkers 0
  qp run domain_pert
  energies="$(ezfio get domain_perturbation domain_energy | tr '[],' '   ')"
  i=1
  for e in $1 ; do
    energy="$(echo $energies | cut -d ' ' -f $i)"
    eq $energy $e $thresh
    i=$((i+1))
  done
}

# Chebyshev series, compared to the full CI energies
function run_cheb() {
  thresh=1.e-5
  test_exe domain_pert || skip
  qp edit --check
  qp set domain_perturbation domain_chebyshev true
  qp set domain_perturbation domain_cheb_auto true
  qp set domain_perturbation n_pert 10
  qp set domain_perturbation n_det_domain 1
  qp set domain_perturbation n_walkers 0
  qp run domain_pert
  qp set domain_perturbation domain_chebyshev false
  energies="$(ezfio get domain_perturbation domain_energy | tr '[],' '   ')"
  i=1
  for e in $1 ; do
    energy="$(echo $energies | cut -d ' ' -f $i)"
    eq $energy $e $thresh
    i=$((i+1))
  done
}

# Stochastic evaluation, compared to the exact one
function run_stoch() {
  test_exe domain_pert || skip
  qp set domain_perturbation domain_chebyshev false
  qp set domain_perturbation n_pert $1
  qp set domain_perturbation n_det_domain $2
  qp set domain_perturbation exact_coupling $3
  qp set domain_perturbation domain_heat_bath $5
  qp set domain_perturbation n_walkers 5000
  qp set domain_perturbation n_replicas 20
  qp set domain_perturbation domain_seed 1
  qp run domain_pert
  energy="$(ezfio get domain_perturbation domain_energy \
    | tr -d '[] ' | cut -d ',' -f 1)"
  error="$(ezfio get domain_perturbation domain_energy_error \
    | tr -d '[] ' | cut -d ',' -f 1)"
  eq "$energy" "$4" "$(LC_ALL=C awk -v e="$error" 'BEGIN { print 6*e }')"
}

@test "H3_2 HF" {
  scf h3_2.ezfio
  qp reset -d
  run "-1.610273153311" 4 1
  run "-1.609959380386" 4 4
}

@test "H3_2 CIS 3 states" {
  scf h3_2.ezfio
  qp reset -d
  qp set determinants n_states 3
  qp run cis
  run "-1.610008725602 -1.311815765541 -1.206677373119" 5 1
  run "-1.610022696198 -1.311852657934 -1.206705811156" 5 6
  run_cheb "-1.610031316672 -1.311908009016"
  qp set determinants n_states 1
}

@test "H4_1 HF" {
  scf h4_1.ezfio
  qp reset -d
  run "-2.017046423901" 4 1
  run "-2.016898059074" 4 4
  run_stoch 4 4 true  "-2.016898059074" false
  run_stoch 4 4 true  "-2.016898059074" true
  run_stoch 4 1 false "-2.017046423901" true
}
