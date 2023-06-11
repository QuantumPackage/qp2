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


function run_stoch() {
  thresh=$2
  test_exe fci || skip
  qp set perturbation do_pt2 True
  qp set determinants n_det_max $3
  qp set determinants n_states  1
  qp set davidson_keywords threshold_davidson 1.e-10
  qp set davidson_keywords n_states_diag 1
  qp run fci
  energy1="$(ezfio get fci energy_pt2 | tr '[]' ' ' | cut -d ',' -f 1)"
  eq $energy1 $1 $thresh
}

@test "H2_1" { # 1s
  qp set_file h2_1.ezfio
  qp set perturbation pt2_max 0.
  run_stoch -1.06415255 1.e-8 10000
}

@test "H2_3" { # 1s
  qp set_file h2_3.ezfio
  qp set perturbation pt2_max 0.
  run_stoch -0.96029881 1.e-8 10000
}

@test "H3_2" { # 3s
  qp set_file h3_2.ezfio
  qp set perturbation pt2_max 0.
  run_stoch -1.61003132 1.e-8 10000
}

@test "H3_4" { # 2s
  qp set_file h3_4.ezfio
  qp set perturbation pt2_max 0.
  run_stoch -1.02434843 1.e-8 10000
}

@test "H4_1" { # 13s
  qp set_file h4_1.ezfio
  qp set perturbation pt2_max 0.
  run_stoch -2.01675062 1.e-8 10000
}

@test "H4_3" { # 10s
  qp set_file h4_3.ezfio
  qp set perturbation pt2_max 0.
  run_stoch -1.95927626 1.e-8 10000
}

@test "H4_5" { # 3s
  qp set_file h4_5.ezfio
  qp set perturbation pt2_max 0.
  run_stoch -1.25852765 1.e-8 10000
}

@test "B-B" { # 10s
  qp set_file b2_stretched.ezfio
  qp set determinants n_det_max 10000
  qp set_frozen_core
  run_stoch -49.14103054419 3.e-4 10000
}

@test "NH3" { # 8s
  qp set_file nh3.ezfio
  qp set_mo_class --core="[1-4]" --act="[5-72]"
  run -56.244753429144986  3.e-4  100000
}

@test "DHNO" { # 8s
  qp set_file dhno.ezfio
  qp set_mo_class --core="[1-7]" --act="[8-64]"
  run -130.466208113547 3.e-4  100000
}

@test "HCO" { # 32s
  qp set_file hco.ezfio
  run -113.395751656985 1.e-3  100000
}

@test "H2O2" { # 21s
  qp set_file h2o2.ezfio
  qp set_mo_class --core="[1-2]" --act="[3-24]" --del="[25-38]"
  run -151.005848404095 1.e-3  100000
}

@test "HBO" { # 18s
  [[ -n $TRAVIS ]] && skip
  qp set_file hbo.ezfio
  run -100.214 1.5e-3  100000
}

@test "H2O" { # 16s
  [[ -n $TRAVIS ]] && skip
  qp set_file h2o.ezfio
  run -76.238051555276  5.e-4  100000
}

@test "ClO" { # 47s
  [[ -n $TRAVIS ]] && skip
  qp set_file clo.ezfio
  run -534.548529710256 1.e-3  100000
}

@test "SO" { # 23s
  [[ -n $TRAVIS ]] && skip
  qp set_file so.ezfio
  run -26.015 3.e-3  100000
}

@test "H2S" { # 37s
  [[ -n $TRAVIS ]] && skip
  qp set_file h2s.ezfio
  run -398.864853669111 5.e-4  100000
}

@test "OH" { # 12s
  [[ -n $TRAVIS ]] && skip
  qp set_file oh.ezfio
  run -75.615 1.5e-3   100000
}

@test "SiH2_3B1" { # 10s
  [[ -n $TRAVIS ]] && skip
  qp set_file sih2_3b1.ezfio
  run -290.0206626734517 3.e-4  100000
}

@test "H3COH" { # 33s
  [[ -n $TRAVIS ]] && skip
  qp set_file h3coh.ezfio
  run -115.206784386204 1.e-3  100000
}

@test "SiH3" { # 15s
  [[ -n $TRAVIS ]] && skip
  qp set_file sih3.ezfio
  run -5.572 1.e-3   100000
}

@test "CH4" { # 16.1612s
  [[ -n $TRAVIS ]] && skip
  qp set_file ch4.ezfio
  qp set_mo_class --core="[1]" --act="[2-30]" --del="[31-59]"
  run -40.2409678239136 3.e-4  100000
}

@test "ClF" { # 16.8864s
  [[ -n $TRAVIS ]] && skip
  qp set_file clf.ezfio
  run -559.174371468224 1.5e-3  100000
}

@test "SO2" { # 17.5645s
  [[ -n $TRAVIS ]] && skip
  qp set_file so2.ezfio
  qp set_mo_class --core="[1-8]" --act="[9-87]"
  run -41.5746738713298 1.5e-3  100000
}

@test "C2H2" { # 17.6827s
  [[ -n $TRAVIS ]] && skip
  qp set_file c2h2.ezfio
  qp set_mo_class --act="[1-30]" --del="[31-36]"
  run -12.367 3.e-3  100000
}

@test "N2" { # 18.0198s
  [[ -n $TRAVIS ]] && skip
  qp set_file n2.ezfio
  qp set_mo_class --core="[1,2]" --act="[3-40]" --del="[41-60]"
  run -109.288   2.e-3  100000
}

@test "N2H4" { # 18.5006s
  [[ -n $TRAVIS ]] && skip
  qp set_file n2h4.ezfio
  qp set_mo_class --core="[1-2]" --act="[3-24]" --del="[25-48]"
  run -111.367332681559 3.e-4  100000
}

@test "CO2" { # 21.1748s
  [[ -n $TRAVIS ]] && skip
  qp set_file co2.ezfio
  qp set_mo_class --core="[1,2]" --act="[3-30]" --del="[31-42]"
  run -187.970184372047 1.6e-3  100000
}

@test "[Cu(NH3)4]2+" { # 25.0417s
  [[ -n $TRAVIS ]] && skip
  qp set_file cu_nh3_4_2plus.ezfio
  qp set_mo_class --core="[1-24]" --act="[25-45]" --del="[46-87]"
  run -1862.98320066637   3.e-04  100000
}

@test "HCN" { # 20.3273s
  [[ -n $TRAVIS ]] && skip
  qp set_file hcn.ezfio
  qp set_mo_class --core="[1,2]" --act="[3-40]" --del="[41-55]"
  run -93.078 2.e-3  100000
}

@test "F2" { # 4.07m
  [[ -n $TRAVIS ]] && skip
  qp set_file f2.ezfio
  qp set_frozen_core
  run_stoch -199.304922384814 3.e-3  100000
}

