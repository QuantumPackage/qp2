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
  qp set davidson threshold_davidson 1.e-10
  qp set davidson n_states_diag 8
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
  qp set davidson threshold_davidson 1.e-10
  qp set davidson n_states_diag 1
  qp run fci
  energy1="$(ezfio get fci energy_pt2 | tr '[]' ' ' | cut -d ',' -f 1)"
  eq $energy1 $1 $thresh
}

@test "B-B" {
  qp set_file b2_stretched.ezfio
  qp set determinants n_det_max 10000
  qp set_frozen_core
  run_stoch -49.14103054419 3.e-4 10000
}

@test "F2" { # 4.07m
  [[ -n $TRAVIS ]] && skip
  qp set_file f2.ezfio
  qp set_frozen_core
  run_stoch -199.304922384814 3.e-4  100000
}

@test "NH3" { # 10.6657s
  qp set_file nh3.ezfio
  qp set_mo_class --core="[1-4]" --act="[5-72]"
  run -56.244753429144986  3.e-4  100000
}

@test "DHNO" { # 11.4721s
  qp set_file dhno.ezfio
  qp set_mo_class --core="[1-7]" --act="[8-64]"
  run -130.459020029816 3.e-4  100000
}

@test "HCO" { # 12.2868s
  qp set_file hco.ezfio
  run -113.389297812482 6.e-4  100000
}

@test "H2O2" { # 12.9214s
  qp set_file h2o2.ezfio
  qp set_mo_class --core="[1-2]" --act="[3-24]" --del="[25-38]"
  run -151.00467 1.e-4  100000
}

@test "HBO" { # 13.3144s
  [[ -n $TRAVIS ]] && skip
  qp set_file hbo.ezfio
  run -100.212560384678 1.e-3  100000
}

@test "H2O" { # 11.3727s
  [[ -n $TRAVIS ]] && skip
  qp set_file h2o.ezfio
  run -76.2361605151999 3.e-4  100000
}

@test "ClO" { # 13.3755s
  [[ -n $TRAVIS ]] && skip
  qp set_file clo.ezfio
  run -534.545616787223 3.e-4  100000
}

@test "SO" { # 13.4952s
  [[ -n $TRAVIS ]] && skip
  qp set_file so.ezfio
  run -26.0096209515081 1.e-3  100000
}

@test "H2S" { # 13.6745s
  [[ -n $TRAVIS ]] && skip
  qp set_file h2s.ezfio
  run -398.859168655255 3.e-4  100000
}

@test "OH" { # 13.865s
  [[ -n $TRAVIS ]] && skip
  qp set_file oh.ezfio
  run -75.6121856748294 3.e-4   100000
}

@test "SiH2_3B1" { # 13.938ss
  [[ -n $TRAVIS ]] && skip
  qp set_file sih2_3b1.ezfio
  run -290.0175411299477 3.e-4  100000
}

@test "H3COH" { # 14.7299s
  [[ -n $TRAVIS ]] && skip
  qp set_file h3coh.ezfio
  run -115.205191406072 3.e-4  100000
}

@test "SiH3" { # 15.99s
  [[ -n $TRAVIS ]] && skip
  qp set_file sih3.ezfio
  run -5.57241217753818 3.e-4   100000
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
  run -559.169313755572 3.e-4  100000
}

@test "SO2" { # 17.5645s
  [[ -n $TRAVIS ]] && skip
  qp set_file so2.ezfio
  qp set_mo_class --core="[1-8]" --act="[9-87]"
  run -41.5746738713298 3.e-4  100000
}

@test "C2H2" { # 17.6827s
  [[ -n $TRAVIS ]] && skip
  qp set_file c2h2.ezfio
  qp set_mo_class --act="[1-30]" --del="[31-36]"
  run -12.3685464085969 3.e-4  100000
}

@test "N2" { # 18.0198s
  [[ -n $TRAVIS ]] && skip
  qp set_file n2.ezfio
  qp set_mo_class --core="[1,2]" --act="[3-40]" --del="[41-60]"
  run -109.28681540699360 3.e-4  100000
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
  run -187.968547952413 3.e-4  100000
}


@test "[Cu(NH3)4]2+" { # 25.0417s
  [[ -n $TRAVIS ]] && skip
  qp set_file cu_nh3_4_2plus.ezfio
  qp set_mo_class --core="[1-24]" --act="[25-45]" --del="[46-87]"
  run -1862.9869374387192  3.e-04  100000
}

@test "HCN" { # 20.3273s
  [[ -n $TRAVIS ]] && skip
  qp set_file hcn.ezfio
  qp set_mo_class --core="[1,2]" --act="[3-40]" --del="[41-55]"
  run -93.0771143355433 3.e-4  100000
}

