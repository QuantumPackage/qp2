#!/usr/bin/env bats

source $QP_ROOT/tests/bats/common.bats.sh
source $QP_ROOT/quantum_package.rc

function run() {
  thresh=$2
  test_exe fci || skip
  qp edit --check
  qp set perturbation do_pt2 False
  qp set determinants n_det_max $3
  qp set determinants n_states  1
  qp set davidson threshold_davidson 1.e-10
  qp set davidson n_states_diag 8
  qp run fci
  energy1="$(qp get fci energy | tr '[]' ' ' | cut -d ',' -f 1)"
  eq $energy1 $1 $thresh
}


function run_stoch() {
  thresh=$2
  test_exe fci || skip
  qp set perturbation do_pt2 True
  qp set perturbation pt2_relative_error 0.005
  qp set determinants n_det_max $3
  qp set determinants n_states  1
  qp set davidson threshold_davidson 1.e-10
  qp set davidson n_states_diag 1
  qp run fci
  energy1="$(qp get fci energy_extrapolated | tr '[]' ' ' | cut -d ',' -f 1)"
  eq $energy1 $1 $thresh
}


@test "B-B" {  # 0:00:10
  qp set_file b2_stretched.ezfio
  qp set determinants n_det_max 10000
  qp set_frozen_core
  run_stoch -49.14097596  0.0001  10000
}

@test "F2" { # 4.07m
  [[ -n $TRAVIS ]] && skip
  qp set_file f2.ezfio
  qp set_frozen_core
  run_stoch -199.307512211742 0.002  100000
}

@test "NH3" { # 10.6657s
  qp set_file nh3.ezfio
  qp set_mo_class --core="[1-4]" --act="[5-72]"
  run -56.24474908 1.e-5  10000
}

@test "DHNO" { # 0:00:10 
  qp set_file dhno.ezfio
  qp set_mo_class --core="[1-7]" --act="[8-64]"
  run -130.45904647  1.e-4  10000
}

@test "HCO" { # 0:01:16 
  qp set_file hco.ezfio
  run_stoch -113.41448940  2.e-3  50000
}

@test "H2O2" { # 0:01:48 
  qp set_file h2o2.ezfio
  qp set_mo_class --core="[1-2]" --act="[3-24]" --del="[25-38]"
  run_stoch -151.02437936 2.e-3  100000
}

@test "HBO" { # 0:00:46 
  [[ -n $TRAVIS ]] && skip
  qp set_file hbo.ezfio
  run_stoch -100.221198108988 2.e-3  50000
}

@test "H2O" { # 0:01:05 
  [[ -n $TRAVIS ]] && skip
  qp set_file h2o.ezfio
  run_stoch -76.241332121813  1.e-3  100000
}

@test "ClO" { # 0:03:07 
  [[ -n $TRAVIS ]] && skip
  qp set_file clo.ezfio
  run_stoch -534.573564655419 1.e-3  100000
}

@test "SO" { # 0:01:49 
  [[ -n $TRAVIS ]] && skip
  qp set_file so.ezfio
  run_stoch -26.04335528  5.e-3  100000
}

@test "H2S" { # 0:01:12
  [[ -n $TRAVIS ]] && skip
  qp set_file h2s.ezfio
  run_stoch -398.865173546866 1.e-3  50000
}

@test "OH" { # 0:00:41
  [[ -n $TRAVIS ]] && skip
  qp set_file oh.ezfio
  run_stoch -75.6193013819546 1.e-3  50000
}

@test "SiH2_3B1" { # 0:00:50
  [[ -n $TRAVIS ]] && skip
  qp set_file sih2_3b1.ezfio
  run_stoch -290.01754869  3.e-5  50000
}

@test "H3COH" { # 0:01:05
  [[ -n $TRAVIS ]] && skip
  qp set_file h3coh.ezfio
  run_stoch -115.224147057725 2.e-3  50000
}

@test "SiH3" { # 0:01:09
  [[ -n $TRAVIS ]] && skip
  qp set_file sih3.ezfio
  run_stoch -5.57812512359276 1.e-3  50000
}

@test "CH4" { # 0:02:06
  [[ -n $TRAVIS ]] && skip
  qp set_file ch4.ezfio
  qp set_mo_class --core="[1]" --act="[2-30]" --del="[31-59]"
  run_stoch -40.2419474611994 1.e-4  100000
}

@test "ClF" { # 0:01:55 
  [[ -n $TRAVIS ]] && skip
  qp set_file clf.ezfio
  run_stoch -559.20666465 1.e-2  50000
}

@test "SO2" { # 0:00:24
  [[ -n $TRAVIS ]] && skip
  qp set_file so2.ezfio
  qp set_mo_class --core="[1-8]" --act="[9-87]"
  run_stoch -41.57468756 1.e-4  50000
}

@test "C2H2" { # 0:00:57 
  [[ -n $TRAVIS ]] && skip
  qp set_file c2h2.ezfio
  qp set_mo_class --act="[1-30]" --del="[31-36]"
  run_stoch -12.3862664765532 1.e-3  50000
}

@test "N2" { # 0:01:15
  [[ -n $TRAVIS ]] && skip
  qp set_file n2.ezfio
  qp set_mo_class --core="[1,2]" --act="[3-40]" --del="[41-60]"
  run_stoch -109.311954243348 2.e-3  50000
}

@test "N2H4" { # 0:00:51
  [[ -n $TRAVIS ]] && skip
  qp set_file n2h4.ezfio
  qp set_mo_class --core="[1-2]" --act="[3-24]" --del="[25-48]"
  run_stoch -111.38119165053 1.e-3 50000
}

@test "CO2" { # 0:01:00
  [[ -n $TRAVIS ]] && skip
  qp set_file co2.ezfio
  qp set_mo_class --core="[1,2]" --act="[3-30]" --del="[31-42]"
  run_stoch -188.002190327443 2.e-3  50000
}

@test "[Cu(NH3)4]2+" { # 0:01:53
  [[ -n $TRAVIS ]] && skip
  qp set_file cu_nh3_4_2plus.ezfio
  qp set_mo_class --core="[1-24]" --act="[25-45]" --del="[46-87]"
  run_stoch -1862.98705340328 1.e-05  50000
}

@test "HCN" { # 0:01:26
  [[ -n $TRAVIS ]] && skip
  qp set_file hcn.ezfio
  qp set_mo_class --core="[1,2]" --act="[3-40]" --del="[41-55]"
  run_stoch -93.0980746734051 5.e-4  50000
}

