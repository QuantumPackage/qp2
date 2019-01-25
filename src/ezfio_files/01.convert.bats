#!/usr/bin/env bats

source $QP_ROOT/tests/bats/common.bats.sh
source $QP_ROOT/quantum_package.rc

function run {
  local INPUT=$1
  local EZ=$2
  cp ${QP_ROOT}/tests/input/$INPUT .
  qp convert_output_to_ezfio $INPUT -o $EZ
  qp set_file $EZ
  qp edit --check
  qp set scf_utils thresh_scf 1.e-12
}

@test "HBO GAMESS" { # 1.14107s
  run hbo.gms.out hbo.ezfio
}

@test "H2O G09" { # 1.19319s
  run h2o.log h2o.ezfio
}

@test "[Cu(NH3)4]2+ GAMESS" { # 1.38541s
  run cu_nh3_4_2plus.gms.out  cu_nh3_4_2plus.ezfio
  qp set scf_utils thresh_scf 1.e-10
}
