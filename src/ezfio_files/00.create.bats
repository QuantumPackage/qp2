#!/usr/bin/env bats

source $QP_ROOT/tests/bats/common.bats.sh
source $QP_ROOT/quantum_package.rc

function run {
  local INPUT=$1
  local EZ=${INPUT/.xyz/.ezfio}
  local EZ=${EZ/.zmt/.ezfio}
  local MULT=$2
  local CHARGE=$3
  local BASIS=$4
  if [[ -n $5 ]] ; then 
    local PSEUDO="-p $5"
  fi
  cp ${QP_ROOT}/tests/input/$INPUT .
  rm -rf $EZ
  qp create_ezfio \
     $INPUT --basis="$BASIS" -m $MULT -c $CHARGE $PSEUDO -o $EZ
  qp edit --check
  qp set scf_utils thresh_scf 1.e-12
  qp set ao_two_e_ints io_ao_two_e_integrals "Write" 
  qp set mo_two_e_ints io_mo_two_e_integrals "Write" 
}

@test "H2_1" {
  run h2_1.xyz 1 0 cc-pvdz
}

@test "H2_3" {
  run h2_3.xyz 3 0 cc-pvdz
}

@test "H3_2" {
  run h3_2.xyz 2 0 cc-pvdz
}

@test "H3_4" {
  run h3_4.xyz 4 0 cc-pvdz
}

@test "H4_1" {
  run h4_1.xyz 1 0 cc-pvdz
}

@test "H4_3" {
  run h4_3.xyz 3 0 cc-pvdz
}

@test "H4_5" {
  run h4_5.xyz 5 0 cc-pvdz
}


@test "B-B" {
  qp set_file b2_stretched.ezfio
  run b2_stretched.zmt 1 0 6-31g
}

@test "C2H2" {
  run c2h2.xyz 1 0 cc-pvdz_ecp_bfd bfd
}

@test "ClO" {
  run clo.xyz 2 0 cc-pvdz
}

@test "DHNO" {
  run dhno.xyz 2 0 "chipman-dzp"
}

@test "H3COH" {
  run h3coh.xyz 1 0 6-31g
}

@test "HCN" {
  run hcn.xyz 1 0 aug-cc-pvdz
}

@test "LiF" {
  run lif.xyz 1 0 cc-pvtz
}

@test "F" {
  run f.xyz 2 0 cc-pvtz
}

@test "Be" {
  run be.xyz 1 0 cc-pvtz
}

@test "N2" {
  run n2.xyz 1 0 cc-pvtz
}

@test "SiH2_3B1" {
  run sih2_3b1.xyz 3 0 6-31g
}

@test "SO" {
  run so.xyz 3 0 cc-pvdz_ecp_bfd bfd
}

@test "CH4" {
  run ch4.xyz 1 0 aug-cc-pvdz
}

@test "CO2" {
  run co2.xyz 1 0 cc-pvdz
}

@test "F2" {
  run f2.zmt 1 0 def2-tzvp
}

@test "HCO" {
  run hco.xyz 2 0 6-31+g
}

@test "NH3" {
  run nh3.xyz 1 0 cc-pvtz
}

@test "SiH3" {
  run sih3.xyz 2 0 cc-pvdz_ecp_bfd bfd
}

@test "ClF" {
  run clf.xyz 1 0 cc-pvdz
}

@test "H2O2" {
  run h2o2.zmt 1 0 cc-pvdz
}

@test "H2S" {
  run h2s.xyz 1 0 cc-pvdz
}

@test "N2H4" {
  run n2h4.zmt 1 0 cc-pvdz
}

@test "OH" {
  run oh.xyz 2 0 aug-ano-pvdz_roos
}

@test "SO2" {
  run so2.xyz 1 0 cc-pvtz_ecp_bfd bfd
}

