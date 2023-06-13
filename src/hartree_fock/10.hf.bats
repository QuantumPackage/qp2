#!/usr/bin/env bats

source $QP_ROOT/tests/bats/common.bats.sh
source $QP_ROOT/quantum_package.rc


function run() {
  thresh=1.e-5
  test_exe scf || skip
  qp set_file $1
  qp edit --check
  qp reset --mos
  qp set scf_utils n_it_scf_max 50
  qp run scf
#  qp set_frozen_core
  energy="$(ezfio get hartree_fock energy)"
  eq $energy $2 $thresh
}


function run_pt_charges() {
  thresh=1.e-5
  cp ${QP_ROOT}/src/nuclei/write_pt_charges.py . 
  cat > hcn.xyz << EOF
3
HCN molecule
C    0.0    0.0    0.0
H    0.0    0.0    1.064
N    0.0    0.0    -1.156
EOF

cat > hcn_charges.xyz << EOF
0.5    2.0   0.0   0.0
0.5   -2.0   0.0   0.0
EOF

EZFIO=hcn_pt_charges
rm -rf $EZFIO
qp create_ezfio -b def2-svp hcn.xyz -o $EZFIO
qp run scf
mv hcn_charges.xyz ${EZFIO}_point_charges.xyz
python write_pt_charges.py ${EZFIO}
qp set nuclei point_charges True
qp run scf | tee ${EZFIO}.pt_charges.out
  energy="$(ezfio get hartree_fock energy)"
good=-92.79920682236470
  eq $energy $good $thresh
rm -rf $EZFIO
}

@test "H2_1" { # 1s
  run h2_1.ezfio -1.005924963288527
}

@test "H2_3" { # 1s
  run h2_3.ezfio -0.9591011604845440
}

@test "H3_2" { # 1s
  run h3_2.ezfio -1.558273529860488
}

@test "H3_4" { # 1s
  run h3_4.ezfio -1.0158684760025190
}

@test "H4_1" { # 1s
  run h4_1.ezfio -1.932022805374405
}

@test "H4_3" { # 1s
  run h4_3.ezfio -1.8948449927787350
}

@test "H4_5" { # 1s
  run h4_5.ezfio -1.2408338805496990
}

@test "point charges" { 
 run_pt_charges
}

@test "HCN" { # 7.792500 8.51926s
  run hcn.ezfio -92.88717500035233
}



@test "B-B" { # 3s
  run b2_stretched.ezfio -48.9950585434279
}

@test "LiF" { # 3 s
  run lif.ezfio -106.9801081911955
}

@test "Be" { # 3 s
  run be.ezfio -14.57287346825270
}

@test "F" { # 3 s
  run f.ezfio -99.40093527229389
}


@test "SiH2_3B1" { # 0.539000  1.51094s
  run sih2_3b1.ezfio -289.9654718453571
}

@test "SO" { # 0.539000  5.70403s
  run so.ezfio -25.7175272905296
}

@test "HCO" { # 0.636700  1.55279s
  run hco.ezfio -113.1841002944744
}

@test "HBO" { # 0.805600 1.4543s
  run  hbo.ezfio  -100.018582259096
}

@test "H2S" { # 1.655600 4.21402s
  run h2s.ezfio -398.6944130421982
}

@test "H3COH" { # 1.751000 2.13527s
  run h3coh.ezfio  -114.9865030596373
}

@test "H2O" { # 1.811100 1.84387s
  run  h2o.ezfio  -0.760270218692179E+02
}

@test "H2O2" { # 2.217000 8.50267s
  run h2o2.ezfio -150.7806608469964
}

@test "ClF" { # 2.797000 6.92182s
  run clf.ezfio -558.8442570663570
}

@test "CO2" { # 2.811100 7.0952s
  run co2.ezfio -187.6507108861204
}

@test "N2H4" { # 4.054600 10.0174s
  run n2h4.ezfio -111.1799916679009
}

@test "ClO" { # 4.927400 7.63417s
  run clo.ezfio -534.2496714154559
}

@test "F2" { # 5.070800 12.6665s
  run f2.ezfio  -198.7643578233773
}

@test "CH4" { # 5.994000 13.3753s
  run ch4.ezfio -40.19961807784367
}


@test "N2" { # 8.648100 13.754s
  run n2.ezfio  -108.9834897852979
}

@test "DHNO" { # 12.856700 16.5908s
  run  dhno.ezfio  -130.427877782432
}

@test "NH3" { # 13.632200 34.7981s
  run nh3.ezfio -56.21783428976567
}

@test "C2H2" { # 19.599000 37.7923s
  run c2h2.ezfio -12.12144044853196
}


@test "SiH3" { # 20.316100 54.0861s
  [[ -n $TRAVIS ]] && skip
  run sih3.ezfio -5.455400439077580
}

@test "OH" { # 32.042200 1.36478m
  [[ -n $TRAVIS ]] && skip
  run oh.ezfio -75.42025413469165
}

@test "[Cu(NH3)4]2+" { # 59.610100 4.18766m
  [[ -n $TRAVIS ]] && skip
  qp set_file cu_nh3_4_2plus.ezfio
  qp set scf_utils thresh_scf 1.e-10
  run  cu_nh3_4_2plus.ezfio -1862.97590358903
}

@test "SO2" { # 71.894900  3.22567m
  [[ -n $TRAVIS ]] && skip
  run so2.ezfio -41.55800401346361
}

