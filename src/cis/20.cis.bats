#!/usr/bin/env bats

source $QP_ROOT/tests/bats/common.bats.sh
source $QP_ROOT/quantum_package.rc

function run() {
  thresh=1.e-5
  test_exe cis || skip
  qp set_file $1
  qp edit --check
  qp set determinants n_states  3
  qp set davidson threshold_davidson 1.e-12
  qp set mo_two_e_ints io_mo_two_e_integrals Write
  qp set_frozen_core 
  qp run cis 
  energy1="$(qp get cis energy | tr '[]' ' ' | cut -d ',' -f 1)"
  energy2="$(qp get cis energy | tr '[]' ' ' | cut -d ',' -f 2)"
  energy3="$(qp get cis energy | tr '[]' ' ' | cut -d ',' -f 3)"
  eq $energy1 $2 $thresh
  eq $energy2 $3 $thresh
  eq $energy3 $4 $thresh
}

@test "B-B" { #  2.0s
  run  b2_stretched.ezfio  -48.995058575280950  -48.974653655601145  -48.974653655601031  

}

@test "SiH2_3B1" { # 1.23281s  1.24958s
  run sih2_3b1.ezfio -289.969297318489 -289.766898643192 -289.737521023380
}

@test "HBO" { # 1.31404s 1.50611s
  run  hbo.ezfio  -100.018582259097  -99.7127500068768  -99.6982683641297
}

@test "HCO" { # 1.33255s 1.32164s
  run hco.ezfio -113.191875527032  -113.10055658881888 -112.99411710845908
}

@test "H2O" { # 1.39318s 1.65722s
  run  h2o.ezfio  -76.02702187043107 -75.6854407466997  -75.61967556334928
}

@test "H3COH" { # 1.40257s 1.54075s
  run h3coh.ezfio  -114.986503059639 -114.649121836046 -114.578365912794
}

@test "H2S" { # 1.44228s 1.61519s
  run h2s.ezfio -398.694413042222 -398.447164835271 -398.412784774083
}



@test "ClF" { # 1.63289s 1.8911s
  [[ -n $TRAVIS ]] && skip
  run clf.ezfio -558.844257066356 -558.664418728406 -558.664418728405
}

@test "ClO" { # 1.65582s 2.06465s
  [[ -n $TRAVIS ]] && skip
  run clo.ezfio -534.263560525680 -534.256601571199 -534.062020844428
}

@test "SO" { # 1.9667s 2.91234s
  [[ -n $TRAVIS ]] && skip
  run so.ezfio -25.750224071640112 -25.586278842164113 -25.582933929660470
}

@test "OH" { # 2.201s 2.65573s
  [[ -n $TRAVIS ]] && skip
  run oh.ezfio -75.4314648243896 -75.4254639668256 -75.2707675632313
}

@test "H2O2" { # 2.27079s 3.07875s
  [[ -n $TRAVIS ]] && skip
  run h2o2.ezfio -150.780660847001 -150.546208866263 -150.483274551717
}

@test "CO2" { # 2.86928s 3.47516s
  [[ -n $TRAVIS ]] && skip
  run co2.ezfio -187.650710886151 -187.300746249524 -187.291641359067
}

@test "C2H2" { # 3.00666s 5.40252s
  [[ -n $TRAVIS ]] && skip
  run c2h2.ezfio -12.1214401949634 -11.8824874421211 -11.8682310791620
}

@test "HCN" { # 4.21678s 6.53796s
  [[ -n $TRAVIS ]] && skip
  run hcn.ezfio -92.8871750003811 -92.6250263755063 -92.6089719143274
}

@test "N2H4" { # 4.81968s 7.439s
  [[ -n $TRAVIS ]] && skip
  run n2h4.ezfio -111.179991667947 -110.894116344878 -110.855788839735
}

@test "SiH3" { # 5.72801s 14.6381s
  [[ -n $TRAVIS ]] && skip
  run sih3.ezfio -5.45916474249436 -5.23512810272682 -5.23512806272007
}

@test "N2" { # 6.11313s 10.555s
  [[ -n $TRAVIS ]] && skip
  run n2.ezfio  -108.983489785305 -108.670192549322 -108.649653940027
}

@test "DHNO" { # 6.42976s 12.9899s
  [[ -n $TRAVIS ]] && skip
  run dhno.ezfio -130.4472288472718 -130.3571808164850 -130.2196257046987
}

@test "CH4" { # 6.4969s 10.9157s
  [[ -n $TRAVIS ]] && skip
  run ch4.ezfio -40.1996180778616 -39.7936150141939 -39.7936150141734
}

@test "F2" { # 10.4758s 13.6221s
  [[ -n $TRAVIS ]] && skip
  run f2.ezfio  -198.764357823385 -198.575548537096 -198.575548537096
}

@test "NH3" { # 14.2066s 29.6974s
  [[ -n ${TRAVIS} ]] && skip
  run nh3.ezfio -56.21783428981829 -55.91997684191139 -55.84753645754046
}

@test "[Cu(NH3)4]2+" { # 29.7711s  3.45478m
  [[ -n ${TRAVIS} ]] && skip
  run  cu_nh3_4_2plus.ezfio -1862.97958885180 -1862.92457657404 -1862.91134959451

}

@test "SO2" { # 32.092s 1.47785m
  [[ -n ${TRAVIS} ]] && skip
  run so2.ezfio -41.5580019075645  -41.38232986913486 -41.35512503680323
}



