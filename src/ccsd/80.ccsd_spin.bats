#!/usr/bin/env bats

source $QP_ROOT/tests/bats/common.bats.sh
source $QP_ROOT/quantum_package.rc


function run() {
  thresh1=1e-6
  thresh2=1e-6
  test_exe scf || skip
  qp set_file $1
  qp edit --check
  #qp run scf
  qp set_frozen_core
  qp set utils_cc cc_par_t true
  qp set utils_cc cc_thresh_conv 1e-12
  file="$(echo $1 | sed 's/.ezfio//g')"
  qp run ccsd_spin_orb | tee $file.ccsd.out
  energy1="$(grep 'E(CCSD)' $file.ccsd.out | tail -n 1 | awk '{printf $3}')"
  energy2="$(grep 'E(T)' $file.ccsd.out | tail -n 1 | awk '{printf $3}')"
  #rm $file.ccsd.out
  eq $energy1 $2 $thresh1
  eq $energy2 $3 $thresh2
}

@test "b2_stretched" {
run  b2_stretched.ezfio -49.136487344382 -0.003497589175
}

@test "be" {
run  be.ezfio -14.623559003577 -0.000230982022
}

@test "c2h2" {
run  c2h2.ezfio -12.394008897618 -0.010790491561
}

@test "ch4" {
run  ch4.ezfio -40.390721785799 -0.004476100282
}

@test "clf" {
run  clf.ezfio -559.186562904081 -0.006577143392
}

@test "clo" {
run  clo.ezfio -534.564874409332 -0.007584571424
}

@test "co2" {
run  co2.ezfio -188.129602527766 -0.018040668885
}

@test "dhno" {
run  dhno.ezfio -130.816650109473 -0.012197331453
}

@test "f2" {
run  f2.ezfio -199.287826338097 -0.017592872692
}

@test "f" {
run  f.ezfio -99.616644511121 -0.003624525307
}

@test "h2o2" {
run  h2o2.ezfio -151.182552729963 -0.009511682086
}

@test "h2o" {
run  h2o.ezfio -76.237710276526 -0.003001800577
}

@test "h2s" {
run  h2s.ezfio -398.861214015390 -0.003300559757
}

@test "h3coh" {
run  h3coh.ezfio -115.221296424969 -0.003566171432
}

@test "hbo" {
run  hbo.ezfio -100.213539770415 -0.006851489212
}

@test "hcn" {
run  hcn.ezfio -93.190247992657 -0.013418135043
}

@test "hco" {
run  hco.ezfio -113.405413962350 -0.007973455337
}

@test "lif" {
run  lif.ezfio -107.270402903250 -0.007742969005
}

@test "n2" {
run  n2.ezfio -109.355358930472 -0.018477744342
}

@test "n2h4" {
run  n2h4.ezfio -111.556885923139 -0.009048077008
}

@test "nh3" {
run  nh3.ezfio -56.465503060954 -0.007638273755
}

@test "oh" {
run  oh.ezfio -75.614606132774 -0.004126661739
}

@test "sih2_3b1" {
run  sih2_3b1.ezfio -290.016780973072 -0.000497825874
}

@test "sih3" {
run  sih3.ezfio -5.575343504534 -0.002094123268
}

@test "so" {
run  so.ezfio -26.035945178665 -0.010594351274
}

#@test "b2_stretched" {
#run  b2_stretched.ezfio -49.136487344382 -49.139984933557
#}
#
#@test "be" {
#run  be.ezfio -14.623559003577 -14.623789985599
#}
#
#@test "c2h2" {
#run  c2h2.ezfio -12.394008897618 -12.404799389179
#}
#
#@test "ch4" {
#run  ch4.ezfio -40.390721784961 -40.395197884406
#}
#
#@test "clf" {
#run  clf.ezfio -559.186562906072 -559.193140046904
#}
#
#@test "clo" {
#run  clo.ezfio -534.564874409333 -534.572458980757
#}
#
#@test "co2" {
#run  co2.ezfio -188.129602511724 -188.147643198675
#}
#
#@test "dhno" {
#run  dhno.ezfio -130.816650109473 -130.828847440925
#}
#
#@test "f2" {
#run  f2.ezfio -199.287826338097 -199.305419210789
#}
#
#@test "f" {
#run  f.ezfio -99.616644511120 -99.620269036428
#}
#
#@test "h2o2" {
#run  h2o2.ezfio -151.182552729963 -151.192064412049
#}
#
#@test "h2o" {
#run  h2o.ezfio -76.237710276526 -76.240712077103
#}
#
#@test "h2s" {
#run  h2s.ezfio -398.861214015416 -398.864514575146
#}
#
#@test "h3coh" {
#run  h3coh.ezfio -115.221296424969 -115.224862596401
#}
#
#@test "hbo" {
#run  hbo.ezfio -100.213539770415 -100.220391259627
#}
#
#@test "hcn" {
#run  hcn.ezfio -93.190247983000 -93.203666131216
#}
#
#@test "hco" {
#run  hco.ezfio -113.405413962350 -113.413387417687
#}
#
#@test "lif" {
#run  lif.ezfio -107.270402903211 -107.278145872216
#}
#
#@test "n2" {
#run  n2.ezfio -109.355358930472 -109.373836674814
#}
#
#@test "n2h4" {
#run  n2h4.ezfio -111.556885922642 -111.565934000556
#}
#
#@test "nh3" {
#run  nh3.ezfio -56.465503060954 -56.473141334709
#}
#
#@test "oh" {
#run  oh.ezfio -75.614606131897 -75.618732794235
#}
#
#@test "sih2_3b1" {
#run  sih2_3b1.ezfio -290.016780973071 -290.017278798946
#}
#
#@test "sih3" {
#run  sih3.ezfio -5.575343504534 -5.577437627802
#}
#
#@test "so" {
#run  so.ezfio -26.035945181998 -26.046539528491
#}

