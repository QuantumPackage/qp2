#!/usr/bin/env bats

source $QP_ROOT/tests/bats/common.bats.sh
source $QP_ROOT/quantum_package.rc

function run() {
  thresh=1.e-5
  test_exe cisd || skip
  qp edit --check
  qp set determinants n_states  2
  qp set davidson threshold_davidson 1.e-12
  qp set davidson n_states_diag 24
  qp run cisd 
  energy1="$(qp get cisd energy | tr '[]' ' ' | cut -d ',' -f 1)"
  energy2="$(qp get cisd energy | tr '[]' ' ' | cut -d ',' -f 2)"
  eq $energy1 $1 $thresh
  eq $energy2 $2 $thresh
}


@test "B-B" { # 
  qp set_file b2_stretched.ezfio
  run  -49.120607088648597 -49.055152453388231
}

@test "SiH2_3B1" { # 1.53842s 3.53856s
  qp set_file sih2_3b1.ezfio
  run -290.015949171697 -289.805036176618
}

@test "HBO" { # 4.42968s  19.6099s
  qp set_file hbo.ezfio
  run -100.2019254455993 -99.79484127741013 
}

@test "HCO" { # 6.6077s 28.6801s
  qp set_file hco.ezfio
  run -113.39088802205114 -113.22204293108558
}

@test "H2O" { # 7.0651s 30.6642s
  qp set_file h2o.ezfio
  run -76.22975602077072 -75.80609108747208 
}





@test "H2S" { # 7.42152s 32.5461s
  [[ -n $TRAVIS ]] && skip
  qp set_file h2s.ezfio
  run -398.853701416768 -398.519020035337
}

@test "N2H4" { # 15.8394s 1.27651m
  [[ -n $TRAVIS ]] && skip
  qp set_file n2h4.ezfio
  qp set_mo_class --core="[1-2]" --act="[3-24]" --del="[25-48]"
  run -111.366247464687 -110.990795989548
}

@test "H2O2" { # 16.3164s 1.46453m
  [[ -n $TRAVIS ]] && skip
  qp set_file h2o2.ezfio
  qp set_mo_class --core="[1-2]" --act="[3-24]" --del="[25-38]"
  run -151.003775695363 -150.650247854914
}

@test "OH" { # 18.2159s 1.28453m
  [[ -n $TRAVIS ]] && skip
  qp set_file oh.ezfio
  run -75.6087472926588 -75.5370393736601
}

@test "CH4" { # 19.821s 1.38648m
  [[ -n $TRAVIS ]] && skip
  qp set_file ch4.ezfio
  qp set_mo_class --core="[1]" --act="[2-30]" --del="[31-59]"
  run -40.2403962667047 -39.8433221754964
}

@test "SiH3" { # 20.2202s 1.38648m
  [[ -n $TRAVIS ]] && skip
  qp set_file sih3.ezfio
  run -5.57096611856522  -5.30950347928823
}

@test "NH3" { # 20.6771s  1.23448m
  [[ -n $TRAVIS ]] && skip
  qp set_file nh3.ezfio
  qp set_mo_class --core="[1-4]" --act="[5-72]"
  run -56.2447484835843 -55.9521689975716
}

@test "DHNO" { # 24.7077s 1.46487m
  [[ -n $TRAVIS ]] && skip
  qp set_file dhno.ezfio
  qp set_mo_class --core="[1-7]" --act="[8-64]" 
  run -130.458814562403 -130.356308303681
}

@test "H3COH" { # 24.7248s 1.85043m
  [[ -n $TRAVIS ]] && skip
  qp set_file h3coh.ezfio
  run -115.204958752377 -114.755913828245
}

@test "[Cu(NH3)4]2+" { # 29.9956s 2.15761m
  [[ -n $TRAVIS ]] && skip
  qp set_file cu_nh3_4_2plus.ezfio
  qp set_mo_class --core="[1-24]" --act="[25-45]" --del="[46-87]"
  run -1862.98689579931  -1862.6883044626563

}

@test "ClF" { # 30.3225s
  [[ -n $TRAVIS ]] && skip
  qp set_file clf.ezfio
  run -559.162476603880  -558.792395927088
}

@test "C2H2" { # 35.3324s
  [[ -n $TRAVIS ]] && skip
  qp set_file c2h2.ezfio
  qp set_mo_class --act="[1-30]" --del="[31-36]"
  run -12.3566731164213 -11.9495394759914 
}

@test "ClO" { # 37.6949s
  [[ -n $TRAVIS ]] && skip
  qp set_file clo.ezfio
  run -534.5404021326773 -534.3818725793897 
}

@test "F2" { # 45.2078s
  [[ -n $TRAVIS ]] && skip
  qp set_file f2.ezfio
  qp set_mo_class --core="[1,2]" --act="[3-30]" --del="[31-62]"
  run -199.056829527539 -198.731828008346
}

@test "SO2" { # 47.6922s
  [[ -n $TRAVIS ]] && skip
  qp set_file so2.ezfio
  qp set_mo_class --core="[1-8]" --act="[9-87]" 
  run -41.5746738710350 -41.3800467740750
}

@test "SO" { # 51.2476s
  [[ -n $TRAVIS ]] && skip
  qp set_file so.ezfio
  run -26.0131812819785 -25.7053111980226
}

@test "CO2" { # 95.3736s
  [[ -n $TRAVIS ]] && skip
  qp set_file co2.ezfio
  qp set_mo_class --core="[1,2]" --act="[3-30]" --del="[31-42]"
  run -187.959378390998  -187.432502050556
}

@test "N2" { # 133.1814
  [[ -n $TRAVIS ]] && skip
  qp set_file n2.ezfio
  qp set_mo_class --core="[1,2]" --act="[3-40]" --del="[41-60]"
  run -109.275693633982 -108.757794570948 
}

@test "HCN" { # 133.8696s
  [[ -n $TRAVIS ]] && skip
  qp set_file hcn.ezfio
  qp set_mo_class --core="[1,2]" --act="[3-40]" --del="[41-55]"
  run -93.0776334511721 -92.6684633795506
}

