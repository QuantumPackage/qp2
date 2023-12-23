#!/usr/bin/env bats

source $QP_ROOT/tests/bats/common.bats.sh
source $QP_ROOT/quantum_package.rc

zero () {
    if [ -z "$1" ]; then echo 0.0; else echo $1; fi
}

function run() {
  thresh1=1e-10
  thresh2=1e-12
  thresh3=1e-4
  test_exe scf || skip
  qp set_file $1
  qp edit --check
  qp reset -d
  qp set_frozen_core
  qp set localization localization_method boys
  file="$(echo $1 | sed 's/.ezfio//g')"
  energy="$(cat $1/hartree_fock/energy)"
  fb_err1="$(qp run debug_gradient_loc | grep 'Max error' | tail -n 1 | awk '{print $3}')"
  fb_err2="$(qp run debug_hessian_loc | grep 'Max error' | tail -n 1 | awk '{print $3}')"
  qp run localization > $file.loc.out
  fb_energy="$(qp run print_energy | grep -A 1 'Nuclear repulsion energy' | tail -n 1 )"
  fb_c="$(cat $file.loc.out | grep 'Criterion:Core' | tail -n 1 | awk '{print $3}')i"
  fb_i="$(cat $file.loc.out | grep 'Criterion:Inactive' | tail -n 1 | awk '{print $3}')"
  fb_a="$(cat $file.loc.out | grep 'Criterion:Active' | tail -n 1 | awk '{print $3}')"
  fb_v="$(cat $file.loc.out | grep 'Criterion:Virtual' | tail -n 1 | awk '{print $3}')"
  qp reset -a
  qp run scf
  qp set_frozen_core
  qp set localization localization_method pipek 
  pm_err1="$(qp run debug_gradient_loc | grep 'Max error' | tail -n 1 | awk '{print $3}')"
  pm_err2="$(qp run debug_hessian_loc | grep 'Max error' | tail -n 1 | awk '{print $3}')"
  qp run localization > $file.loc.out
  pm_c="$(cat $file.loc.out | grep 'Criterion:Core' | tail -n 1 | awk '{print $3}')i"
  pm_i="$(cat $file.loc.out | grep 'Criterion:Inactive' | tail -n 1 | awk '{print $3}')"
  pm_a="$(cat $file.loc.out | grep 'Criterion:Active' | tail -n 1 | awk '{print $3}')"
  pm_v="$(cat $file.loc.out | grep 'Criterion:Virtual' | tail -n 1 | awk '{print $3}')"
  pm_energy="$(qp run print_energy | grep -A 1 'Nuclear repulsion energy' | tail -n 1 )"
  qp set localization localization_method boys
  qp reset -a
  qp run scf
  qp set_frozen_core
  eq $energy $fb_energy $thresh1 
  eq $fb_err1 0.0 $thresh2
  eq $fb_err2 0.0 $thresh2
  eq $energy $pm_energy $thresh1
  eq $pm_err1 0.0 $thresh2
  eq $pm_err2 0.0 $thresh2
  fb_c=$(zero $fb_c)
  fb_i=$(zero $fb_i)
  fb_a=$(zero $fb_a)
  fb_v=$(zero $fb_v)
  pm_c=$(zero $pm_c)
  pm_i=$(zero $pm_i)
  pm_a=$(zero $pm_a)
  pm_v=$(zero $pm_v)
  eq $fb_c $2 $thresh3
  eq $fb_i $3 $thresh3
  eq $fb_a $4 $thresh3
  eq $fb_v $5 $thresh3
  eq $pm_c $6 $thresh3
  eq $pm_i $7 $thresh3
  eq $pm_a $8 $thresh3
  eq $pm_v $9 $thresh3
}

@test "b2_stretched" {
run  b2_stretched.ezfio -32.1357551678876 -47.0041982094667 0.0 -223.470015856259 -1.99990778964451 -2.51376723927071 0.0 -12.8490602539275
}

@test "clo" {
run  clo.ezfio -44.1624001765291 -32.4386660941387 0.0 -103.666309287187 -5.99985418946811 -5.46871580225222 0.0 -20.2480064922275
}

@test "clf" {
run  clf.ezfio -47.5143398826967 -35.7206886315104 0.0 -107.043029033468 -5.99994222062230 -6.63916513458470 0.0 -19.7035159913484
}

@test "h2o2" {
run  h2o2.ezfio -7.76848143170524 -30.9694344369829 0.0 -175.898343829453 -1.99990497554575 -5.62980322957485 0.0 -33.5699813186666
}

@test "h2o" {
run  h2o.ezfio 0.0 -2.52317434969591 0.0 -45.3136377925359 0.0 -3.01248365356981 0.0 -22.4470831240924
}

@test "h3coh" {
run  h3coh.ezfio -3.66763692804590 -24.0463089480870 0.0 -111.485948435075 -1.99714061342078 -4.89242181322988 0.0 -23.6405412057679
}

@test "n2h4" {
run  n2h4.ezfio -7.46608163002070 -35.7632174051822 0.0 -305.913449004632 -1.99989326143356 -4.62496615892268 0.0 -51.5171904685553
}

