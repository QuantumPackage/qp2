#!/usr/bin/env bats

source $QP_ROOT/tests/bats/common.bats.sh
source $QP_ROOT/quantum_package.rc

# Run test_libcint with both integral engines and compare the tagged
# integral lines. A line passes when |a-b| <= tol(tag) or
# |a-b| <= 1e-12*|b|; tol is 1e-10 for the overlaps because QP's
# overlap_gaussian_xyz returns a hard-coded 1e-10 when a primitive pair
# underflows (src/utils/one_e_integration.irp.f), and 1e-11 for the
# other operators.
#
# Extra arguments are passed to create_ezfio (e.g. -m 2 for an odd
# electron count).
function compare_engines() {
  local INPUT=$1
  local BASIS=$2
  local LABEL=$3
  local EZFIO=libcint_${LABEL}.ezfio
  shift 3

  test_exe test_libcint || skip

  cp ${QP_ROOT}/tests/input/$INPUT .
  rm -rf $EZFIO
  qp create_ezfio -b $BASIS "$@" $INPUT -o $EZFIO
  qp set_file $EZFIO

  qp set ao_basis use_libcint true
  qp run test_libcint > libcint_${LABEL}_true.out

  # Without this, a silent fallback would compare QP against itself
  grep -qE '^INFO use_libcint +T$' libcint_${LABEL}_true.out
  grep -qE '^INFO libcint_avail +T$' libcint_${LABEL}_true.out
  grep -qE '^INFO libcint_nuc_avail +T$' libcint_${LABEL}_true.out

  grep -E '^(OVLP|KIN|NAI|ERI|ERF) ' libcint_${LABEL}_true.out \
    > libcint_${LABEL}_true.dat

  qp set ao_basis use_libcint false
  qp run test_libcint > libcint_${LABEL}_false.out
  grep -E '^(OVLP|KIN|NAI|ERI|ERF) ' libcint_${LABEL}_false.out \
    > libcint_${LABEL}_false.dat

  run awk '
    function key(   i, s) {
      s = $1
      for (i = 2; i < NF; i++) s = s " " $i
      return s
    }
    NR == FNR { ref[key()] = $NF; next }
    {
      k = key()
      d = ref[k] - $NF
      if (d < 0) d = -d
      a = ($NF < 0 ? -$NF : $NF)
      tol = ($1 == "OVLP" ? 1.e-10 : 1.e-11)
      if (d > tol && d > 1.e-12 * a) {
        if (bad < 20) {
          printf "MISMATCH %s: %.16e vs %.16e\n", k, ref[k], $NF
        }
        bad++
      }
    }
    END {
      if (bad > 0) {
        printf "%d mismatched integrals\n", bad > "/dev/stderr"
      }
      exit (bad > 0)
    }
  ' libcint_${LABEL}_false.dat libcint_${LABEL}_true.dat
  [ "$status" -eq 0 ]

  rm -rf $EZFIO
  rm -f libcint_${LABEL}_true.out libcint_${LABEL}_false.out \
        libcint_${LABEL}_true.dat libcint_${LABEL}_false.dat
}

@test "libcint p-d" {
  compare_engines ch4.xyz cc-pvdz ch4
}

@test "libcint f" {
  compare_engines f.xyz cc-pvtz f -m 2
}
