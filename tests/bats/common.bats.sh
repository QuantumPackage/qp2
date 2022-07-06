#!/usr/bin/env bats

# floating point number comparison
# Compare two numbers ($1, $2) with a given precision ($3)
# If the numbers are not equal, the exit code is 1 else it is 0
# So we strip the "-", is the abs value of the poor
function eq() {
    declare -a diff
    diff=($(awk -v d1=$1 -v d2=$2 -v n1=${1#-} -v n2=${2#-} -v p=$3 'BEGIN{ if ((n1-n2)^2 < p^2) print 0; print 1 " " (d1-d2) " " d1 " " d2 }'))
    if [[ "${diff[0]}" == "0" ]]
    then
       return 0
    else
       echo "#~-~-~-~-~- Test Failed -~-~-~-~-~-#"
       echo "Test      : " ${BATS_TEST_DESCRIPTION}
       echo "Error     : " ${diff[1]}
       echo "Reference : " ${diff[3]}
       echo "Computed  : " ${diff[2]}
       echo "#~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-#"
       exit 1
    fi
}


#   ___           
#    |  ._  o _|_ 
#   _|_ | | |  |_ 
#                 
source ${QP_EZFIO}/Bash/ezfio.sh
TEST_DIR=${QP_ROOT}/tests/work/

mkdir -p "${TEST_DIR}"

cd "${TEST_DIR}" || exit 1

function test_exe() {
  l_EXE=$(awk "/^$1 / { print \$2 }" < "${QP_ROOT}"/data/executables)
  l_EXE=$(echo $l_EXE | sed "s|\$QP_ROOT|$QP_ROOT|")
  if [[ -x "$l_EXE" ]]
  then
    return 0
  else
    return 127
  fi
}

run_only_test() {
  if [[ "$BATS_TEST_DESCRIPTION" != "$1" ]] && [[ "$BATS_TEST_NUMBER" != "$1" ]]; then
    if [[ -z "$BATS_TEST_FILENAME" ]] ; then
      exit 0
    else
      skip
    fi
  fi
#  sleep 1
}

setup() {
  if [[ -n $TEST ]] ; then
    run_only_test $TEST || exit 0
  fi
}
