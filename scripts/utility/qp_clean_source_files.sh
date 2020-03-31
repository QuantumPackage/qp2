#!/bin/bash
#
# Cleans the source files from non-ascii characters
#
# Tue Mar 31 18:28:42 CEST 2020
#

function help() {
  cat << EOF
Cleans the source files of QP from non-ascii characters.

Usage:

    $(basename $0) [-h|--help]

Options:

    -h    --help    Prints the help message

EOF
  exit 0
}

# Check the QP_ROOT directory
if [[ -z ${QP_ROOT} ]] ; then
  echo "The QP_ROOT environment variable is not set."
  echo "Please reload the quantum_package.rc file."
  exit 1
fi


FILES=$(grep -P "\xA0" ${QP_ROOT}/src/*/*.f | cut -d ':' -f 1 | sort | uniq)
for F in $FILES ; do
  echo "Cleaning $F"
  vim -c "% s/\%xA0/ /g" -c ":wq" $F
done

