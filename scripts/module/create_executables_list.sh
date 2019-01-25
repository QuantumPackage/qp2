#!/bin/bash 
#
# Thu Mar 26 01:27:14 CET 2015

if [[ -z ${QP_ROOT} ]]
then
  print "The QP_ROOT environment variable is not set."
  print "Please reload the quantum_package.rc file."
  exit -1
fi

cd ${QP_ROOT}/data
rm -f executables
EXES=$(find -L ${QP_ROOT}/src -maxdepth 2 -depth -executable -type f | grep -e "${QP_ROOT}/src/[^/]*/[^/]*$" |sort ) 

for EXE in $EXES
do
   printf "%-30s %s\n" $(basename $EXE) $EXE | sed "s|${QP_ROOT}|\$QP_ROOT|g" >> executables 
done
