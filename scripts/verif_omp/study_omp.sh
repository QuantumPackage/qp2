#!/bin/sh

list_comp="ifort gfortran-7 gfortran-8 gfortran-9"

FILE=results.dat

touch $FILE
rm $FILE

echo "1: omp_set_max_active_levels(5)" >> $FILE
echo "2: omp_set_nested(.True.)" >> $FILE
echo "3: 1 + 2" >> $FILE
echo "" >> $FILE
echo "1 2 3" >> $FILE
for comp in $list_comp
do
	$comp --version > /dev/null \
        && $comp -O0 -fopenmp check_omp_v2.f90 \
	&& echo $(./a.out | grep "Tests:" | cut -d ":" -f2- ) $(echo " : ") $($comp --version | head -n 1) >> $FILE

done

cat $FILE

