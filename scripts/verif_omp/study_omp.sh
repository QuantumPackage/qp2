#!/bin/sh

# list of compilers
list_comp="ifort gfortran-7 gfortran-8 gfortran-9"

# file to store the results
FILE=results.dat

touch $FILE
rm $FILE

# Comments
echo "1: omp_set_max_active_levels(5)" >> $FILE
echo "2: omp_set_nested(.True.)" >> $FILE
echo "3: 1 + 2" >> $FILE
echo "" >> $FILE
echo "1 2 3" >> $FILE

# loop on the comp
for comp in $list_comp
do
	$comp --version > /dev/null \
        && $comp -O0 -fopenmp check_omp.f90 \
	&& echo $(./a.out | grep "Tests:" | cut -d ":" -f2- ) $(echo " : ") $($comp --version | head -n 1) >> $FILE

done

# Display
cat $FILE

