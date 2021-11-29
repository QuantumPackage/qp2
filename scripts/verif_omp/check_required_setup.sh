#!/bin/sh

# take one argument which is the compiler used
# return the required IRPF90_FLAGS for the $1 compiler

if [ -z "$1" ]
then
    echo "Give the compiler in argument"
else

$1 --version > /dev/null \
&& $1 -O0 -fopenmp check_omp.f90 \
&& ./a.out | tail -n 1


# if there is an error or if the compiler is not found
$1 --version > /dev/null || echo 'compiler not found'

fi
