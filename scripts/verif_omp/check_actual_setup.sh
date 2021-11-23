#!/bin/sh

# go in qp2/src/fci to run check_omp_actual_setup
# to see if we can run in parallel an omp section in another one
echo ""
echo "Please wait..."
echo ""
cd ../../src/fci
ninja || echo "Please recompile from the root"  
echo ""
./check_omp_actual_setup
cd ../../scripts/verif_omp
