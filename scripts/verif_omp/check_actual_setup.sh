#!/bin/sh

echo ""
echo "Please wait..."
echo ""
cd ../../src/fci
ninja || echo "Please recompile from the root"  
echo ""
./check_omp_actual_setup
cd ../../scripts/verif_omp
