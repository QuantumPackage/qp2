#!/bin/bash
# Stage 3

# Extract cache from compile stage
cd $HOME
tar -zxf ./cache/compil.tgz
rm ./cache/compil.tgz

# Configure QP2
cd $HOME/QuantumPackage/qp2
source ./quantum_package.rc
qp_test -a


