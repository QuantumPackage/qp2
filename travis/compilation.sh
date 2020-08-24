#!/bin/bash
# Stage 2

# Extract cache from config stage
cd $HOME
tar -zxf ./cache/config.tgz
rm ./cache/config.tgz

# Configure QP2
cd $HOME/QuantumPackage/qp2
source ./quantum_package.rc
ninja -j 1 -v

# Create cache
cd $HOME
tar -zcf ./cache/compil.tgz QuantumPackage

