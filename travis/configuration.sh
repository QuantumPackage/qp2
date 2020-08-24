#!/bin/bash
# Stage 1

# Configure QP2
cd $HOME/QuantumPackage/qp2
./configure --install all --config ./config/travis.cfg

# Create cache
cd $HOME
tar -zcf ./cache/config.tgz QuantumPackage
