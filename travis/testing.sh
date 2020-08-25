#!/bin/bash
# Stage 3

# Extract cache from compile stage
cd ../
tar -zxf $HOME/cache/compil.tgz

# Configure QP2
cd qp2
source ./quantum_package.rc
qp_test -a && rm $HOME/cache/compil.tgz





