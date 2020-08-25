#!/bin/bash
# Stage 2

# Extract cache from config stage
cd ../
tar -zxf $HOME/cache/config.tgz

# Configure QP2
cd qp2
source ./quantum_package.rc
ninja -j 1 -v

# Create cache
cd ..
tar -zcf $HOME/cache/compil.tgz qp2 && rm $HOME/cache/config.tgz

