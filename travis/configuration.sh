#!/bin/bash
# Stage 1

# Configure QP2
./configure --install all --config ./config/travis.cfg || exit -1

# Create cache
cd ../
tar -zcf $HOME/cache/config.tgz qp2

