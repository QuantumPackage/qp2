#!/bin/bash
# Stage 1

# Configure QP2
./configure --install all --config ./config/travis.cfg

# Create cache
cd ../
tar -zcf $HOME/cache/config.tgz qp2

