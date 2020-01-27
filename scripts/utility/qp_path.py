#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import print_function
import os
import sys

try:
    QP_ROOT = os.environ['QP_ROOT']
except:
    print("source quantum_package.rc")
    sys.exit(1)
else:
    QP_EZFIO   = os.environ["QP_EZFIO"]
    QP_SRC     = os.path.join(QP_ROOT, "src")
    QP_PLUGINS = os.path.join(QP_ROOT, "plugins")
    QP_OCAML   = os.path.join(QP_ROOT, "ocaml")
    QP_TESTS   = os.path.join(QP_ROOT, "tests")
    QP_DATA    = os.path.join(QP_ROOT, "data")
