#!/bin/bash

ezfio_name=h2_cc-pvdz
cat <<EOF > convert.py
#!/usr/bin/env python
from pyscf import gto, scf
mol = gto.M(atom='H 0 0 0; H 0 0 1.2', basis='ccpvdz')
mf = scf.RHF(mol)
mf.kernel()

import importlib
import MolPyscfToQPkpts
from MolPyscfToQPkpts import *
pyscf2QP2_mol(mf)
EOF
python convert.py 
mv qpdat.h5  $ezfio_name.h5
qp_convert_h5_to_ezfio ${ezfio_name}.h5
