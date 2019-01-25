BEGIN_SHELL [ /usr/bin/env python2 ]
from perturbation import perturbations
import os

filename = os.environ["QP_ROOT"]+"/src/perturbation/perturbation.template.f"
file = open(filename,'r')
template = file.read()
file.close()

for p in perturbations:
  print template.replace("$PERT",p)

END_SHELL
