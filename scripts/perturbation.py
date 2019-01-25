#!/usr/bin/env python2

import os

Pert_dir = os.environ["QP_ROOT"]+"/src/perturbation/"

perturbations = []

for filename in filter(lambda x: x.endswith(".irp.f"), os.listdir(Pert_dir)):

  filename = Pert_dir+filename
  file = open(filename,'r')
  lines = file.readlines()
  file.close()
  for line in lines:
      buffer = line.lower().lstrip().split()
      if len(buffer) > 1:
        if buffer[0] == "subroutine" and buffer[1].startswith("pt2_"):
          p = (buffer[1].split('(')[0])[4:]
          perturbations.append( p )


if __name__ == '__main__':
  print 'Perturbations:'
  for k in perturbations:
    print '* ', k
