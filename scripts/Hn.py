#!/usr/bin/env python3

import sys
from math import *
arg = sys.argv
#f = open('data_dft','r')
n = int(sys.argv[1])
r = float(sys.argv[2])
f = open('H'+str(n)+'_'+str(r)+'.xyz','w')
string=str(n)+"\n"
f.write(string)
string="\n"
f.write(string)
for i in range(n):
 x = r * cos(2.* i* pi/n)
 y = r * sin(2.* i* pi/n)
 z = 0.
 string="H   "+str(x)+"  "+str(y)+"  "+str(z)+"\n"
 f.write(string)
 
#lines = f.readlines()
#cipsi_dft= []
#
#dissoc = []
#dissoc.append(float(-76.0179223470363))
#dissoc.append(float(-76.0592367866993))
#dissoc.append(float(-76.0678739715659))
#delta_e = []
#
#for line in lines:
#  data = line.split()
#  if(len(data)>0):
#   dft=float(data[1])
#   fci=float(data[2])
#   e=fci+dft
#   cipsi_dft.append(e)
#
#print(*cipsi_dft,sep="  &  ")
#
#for i in 0,1,2:
# delta_e.append(1000.*(dissoc[i] - cipsi_dft[i]))
#
#print(*delta_e,sep="  &  ")
#
