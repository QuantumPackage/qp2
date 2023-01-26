#!/usr/bin/env python  
import os
import sys

# First argument is the EZFIO file
# It reads a file EZFIO_point_charges.xyz written in this way:
# charge x y z (Angstrom)
# for all charges


def zip_in_ezfio(ezfio,tmp):
  tmpzip=tmp+".gz"
  cmdzip="gzip -c "+tmp+" > "+tmpzip
  os.system(cmdzip)
  os.system("rm "+tmp)
  cmdmv="mv "+tmpzip+" "+EZFIO+"/nuclei/"+tmpzip
  os.system(cmdmv)

def mv_in_ezfio(ezfio,tmp):
  cmdmv="mv "+tmp+" "+EZFIO+"/nuclei/"+tmp
  os.system(cmdmv)


# Getting the EZFIO
EZFIO=sys.argv[1]
EZFIO=EZFIO.replace("/", "")
print(EZFIO)

# Reading the point charges and convert the Angstrom geometry in Bohr for QP
f = open(EZFIO+'_point_charges.xyz','r')
lines = f.readlines()
convert_angs_to_bohr=1.8897259885789233

n_charges=0
coord_x=[]
coord_y=[]
coord_z=[]
charges=[]
for line in lines:
  data = line.split()
  if(len(data)>0):
   n_charges += 1
   charges.append(str(data[0]))
   coord_x.append(str(convert_angs_to_bohr*float(data[1])))
   coord_y.append(str(convert_angs_to_bohr*float(data[2])))
   coord_z.append(str(convert_angs_to_bohr*float(data[3])))

# Write the file containing the number of charges and set in EZFIO folder
tmp="n_pts_charge"
fncharges = open(tmp,'w')
fncharges.write(" "+str(n_charges)+'\n')
fncharges.close()
mv_in_ezfio(EZFIO,tmp)

# Write the file containing the charges and set in EZFIO folder 
tmp="pts_charge_z"
fcharges = open(tmp,'w')
fcharges.write(" 1\n")
fcharges.write(" "+str(n_charges)+'\n')
for i in range(n_charges):
 fcharges.write(charges[i]+'\n')
fcharges.close()
zip_in_ezfio(EZFIO,tmp)

# Write the file containing the charge coordinates and set in EZFIO folder
tmp="pts_charge_coord"
fcoord = open(tmp,'w')
fcoord.write("  2\n")
fcoord.write("                   "+str(n_charges)+'                    3\n')
#fcoord.write(" "+'   3 '+str(n_charges)+' \n')
for i in range(n_charges):
 fcoord.write('   '+coord_x[i]+'\n')
for i in range(n_charges):
 fcoord.write('   '+coord_y[i]+'\n')
for i in range(n_charges):
 fcoord.write('   '+coord_z[i]+'\n')
fcoord.close()
zip_in_ezfio(EZFIO,tmp)

