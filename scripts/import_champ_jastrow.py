#!/usr/bin/env python3

conv = [ 0, 0, 2 , 6 , 13 , 23 , 37 , 55 , 78 , 106 , 140 ]


def import_jastrow(jastrow_filename):
  with open(jastrow_filename,'r') as jastrow_file:
    lines = [ line.strip() for line in jastrow_file.readlines() ]
  lines = [ line for line in lines if line != "" ]
  start = 0
  end = len(lines)
  for i,line in enumerate(lines):
    if line.startswith("jastrow_parameter"):
      start = i
    elif line.startswith("end"):
      end = i
  lines = lines[start:end]
  type_num = (len(lines)-4)//2
  nord_a,nord_b,nord_c = [ int(i) for i in lines[1].split()[:3] ]
  scale_k = float(lines[2].split()[0])
  vec_a = []
  for j in range(type_num):
    vec_a += [ float(i) for i in lines[3+j].split()[:nord_a+1] ]
  vec_b = [ float(i) for i in lines[3+type_num].split()[:nord_b+1] ]
  vec_c = []
  for j in range(type_num):
    vec_c += [ float(i) for i in lines[4+type_num+j].split()[:conv[nord_c]] ]

  return {
    'type_num' : type_num,
    'scale_k' : scale_k,
    'nord_a' : nord_a,
    'nord_b' : nord_b,
    'nord_c' : nord_c,
    'vec_a' : vec_a,
    'vec_b' : vec_b,
    'vec_c' : vec_c,
  }


if __name__ == '__main__':
  import sys
  from ezfio import ezfio
  ezfio.set_file(sys.argv[1])
  jastrow_file = sys.argv[2]
  jastrow = import_jastrow(jastrow_file)
  print (jastrow)
  ezfio.set_jastrow_jast_type("Qmckl")
  ezfio.set_jastrow_jast_qmckl_type_nucl_num(jastrow['type_num'])
  charges = ezfio.get_nuclei_nucl_charge()
  types = {}
  k = 1
  for c in charges:
    if c not in types:
      types[c] = k
      k += 1
  type_nucl_vector = [types[c] for c in charges]
  print(type_nucl_vector)
  ezfio.set_jastrow_jast_qmckl_type_nucl_vector(type_nucl_vector)
  ezfio.set_jastrow_jast_qmckl_rescale_ee(jastrow['scale_k'])
  ezfio.set_jastrow_jast_qmckl_rescale_en([jastrow['scale_k'] for i in type_nucl_vector])
  ezfio.set_jastrow_jast_qmckl_aord_num(jastrow['nord_a'])
  ezfio.set_jastrow_jast_qmckl_bord_num(jastrow['nord_b'])
  ezfio.set_jastrow_jast_qmckl_cord_num(jastrow['nord_c'])
  ezfio.set_jastrow_jast_qmckl_c_vector_size(len(jastrow['vec_c']))
  ezfio.set_jastrow_jast_qmckl_a_vector(jastrow['vec_a'])
  ezfio.set_jastrow_jast_qmckl_b_vector(jastrow['vec_b'])
  ezfio.set_jastrow_jast_qmckl_c_vector(jastrow['vec_c'])

