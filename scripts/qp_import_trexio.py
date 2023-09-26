#!/usr/bin/env python3
"""
convert TREXIO file to EZFIO

Usage:
    qp_import_trexio [-o EZFIO_DIR] FILE

Options:
    -o --output=EZFIO_DIR    Produced directory
                             by default is FILE.ezfio

"""

import sys
import os
import numpy as np
from functools import reduce
from ezfio import ezfio
from docopt import docopt
import qp_bitmasks

try:
  import trexio
except ImportError:
    print("Error: trexio python module is not found. Try python3 -m pip install trexio")
    sys.exit(1)


try:
    QP_ROOT = os.environ["QP_ROOT"]
    QP_EZFIO = os.environ["QP_EZFIO"]
except KeyError:
    print("Error: QP_ROOT environment variable not found.")
    sys.exit(1)
else:
    sys.path = [QP_EZFIO + "/Python",
                QP_ROOT + "/install/resultsFile",
                QP_ROOT + "/install",
                QP_ROOT + "/scripts"] + sys.path

def uint64_to_int64(u):
    # Check if the most significant bit is set
    if u & (1 << 63):
        # Calculate the two's complement
        result = -int(np.bitwise_not(np.uint64(u))+1)
    else:
        # The number is already positive
        result = u
    return result

def generate_xyz(l):

    def create_z(x,y,z):
       return (x, y, l-(x+y))

    def create_y(accu,x,y,z):
       if y == 0:
          result = [create_z(x,y,z)] + accu
       else:
          result = create_y([create_z(x,y,z)] + accu , x, y-1, z)
       return result

    def create_x(accu,x,y,z):
       if x == 0:
          result = create_y([], x,y,z) + accu
       else:
          xnew = x-1
          ynew = l-xnew
          result = create_x(create_y([],x,y,z) + accu , xnew, ynew, z)
       return result

    result = create_x([], l, 0, 0)
    result.reverse()
    return result



def write_ezfio(trexio_filename, filename):

    try:
        trexio_file = trexio.File(trexio_filename,mode='r',back_end=trexio.TREXIO_TEXT)
    except:
        trexio_file = trexio.File(trexio_filename,mode='r',back_end=trexio.TREXIO_HDF5)

    ezfio.set_file(filename)
    ezfio.set_trexio_trexio_file(trexio_filename)

    print("Nuclei\t\t...\t", end=' ')

    charge = [0.]
    if trexio.has_nucleus(trexio_file):
        charge = trexio.read_nucleus_charge(trexio_file)
        ezfio.set_nuclei_nucl_num(len(charge))
        ezfio.set_nuclei_nucl_charge(charge)

        coord = trexio.read_nucleus_coord(trexio_file)
        coord = np.transpose(coord)
        ezfio.set_nuclei_nucl_coord(coord)

        label = trexio.read_nucleus_label(trexio_file)
        nucl_num = trexio.read_nucleus_num(trexio_file)

        # Transformt H1 into H
        import re
        p = re.compile(r'(\d*)$')
        label = [p.sub("", x).capitalize() for x in label]
        ezfio.set_nuclei_nucl_label(label)
        print("OK")

    else:
        ezfio.set_nuclei_nucl_num(1)
        ezfio.set_nuclei_nucl_charge([0.])
        ezfio.set_nuclei_nucl_coord([0.,0.,0.])
        ezfio.set_nuclei_nucl_label(["X"])
        print("None")



    print("Electrons\t...\t", end=' ')

    try:
        num_beta = trexio.read_electron_dn_num(trexio_file)
    except:
        num_beta = int(sum(charge))//2

    try:
        num_alpha = trexio.read_electron_up_num(trexio_file)
    except:
        num_alpha = int(sum(charge)) - num_beta

    if num_alpha == 0:
        print("\n\nError: There are zero electrons in the TREXIO file.\n\n")
        sys.exit(1)
    ezfio.set_electrons_elec_alpha_num(num_alpha)
    ezfio.set_electrons_elec_beta_num(num_beta)

    print(f"{num_alpha} {num_beta}")

    print("Basis\t\t...\t", end=' ')

    shell_num = 0
    try:
        basis_type = trexio.read_basis_type(trexio_file)

        if basis_type.lower() in ["gaussian", "slater"]:
            shell_num   = trexio.read_basis_shell_num(trexio_file)
            prim_num    = trexio.read_basis_prim_num(trexio_file)
            ang_mom     = trexio.read_basis_shell_ang_mom(trexio_file)
            nucl_index  = trexio.read_basis_nucleus_index(trexio_file)
            exponent    = trexio.read_basis_exponent(trexio_file)
            coefficient = trexio.read_basis_coefficient(trexio_file)
            shell_index = trexio.read_basis_shell_index(trexio_file)
            ao_shell    = trexio.read_ao_shell(trexio_file)

            ezfio.set_basis_basis("Read from TREXIO")
            ezfio.set_ao_basis_ao_basis("Read from TREXIO")
            ezfio.set_basis_shell_num(shell_num)
            ezfio.set_basis_prim_num(prim_num)
            ezfio.set_basis_shell_ang_mom(ang_mom)
            ezfio.set_basis_basis_nucleus_index([ x+1 for x in nucl_index ])
            ezfio.set_basis_prim_expo(exponent)
            ezfio.set_basis_prim_coef(coefficient)

            nucl_shell_num = []
            prev = None
            m = 0
            for i in ao_shell:
                if i != prev:
                   m += 1
                   if prev is None or nucl_index[i] != nucl_index[prev]:
                        nucl_shell_num.append(m)
                        m = 0
                prev = i
            assert (len(nucl_shell_num) == nucl_num)

            shell_prim_num = []
            prev = shell_index[0]
            count = 0
            for i in shell_index:
                if i != prev:
                   shell_prim_num.append(count)
                   count = 0
                count += 1
                prev = i
            shell_prim_num.append(count)

            assert (len(shell_prim_num) == shell_num)

            ezfio.set_basis_shell_prim_num(shell_prim_num)
            ezfio.set_basis_shell_index([x+1 for x in shell_index])
            ezfio.set_basis_nucleus_shell_num(nucl_shell_num)


            shell_factor = trexio.read_basis_shell_factor(trexio_file)
            prim_factor  = trexio.read_basis_prim_factor(trexio_file)

        elif basis_type.lower() == "numerical":

            shell_num   = trexio.read_basis_shell_num(trexio_file)
            prim_num    = shell_num
            ang_mom     = trexio.read_basis_shell_ang_mom(trexio_file)
            nucl_index  = trexio.read_basis_nucleus_index(trexio_file)
            exponent    = [1.]*prim_num
            coefficient = [1.]*prim_num
            shell_index = [i for i in range(shell_num)]
            ao_shell    = trexio.read_ao_shell(trexio_file)

            ezfio.set_basis_basis("None")
            ezfio.set_ao_basis_ao_basis("None")
            ezfio.set_basis_shell_num(shell_num)
            ezfio.set_basis_prim_num(prim_num)
            ezfio.set_basis_shell_ang_mom(ang_mom)
            ezfio.set_basis_basis_nucleus_index([ x+1 for x in nucl_index ])
            ezfio.set_basis_prim_expo(exponent)
            ezfio.set_basis_prim_coef(coefficient)

            nucl_shell_num = []
            prev = None
            m = 0
            for i in ao_shell:
                if i != prev:
                   m += 1
                   if prev is None or nucl_index[i] != nucl_index[prev]:
                        nucl_shell_num.append(m)
                        m = 0
                prev = i
            assert (len(nucl_shell_num) == nucl_num)

            shell_prim_num = []
            prev = shell_index[0]
            count = 0
            for i in shell_index:
                if i != prev:
                   shell_prim_num.append(count)
                   count = 0
                count += 1
                prev = i
            shell_prim_num.append(count)

            assert (len(shell_prim_num) == shell_num)

            ezfio.set_basis_shell_prim_num(shell_prim_num)
            ezfio.set_basis_shell_index([x+1 for x in shell_index])
            ezfio.set_basis_nucleus_shell_num(nucl_shell_num)

            shell_factor = trexio.read_basis_shell_factor(trexio_file)
            prim_factor  = [1.]*prim_num
        else:
           raise TypeError

        print(basis_type)
    except:
        print("None")
        ezfio.set_ao_basis_ao_cartesian(True)

    print("AOS\t\t...\t", end=' ')

    try:
        cartesian = trexio.read_ao_cartesian(trexio_file)
    except:
        cartesian = True

    if not cartesian:
        raise TypeError('Only cartesian TREXIO files can be converted')

    ao_num = trexio.read_ao_num(trexio_file)
    ezfio.set_ao_basis_ao_num(ao_num)

    if shell_num > 0:
        ao_shell    = trexio.read_ao_shell(trexio_file)
        at = [ nucl_index[i]+1 for i in ao_shell ]
        ezfio.set_ao_basis_ao_nucl(at)

        num_prim0 = [ 0 for i in range(shell_num) ]
        for i in shell_index:
            num_prim0[i] += 1

        coef = {}
        expo = {}
        for i,c in enumerate(coefficient):
            idx = shell_index[i]
            if idx in coef:
              coef[idx].append(c)
              expo[idx].append(exponent[i])
            else:
              coef[idx] = [c]
              expo[idx] = [exponent[i]]

        coefficient = []
        exponent    = []
        power_x     = []
        power_y     = []
        power_z     = []
        num_prim    = []

        for i in range(shell_num):
            for x,y,z in generate_xyz(ang_mom[i]):
                power_x.append(x)
                power_y.append(y)
                power_z.append(z)
                coefficient.append(coef[i])
                exponent.append(expo[i])
                num_prim.append(num_prim0[i])

        assert (len(coefficient) == ao_num)
        ezfio.set_ao_basis_ao_power(power_x + power_y + power_z)
        ezfio.set_ao_basis_ao_prim_num(num_prim)

        prim_num_max = max( [ len(x) for x in coefficient ] )

        for i in range(ao_num):
            coefficient[i] += [0. for j in range(len(coefficient[i]), prim_num_max)]
            exponent   [i] += [0. for j in range(len(exponent[i]), prim_num_max)]

        coefficient = reduce(lambda x, y: x + y, coefficient, [])
        exponent    = reduce(lambda x, y: x + y, exponent   , [])

        coef = []
        expo = []
        for i in range(prim_num_max):
            for j in range(i, len(coefficient), prim_num_max):
                coef.append(coefficient[j])
                expo.append(exponent[j])

#        ezfio.set_ao_basis_ao_prim_num_max(prim_num_max)
        ezfio.set_ao_basis_ao_coef(coef)
        ezfio.set_ao_basis_ao_expo(expo)

        print("OK")

    else:
        print("None")


    #                _
    # |\/|  _   _   |_)  _.  _ o  _
    # |  | (_) _>   |_) (_| _> | _>
    #

    print("MOS\t\t...\t", end=' ')

    labels = { "Canonical" : "Canonical",
               "RHF" : "Canonical",
               "BOYS" : "Localized",
               "ROHF" : "Canonical",
               "UHF" : "Canonical",
               "Natural": "Natural" }
    try:
      label = labels[trexio.read_mo_type(trexio_file)]
    except:
      label = "None"
    ezfio.set_mo_basis_mo_label(label)
    ezfio.set_determinants_mo_label(label)

    try:
      clss = trexio.read_mo_class(trexio_file)
      core     = [ i for i in clss if i.lower() == "core" ]
      inactive = [ i for i in clss if i.lower() == "inactive" ]
      active   = [ i for i in clss if i.lower() == "active" ]
      virtual  = [ i for i in clss if i.lower() == "virtual" ]
      deleted  = [ i for i in clss if i.lower() == "deleted" ]
    except trexio.Error:
      pass

    try:
      mo_num = trexio.read_mo_num(trexio_file)
      ezfio.set_mo_basis_mo_num(mo_num)

      MoMatrix = trexio.read_mo_coefficient(trexio_file)
      ezfio.set_mo_basis_mo_coef(MoMatrix)

      mo_occ = [ 0. for i in range(mo_num) ]
      for i in range(num_alpha):
         mo_occ[i] += 1.
      for i in range(num_beta):
         mo_occ[i] += 1.
      ezfio.set_mo_basis_mo_occ(mo_occ)
      print("OK")
    except:
      print("None")



    print("Pseudos\t\t...\t", end=' ')

    ezfio.set_pseudo_do_pseudo(False)

    if trexio.has_ecp_ang_mom(trexio_file):
        ezfio.set_pseudo_do_pseudo(True)
        max_ang_mom_plus_1 = trexio.read_ecp_max_ang_mom_plus_1(trexio_file)
        z_core = trexio.read_ecp_z_core(trexio_file)
        ang_mom = trexio.read_ecp_ang_mom(trexio_file)
        nucleus_index = trexio.read_ecp_nucleus_index(trexio_file)
        exponent = trexio.read_ecp_exponent(trexio_file)
        coefficient = trexio.read_ecp_coefficient(trexio_file)
        power = trexio.read_ecp_power(trexio_file)

        lmax = max( max_ang_mom_plus_1 ) - 1
        ezfio.set_pseudo_pseudo_lmax(lmax)
        ezfio.set_pseudo_nucl_charge_remove(z_core)

        prev_center = None
        ecp = {}
        for i in range(len(ang_mom)):
            center = nucleus_index[i]
            if center != prev_center:
               ecp[center] = { "lmax": max_ang_mom_plus_1[center],
                               "zcore": z_core[center],
                               "contr": {} }
               for j in range(max_ang_mom_plus_1[center]+1):
                    ecp[center]["contr"][j] = []

            ecp[center]["contr"][ang_mom[i]].append( (coefficient[i], power[i], exponent[i]) )
            prev_center = center

        ecp_loc = {}
        ecp_nl  = {}
        kmax    = 0
        klocmax    = 0
        for center in ecp:
            ecp_nl [center] = {}
            for k in ecp[center]["contr"]:
                if k == ecp[center]["lmax"]:
                    ecp_loc[center] = ecp[center]["contr"][k]
                    klocmax = max(len(ecp_loc[center]), klocmax)
                else:
                    ecp_nl [center][k] = ecp[center]["contr"][k]
                    kmax = max(len(ecp_nl [center][k]), kmax)

        ezfio.set_pseudo_pseudo_klocmax(klocmax)
        ezfio.set_pseudo_pseudo_kmax(kmax)

        pseudo_n_k = [[0  for _ in range(nucl_num)] for _ in range(klocmax)]
        pseudo_v_k = [[0. for _ in range(nucl_num)] for _ in range(klocmax)]
        pseudo_dz_k = [[0. for _ in range(nucl_num)] for _ in range(klocmax)]
        pseudo_n_kl = [[[0  for _ in range(nucl_num)] for _ in range(kmax)] for _ in range(lmax+1)]
        pseudo_v_kl = [[[0. for _ in range(nucl_num)] for _ in range(kmax)] for _ in range(lmax+1)]
        pseudo_dz_kl = [[[0. for _ in range(nucl_num)] for _ in range(kmax)] for _ in range(lmax+1)]
        for center in ecp_loc:
            for k in range( len(ecp_loc[center]) ):
                v, n, dz = ecp_loc[center][k]
                pseudo_n_k[k][center] = n
                pseudo_v_k[k][center] = v
                pseudo_dz_k[k][center] = dz

        ezfio.set_pseudo_pseudo_n_k(pseudo_n_k)
        ezfio.set_pseudo_pseudo_v_k(pseudo_v_k)
        ezfio.set_pseudo_pseudo_dz_k(pseudo_dz_k)

        for center in ecp_nl:
            for l in range( len(ecp_nl[center]) ):
                for k in range( len(ecp_nl[center][l]) ):
                    v, n, dz = ecp_nl[center][l][k]
                    pseudo_n_kl[l][k][center] = n
                    pseudo_v_kl[l][k][center] = v
                    pseudo_dz_kl[l][k][center] = dz

        ezfio.set_pseudo_pseudo_n_kl(pseudo_n_kl)
        ezfio.set_pseudo_pseudo_v_kl(pseudo_v_kl)
        ezfio.set_pseudo_pseudo_dz_kl(pseudo_dz_kl)
        print("OK")

    else:
        print("None")

    print("Determinant\t...\t", end=' ')
    alpha = [ i for i in range(num_alpha) ]
    beta  = [ i for i in range(num_beta) ]
    if trexio.has_mo_spin(trexio_file):
       spin = trexio.read_mo_spin(trexio_file)
       if max(spin) == 1:
         alpha = [ i for i in range(len(spin)) if spin[i] == 0 ]
         alpha = [ alpha[i] for i in range(num_alpha) ]
         beta  = [ i for i in range(len(spin)) if spin[i] == 1 ]
         beta  = [ beta[i] for i in range(num_beta) ]
         print("Warning -- UHF orbitals --", end=' ')
    alpha_s = ['0']*mo_num
    beta_s  = ['0']*mo_num
    for i in alpha:
      alpha_s[i] = '1'
    for i in beta:
      beta_s[i] = '1'
    alpha_s = ''.join(alpha_s)[::-1]
    beta_s = ''.join(beta_s)[::-1]
    def conv(i):
      try:
        result = np.int64(i)
      except:
        result = np.int64(i-2**63-1)
      return result

    alpha = [ uint64_to_int64(int(i,2)) for i in qp_bitmasks.string_to_bitmask(alpha_s) ][::-1]
    beta  = [ uint64_to_int64(int(i,2)) for i in qp_bitmasks.string_to_bitmask(beta_s ) ][::-1]
    ezfio.set_determinants_bit_kind(8)
    ezfio.set_determinants_n_int(1+mo_num//64)
    ezfio.set_determinants_n_det(1)
    ezfio.set_determinants_n_states(1)
    ezfio.set_determinants_psi_det(alpha+beta)
    ezfio.set_determinants_psi_coef([[1.0]])
    print("OK")




def get_full_path(file_path):
    file_path = os.path.expanduser(file_path)
    file_path = os.path.expandvars(file_path)
    return file_path


if __name__ == '__main__':
    ARGUMENTS = docopt(__doc__)

    FILE = get_full_path(ARGUMENTS['FILE'])
    trexio_filename = FILE

    if ARGUMENTS["--output"]:
        EZFIO_FILE = get_full_path(ARGUMENTS["--output"])
    else:
        EZFIO_FILE = "{0}.ezfio".format(FILE)

    write_ezfio(trexio_filename, EZFIO_FILE)
    sys.stdout.flush()

