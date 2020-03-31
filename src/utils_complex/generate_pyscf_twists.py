#! /usr/bin/env python

from nexus import settings,job,run_project,obj
from nexus import generate_physical_system
from nexus import generate_pyscf
import spglib
import numpy as np
import os

settings(
    results = '',
    sleep   = 3,
    machine = 'ws8',
    generate_only=1,
    )

solid_tmp_file = './.tmp_solid_template.py'
jobparams={}
jobparams['dfname'] = 'df_ints.h5'
jobparams['chkname'] = 'diamond.chk'
jobparams['dftype'] = 'GDF'
jobparams['auxbasis'] = 'weigend'
jobparams['xc'] = 'b3lyp'

show_kmap = True

pyscf_job = job(cores=1,serial=True)

cell_types = [
    'diamond_8_real',
#    'diamond_8_comp',
    ]

cell_info = obj(
    diamond_8_real = obj(
#        tiling = [[ 1, -1,  1],
#                  [ 1,  1, -1],
#                  [-1,  1,  1]],
        tiling = (2,2,2),
        kgrid  = (6,6,6),
        ),
    )

a = 3.37316115
axes = np.array([[a,a,0],[0,a,a],[a,0,a]])
elem = ['C','C']
pos = [[0,0,0],[a/2,a/2,a/2]]

scf_info = obj(
        basis = 'bfd-vdz',
        ecp = 'bfd',
        drop_exponent = 0.1,
        verbose = 5,
        )

tempstr = """
#!/usr/bin/env python

'''
Gamma point post-HF calculation needs only real integrals.
Methods implemented in finite-size system can be directly used here without
any modification.
'''

import numpy as np
from pyscf import lib
from pyscf.pbc import gto, scf, dft
from pyscf import gto as Mgto
from pyscf.pbc import df 
from pyscf.pbc import ao2mo
from pyscf.pbc import tools
from pyscf.pbc.tools.pbc import super_cell
from functools import reduce
import scipy.linalg as la
import os

restart = False

$system

$twistinfo

pwd_top = os.path.dirname(os.path.realpath(__file__))
for i in range(len(allkpts)):
    jobdir=pwd_top + '/twist-{:02d}/'.format(i)
    if not restart:
        os.mkdir(jobdir)
    os.chdir(jobdir) 
    sp_twist=supTwist[i]
    kpts_i=allkpts[i]  
    print("i ",i, kpts_i)
    supcell=cell
    mydf = df.$dftype(supcell,kpts_i)
    mydf.auxbasis = '$auxbasis'
    dfpath = jobdir+'$dfname'
    if not restart:
        mydf._cderi_to_save = dfpath   # new
        mydf.build()                         # new
    #end if
    mf = scf.KROKS(supcell,kpts_i).density_fit()
    mf.xc='$xc'
   
    mf.tol           = 1e-10 
   
    mf.exxdiv = 'ewald'
    mf.with_df = mydf
    chkpath = jobdir + '$chkname'
    mf.chkfile = chkpath
    mf.with_df._cderi = dfpath
    if restart:
        dm = mf.from_chk(chkpath)    # restart
        e_scf=mf.kernel(dm)                    # restart
    else:
        e_scf=mf.kernel()                      # new
    #end if
   
    with open('e_scf','w') as ener:
        ener.write('%s\\n' % (e_scf))
    print('e_scf',e_scf)
   
    #title="S8-twist%s"%i
    #### generated conversion text ###
    #from PyscfToQmcpack import savetoqmcpack
    #savetoqmcpack(cell,mf,title=title,kmesh=kmesh,kpts=kpts_i,sp_twist=sp_twist)
   
    mycas = list(range(0,30))
    #title="S8-Cas30-twist%s"%i
    #### generated conversion text ###
    #from PyscfToQmcpack import savetoqmcpack
    #savetoqmcpack(cell,mf,title=title,kmesh=kmesh,kpts=kpts_i,sp_twist=sp_twist, cas_idx=mycas)
    #### end generated conversion text ###
   
    from MolPyscfToQPkpts import pyscf2QP2
    pyscf2QP2(supcell,mf,kpts=kpts_i,int_threshold = 1E-15,cas_idx=mycas)
    print('Done for Tw%s'%i)
    os.chdir(pwd_top) 

"""

#replace $dfname with df_ints.h5
#replace $chkname with checkpoint file name
#$dftype is GDF, MDF, etc.
#$auxbasis weigend?
#$xc b3lyp
#$twistinfo


for cell_type in cell_types:
    cell_tiling = cell_info[cell_type].tiling
    cell_kgrid  = cell_info[cell_type].kgrid

    diamond = generate_physical_system(
#        axes     = '''
#                   3.37316115  3.37316115  0.00000000
#                   0.00000000  3.37316115  3.37316115
#                   3.37316115  0.00000000  3.37316115''',
#        elem_pos = '''  
#                   C  0.00000000   0.00000000   0.00000000
#                   C  1.686580575  1.686580575  1.686580575''',
        axes = axes,
        elem = elem,
        pos = pos,
        units    = 'B',
        tiling   = cell_tiling,
        kgrid    = cell_kgrid,
        kshift   = (0,0,0),
        C        = 4,
        symm_kgrid=True,
        )


    jobparams['twistinfo'] = ''
    if show_kmap:
        print (cell_type)
        print ('===============================')
        s = diamond.structure.copy()
        kmap = s.kmap()

        print ('supercell kpoints/twists')
        jobparams['twistinfo']+='# supercell kpoints/twists\n'
        jobparams['twistinfo']+='supTwist=array([\n'
        for i,k in enumerate(s.kpoints):
            print ('  ',i,k)
            jobparams['twistinfo']+=(str(list(k))+',\n')
        #end for
        jobparams['twistinfo']+='])\n'

        print ('primitive cell kpoints')
        # this should already be written by nexus as part of $system
        #jobparams['twistinfo']+='# primitive cell kpoints\n')
        #jobparams['twistinfo']+='orig=array([\n'
        for i,k in enumerate(s.folded_structure.kpoints):
            print ('  ',i,k)
            #jobparams['twistinfo']+=(str(list(k))+',\n')
        #end for
        #jobparams['twistinfo']+='])\n'

        jobparams['twistinfo']+='# mapping from supercell to primitive cell k-points\n'
        jobparams['twistinfo']+='mymap=array([\n'
        for kmapkey in kmap.sorted_keys():
            jobparams['twistinfo']+=(str(list(kmap[kmapkey]))+',\n')
        jobparams['twistinfo']+='])\n'
        print ('mapping from supercell to primitive cell k-points')
        print (kmap)

        #jobparams['twistinfo']+=('allkpts=array(list(map(lambda xs: list(map(lambda x: orig[x], xs)), mymap)))\n')
        jobparams['twistinfo']+=('allkpts=array(list(map(lambda xs: list(map(lambda x: kpts[x], xs)), mymap)))\n')
        jobparams['twistinfo']+=('kweights=array('+str(list(s.kweights))+')\n')
    #end if

    tmp_template = tempstr
    for key in jobparams.keys():
        tmp_template = tmp_template.replace('$'+key,jobparams[key])
    
    if os.path.isfile(solid_tmp_file):
        raise Exception(solid_tmp_file,'solid_tmp_file already exists: delete file and try again')
    with open(solid_tmp_file,'w') as stmp:
        stmp.write(tmp_template)

    scf = generate_pyscf(
        identifier = 'scf',
        path       = cell_type,
        job        = pyscf_job,
        template   = solid_tmp_file,
        system     = diamond,
        cell = scf_info,
#        cell       = obj(
#            basis         = 'bfd-vtz',
#            ecp           = 'bfd',
#            drop_exponent = 0.1,
#            verbose       = 5,
#            ),
        save_qmc   = True ,
        )
    if os.path.isfile(solid_tmp_file):
        os.remove(solid_tmp_file)

#end for

#if show_kmap:
#    exit()
#end if

run_project()
