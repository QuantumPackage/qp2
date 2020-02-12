import numpy as np
from functools import reduce


def memoize(f):
    memo = {}
    def helper(x):
        if x not in memo:
            memo[x] = f(x)
        return memo[x]
    return helper

@memoize
def idx2_tri(iijj):
    '''
    iijj should be a 2-tuple
    return triangular compound index for (0-indexed counting)
    '''
    ij1=min(iijj)
    ij2=max(iijj)
    return ij1+(ij2*(ij2+1))//2
#    return ij1+(ij2*(ij2-1))//2

def pad(arr_in,outshape):
    arr_out = np.zeros(outshape,dtype=np.complex128)
    dataslice = tuple(slice(0,arr_in.shape[dim]) for dim in range(len(outshape)))
    arr_out[dataslice] = arr_in
    return arr_out

def makesq(vlist,n1,n2):
    '''
    make hermitian matrices of size (n2 x n2) from from lower triangles
    vlist is n1 lower triangles in flattened form
    given: ([a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t],2,4)
           output a 2x4x4 array, where each 4x4 is the square constructed from the lower triangle
    [
     [
      [a  b* d* g*]
      [b  c  e* h*]
      [d  e  f  i*]
      [g  h  i  j ]
     ],
     [
      [k  l* n* q*]
      [l  m  o* r*]
      [n  o  p  s*]
      [q  r  s  t ]
     ]
    ]
    '''
    out=np.zeros([n1,n2,n2],dtype=np.complex128)
    n0 = vlist.shape[0]
    lmask=np.tri(n2,dtype=bool)
    for i in range(n0):
        out[i][lmask] = vlist[i].conj()
    out2=out.transpose([0,2,1])
    for i in range(n0):
        out2[i][lmask] = vlist[i]
    return out2


def makesq3(vlist,n2):
    '''
    make hermitian matrices of size (n2 x n2) from from lower triangles
    vlist is n1 lower triangles in flattened form
    given: ([a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t],2,4)
           output a 2x4x4 array, where each 4x4 is the square constructed from the lower triangle
    [
     [
      [a  b* d* g*]
      [b  c  e* h*]
      [d  e  f  i*]
      [g  h  i  j ]
     ],
     [
      [k  l* n* q*]
      [l  m  o* r*]
      [n  o  p  s*]
      [q  r  s  t ]
     ]
    ]
    '''
    n0 = vlist.shape[0]
    out=np.zeros([n0,n2,n2],dtype=np.complex128)
    lmask=np.tri(n2,dtype=bool)
    for i in range(n0):
        out[i][lmask] = vlist[i].conj()
    out2=out.transpose([0,2,1])
    for i in range(n0):
        out2[i][lmask] = vlist[i]
    return out2

def makesq2(vlist,n1,n2):
    out=np.zeros([n1,n2,n2],dtype=np.complex128)
    lmask=np.tri(n2,dtype=bool)
    tmp=np.zeros([n2,n2],dtype=np.complex128)
    tmp2=np.zeros([n2,n2],dtype=np.complex128)
    for i in range(n1):
        tmp[lmask] = vlist[i].conj()
        tmp2=tmp.T
        tmp2[lmask] = vlist[i]
        out[i] = tmp2.copy()
    return out


def get_phase(cell, kpts, kmesh=None):
    '''
    The unitary transformation that transforms the supercell basis k-mesh
    adapted basis.
    '''
    from pyscf.pbc import tools
    from pyscf import lib

    latt_vec = cell.lattice_vectors()
    if kmesh is None:
        # Guess kmesh
        scaled_k = cell.get_scaled_kpts(kpts).round(8)
        kmesh = (len(np.unique(scaled_k[:,0])),
                 len(np.unique(scaled_k[:,1])),
                 len(np.unique(scaled_k[:,2])))

    R_rel_a = np.arange(kmesh[0])
    R_rel_b = np.arange(kmesh[1])
    R_rel_c = np.arange(kmesh[2])
    R_vec_rel = lib.cartesian_prod((R_rel_a, R_rel_b, R_rel_c))
    R_vec_abs = np.einsum('nu, uv -> nv', R_vec_rel, latt_vec)

    NR = len(R_vec_abs)
    phase = np.exp(1j*np.einsum('Ru, ku -> Rk', R_vec_abs, kpts))
    phase /= np.sqrt(NR)  # normalization in supercell

    # R_rel_mesh has to be construct exactly same to the Ts in super_cell function
    scell = tools.super_cell(cell, kmesh)
    return scell, phase

def mo_k2gamma(cell, mo_energy, mo_coeff, kpts, kmesh=None):
    '''
    Transform MOs in Kpoints to the equivalents supercell
    '''
    from pyscf import lib
    import scipy.linalg as la
    scell, phase = get_phase(cell, kpts, kmesh)

    E_g = np.hstack(mo_energy)
    C_k = np.asarray(mo_coeff)
    Nk, Nao, Nmo = C_k.shape
    NR = phase.shape[0]

    # Transform AO indices
    C_gamma = np.einsum('Rk, kum -> Rukm', phase, C_k)
    C_gamma = C_gamma.reshape(Nao*NR, Nk*Nmo)

    E_sort_idx = np.argsort(E_g)
    E_g = E_g[E_sort_idx]
    C_gamma = C_gamma[:,E_sort_idx]
    s = scell.pbc_intor('int1e_ovlp')
    assert(abs(reduce(np.dot, (C_gamma.conj().T, s, C_gamma))
               - np.eye(Nmo*Nk)).max() < 1e-7)

    # Transform MO indices
    E_k_degen = abs(E_g[1:] - E_g[:-1]).max() < 1e-5
    if np.any(E_k_degen):
        degen_mask = np.append(False, E_k_degen) | np.append(E_k_degen, False)
        shift = min(E_g[degen_mask]) - .1
        f = np.dot(C_gamma[:,degen_mask] * (E_g[degen_mask] - shift),
                   C_gamma[:,degen_mask].conj().T)
        assert(abs(f.imag).max() < 1e-5)

        e, na_orb = la.eigh(f.real, s, type=2)
        C_gamma[:,degen_mask] = na_orb[:, e>0]

    if abs(C_gamma.imag).max() < 1e-7:
        print('!Warning  Some complexe pollutions in MOs are present')

    C_gamma = C_gamma.real
    if  abs(reduce(np.dot, (C_gamma.conj().T, s, C_gamma)) - np.eye(Nmo*Nk)).max() < 1e-7:
        print('!Warning  Some complexe pollutions in MOs are present')

    s_k = cell.pbc_intor('int1e_ovlp', kpts=kpts)
    # overlap between k-point unitcell and gamma-point supercell
    s_k_g = np.einsum('kuv,Rk->kuRv', s_k, phase.conj()).reshape(Nk,Nao,NR*Nao)
    # The unitary transformation from k-adapted orbitals to gamma-point orbitals
    mo_phase = lib.einsum('kum,kuv,vi->kmi', C_k.conj(), s_k_g, C_gamma)

    return mo_phase

def qp2rename():
    import shutil
    qp2names={}
    qp2names['mo_coef_complex'] = 'C.qp'
    qp2names['bielec_ao_complex'] = 'W.qp'

    qp2names['kinetic_ao_complex'] = 'T.qp'
    qp2names['ne_ao_complex'] = 'V.qp'
    qp2names['overlap_ao_complex'] = 'S.qp'


    for old,new in qp2names.items():
        shutil.move(old,new)
    shutil.copy('e_nuc','E.qp')

def pyscf2QP(cell,mf, kpts, kmesh=None, cas_idx=None, int_threshold = 1E-8, 
        print_ao_ints_bi=False, 
        print_mo_ints_bi=False, 
        print_ao_ints_df=True, 
        print_mo_ints_df=False, 
        print_ao_ints_mono=True, 
        print_mo_ints_mono=False):
    '''
    kpts = List of kpoints coordinates. Cannot be null, for gamma is other script
    kmesh = Mesh of kpoints (optional)
    cas_idx = List of active MOs. If not specified all MOs are actives
    int_threshold = The integral will be not printed in they are bellow that
    '''
  
    from pyscf.pbc import ao2mo
    from pyscf.pbc import tools
    from pyscf.pbc.gto import ecp
    import h5py
    
    mo_coef_threshold = int_threshold
    ovlp_threshold = int_threshold
    kin_threshold = int_threshold
    ne_threshold = int_threshold
    bielec_int_threshold = int_threshold

    natom = len(cell.atom_coords())
    print('n_atom per kpt',   natom)
    print('num_elec per kpt', cell.nelectron)
  
    mo_coeff = mf.mo_coeff
    # Mo_coeff actif
    mo_k = np.array([c[:,cas_idx] for c in mo_coeff] if cas_idx is not None else mo_coeff)
    e_k =  np.array([e[cas_idx] for e in mf.mo_energy] if cas_idx is not None else mf.mo_energy)
  
    Nk, nao, nmo = mo_k.shape
    print("n Kpts", Nk)
    print("n active Mos per kpt", nmo)
    print("n AOs per kpt", nao)

    naux = mf.with_df.auxcell.nao
    print("n df fitting functions", naux)
    with open('num_df','w') as f:
        f.write(str(naux))
  
    # Write all the parameter need to creat a dummy EZFIO folder who will containt the integral after.
    # More an implentation detail than a real thing
    with open('param','w') as f:
    # Note the use of nmo_tot
        f.write(' '.join(map(str,(cell.nelectron*Nk, Nk*nmo, natom*Nk))))
  
    with open('num_ao','w') as f:
        f.write(str(nao*Nk))
    with open('num_kpts','w') as f:
        f.write(str(Nk))
    #                             _                             
    # |\ |      _ |  _   _. ._   |_)  _  ._      |  _ o  _  ._  
    # | \| |_| (_ | (/_ (_| |    | \ (/_ |_) |_| | _> | (_) | | 
    #                                    |                      
    
    #Total energy shift due to Ewald probe charge = -1/2 * Nelec*madelung/cell.vol =
    shift = tools.pbc.madelung(cell, kpts)*cell.nelectron * -.5 
    e_nuc = (cell.energy_nuc() + shift)*Nk
  
    print('nucl_repul', e_nuc)
    with open('e_nuc','w') as f:
        f.write(str(e_nuc))
  
  
  
    #       __    __          _                                 
    # |\/| |  |  |    _   _  |_  _ 
    # |  | |__|  |__ (_) (/_ |  _> 
    #                                               
    with open('mo_coef_complex','w') as outfile:
        c_kpts = np.reshape(mo_k,(Nk,nao,nmo))

        for ik in range(Nk):
            shift1=ik*nao+1
            shift2=ik*nmo+1
            for i in range(nao):
                for j in range(nmo):
                    cij = c_kpts[ik,i,j]
                    if abs(cij) > mo_coef_threshold:
                        outfile.write('%s %s %s %s\n' % (i+shift1, j+shift2, cij.real, cij.imag))
    
    # ___                                              
    #  |  ._ _|_  _   _  ._ _. |  _   |\/|  _  ._   _  
    # _|_ | | |_ (/_ (_| | (_| | _>   |  | (_) | | (_) 
    #                 _|                              
   
    if mf.cell.pseudo:
        v_kpts_ao = np.reshape(mf.with_df.get_pp(kpts=kpts),(Nk,nao,nao))
    else:
        v_kpts_ao = np.reshape(mf.with_df.get_nuc(kpts=kpts),(Nk,nao,nao))
    if len(cell._ecpbas) > 0:
        v_kpts_ao += np.reshape(ecp.ecp_int(cell, kpts),(Nk,nao,nao))

    ne_ao = ('ne',v_kpts_ao,ne_threshold)
    ovlp_ao = ('overlap',np.reshape(mf.get_ovlp(cell=cell,kpts=kpts),(Nk,nao,nao)),ovlp_threshold)
    kin_ao = ('kinetic',np.reshape(cell.pbc_intor('int1e_kin',1,1,kpts=kpts),(Nk,nao,nao)),kin_threshold)

    for name, intval_kpts_ao, thresh in (ne_ao, ovlp_ao, kin_ao):
        if print_ao_ints_mono:
            with open('%s_ao_complex' % name,'w') as outfile:
                for ik in range(Nk):
                    shift=ik*nao+1
                    for i in range(nao):
                        for j in range(i,nao):
                            int_ij = intval_kpts_ao[ik,i,j]
                            if abs(int_ij) > thresh:
                                outfile.write('%s %s %s %s\n' % (i+shift, j+shift, int_ij.real, int_ij.imag))
        if print_mo_ints_mono:
            intval_kpts_mo = np.einsum('kim,kij,kjn->kmn',mo_k.conj(),intval_kpts_ao,mo_k)
            with open('%s_mo_complex' % name,'w') as outfile:
                for ik in range(Nk):
                    shift=ik*nmo+1
                    for i in range(nmo):
                        for j in range(i,nmo):
                            int_ij = intval_kpts_mo[ik,i,j]
                            if abs(int_ij) > thresh:
                                outfile.write('%s %s %s %s\n' % (i+shift, j+shift, int_ij.real, int_ij.imag))

  
    # ___                              _    
    #  |  ._ _|_  _   _  ._ _. |  _   |_) o 
    # _|_ | | |_ (/_ (_| | (_| | _>   |_) | 
    #                 _|                    
    #
    kconserv = tools.get_kconserv(cell, kpts)

    with open('kconserv_complex','w') as outfile:
        for a in range(Nk):
            for b in range(Nk):
                for c in range(Nk):
                    d = kconserv[a,b,c]
                    outfile.write('%s %s %s %s\n' % (a+1,c+1,b+1,d+1))
    

    intfile=h5py.File(mf.with_df._cderi,'r')

    j3c = intfile.get('j3c')
    naosq = nao*nao
    naotri = (nao*(nao+1))//2
    j3ckeys = list(j3c.keys())
    j3ckeys.sort(key=lambda strkey:int(strkey))

    # in new(?) version of PySCF, there is an extra layer of groups before the datasets
    # datasets used to be [/j3c/0,   /j3c/1,   /j3c/2,   ...]
    # datasets now are    [/j3c/0/0, /j3c/1/0, /j3c/2/0, ...]
    j3clist = [j3c.get(i+'/0') for i in j3ckeys]
    if j3clist==[None]*len(j3clist):
    # if using older version, stop before last level
        j3clist = [j3c.get(i) for i in j3ckeys]

    nkinvsq = 1./np.sqrt(Nk)

    # dimensions are (kikj,iaux,jao,kao), where kikj is compound index of kpts i and j
    # output dimensions should be reversed (nao, nao, naux, nkptpairs)
    j3arr=np.array([(i.value.reshape([-1,nao,nao]) if (i.shape[1] == naosq) else makesq3(i.value,nao)) * nkinvsq for i in j3clist])

    nkpt_pairs = j3arr.shape[0]

    if print_ao_ints_df:
        with open('df_ao_integral_array','w') as outfile:
            pass
        with open('df_ao_integral_array','a') as outfile:
            for k,kpt_pair in enumerate(j3arr):
                for iaux,dfbasfunc in enumerate(kpt_pair):
                    for i,i0 in enumerate(dfbasfunc):
                        for j,v in enumerate(i0):
                            if (abs(v) > bielec_int_threshold):
                                outfile.write('%s %s %s %s %s %s\n' % (i+1,j+1,iaux+1,k+1,v.real,v.imag))

    if print_mo_ints_df:
        kpair_list=[]
        for i in range(Nk):
            for j in range(Nk):
                if(i>=j):
                    kpair_list.append((i,j,idx2_tri((i,j))))
        j3mo = np.array([np.einsum('mij,ik,jl->mkl',j3arr[kij],mo_k[ki].conj(),mo_k[kj]) for ki,kj,kij in kpair_list])
        with open('df_mo_integral_array','w') as outfile:
            pass
        with open('df_mo_integral_array','a') as outfile:
            for k,kpt_pair in enumerate(j3mo):
                for iaux,dfbasfunc in enumerate(kpt_pair):
                    for i,i0 in enumerate(dfbasfunc):
                        for j,v in enumerate(i0):
                            if (abs(v) > bielec_int_threshold):
                                outfile.write('%s %s %s %s %s %s\n' % (i+1,j+1,iaux+1,k+1,v.real,v.imag))




#    eri_4d_ao = np.zeros((Nk,nao,Nk,nao,Nk,nao,Nk,nao), dtype=np.complex)
#    for d, kd in enumerate(kpts):
#        for c, kc in enumerate(kpts):
#            if c > d: break
#            idx2_cd = idx2_tri(c,d)
#            for b, kb in enumerate(kpts):
#                if b > d: break
#                a = kconserv[b,c,d]
#                if idx2_tri(a,b) > idx2_cd: continue
#                if ((c==d) and (a>b)): continue
#                ka = kpts[a]
#                v = mf.with_df.get_ao_eri(kpts=[ka,kb,kc,kd],compact=False).reshape((nao,)*4)
#                v *= 1./Nk
#                eri_4d_ao[a,:,b,:,c,:,d] = v
#    
#    eri_4d_ao = eri_4d_ao.reshape([Nk*nao]*4)


    if (print_ao_ints_bi or print_mo_ints_bi):
        if print_ao_ints_bi:
            with open('bielec_ao_complex','w') as outfile: 
                pass
        if print_mo_ints_bi:
            with open('bielec_mo_complex','w') as outfile: 
                pass
        for d, kd in enumerate(kpts):
            for c, kc in enumerate(kpts):
                if c > d: break
                idx2_cd = idx2_tri((c,d))
                for b, kb in enumerate(kpts):
                    if b > d: break
                    a = kconserv[b,c,d]
                    if idx2_tri((a,b)) > idx2_cd: continue
                    if ((c==d) and (a>b)): continue
                    ka = kpts[a]

                    if print_ao_ints_bi:
                        with open('bielec_ao_complex','a') as outfile:
                            eri_4d_ao_kpt = mf.with_df.get_ao_eri(kpts=[ka,kb,kc,kd],compact=False).reshape((nao,)*4)
                            eri_4d_ao_kpt *= 1./Nk
                            for l in range(nao):
                                ll=l+d*nao
                                for j in range(nao):
                                    jj=j+c*nao
                                    if jj>ll: break
                                    idx2_jjll = idx2_tri((jj,ll))
                                    for k in range(nao):
                                        kk=k+b*nao
                                        if kk>ll: break
                                        for i in range(nao):
                                            ii=i+a*nao
                                            if idx2_tri((ii,kk)) > idx2_jjll: break
                                            if ((jj==ll) and (ii>kk)): break
                                            v=eri_4d_ao_kpt[i,k,j,l]
                                            if (abs(v) > bielec_int_threshold):
                                                outfile.write('%s %s %s %s %s %s\n' % (ii+1,jj+1,kk+1,ll+1,v.real,v.imag))
            
                    if print_mo_ints_bi:
                        with open('bielec_mo_complex','a') as outfile:
                            eri_4d_mo_kpt = mf.with_df.ao2mo([mo_k[a], mo_k[b], mo_k[c], mo_k[d]],
                                                              [ka,kb,kc,kd],compact=False).reshape((nmo,)*4)
                            eri_4d_mo_kpt *= 1./Nk
                            for l in range(nmo):
                                ll=l+d*nmo
                                for j in range(nmo):
                                    jj=j+c*nmo
                                    if jj>ll: break
                                    idx2_jjll = idx2_tri((jj,ll))
                                    for k in range(nmo):
                                        kk=k+b*nmo
                                        if kk>ll: break
                                        for i in range(nmo):
                                            ii=i+a*nmo
                                            if idx2_tri((ii,kk)) > idx2_jjll: break
                                            if ((jj==ll) and (ii>kk)): break
                                            v=eri_4d_mo_kpt[i,k,j,l]
                                            if (abs(v) > bielec_int_threshold):
                                                outfile.write('%s %s %s %s %s %s\n' % (ii+1,jj+1,kk+1,ll+1,v.real,v.imag))


def pyscf2QP2(cell,mf, kpts, kmesh=None, cas_idx=None, int_threshold = 1E-8, 
        print_ao_ints_bi=False, 
        print_mo_ints_bi=False, 
        print_ao_ints_df=True, 
        print_mo_ints_df=False, 
        print_ao_ints_mono=True, 
        print_mo_ints_mono=False):
    '''
    kpts = List of kpoints coordinates. Cannot be null, for gamma is other script
    kmesh = Mesh of kpoints (optional)
    cas_idx = List of active MOs. If not specified all MOs are actives
    int_threshold = The integral will be not printed in they are bellow that
    '''
  
    from pyscf.pbc import ao2mo
    from pyscf.pbc import tools
    from pyscf.pbc.gto import ecp
    from pyscf.data import nist
    import h5py
    import scipy

    qph5=h5py.File('qpdat.h5')
    qph5.create_group('nuclei')
    qph5.create_group('electrons')
    qph5.create_group('ao_basis')
    qph5.create_group('mo_basis')

    
    mo_coef_threshold = int_threshold
    ovlp_threshold = int_threshold
    kin_threshold = int_threshold
    ne_threshold = int_threshold
    bielec_int_threshold = int_threshold

    natom = cell.natm
    nelec = cell.nelectron
    atom_xyz = mf.cell.atom_coords()
    if not(mf.cell.unit.startswith(('B','b','au','AU'))):
        atom_xyz /= nist.BOHR # always convert to au
    
    strtype=h5py.special_dtype(vlen=str)
    atom_dset=qph5.create_dataset('nuclei/nucl_label',(natom,),dtype=strtype)
    for i in range(natom):
        atom_dset[i] = mf.cell.atom_pure_symbol(i)
    qph5.create_dataset('nuclei/nucl_coord',data=atom_xyz)
    qph5.create_dataset('nuclei/nucl_charge',data=mf.cell.atom_charges())


    print('n_atom per kpt',   natom)
    print('num_elec per kpt', nelec)

    mo_coeff = mf.mo_coeff
    # Mo_coeff actif
    mo_k = np.array([c[:,cas_idx] for c in mo_coeff] if cas_idx is not None else mo_coeff)
    e_k =  np.array([e[cas_idx] for e in mf.mo_energy] if cas_idx is not None else mf.mo_energy)
  
    Nk, nao, nmo = mo_k.shape
    print("n Kpts", Nk)
    print("n active Mos per kpt", nmo)
    print("n AOs per kpt", nao)

    naux = mf.with_df.auxcell.nao
    print("n df fitting functions", naux)
  
    #in old version: param << nelec*Nk, nmo*Nk, natom*Nk
    qph5['electrons'].attrs['elec_alpha_num']=nelec*Nk 
    qph5['electrons'].attrs['elec_beta_num']=nelec*Nk
    qph5['mo_basis'].attrs['mo_num']=Nk*nmo
    qph5['ao_basis'].attrs['ao_num']=Nk*nao
    qph5['nuclei'].attrs['nucl_num']=Nk*natom
    qph5['nuclei'].attrs['kpt_num']=Nk
    qph5.create_group('ao_two_e_ints')
    qph5['ao_two_e_ints'].attrs['df_num']=naux

    qph5['ao_basis'].attrs['ao_basis']=mf.cell.basis
    ao_nucl=[mf.cell.bas_atom(i)+1 for i in range(nao)]
    qph5.create_dataset('ao_basis/ao_nucl',data=Nk*ao_nucl)


    #                             _                             
    # |\ |      _ |  _   _. ._   |_)  _  ._      |  _ o  _  ._  
    # | \| |_| (_ | (/_ (_| |    | \ (/_ |_) |_| | _> | (_) | | 
    #                                    |                      
    
    #Total energy shift due to Ewald probe charge = -1/2 * Nelec*madelung/cell.vol =
    shift = tools.pbc.madelung(cell, kpts)*cell.nelectron * -.5 
    e_nuc = (cell.energy_nuc() + shift)*Nk
  
    print('nucl_repul', e_nuc)
    qph5['nuclei'].attrs['nuclear_repulsion']=e_nuc
  
    #       __    __          _                                 
    # |\/| |  |  |    _   _  |_  _ 
    # |  | |__|  |__ (_) (/_ |  _> 
    #                                               
    mo_coef_blocked=scipy.linalg.block_diag(*mo_k)
    qph5.create_dataset('mo_basis/mo_coef_real',data=mo_coef_blocked.real)
    qph5.create_dataset('mo_basis/mo_coef_imag',data=mo_coef_blocked.imag)
    qph5.create_dataset('mo_basis/mo_coef_kpts_real',data=mo_k.real)
    qph5.create_dataset('mo_basis/mo_coef_kpts_imag',data=mo_k.imag)

    # ___                                              
    #  |  ._ _|_  _   _  ._ _. |  _   |\/|  _  ._   _  
    # _|_ | | |_ (/_ (_| | (_| | _>   |  | (_) | | (_) 
    #                 _|                              
   
    if mf.cell.pseudo:
        v_kpts_ao = np.reshape(mf.with_df.get_pp(kpts=kpts),(Nk,nao,nao))
    else:
        v_kpts_ao = np.reshape(mf.with_df.get_nuc(kpts=kpts),(Nk,nao,nao))
    if len(cell._ecpbas) > 0:
        v_kpts_ao += np.reshape(ecp.ecp_int(cell, kpts),(Nk,nao,nao))
      
    ne_ao = ('V',v_kpts_ao,ne_threshold)
    ovlp_ao = ('S',np.reshape(mf.get_ovlp(cell=cell,kpts=kpts),(Nk,nao,nao)),ovlp_threshold)
    kin_ao = ('T',np.reshape(cell.pbc_intor('int1e_kin',1,1,kpts=kpts),(Nk,nao,nao)),kin_threshold)
    
    kin_ao_blocked=scipy.linalg.block_diag(*kin_ao[1])
    ovlp_ao_blocked=scipy.linalg.block_diag(*ovlp_ao[1])
    ne_ao_blocked=scipy.linalg.block_diag(*v_kpts_ao)
    qph5.create_dataset('ao_one_e_ints/ao_integrals_kinetic_real',data=kin_ao_blocked.real)
    qph5.create_dataset('ao_one_e_ints/ao_integrals_kinetic_imag',data=kin_ao_blocked.imag)
    qph5.create_dataset('ao_one_e_ints/ao_integrals_overlap_real',data=ovlp_ao_blocked.real)
    qph5.create_dataset('ao_one_e_ints/ao_integrals_overlap_imag',data=ovlp_ao_blocked.imag)
    qph5.create_dataset('ao_one_e_ints/ao_integrals_n_e_real',    data=ne_ao_blocked.real)
    qph5.create_dataset('ao_one_e_ints/ao_integrals_n_e_imag',    data=ne_ao_blocked.imag)
    


    for name, intval_kpts_ao, thresh in (ne_ao, ovlp_ao, kin_ao):
        if print_ao_ints_mono:
            with open('%s.qp' % name,'w') as outfile:
                for ik in range(Nk):
                    shift=ik*nao+1
                    for i in range(nao):
                        for j in range(i,nao):
                            int_ij = intval_kpts_ao[ik,i,j]
                            if abs(int_ij) > thresh:
                                outfile.write('%s %s %s %s\n' % (i+shift, j+shift, int_ij.real, int_ij.imag))
        if print_mo_ints_mono:
            intval_kpts_mo = np.einsum('kim,kij,kjn->kmn',mo_k.conj(),intval_kpts_ao,mo_k)
            with open('%s_mo.qp' % name,'w') as outfile:
                for ik in range(Nk):
                    shift=ik*nmo+1
                    for i in range(nmo):
                        for j in range(i,nmo):
                            int_ij = intval_kpts_mo[ik,i,j]
                            if abs(int_ij) > thresh:
                                outfile.write('%s %s %s %s\n' % (i+shift, j+shift, int_ij.real, int_ij.imag))

  
    # ___                              _    
    #  |  ._ _|_  _   _  ._ _. |  _   |_) o 
    # _|_ | | |_ (/_ (_| | (_| | _>   |_) | 
    #                 _|                    
    #
    kconserv = tools.get_kconserv(cell, kpts)
    qph5.create_dataset('nuclei/kconserv',data=np.transpose(kconserv+1,(0,2,1)))
    kcon_test = np.zeros((Nk,Nk,Nk),dtype=int)
    for a in range(Nk):
        for b in range(Nk):
            for c in range(Nk):
                kcon_test[a,c,b] = kconserv[a,b,c]+1
    qph5.create_dataset('nuclei/kconserv_test',data=kcon_test)

    
    with open('K.qp','w') as outfile:
        for a in range(Nk):
            for b in range(Nk):
                for c in range(Nk):
                    d = kconserv[a,b,c]
                    outfile.write('%s %s %s %s\n' % (a+1,c+1,b+1,d+1))
    

    intfile=h5py.File(mf.with_df._cderi,'r')

    j3c = intfile.get('j3c')
    naosq = nao*nao
    naotri = (nao*(nao+1))//2
    j3ckeys = list(j3c.keys())
    j3ckeys.sort(key=lambda strkey:int(strkey))

    # in new(?) version of PySCF, there is an extra layer of groups before the datasets
    # datasets used to be [/j3c/0,   /j3c/1,   /j3c/2,   ...]
    # datasets now are    [/j3c/0/0, /j3c/1/0, /j3c/2/0, ...]
    j3clist = [j3c.get(i+'/0') for i in j3ckeys]
    if j3clist==[None]*len(j3clist):
    # if using older version, stop before last level
        j3clist = [j3c.get(i) for i in j3ckeys]

    nkinvsq = 1./np.sqrt(Nk)

    # dimensions are (kikj,iaux,jao,kao), where kikj is compound index of kpts i and j
    # output dimensions should be reversed (nao, nao, naux, nkptpairs)
    j3arr=np.array([(i.value.reshape([-1,nao,nao]) if (i.shape[1] == naosq) else makesq3(i.value,nao)) * nkinvsq for i in j3clist])

    nkpt_pairs = j3arr.shape[0]
    df_ao_tmp = np.zeros((nao,nao,naux,nkpt_pairs),dtype=np.complex128)

    if print_ao_ints_df:
        with open('D.qp','w') as outfile:
            pass
        with open('D.qp','a') as outfile:
            for k,kpt_pair in enumerate(j3arr):
                for iaux,dfbasfunc in enumerate(kpt_pair):
                    for i,i0 in enumerate(dfbasfunc):
                        for j,v in enumerate(i0):
                            if (abs(v) > bielec_int_threshold):
                                outfile.write('%s %s %s %s %s %s\n' % (i+1,j+1,iaux+1,k+1,v.real,v.imag))
                                df_ao_tmp[i,j,iaux,k]=v
        
        qph5.create_dataset('ao_two_e_ints/df_ao_integrals_real',data=df_ao_tmp.real)
        qph5.create_dataset('ao_two_e_ints/df_ao_integrals_imag',data=df_ao_tmp.imag)

    if print_mo_ints_df:
        kpair_list=[]
        for i in range(Nk):
            for j in range(Nk):
                if(i>=j):
                    kpair_list.append((i,j,idx2_tri((i,j))))
        j3mo = np.array([np.einsum('mij,ik,jl->mkl',j3arr[kij],mo_k[ki].conj(),mo_k[kj]) for ki,kj,kij in kpair_list])
        df_mo_tmp = np.zeros((nmo,nmo,naux,nkpt_pairs),dtype=np.complex128)
        with open('D_mo.qp','w') as outfile:
            pass
        with open('D_mo.qp','a') as outfile:
            for k,kpt_pair in enumerate(j3mo):
                for iaux,dfbasfunc in enumerate(kpt_pair):
                    for i,i0 in enumerate(dfbasfunc):
                        for j,v in enumerate(i0):
                            if (abs(v) > bielec_int_threshold):
                                outfile.write('%s %s %s %s %s %s\n' % (i+1,j+1,iaux+1,k+1,v.real,v.imag))
                                df_mo_tmp[i,j,iaux,k]=v
        qph5.create_dataset('mo_two_e_ints/df_mo_integrals_real',data=df_mo_tmp.real)
        qph5.create_dataset('mo_two_e_ints/df_mo_integrals_imag',data=df_mo_tmp.imag)



#    eri_4d_ao = np.zeros((Nk,nao,Nk,nao,Nk,nao,Nk,nao), dtype=np.complex)
#    for d, kd in enumerate(kpts):
#        for c, kc in enumerate(kpts):
#            if c > d: break
#            idx2_cd = idx2_tri(c,d)
#            for b, kb in enumerate(kpts):
#                if b > d: break
#                a = kconserv[b,c,d]
#                if idx2_tri(a,b) > idx2_cd: continue
#                if ((c==d) and (a>b)): continue
#                ka = kpts[a]
#                v = mf.with_df.get_ao_eri(kpts=[ka,kb,kc,kd],compact=False).reshape((nao,)*4)
#                v *= 1./Nk
#                eri_4d_ao[a,:,b,:,c,:,d] = v
#    
#    eri_4d_ao = eri_4d_ao.reshape([Nk*nao]*4)


    if (print_ao_ints_bi or print_mo_ints_bi):
        if print_ao_ints_bi:
            with open('W.qp','w') as outfile: 
                pass
        if print_mo_ints_bi:
            with open('W_mo.qp','w') as outfile: 
                pass
        for d, kd in enumerate(kpts):
            for c, kc in enumerate(kpts):
                if c > d: break
                idx2_cd = idx2_tri((c,d))
                for b, kb in enumerate(kpts):
                    if b > d: break
                    a = kconserv[b,c,d]
                    if idx2_tri((a,b)) > idx2_cd: continue
                    if ((c==d) and (a>b)): continue
                    ka = kpts[a]

                    if print_ao_ints_bi:
                        with open('W.qp','a') as outfile:
                            eri_4d_ao_kpt = mf.with_df.get_ao_eri(kpts=[ka,kb,kc,kd],compact=False).reshape((nao,)*4)
                            eri_4d_ao_kpt *= 1./Nk
                            for l in range(nao):
                                ll=l+d*nao
                                for j in range(nao):
                                    jj=j+c*nao
                                    if jj>ll: break
                                    idx2_jjll = idx2_tri((jj,ll))
                                    for k in range(nao):
                                        kk=k+b*nao
                                        if kk>ll: break
                                        for i in range(nao):
                                            ii=i+a*nao
                                            if idx2_tri((ii,kk)) > idx2_jjll: break
                                            if ((jj==ll) and (ii>kk)): break
                                            v=eri_4d_ao_kpt[i,k,j,l]
                                            if (abs(v) > bielec_int_threshold):
                                                outfile.write('%s %s %s %s %s %s\n' % (ii+1,jj+1,kk+1,ll+1,v.real,v.imag))
            
                    if print_mo_ints_bi:
                        with open('W_mo.qp','a') as outfile:
                            eri_4d_mo_kpt = mf.with_df.ao2mo([mo_k[a], mo_k[b], mo_k[c], mo_k[d]],
                                                              [ka,kb,kc,kd],compact=False).reshape((nmo,)*4)
                            eri_4d_mo_kpt *= 1./Nk
                            for l in range(nmo):
                                ll=l+d*nmo
                                for j in range(nmo):
                                    jj=j+c*nmo
                                    if jj>ll: break
                                    idx2_jjll = idx2_tri((jj,ll))
                                    for k in range(nmo):
                                        kk=k+b*nmo
                                        if kk>ll: break
                                        for i in range(nmo):
                                            ii=i+a*nmo
                                            if idx2_tri((ii,kk)) > idx2_jjll: break
                                            if ((jj==ll) and (ii>kk)): break
                                            v=eri_4d_mo_kpt[i,k,j,l]
                                            if (abs(v) > bielec_int_threshold):
                                                outfile.write('%s %s %s %s %s %s\n' % (ii+1,jj+1,kk+1,ll+1,v.real,v.imag))

    
  
#def testpyscf2QP(cell,mf, kpts, kmesh=None, cas_idx=None, int_threshold = 1E-8):
#    '''
#    kpts = List of kpoints coordinates. Cannot be null, for gamma is other script
#    kmesh = Mesh of kpoints (optional)
#    cas_idx = List of active MOs. If not specified all MOs are actives
#    int_threshold = The integral will be not printed in they are bellow that
#    '''
#
#    from pyscf.pbc import ao2mo
#    from pyscf.pbc import tools
#    from pyscf.pbc.gto import ecp
#
#    mo_coef_threshold = int_threshold
#    ovlp_threshold = int_threshold
#    kin_threshold = int_threshold
#    ne_threshold = int_threshold
#    bielec_int_threshold = int_threshold
#
#    natom = len(cell.atom_coords())
#    print('n_atom per kpt',   natom)
#    print('num_elec per kpt', cell.nelectron)
#
#    mo_coeff = mf.mo_coeff
#    # Mo_coeff actif
#    mo_k = np.array([c[:,cas_idx] for c in mo_coeff] if cas_idx is not None else mo_coeff)
#    e_k =  np.array([e[cas_idx] for e in mf.mo_energy] if cas_idx is not None else mf.mo_energy)
#
#    Nk, nao, nmo = mo_k.shape
#    print("n Kpts", Nk)
#    print("n active Mos per kpt", nmo)
#    print("n AOs per kpt", nao)
#
#    naux = mf.with_df.get_naoaux()
#    print("n df fitting functions", naux)
#
#    #                             _
#    # |\ |      _ |  _   _. ._   |_)  _  ._      |  _ o  _  ._
#    # | \| |_| (_ | (/_ (_| |    | \ (/_ |_) |_| | _> | (_) | |
#    #                                    |
#
#    #Total energy shift due to Ewald probe charge = -1/2 * Nelec*madelung/cell.vol =
#    shift = tools.pbc.madelung(cell, kpts)*cell.nelectron * -.5
#    e_nuc = (cell.energy_nuc() + shift)*Nk
#
#    print('nucl_repul', e_nuc)
#
#
#    # ___
#    #  |  ._ _|_  _   _  ._ _. |  _   |\/|  _  ._   _  
#    # _|_ | | |_ (/_ (_| | (_| | _>   |  | (_) | | (_) 
#    #                 _|                              
#   
#    if mf.cell.pseudo:
#        v_kpts_ao = np.reshape(mf.with_df.get_pp(kpts=kpts),(Nk,nao,nao))
#    else:
#        v_kpts_ao = np.reshape(mf.with_df.get_nuc(kpts=kpts),(Nk,nao,nao))
#    if len(cell._ecpbas) > 0:
#        v_kpts_ao += np.reshape(ecp.ecp_int(cell, kpts),(Nk,nao,nao))
#
#    ne_ao = ('ne',v_kpts_ao,ne_threshold)
#    ovlp_ao = ('overlap',np.reshape(mf.get_ovlp(cell=cell,kpts=kpts),(Nk,nao,nao)),ovlp_threshold)
#    kin_ao = ('kinetic',np.reshape(cell.pbc_intor('int1e_kin',1,1,kpts=kpts),(Nk,nao,nao)),kin_threshold)
#
#  
#    # ___                              _    
#    #  |  ._ _|_  _   _  ._ _. |  _   |_) o 
#    # _|_ | | |_ (/_ (_| | (_| | _>   |_) | 
#    #                 _|                    
#    #
#    kconserv = tools.get_kconserv(cell, kpts)
#
#    
#    import h5py
#
#    intfile=h5py.File(mf.with_df._cderi,'r')
#
#    j3c = intfile.get('j3c')
#    naosq = nao*nao
#    naotri = (nao*(nao+1))//2
#    j3keys = list(j3c.keys())
#    j3keys.sort(key=lambda x:int(x))
#    j3clist = [j3c.get(i) for i in j3keys]
#    nkinvsq = 1./np.sqrt(Nk)
#    
#    # dimensions are (kikj,iaux,jao,kao), where kikj is compound index of kpts i and j
#    # output dimensions should be reversed (nao, nao, naux, nkptpairs)
#    j3arr=np.array([(pad(i.value.reshape([-1,nao,nao]),[naux,nao,nao]) if (i.shape[1] == naosq) else makesq(i.value,naux,nao)) * nkinvsq for i in j3clist])
#
#    nkpt_pairs = j3arr.shape[0]
#
#    kpair_list=[]
#    for i in range(Nk):
#        for j in range(Nk):
#            if(i>=j):
#                kpair_list.append((i,j,idx2_tri((i,j))))
#    j3mo = np.array([np.einsum('mij,ik,jl->mkl',j3arr[kij,:,:,:],mo_k[ki,:,:].conj(),mo_k[kj,:,:]) for ki,kj,kij in kpair_list])
#
#
#
#    eri_mo = np.zeros(4*[nmo*Nk],dtype=np.complex128)
#    eri_ao = np.zeros(4*[nao*Nk],dtype=np.complex128)
#
#    for d, kd in enumerate(kpts):
#        for c, kc in enumerate(kpts):
#            for b, kb in enumerate(kpts):
#                a = kconserv[b,c,d]
#                ka = kpts[a]
#                eri_4d_ao_kpt = mf.with_df.get_ao_eri(kpts=[ka,kb,kc,kd],compact=False).reshape((nao,)*4)
#                eri_4d_ao_kpt *= 1./Nk
#                for l in range(nao):
#                    ll=l+d*nao
#                    for j in range(nao):
#                        jj=j+c*nao
#                        for k in range(nao):
#                            kk=k+b*nao
#                            for i in range(nao):
#                                ii=i+a*nao
#                                v=eri_4d_ao_kpt[i,k,j,l]
#                                eri_ao[ii,kk,jj,ll]=v
#        
#                eri_4d_mo_kpt = mf.with_df.ao2mo([mo_k[a], mo_k[b], mo_k[c], mo_k[d]],
#                                                  [ka,kb,kc,kd],compact=False).reshape((nmo,)*4)
#                eri_4d_mo_kpt *= 1./Nk
#                for l in range(nmo):
#                    ll=l+d*nmo
#                    for j in range(nmo):
#                        jj=j+c*nmo
#                        for k in range(nmo):
#                            kk=k+b*nmo
#                            for i in range(nmo):
#                                ii=i+a*nmo
#                                v=eri_4d_mo_kpt[i,k,j,l]
#                                eri_mo[ii,kk,jj,ll]=v
#        
#    return (mo_k,j3arr,j3mo,eri_ao,eri_mo,kpair_list)
    
  
