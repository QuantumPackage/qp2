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

def xyzcount(s):
    return list(map(s.count,['x','y','z']))

def idx40(i,j,k,l):
    return idx2_tri((idx2_tri((i,k)),idx2_tri((j,l))))

def idx4(i,j,k,l):
    return idx2_tri((idx2_tri((i-1,k-1)),idx2_tri((j-1,l-1))))+1

def stri4(i,j,k,l):
    return (4*'{:5d}').format(i,j,k,l)

def stri4z(i,j,k,l,zr,zi):
    return (4*'{:5d}'+2*'{:25.16e}').format(i,j,k,l,zr,zi)

def stri2z(i,j,zr,zi):
    return (2*'{:5d}'+2*'{:25.16e}').format(i,j,zr,zi)

def strijklikjli4z(i,j,k,l,zr,zi):
    return ('{:10d}'+ 2*'{:8d}'+4*'{:5d}'+2*'{:25.16e}').format(idx4(i,j,k,l),idx2_tri((i-1,k-1))+1,idx2_tri((j-1,l-1))+1,i,j,k,l,zr,zi)


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

def print_mo_bi(mf,kconserv=None,outfilename='W.mo.qp',cas_idx=None,bielec_int_threshold = 1E-8):

    cell = mf.cell
    kpts = mf.kpts
    #nao = mf.cell.nao
    #Nk = kpts.shape[0]
    
    mo_coeff = mf.mo_coeff
    # Mo_coeff actif
    mo_k = np.array([c[:,cas_idx] for c in mo_coeff] if cas_idx is not None else mo_coeff)
  
    Nk, nao, nmo = mo_k.shape

    if (kconserv is None):
        from pyscf.pbc import tools
        kconserv = tools.get_kconserv(cell, kpts)

    with open(outfilename,'w') as outfile:
        for d, kd in enumerate(kpts):
            for c, kc in enumerate(kpts):
                if c > d: break
                #idx2_cd = idx2_tri((c,d))
                for b, kb in enumerate(kpts):
                    if b > d: break
                    a = kconserv[b,c,d]
                    if a > d: continue
                    #if idx2_tri((a,b)) > idx2_cd: continue
                    #if ((c==d) and (a>b)): continue
                    ka = kpts[a]
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
                                        outfile.write(stri4z(ii+1,jj+1,kk+1,ll+1,
                                                             v.real,v.imag)+'\n')


def print_ao_bi(mf,kconserv=None,outfilename='W.ao.qp',bielec_int_threshold = 1E-8):

    cell = mf.cell
    kpts = mf.kpts
    nao = mf.cell.nao
    Nk = kpts.shape[0]

    if (kconserv is None):
        from pyscf.pbc.tools import get_kconserv
        kconserv = get_kconserv(cell, kpts)

    with open(outfilename,'w') as outfile:
        for d, kd in enumerate(kpts):
            for c, kc in enumerate(kpts):
                if c > d: break
                #idx2_cd = idx2_tri((c,d))
                for b, kb in enumerate(kpts):
                    if b > d: break
                    a = kconserv[b,c,d]
                    if a > d: continue
                    #if idx2_tri((a,b)) > idx2_cd: continue
                    #if ((c==d) and (a>b)): continue
                    ka = kpts[a]

                    eri_4d_ao_kpt = mf.with_df.get_ao_eri(kpts=[ka,kb,kc,kd],
                                                          compact=False).reshape((nao,)*4)
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
                                        outfile.write(stri4z(ii+1,jj+1,kk+1,ll+1,
                                                             v.real,v.imag)+'\n')


def print_kcon_chem_to_phys(kcon,fname):
    '''
    input: kconserv in chem notation kcon_c[a,b,c] = d
           where (ab|cd) is allowed by symmetry
    output: kconserv in phys notation kcon_p[i,j,k] = l
            where <ij|kl> is allowed by symmetry
            (printed to file)
    '''
    Nk,n2,n3 = kcon.shape
    if (n2!=n3 or Nk!=n2):
        raise Exception('print_kcon_chem_to_phys called with non-cubic array')

    with open(fname,'w') as outfile:
        for a in range(Nk):
            for b in range(Nk):
                for c in range(Nk):
                    d = kcon[a,b,c]
                    outfile.write(stri4(a+1,c+1,b+1,d+1)+'\n')
    
def print_kpts_unblocked(ints_k,outfilename,thresh):
    '''
    for ints_k of shape (Nk,n1,n2),
    print the elements of the corresponding block-diagonal matrix of shape (Nk*n1,Nk*n2) in file
    '''
    Nk,n1,n2 = ints_k.shape
    with open(outfilename,'w') as outfile:
        for ik in range(Nk):
            shift1 = ik*n1+1
            shift2 = ik*n2+1
            for i1 in range(n1):
                for i2 in range(n2):
                    int_ij = ints_k[ik,i1,i2]
                    if abs(int_ij) > thresh:
                        outfile.write(stri2z(i1+shift1, i2+shift2, int_ij.real, int_ij.imag)+'\n')
    return

def print_kpts_unblocked_upper(ints_k,outfilename,thresh):
    '''
    for hermitian ints_k of shape (Nk,n1,n1),
    print the elements of the corresponding block-diagonal matrix of shape (Nk*n1,Nk*n1) in file
    (only upper triangle is printed)
    '''
    Nk,n1,n2 = ints_k.shape
    if (n1!=n2):
        raise Exception('print_kpts_unblocked_upper called with non-square matrix')

    with open(outfilename,'w') as outfile:
        for ik in range(Nk):
            shift = ik*n1+1
            for i1 in range(n1):
                for i2 in range(i1,n1):
                    int_ij = ints_k[ik,i1,i2]
                    if abs(int_ij) > thresh:
                        outfile.write(stri2z(i1+shift, i2+shift, int_ij.real, int_ij.imag)+'\n')
    return



def get_kin_ao(mf):
    nao = mf.cell.nao_nr()
    Nk = len(mf.kpts)
    return np.reshape(mf.cell.pbc_intor('int1e_kin',1,1,kpts=mf.kpts),(Nk,nao,nao))

def get_ovlp_ao(mf):
    nao = mf.cell.nao_nr()
    Nk = len(mf.kpts)
    return np.reshape(mf.get_ovlp(cell=mf.cell,kpts=mf.kpts),(Nk,nao,nao))

def get_pot_ao(mf):
    nao = mf.cell.nao_nr()
    Nk = len(mf.kpts)

    if mf.cell.pseudo:
        v_kpts_ao = np.reshape(mf.with_df.get_pp(kpts=mf.kpts),(Nk,nao,nao))
    else:
        v_kpts_ao = np.reshape(mf.with_df.get_nuc(kpts=mf.kpts),(Nk,nao,nao))

    if len(mf.cell._ecpbas) > 0:
        from pyscf.pbc.gto import ecp
        v_kpts_ao += np.reshape(ecp.ecp_int(mf.cell, mf.kpts),(Nk,nao,nao))

    return v_kpts_ao

def ao_to_mo_1e(ao_kpts,mo_coef):
    return np.einsum('kim,kij,kjn->kmn',mo_coef.conj(),ao_kpts,mo_coef)

def get_j3ao_old(fname,nao,Nk):
    '''
    returns list of Nk_pair arrays of shape (naux,nao,nao)
    if naux is the same for each pair, returns numpy array
    if naux is not the same for each pair, returns array of arrays
    '''
    import h5py
    with h5py.File(fname,'r') as intfile:
        j3c = intfile.get('j3c')
        j3ckeys = list(j3c.keys())
        j3ckeys.sort(key=lambda strkey:int(strkey))
    
        # in new(?) version of PySCF, there is an extra layer of groups before the datasets
        # datasets used to be [/j3c/0,   /j3c/1,   /j3c/2,   ...]
        # datasets now are    [/j3c/0/0, /j3c/1/0, /j3c/2/0, ...]
        j3clist = [j3c.get(i+'/0') for i in j3ckeys]
        #if j3clist==[None]*len(j3clist):
        if not(any(j3clist)):
        # if using older version, stop before last level
            j3clist = [j3c.get(i) for i in j3ckeys]
    
        naosq = nao*nao
        naotri = (nao*(nao+1))//2
        nkinvsq = 1./np.sqrt(Nk)
    
        # dimensions are (kikj,iaux,jao,kao), where kikj is compound index of kpts i and j
        # output dimensions should be reversed (nao, nao, naux, nkptpairs)
        return np.array([(i.value.reshape([-1,nao,nao]) if (i.shape[1] == naosq) else makesq3(i.value,nao)) * nkinvsq for i in j3clist])

def get_j3ao(fname,nao,Nk):
    '''
    returns padded df AO array
    fills in zeros when functions are dropped due to linear dependency
    last AO index corresponds to smallest kpt index?
    (k, mu, i, j) where i.kpt >= j.kpt
    '''
    import h5py
    with h5py.File(fname,'r') as intfile:
        j3c = intfile.get('j3c')
        j3ckeys = list(j3c.keys())
        nkpairs = len(j3ckeys)

        # get num order instead of lex order
        j3ckeys.sort(key=lambda strkey:int(strkey))

        # in new(?) version of PySCF, there is an extra layer of groups before the datasets
        # datasets used to be [/j3c/0,   /j3c/1,   /j3c/2,   ...]
        # datasets now are    [/j3c/0/0, /j3c/1/0, /j3c/2/0, ...]
        keysub = '/0' if bool(j3c.get('0/0',getclass=True)) else ''

        naux = max(map(lambda k: j3c[k+keysub].shape[0],j3c.keys()))

        naosq = nao*nao
        naotri = (nao*(nao+1))//2
        nkinvsq = 1./np.sqrt(Nk)

        j3arr = np.zeros((nkpairs,naux,nao,nao),dtype=np.complex128)

        for i,kpair in enumerate(j3ckeys):
            iaux,dim2 = j3c[kpair+keysub].shape
            if (dim2==naosq):
                j3arr[i,:iaux,:,:] = j3c[kpair+keysub][()].reshape([iaux,nao,nao]) * nkinvsq
                #j3arr[i,:iaux,:,:] = j3c[kpair+keysub][()].reshape([iaux,nao,nao]).transpose((0,2,1)) * nkinvsq
            else:
                j3arr[i,:iaux,:,:] = makesq3(j3c[kpair+keysub][()],nao) * nkinvsq
                #j3arr[i,:iaux,:,:] = makesq3(j3c[kpair+keysub][()].conj(),nao) * nkinvsq

        return j3arr

def get_j3ao_new(fname,nao,Nk):
    '''
    returns padded df AO array
    fills in zeros when functions are dropped due to linear dependency
    last AO index corresponds to largest kpt index?
    (k, mu, j, i) where i.kpt >= j.kpt
    '''
    import h5py
    with h5py.File(fname,'r') as intfile:
        j3c = intfile.get('j3c')
        j3ckeys = list(j3c.keys())
        nkpairs = len(j3ckeys)

        # get num order instead of lex order
        j3ckeys.sort(key=lambda strkey:int(strkey))

        # in new(?) version of PySCF, there is an extra layer of groups before the datasets
        # datasets used to be [/j3c/0,   /j3c/1,   /j3c/2,   ...]
        # datasets now are    [/j3c/0/0, /j3c/1/0, /j3c/2/0, ...]
        keysub = '/0' if bool(j3c.get('0/0',getclass=True)) else ''

        naux = max(map(lambda k: j3c[k+keysub].shape[0],j3c.keys()))

        naosq = nao*nao
        naotri = (nao*(nao+1))//2
        nkinvsq = 1./np.sqrt(Nk)

        j3arr = np.zeros((nkpairs,naux,nao,nao),dtype=np.complex128)

        for i,kpair in enumerate(j3ckeys):
            iaux,dim2 = j3c[kpair+keysub].shape
            if (dim2==naosq):
                j3arr[i,:iaux,:,:] = j3c[kpair+keysub][()].reshape([iaux,nao,nao]).transpose((0,2,1)) * nkinvsq
            else:
                j3arr[i,:iaux,:,:] = makesq3(j3c[kpair+keysub][()].conj(),nao) * nkinvsq

        return j3arr

def print_df(j3arr,fname,thresh):
    with open(fname,'w') as outfile:
        for k,kpt_pair in enumerate(j3arr):
            for iaux,dfbasfunc in enumerate(kpt_pair):
                for i,i0 in enumerate(dfbasfunc):
                    for j,v in enumerate(i0):
                        if (abs(v) > thresh):
                            outfile.write(stri4z(i+1,j+1,iaux+1,k+1,v.real,v.imag)+'\n')
    return

def df_pad_ref_test(j3arr,nao,naux,nkpt_pairs):
    df_ao_tmp = np.zeros((nao,nao,naux,nkpt_pairs),dtype=np.complex128)
    for k,kpt_pair in enumerate(j3arr):
        for iaux,dfbasfunc in enumerate(kpt_pair):
            for i,i0 in enumerate(dfbasfunc):
                for j,v in enumerate(i0):
                    df_ao_tmp[i,j,iaux,k]=v
    return df_ao_tmp


def df_ao_to_mo(j3ao,mo_coef):
    from itertools import product
    Nk = mo_coef.shape[0]
    kpair_list = ((i,j,idx2_tri((i,j))) for (i,j) in product(range(Nk),repeat=2) if (i>=j))
    return np.array([
        np.einsum('mij,ik,jl->mkl',j3ao[kij],mo_coef[ki].conj(),mo_coef[kj])
        for ki,kj,kij in kpair_list])

    
def df_ao_to_mo_new(j3ao,mo_coef):
    #TODO: fix this (C/F ordering, conj, transpose, view cmplx->float)

    from itertools import product
    Nk = mo_coef.shape[0]
    return np.array([
        np.einsum('mji,ik,jl->mlk',j3ao[idx2_tri((ki,kj))],mo_coef[ki].conj(),mo_coef[kj])
        for ki,kj in product(range(Nk),repeat=2) if (ki>=kj)],dtype=np.complex128)

def df_ao_to_mo_test(j3ao,mo_coef):
    from itertools import product
    Nk = mo_coef.shape[0]
    return np.array([
        np.einsum('mij,ik,jl->mkl',j3ao[idx2_tri((ki,kj))],mo_coef[ki].conj(),mo_coef[kj])
        for ki,kj in product(range(Nk),repeat=2) if (ki>=kj)])

def pyscf2QP2_mo(cell,mf,kpts,kmesh=None,cas_idx=None, int_threshold = 1E-8,qph5path='qpdat.h5'):
    pyscf2QP2(cell,mf,kpts,kmesh,cas_idx,int_threshold,qph5path,
            print_ao_ints_df=False,
            print_mo_ints_df=True,
            print_ao_ints_mono=False,
            print_mo_ints_mono=True)
    return



def pyscf2QP2(cell,mf, kpts, kmesh=None, cas_idx=None, int_threshold = 1E-8, 
        qph5path = 'qpdat.h5', sp_twist=None,
        print_ao_ints_bi=False, 
        print_mo_ints_bi=False, 
        print_ao_ints_df=True, 
        print_mo_ints_df=False, 
        print_ao_ints_mono=True, 
        print_mo_ints_mono=False,
        print_debug=False):
    '''
    kpts = List of kpoints coordinates. Cannot be null, for gamma is other script
    kmesh = Mesh of kpoints (optional)
    cas_idx = List of active MOs. If not specified all MOs are actives
    int_threshold = The integral will be not printed in they are bellow that
    '''
  
#    from pyscf.pbc import ao2mo
    from pyscf.pbc import tools
    import h5py
#    import scipy
    from scipy.linalg import block_diag

    mo_coef_threshold = int_threshold
    ovlp_threshold = int_threshold
    kin_threshold = int_threshold
    ne_threshold = int_threshold
    bielec_int_threshold = int_threshold
    thresh_mono = int_threshold
    
    
#    qph5path = 'qpdat.h5'
    # create hdf5 file, delete old data if exists
    with h5py.File(qph5path,'w') as qph5:
        qph5.create_group('nuclei')
        qph5.create_group('electrons')
        qph5.create_group('ao_basis')
        qph5.create_group('mo_basis')
        qph5.create_group('qmcpack')
        #qph5.create_group('pseudo')
        #qph5['pseudo'].attrs['do_pseudo']=False

    #if mf.cell.cart:
    #    mo_coeff = mf.mo_coeff.copy()
    #else:
    #    # normalized can be one of ['all','sp',None]
    #    # we can either normalize here or after qp
    #    c2s = mf.cell.cart2sph_coeff(normalized='sp')
    #    mo_coeff = list(map(lambda x: np.dot(c2s,x),mf.mo_coeff))
    mo_coeff = mf.mo_coeff
    # Mo_coeff actif
    mo_k = np.array([c[:,cas_idx] for c in mo_coeff] if cas_idx is not None else mo_coeff)
    e_k =  np.array([e[cas_idx] for e in mf.mo_energy] if cas_idx is not None else mf.mo_energy)
  
    Nk, nao, nmo = mo_k.shape

    print("n Kpts", Nk)
    print("n active Mos per kpt", nmo)
    print("n AOs per kpt", nao)

    ##########################################
    #                                        #
    #                Nuclei                  #
    #                                        #
    ##########################################

    natom = cell.natm
    print('n_atom per kpt',   natom)

    atom_xyz = mf.cell.atom_coords()
    if not(mf.cell.unit.startswith(('B','b','au','AU'))):
        from pyscf.data.nist import BOHR
        atom_xyz /= BOHR # always convert to au

    with h5py.File(qph5path,'a') as qph5:
        qph5['nuclei'].attrs['kpt_num']=Nk
        qph5['nuclei'].attrs['nucl_num']=natom
        qph5.create_dataset('nuclei/nucl_coord',data=atom_xyz)
        qph5.create_dataset('nuclei/nucl_charge',data=mf.cell.atom_charges())

        strtype=h5py.special_dtype(vlen=str)
        atom_dset=qph5.create_dataset('nuclei/nucl_label',(natom,),dtype=strtype)
        for i in range(natom):
            atom_dset[i] = mf.cell.atom_pure_symbol(i)

    ##########################################
    #                                        #
    #                Basis                   #
    #                                        #
    ##########################################

    # nucleus on which each AO is centered
    ao_nucl=[i[0] for i in mf.cell.ao_labels(fmt=False,base=1)]


    nprim_max = 0
    nfunc_tot = 0
    for iatom, (sh0,sh1,ao0,ao1) in enumerate(cell.aoslice_by_atom()):
        for ib in range(sh0,sh1): # sets of contracted exponents
            nprim = cell.bas_nprim(ib)
            nfunc_tot += cell.bas_nctr(ib)
            if (nprim > nprim_max):
                nprim_max = nprim

    qp_prim_num = np.zeros((nfunc_tot),dtype=int)
    qp_coef = np.zeros((nfunc_tot,nprim_max))
    qp_expo = np.zeros((nfunc_tot,nprim_max))
    qp_nucl = np.zeros((nfunc_tot),dtype=int)
    qp_lbas = np.zeros((nfunc_tot),dtype=int)

    clabels = cell.cart_labels(fmt=False)

    tmp_idx=0
    for iatom, (sh0,sh1,ao0,ao1) in enumerate(cell.aoslice_by_atom()):
        # shell start,end; AO start,end (sph or cart) for each atom
        for ib in range(sh0,sh1): # sets of contracted exponents
            l = cell.bas_angular(ib)    # angular momentum
            nprim = cell.bas_nprim(ib)  # numer of primitives
            es = cell.bas_exp(ib)       # exponents
            cs = cell.bas_ctr_coeff(ib) # coeffs
            nctr = cell.bas_nctr(ib)    # number of contractions
            print(iatom,ib,l,nprim,nctr,tmp_idx)
            #if cell.cart:
            #    nfuncmax = ((l+1)*(l+2))//2
            #else:
            #    nfuncmax = 2*l+1
            for ic in range(nctr): # sets of contraction coeffs
                qp_expo[tmp_idx,:nprim] = es[:]
                qp_coef[tmp_idx,:nprim] = cs[:,ic]
                qp_nucl[tmp_idx] = iatom 
                qp_lbas[tmp_idx] = l
                qp_prim_num[tmp_idx] = nprim
                tmp_idx += 1
                #for nfunc in range(nfuncmax):
                #for nfunc in range(((l+1)*(l+2))//2): # always use cart for qp ao basis?
                #    qp_expo[tmp_idx,:nprim] = es[:]
                #    qp_coef[tmp_idx,:nprim] = cs[:,ic]
                #    qp_nucl[tmp_idx] = iatom + 1
                #    qp_pwr[tmp_idx,:] = xyzcount(clabels[tmp_idx][3])
                #    qp_prim_num[tmp_idx] = nprim
                #    tmp_idx += 1

    with h5py.File(qph5path,'a') as qph5:
        qph5['mo_basis'].attrs['mo_num']=Nk*nmo
        qph5['ao_basis'].attrs['ao_num']=Nk*nao

        #qph5['ao_basis'].attrs['ao_basis']=mf.cell.basis
        qph5['ao_basis'].attrs['ao_basis']="dummy basis"

        qph5.create_dataset('ao_basis/ao_nucl',data=Nk*ao_nucl)

        qph5['qmcpack'].attrs['qmc_nshell']=nfunc_tot
        qph5['qmcpack'].attrs['qmc_prim_num_max']=nprim_max
        qph5.create_dataset('qmcpack/qmc_nucl',data=qp_nucl)
        qph5.create_dataset('qmcpack/qmc_prim_num',data=qp_prim_num)
        qph5.create_dataset('qmcpack/qmc_expo',data=qp_expo.T)
        qph5.create_dataset('qmcpack/qmc_coef',data=qp_coef.T)
        qph5.create_dataset('qmcpack/qmc_lbas',data=qp_lbas)


#    with h5py.File(qph5path,'a') as qph5:
#        qph5['mo_basis'].attrs['mo_num']=Nk*nmo
#        qph5['ao_basis'].attrs['ao_num']=Nk*nao
#
#        #qph5['ao_basis'].attrs['ao_basis']=mf.cell.basis
#        qph5['ao_basis'].attrs['ao_basis']="dummy basis"
#
#        qph5.create_dataset('ao_basis/ao_nucl',data=Nk*ao_nucl)

    ##########################################
    #                                        #
    #              Electrons                 #
    #                                        #
    ##########################################

    nelec = cell.nelectron
    neleca,nelecb = cell.nelec

    print('num_elec per kpt', nelec)

    with h5py.File(qph5path,'a') as qph5:
        #in old version: param << nelec*Nk, nmo*Nk, natom*Nk
        qph5['electrons'].attrs['elec_alpha_num']=neleca*Nk 
        qph5['electrons'].attrs['elec_beta_num']=nelecb*Nk

    ##########################################
    #                                        #
    #           Nuclear Repulsion            #
    #                                        #
    ##########################################

    #Total energy shift due to Ewald probe charge = -1/2 * Nelec*madelung/cell.vol =
    shift = tools.pbc.madelung(cell, kpts)*cell.nelectron * -.5 
    e_nuc = (cell.energy_nuc() + shift)*Nk
  
    print('nucl_repul', e_nuc)

    with h5py.File(qph5path,'a') as qph5:
        qph5['nuclei'].attrs['nuclear_repulsion']=e_nuc
  
    ##########################################
    #                                        #
    #               MO Coef                  #
    #                                        #
    ##########################################
    
    with h5py.File(qph5path,'a') as qph5:
        # k,mo,ao(,2)
        mo_coef_f = np.array(mo_k.transpose((0,2,1)),order='c',dtype=np.complex128)
        #mo_coef_blocked=block_diag(*mo_k)
        mo_coef_blocked_f = block_diag(*mo_coef_f)
        qph5.create_dataset('mo_basis/mo_coef_complex',data=mo_coef_blocked_f.view(dtype=np.float64).reshape((Nk*nmo,Nk*nao,2)))
        qph5.create_dataset('mo_basis/mo_coef_kpts',data=mo_coef_f.view(dtype=np.float64).reshape((Nk,nmo,nao,2)))
   
    if print_debug:
        print_kpts_unblocked(mo_k,'C.qp',mo_coef_threshold)   

    ##########################################
    #                                        #
    #            Integrals Mono              #
    #                                        #
    ##########################################
    
    ne_ao = get_pot_ao(mf)
    kin_ao = get_kin_ao(mf)
    ovlp_ao = get_ovlp_ao(mf)

    if print_ao_ints_mono:

        with h5py.File(qph5path,'a') as qph5:
            kin_ao_f =  np.array(kin_ao.transpose((0,2,1)),order='c',dtype=np.complex128)
            ovlp_ao_f = np.array(ovlp_ao.transpose((0,2,1)),order='c',dtype=np.complex128)
            ne_ao_f =   np.array(ne_ao.transpose((0,2,1)),order='c',dtype=np.complex128)

            qph5.create_dataset('ao_one_e_ints/ao_integrals_kinetic_kpts',data=kin_ao_f.view(dtype=np.float64).reshape((Nk,nao,nao,2)))
            qph5.create_dataset('ao_one_e_ints/ao_integrals_overlap_kpts',data=ovlp_ao_f.view(dtype=np.float64).reshape((Nk,nao,nao,2)))
            qph5.create_dataset('ao_one_e_ints/ao_integrals_n_e_kpts',    data=ne_ao_f.view(dtype=np.float64).reshape((Nk,nao,nao,2)))

        if print_debug:
            for fname,ints in zip(('S.qp','V.qp','T.qp'),
                                  (ovlp_ao, ne_ao, kin_ao)):
                print_kpts_unblocked_upper(ints,fname,thresh_mono)

    if print_mo_ints_mono:
        kin_mo = ao_to_mo_1e(kin_ao,mo_k)
        ovlp_mo = ao_to_mo_1e(ovlp_ao,mo_k)
        ne_mo = ao_to_mo_1e(ne_ao,mo_k)

        with h5py.File(qph5path,'a') as qph5:
            kin_mo_f =  np.array(kin_mo.transpose((0,2,1)),order='c',dtype=np.complex128)
            ovlp_mo_f = np.array(ovlp_mo.transpose((0,2,1)),order='c',dtype=np.complex128)
            ne_mo_f =   np.array(ne_mo.transpose((0,2,1)),order='c',dtype=np.complex128)

            qph5.create_dataset('mo_one_e_ints/mo_integrals_kinetic_kpts',data=kin_mo_f.view(dtype=np.float64).reshape((Nk,nmo,nmo,2)))
            qph5.create_dataset('mo_one_e_ints/mo_integrals_overlap_kpts',data=ovlp_mo_f.view(dtype=np.float64).reshape((Nk,nmo,nmo,2)))
            qph5.create_dataset('mo_one_e_ints/mo_integrals_n_e_kpts',    data=ne_mo_f.view(dtype=np.float64).reshape((Nk,nmo,nmo,2)))
        if print_debug:
            for fname,ints in zip(('S.mo.qp','V.mo.qp','T.mo.qp'),
                                  (ovlp_mo, ne_mo, kin_mo)):
                print_kpts_unblocked_upper(ints,fname,thresh_mono)

  
    ##########################################
    #                                        #
    #               k-points                 #
    #                                        #
    ##########################################

    kconserv = tools.get_kconserv(cell, kpts)

    with h5py.File(qph5path,'a') as qph5:
        kcon_f_phys = np.array(kconserv.transpose((1,2,0)),order='c')
        qph5.create_dataset('nuclei/kconserv',data=kcon_f_phys+1)
    
    if print_debug:
        print_kcon_chem_to_phys(kconserv,'K.qp')
  
    ##########################################
    #                                        #
    #             Integrals Bi               #
    #                                        #
    ##########################################

    j3ao_new = get_j3ao_new(mf.with_df._cderi,nao,Nk)

    # test? nkpt_pairs should be (Nk*(Nk+1))//2
    nkpt_pairs, naux, _, _ = j3ao_new.shape

    print("n df fitting functions", naux)
    with h5py.File(qph5path,'a') as qph5:
        qph5.create_group('ao_two_e_ints')
        qph5['ao_two_e_ints'].attrs['df_num']=naux

    if print_ao_ints_df:
        if print_debug:
            print_df(j3ao_new,'D.qp',bielec_int_threshold)

        with h5py.File(qph5path,'a') as qph5:
            qph5.create_dataset('ao_two_e_ints/df_ao_integrals',data=j3ao_new.view(dtype=np.float64).reshape((nkpt_pairs,naux,nao,nao,2)))

    if print_mo_ints_df:

        j3mo_new = df_ao_to_mo_new(j3ao_new,mo_k)

        if print_debug:
            print_df(j3mo_new,'D.mo.qp',bielec_int_threshold)

        with h5py.File(qph5path,'a') as qph5:
            qph5.create_dataset('mo_two_e_ints/df_mo_integrals',data=j3mo_new.view(dtype=np.float64).reshape((nkpt_pairs,naux,nmo,nmo,2)))

    if (print_ao_ints_bi):
        print_ao_bi(mf,kconserv,'W.qp',bielec_int_threshold)
    if (print_mo_ints_bi):
        print_mo_bi(mf,kconserv,'W.mo.qp',cas_idx,bielec_int_threshold)


    ##########################################
    #                                        #
    #       PBC DATA To transfer to qmc      #
    #                                        #
    ##########################################

    if len(kpts)== 0:
        sp_twist=[0.0,0.0,0.0]


    with h5py.File(qph5path,'a') as qph5:
        qph5['qmcpack'].attrs['PBC']=True 
        qph5.create_dataset('qmcpack/Super_Twist',(1,3),dtype="f8",data=sp_twist)
        qph5.create_dataset('qmcpack/LatticeVectors',(3,3),dtype="f8",data=cell.lattice_vectors().T)
        qph5.create_dataset('qmcpack/eigenval',(1,Nk*nmo),dtype="f8",data=mf.mo_energy)
                                                                                                                      

#    ##########################################
#    #                                        #
#    #                 ECP                    #
#    #                                        #
#    ##########################################
#
#    if (cell.has_ecp()):
#        #atsymb = [mol.atom_pure_symbol(i) for i in range(natom)]
#        #pyecp = mol._ecp
#        ## nelec to remove for each atom
#        #nuc_z_remov = [pyecp[i][0] for i in atsymb]
#        #nl_per_atom = [len(pyecp[i][1]) for i in atsymb]
#        ## list of l-values for channels of each atom
#        #ecp_l = [[pyecp[i][1][j][0] for j in range(len(pyecp[i][1]))] for i in atsymb]
#        ## list of [exp,coef] for each channel (r**0,1,2,3,4,5,)
#        #ecp_ac = [[pyecp[i][1][j][1] for j in range(len(pyecp[i][1]))] for i in atsymb]
#        pyecp = [cell._ecp[cell.atom_pure_symbol(i)] for i in range(natom)]
#        nzrmv=[0]*natom
#        lmax=0
#        klocmax=0
#        knlmax=0
#        for i,(nz,dat) in enumerate(pyecp):
#            nzrmv[i]=nz
#            for lval,ac in dat:
#                if (lval==-1):
#                    klocmax=max(sum(len(j) for j in ac),klocmax)
#                else:
#                    lmax=max(lval,lmax)
#                    knlmax=max(sum(len(j) for j in ac),knlmax)
#        #psd_nk = np.zeros((natom,klocmax),dtype=int)
#        #psd_vk = np.zeros((natom,klocmax),dtype=float)
#        #psd_dzk = np.zeros((natom,klocmax),dtype=float)
#        #psd_nkl = np.zeros((natom,knlmax,lmax+1),dtype=int)
#        #psd_vkl = np.zeros((natom,knlmax,lmax+1),dtype=float)
#        #psd_dzkl = np.zeros((natom,knlmax,lmax+1),dtype=float)
#        klnlmax=max(klocmax,knlmax)
#        psd_n = np.zeros((lmax+2,klnlmax,natom),dtype=int)
#        psd_v = np.zeros((lmax+2,klnlmax,natom),dtype=float)
#        psd_dz = np.zeros((lmax+2,klnlmax,natom),dtype=float)
#        for i,(_,dat) in enumerate(pyecp):
#            for lval,ac in dat:
#                count=0
#                for ri,aici in enumerate(ac):
#                    for ai,ci in aici:
#                        psd_n[lval+1,count,i] = ri-2
#                        psd_v[lval+1,count,i] = ci
#                        psd_dz[lval+1,count,i] = ai
#                        count += 1
#        psd_nk = psd_n[0,:klocmax]
#        psd_vk = psd_v[0,:klocmax]
#        psd_dzk = psd_dz[0,:klocmax]
#        psd_nkl = psd_n[1:,:knlmax]
#        psd_vkl = psd_v[1:,:knlmax]
#        psd_dzkl = psd_dz[1:,:knlmax]
#        with h5py.File(qph5path,'a') as qph5:
#            qph5['pseudo'].attrs['do_pseudo']=True
#            qph5['pseudo'].attrs['pseudo_lmax']=lmax
#            qph5['pseudo'].attrs['pseudo_klocmax']=klocmax
#            qph5['pseudo'].attrs['pseudo_kmax']=knlmax
#            qph5.create_dataset('pseudo/nucl_charge_remove',data=nzrmv)
#            qph5.create_dataset('pseudo/pseudo_n_k',data=psd_nk)
#            qph5.create_dataset('pseudo/pseudo_n_kl',data=psd_nkl)
#            qph5.create_dataset('pseudo/pseudo_v_k',data=psd_vk)
#            qph5.create_dataset('pseudo/pseudo_v_kl',data=psd_vkl)
#            qph5.create_dataset('pseudo/pseudo_dz_k',data=psd_dzk)
#            qph5.create_dataset('pseudo/pseudo_dz_kl',data=psd_dzkl)
#
#        ## nelec to remove for each atom
#        #nuc_z_remov = [i[0] for i in pyecp]
#        #nl_per_atom = [len(i[1]) for i in pyecp]
#        ## list of l-values for channels of each atom
#        #ecp_l = [[ j[0] for j in i[1] ] for i in pyecp]
#        #lmax = max(map(max,ecp_l))
#        ## list of [exp,coef] for each channel (r**0,1,2,3,4,5,)
#        #ecp_ac = [[ j[1] for j in i[1] ] for i in pyecp]

    return

def pyscf2QP2_mol(mf, cas_idx=None, int_threshold = 1E-8, 
        qph5path = 'qpdat.h5',
        print_debug=False):
    '''
    cas_idx = List of active MOs. If not specified all MOs are actives
    int_threshold = The integral will be not printed in they are bellow that
    norm should be one of 'sp', 'all', or None
    '''
  
    import h5py

    norm='sp'
    mol = mf.mol
    nao_c = mol.nao_cart()

    mo_coef_threshold = int_threshold
    ovlp_threshold = int_threshold
    kin_threshold = int_threshold
    ne_threshold = int_threshold
    bielec_int_threshold = int_threshold
    thresh_mono = int_threshold
    
    
#    qph5path = 'qpdat.h5'
    # create hdf5 file, delete old data if exists
    with h5py.File(qph5path,'w') as qph5:
        qph5.create_group('nuclei')
        qph5.create_group('electrons')
        qph5.create_group('ao_basis')
        qph5.create_group('mo_basis')
        qph5.create_group('pseudo')
        qph5['pseudo'].attrs['do_pseudo']=False

    if mf.mol.cart:
        mo_coeff = mf.mo_coeff.copy()
    else:
        #c2s = mol.cart2sph_coeff(normalized=norm)
        c2s = mol.cart2sph_coeff(normalized='sp')
        #c2s = mol.cart2sph_coeff(normalized='all')
        #c2s = mol.cart2sph_coeff(normalized=None)
        mo_coeff = np.dot(c2s,mf.mo_coeff)
    #TODO: clean this up; use mol.cart_labels(fmt=False)
    dnormlbl1=["dxx","dyy","dzz"]
    dnormfac1 = 2.0*np.sqrt(np.pi/5)

    dnormlbl2=["dxy","dxz","dyz"]
    dnormfac2 = 2.0*np.sqrt(np.pi/15)

    fnormlbl1=["fxxx","fyyy","fzzz"]
    fnormfac1 = 2.0*np.sqrt(np.pi/7)

    fnormlbl2=["fxxy","fxxz","fxyy","fxzz","fyyz","fyzz"]
    fnormfac2 = 2.0*np.sqrt(np.pi/35)

    fnormlbl3=["fxyz"]
    fnormfac3 = 2.0*np.sqrt(np.pi/105)

    gnormlbl1=["gxxxx","gyyyy","gzzzz"]
    gnormfac1 = 2.0*np.sqrt(np.pi/9)

    gnormlbl2=["gxxxy","gxxxz","gxyyy","gxzzz","gyyyz","gyzzz"]
    gnormfac2 = 2.0*np.sqrt(np.pi/63)

    gnormlbl3=["gxxyy","gxxzz","gyyzz"]
    gnormfac3 = 2.0*np.sqrt(np.pi/105)

    gnormlbl4=["gxxyz","gxyyz","gxyzz"]
    gnormfac4 = 2.0*np.sqrt(np.pi/315)

    hnormlbl1=["hxxxxx","hyyyyy","hzzzzz"]
    hnormfac1 = 2.0*np.sqrt(np.pi/11)

    hnormlbl2=["hxxxxy","hxxxxz","hxyyyy","hxzzzz","hyyyyz","hyzzzz"]
    hnormfac2 = 2.0*np.sqrt(np.pi/99)

    hnormlbl3=["hxxxyy","hxxxzz","hxxyyy","hxxzzz","hyyyzz","hyyzzz"]
    hnormfac3 = 2.0*np.sqrt(np.pi/231)

    hnormlbl4=["hxxxyz","hxyyyz","hxyzzz"]
    hnormfac4 = 2.0*np.sqrt(np.pi/693)

    hnormlbl5=["hxxyyz","hxxyzz","hxyyzz"]
    hnormfac5 = 2.0*np.sqrt(np.pi/1155)

    for i_lbl,mo_lbl in enumerate(mol.cart_labels()):
        if any(i in mo_lbl for i in dnormlbl1):
            mo_coeff[i_lbl,:] *= dnormfac1
        elif any(i in mo_lbl for i in dnormlbl2):
            mo_coeff[i_lbl,:] *= dnormfac2
        elif any(i in mo_lbl for i in fnormlbl1):
            mo_coeff[i_lbl,:] *= fnormfac1
        elif any(i in mo_lbl for i in fnormlbl2):
            mo_coeff[i_lbl,:] *= fnormfac2
        elif any(i in mo_lbl for i in fnormlbl3):
            mo_coeff[i_lbl,:] *= fnormfac3
        elif any(i in mo_lbl for i in gnormlbl1):
            mo_coeff[i_lbl,:] *= gnormfac1
        elif any(i in mo_lbl for i in gnormlbl2):
            mo_coeff[i_lbl,:] *= gnormfac2
        elif any(i in mo_lbl for i in gnormlbl3):
            mo_coeff[i_lbl,:] *= gnormfac3
        elif any(i in mo_lbl for i in gnormlbl4):
            mo_coeff[i_lbl,:] *= gnormfac4
        elif any(i in mo_lbl for i in hnormlbl1):
            mo_coeff[i_lbl,:] *= hnormfac1
        elif any(i in mo_lbl for i in hnormlbl2):
            mo_coeff[i_lbl,:] *= hnormfac2
        elif any(i in mo_lbl for i in hnormlbl3):
            mo_coeff[i_lbl,:] *= hnormfac3
        elif any(i in mo_lbl for i in hnormlbl4):
            mo_coeff[i_lbl,:] *= hnormfac4
        elif any(i in mo_lbl for i in hnormlbl5):
            mo_coeff[i_lbl,:] *= hnormfac5

    # Mo_coeff actif
    mo_c = np.array([c[:,cas_idx] for c in mo_coeff] if cas_idx is not None else mo_coeff)
    e_c =  np.array([e[cas_idx] for e in mf.mo_energy] if cas_idx is not None else mf.mo_energy)
  
    nao, nmo = mo_c.shape

    print("n active MOs", nmo)
    print("n AOs", nao)
    assert nao==nao_c, "wrong number of AOs"

    ##########################################
    #                                        #
    #                Nuclei                  #
    #                                        #
    ##########################################

    natom = mol.natm
    print('n_atom',   natom)

    atom_xyz = mol.atom_coords(unit='Bohr')
    #if not(mol.unit.startswith(('B','b','au','AU'))):
    #    from pyscf.data.nist import BOHR
    #    atom_xyz /= BOHR # always convert to au

    with h5py.File(qph5path,'a') as qph5:
        qph5['nuclei'].attrs['nucl_num']=natom
        qph5.create_dataset('nuclei/nucl_coord',data=atom_xyz)
        qph5.create_dataset('nuclei/nucl_charge',data=mol.atom_charges())

        strtype=h5py.special_dtype(vlen=str)
        atom_dset=qph5.create_dataset('nuclei/nucl_label',(natom,),dtype=strtype)
        for i in range(natom):
            atom_dset[i] = mol.atom_pure_symbol(i)

    ##########################################
    #                                        #
    #                 ECP                    #
    #                                        #
    ##########################################

    if (mol.has_ecp()):
        #atsymb = [mol.atom_pure_symbol(i) for i in range(natom)]
        #pyecp = mol._ecp
        ## nelec to remove for each atom
        #nuc_z_remov = [pyecp[i][0] for i in atsymb]
        #nl_per_atom = [len(pyecp[i][1]) for i in atsymb]
        ## list of l-values for channels of each atom
        #ecp_l = [[pyecp[i][1][j][0] for j in range(len(pyecp[i][1]))] for i in atsymb]
        ## list of [exp,coef] for each channel (r**0,1,2,3,4,5,)
        #ecp_ac = [[pyecp[i][1][j][1] for j in range(len(pyecp[i][1]))] for i in atsymb]
        pyecp = [mol._ecp[mol.atom_pure_symbol(i)] for i in range(natom)]
        nzrmv=[0]*natom
        lmax=0
        klocmax=0
        knlmax=0
        for i,(nz,dat) in enumerate(pyecp):
            nzrmv[i]=nz
            for lval,ac in dat:
                if (lval==-1):
                    klocmax=max(sum(len(j) for j in ac),klocmax)
                else:
                    lmax=max(lval,lmax)
                    knlmax=max(sum(len(j) for j in ac),knlmax)
        #psd_nk = np.zeros((natom,klocmax),dtype=int)
        #psd_vk = np.zeros((natom,klocmax),dtype=float)
        #psd_dzk = np.zeros((natom,klocmax),dtype=float)
        #psd_nkl = np.zeros((natom,knlmax,lmax+1),dtype=int)
        #psd_vkl = np.zeros((natom,knlmax,lmax+1),dtype=float)
        #psd_dzkl = np.zeros((natom,knlmax,lmax+1),dtype=float)
        klnlmax=max(klocmax,knlmax)
        psd_n = np.zeros((lmax+2,klnlmax,natom),dtype=int)
        psd_v = np.zeros((lmax+2,klnlmax,natom),dtype=float)
        psd_dz = np.zeros((lmax+2,klnlmax,natom),dtype=float)
        for i,(_,dat) in enumerate(pyecp):
            for lval,ac in dat:
                count=0
                for ri,aici in enumerate(ac):
                    for ai,ci in aici:
                        psd_n[lval+1,count,i] = ri-2
                        psd_v[lval+1,count,i] = ci
                        psd_dz[lval+1,count,i] = ai
                        count += 1
        psd_nk = psd_n[0,:klocmax]
        psd_vk = psd_v[0,:klocmax]
        psd_dzk = psd_dz[0,:klocmax]
        psd_nkl = psd_n[1:,:knlmax]
        psd_vkl = psd_v[1:,:knlmax]
        psd_dzkl = psd_dz[1:,:knlmax]
        with h5py.File(qph5path,'a') as qph5:
            qph5['pseudo'].attrs['do_pseudo']=True
            qph5['pseudo'].attrs['pseudo_lmax']=lmax
            qph5['pseudo'].attrs['pseudo_klocmax']=klocmax
            qph5['pseudo'].attrs['pseudo_kmax']=knlmax
            qph5.create_dataset('pseudo/nucl_charge_remove',data=nzrmv)
            qph5.create_dataset('pseudo/pseudo_n_k',data=psd_nk)
            qph5.create_dataset('pseudo/pseudo_n_kl',data=psd_nkl)
            qph5.create_dataset('pseudo/pseudo_v_k',data=psd_vk)
            qph5.create_dataset('pseudo/pseudo_v_kl',data=psd_vkl)
            qph5.create_dataset('pseudo/pseudo_dz_k',data=psd_dzk)
            qph5.create_dataset('pseudo/pseudo_dz_kl',data=psd_dzkl)

        ## nelec to remove for each atom
        #nuc_z_remov = [i[0] for i in pyecp]
        #nl_per_atom = [len(i[1]) for i in pyecp]
        ## list of l-values for channels of each atom
        #ecp_l = [[ j[0] for j in i[1] ] for i in pyecp]
        #lmax = max(map(max,ecp_l))
        ## list of [exp,coef] for each channel (r**0,1,2,3,4,5,)
        #ecp_ac = [[ j[1] for j in i[1] ] for i in pyecp]


    ##########################################
    #                                        #
    #                Basis                   #
    #                                        #
    ##########################################

    # nucleus on which each AO is centered
    ao_nucl=[i[0] for i in mf.mol.ao_labels(fmt=False,base=1)]


    nprim_max = 0
    for iatom, (sh0,sh1,ao0,ao1) in enumerate(mol.aoslice_by_atom()):
        for ib in range(sh0,sh1): # sets of contracted exponents
            nprim = mol.bas_nprim(ib)
            if (nprim > nprim_max):
                nprim_max = nprim
    
    qp_prim_num = np.zeros((nao),dtype=int)
    qp_coef = np.zeros((nao,nprim_max))
    qp_expo = np.zeros((nao,nprim_max))
    qp_nucl = np.zeros((nao),dtype=int)
    qp_pwr = np.zeros((nao,3),dtype=int)
    
    clabels = mol.cart_labels(fmt=False)
    
    tmp_idx=0
    for iatom, (sh0,sh1,ao0,ao1) in enumerate(mol.aoslice_by_atom()):
        # shell start,end; AO start,end (sph or cart) for each atom
        for ib in range(sh0,sh1): # sets of contracted exponents
            l = mol.bas_angular(ib)    # angular momentum
            nprim = mol.bas_nprim(ib)  # numer of primitives
            es = mol.bas_exp(ib)       # exponents
            cs = mol.bas_ctr_coeff(ib) # coeffs
            nctr = mol.bas_nctr(ib)    # number of contractions
            print(iatom,ib,l,nprim,nctr,tmp_idx)
            for ic in range(nctr): # sets of contraction coeffs
                for nfunc in range(((l+1)*(l+2))//2): # always use cart for qp ao basis?
                    qp_expo[tmp_idx,:nprim] = es[:]
                    qp_coef[tmp_idx,:nprim] = cs[:,ic]
                    qp_nucl[tmp_idx] = iatom + 1
                    qp_pwr[tmp_idx,:] = xyzcount(clabels[tmp_idx][3])
                    qp_prim_num[tmp_idx] = nprim
                    tmp_idx += 1

    with h5py.File(qph5path,'a') as qph5:
        qph5['mo_basis'].attrs['mo_num']=nmo
        qph5['ao_basis'].attrs['ao_num']=nao

        #qph5['ao_basis'].attrs['ao_basis']=mf.cell.basis
        qph5['ao_basis'].attrs['ao_basis']="dummy basis"

        qph5.create_dataset('ao_basis/ao_nucl',data=qp_nucl)
        qph5.create_dataset('ao_basis/ao_prim_num',data=qp_prim_num)
        qph5.create_dataset('ao_basis/ao_expo',data=qp_expo.T)
        qph5.create_dataset('ao_basis/ao_coef',data=qp_coef.T)
        qph5.create_dataset('ao_basis/ao_power',data=qp_pwr.T)

    ##########################################
    #                                        #
    #              Electrons                 #
    #                                        #
    ##########################################

    nelec = mol.nelectron
    neleca,nelecb = mol.nelec

    print('num_elec', nelec)

    with h5py.File(qph5path,'a') as qph5:
        qph5['electrons'].attrs['elec_alpha_num']=neleca 
        qph5['electrons'].attrs['elec_beta_num']=nelecb

    ##########################################
    #                                        #
    #           Nuclear Repulsion            #
    #                                        #
    ##########################################

    e_nuc = mol.energy_nuc()
  
    print('nucl_repul', e_nuc)

    with h5py.File(qph5path,'a') as qph5:
        qph5['nuclei'].attrs['nuclear_repulsion']=e_nuc
  
    ##########################################
    #                                        #
    #               MO Coef                  #
    #                                        #
    ##########################################
    
    with h5py.File(qph5path,'a') as qph5:
        qph5.create_dataset('mo_basis/mo_coef',data=mo_c.T)
   
    return
