import numpy as np
from functools import reduce

def pyscf2QP(cell,mf, kpts, kmesh=None, cas_idx=None, int_threshold = 1E-8):
    '''
    kpts = List of kpoints coordinates. Cannot be null, for gamma is other script
    kmesh = Mesh of kpoints (optional)
    cas_idx = List of active MOs. If not specified all MOs are actives
    int_threshold = The integral will be not printed in they are bellow that
    '''
  
    from pyscf.pbc import ao2mo
    from pyscf.pbc import tools
    from pyscf import lib
    from pyscf.pbc.gto import ecp
    
    mo_coef_threshold = int_threshold
    ovlp_threshold = int_threshold
    kin_threshold = int_threshold
    ne_threshold = int_threshold
    bielec_int_threshold = int_threshold

    natom = len(cell.atom_coords())
    print('n_atom',   natom)
    print('num_elec', cell.nelectron)
  
    mo_coeff = mf.mo_coeff
    # Mo_coeff actif
    mo_k = np.array([c[:,cas_idx] for c in mo_coeff] if cas_idx is not None else mo_coeff)
    e_k =  np.array([e[cas_idx] for e in mf.mo_energy] if cas_idx is not None else mf.mo_energy)
  
    Nk, nao, nmo = mo_k.shape
    print("n Kpts", Nk)
    print("n active Mos", nmo)
  
    # Write all the parameter need to creat a dummy EZFIO folder who will containt the integral after.
    # More an implentation detail than a real thing
    with open('param','w') as f:
    # Note the use of nmo_tot
        f.write(' '.join(map(str,(cell.nelectron*Nk, Nk*nmo, natom*Nk))))
  
    with open('num_ao','w') as f:
        f.write(str(nao*Nk))
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
  
    def get_phase(cell, kpts, kmesh=None):
        '''
        The unitary transformation that transforms the supercell basis k-mesh
        adapted basis.
        '''
  
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
  
    #       __    __          _                                 
    # |\/| |  |  |    _   _  |_  _ 
    # |  | |__|  |__ (_) (/_ |  _> 
    #                                               
    with open('mo_coef_complex','w') as outfile:
        c_kpts = np.reshape(mf.mo_coeff_kpts,(Nk,nao,nmo))

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
    
    with open('overlap_ao_complex','w') as outfile:
        #s_kpts_ao = np.reshape(mf.cell.pbc_intor('int1e_ovlp_sph',kpts=kpts),(Nk,nao,nao))
        s_kpts_ao = np.reshape(mf.get_ovlp(cell=cell,kpts=kpts),(Nk,nao,nao))

        for ik in range(Nk):
            shift=ik*nao+1
            for i in range(nao):
                for j in range(i,nao):
                    sij = s_kpts_ao[ik,i,j]
                    if abs(sij) > ovlp_threshold:
                        outfile.write('%s %s %s %s\n' % (i+shift, j+shift, sij.real, sij.imag))
    
    with open('kinetic_ao_complex','w') as outfile:
        #t_kpts_ao = np.reshape(mf.cell.pbc_intor('int1e_kin_sph',kpts=kpts),(Nk,nao,nao))
        t_kpts_ao = np.reshape(cell.pbc_intor('int1e_kin',1,1,kpts=kpts),(Nk,nao,nao))

        for ik in range(Nk):
            shift=ik*nao+1
            for i in range(nao):
                for j in range(i,nao):
                    tij = t_kpts_ao[ik,i,j]
                    if abs(tij) > kin_threshold:
                        outfile.write('%s %s %s %s\n' % (i+shift, j+shift, tij.real, tij.imag))
    
    with open('ne_ao_complex','w') as outfile:
        if mf.cell.pseudo:
            v_kpts_ao = np.reshape(mf.with_df.get_pp(kpts=kpts),(Nk,nao,nao))
        else:
            v_kpts_ao = np.reshape(mf.with_df.get_nuc(kpts=kpts),(Nk,nao,nao))
        if len(cell._ecpbas) > 0:
            v_kpts_ao += np.reshape(ecp.ecp_int(cell, kpts),(Nk,nao,nao))

        for ik in range(Nk):
            shift=ik*nao+1
            for i in range(nao):
                for j in range(i,nao):
                    vij = v_kpts_ao[ik,i,j]
                    if abs(vij) > ne_threshold:
                        outfile.write('%s %s %s %s\n' % (i+shift, j+shift, vij.real, vij.imag))

  
    # ___                              _    
    #  |  ._ _|_  _   _  ._ _. |  _   |_) o 
    # _|_ | | |_ (/_ (_| | (_| | _>   |_) | 
    #                 _|                    
    #
    kconserv = tools.get_kconserv(cell, kpts)
    def idx2_tri(i,j):
        ij1=min(i,j)
        ij2=max(i,j)
        return ij1+(ij2*(ij2-1))//2
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
#
#    with open('bielec_ao_complex','w') as outfile: 
#        for d in range(Nk):
#            for c in range(Nk):
#                if c > d: break
#                idx2_cd = idx2_tri(c,d)
#                for b in range(Nk):
#                    if b > d: break
#                    a = kconserv[b,c,d]
#                    if idx2_tri(a,b) > idx2_cd: continue
#                    if ((c==d) and (a>b)): continue
#                    for l in range(nao):
#                        ll=l+d*nao
#                        for j in range(nao):
#                            jj=j+c*nao
#                            if jj>ll: break
#                            idx2_jjll = idx2_tri(jj,ll)
#                            for k in range(nao):
#                                kk=k+b*nao
#                                if kk>ll: break
#                                for i in range(nao):
#                                    ii=k+a*nao
#                                    if idx2_tri(ii,kk) > idx2_jjll: break
#                                    if ((jj==ll) and (ii>kk)): break
#                                        v=eri_4d_ao[ii,kk,jj,ll]
#                                        if (abs(v) > bielec_int_threshold):
#                                            outfile.write('%s %s %s %s %s %s\n' % (ii+1,jj+1,kk+1,ll+1,v.real,v.imag))

    with open('bielec_ao_complex','w') as outfile: 
        for d, kd in enumerate(kpts):
            for c, kc in enumerate(kpts):
                if c > d: break
                idx2_cd = idx2_tri(c,d)
                for b, kb in enumerate(kpts):
                    if b > d: break
                    a = kconserv[b,c,d]
                    if idx2_tri(a,b) > idx2_cd: continue
                    if ((c==d) and (a>b)): continue
                    ka = kpts[a]
                    eri_4d_ao_kpt = mf.with_df.get_ao_eri(kpts=[ka,kb,kc,kd],compact=False).reshape((nao,)*4)
                    eri_4d_ao_kpt *= 1./Nk
                    for l in range(nao):
                        ll=l+d*nao
                        for j in range(nao):
                            jj=j+c*nao
                            if jj>ll: break
                            idx2_jjll = idx2_tri(jj,ll)
                            for k in range(nao):
                                kk=k+b*nao
                                if kk>ll: break
                                for i in range(nao):
                                    ii=i+a*nao
                                    if idx2_tri(ii,kk) > idx2_jjll: break
                                    if ((jj==ll) and (ii>kk)): break
                                    v=eri_4d_ao_kpt[i,k,j,l]
                                    if (abs(v) > bielec_int_threshold):
                                        outfile.write('%s %s %s %s %s %s\n' % (ii+1,jj+1,kk+1,ll+1,v.real,v.imag))

    
    
  
