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
    print("n active Mos per kpt", nmo)
  
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

    with open('kconserv_complex','w') as outfile:
        for a in range(Nk):
            for b in range(Nk):
                for c in range(Nk):
                    d = kconserv[a,b,c]
                    outfile.write('%s %s %s %s\n' % (a+1,c+1,b+1,d+1))

    
    
  
