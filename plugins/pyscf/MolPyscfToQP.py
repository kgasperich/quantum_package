
import numpy,re,sys
#from functools import reduce

def pyscf2QP(cell,mf, kpts=[], int_threshold = 1E-15):
  # The integral will be not printed in they are bellow that


  PBC=False
  Gamma=False 
  nkpts=len(kpts)
  print('nkpts=',nkpts)
  ComputeMode= re.split('[. ]', str(mf))
  print('ComputeMode=',ComputeMode)

  for n in ComputeMode:
    if n in ("UHF","KUHF","UKS"):
      sys.exit('Unrestricted calculation unsupported in Quantum Package')
    if n == "pbc":
      PBC=True

  if PBC and len(kpts) == 0:
    #sys.exit("ERROR (read!): You need to specify explicit the list of K-point (including gamma)")
    Gamma=True

  print('Performing PBC?:',PBC)
  if PBC:
    from pyscf.pbc import ao2mo
    from pyscf.pbc import tools
  else:
    from pyscf import ao2mo

  if not Gamma:
    natom = len(cell.atom_coords())*nkpts
    nelectron=cell.nelectron*nkpts
    print('n_atom',   natom)
    print('num_elec', nelectron)
    print('nucl_num', len(cell.atom_coords())*nkpts)
    #mo_coeff = numpy.concatenate((mf.mo_coeff[0],mf.mo_coeff[1]),axis=0)
    l_mo_coeff = mf.mo_coeff
  else:
    natom = len(cell.atom_coords())
    nelectron=cell.nelectron
    print('n_atom',   natom)
    print('num_elec', nelectron)
    print('nucl_num', len(cell.atom_coords()))
    mo_coeff = mf.mo_coeff
    l_mo_coeff = [mo_coeff]

  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print('Shape mo_coeff  [ (mo,ao) ]')
  l_mo_coeff_shape = [mo_coeff.shape for mo_coeff in l_mo_coeff]
  print(l_mo_coeff_shape)

  nmo = sum(shape[-1] for shape in l_mo_coeff_shape)
  nmo_kpts = l_mo_coeff[0].shape[-1]

  print('nmo', nmo)
  print('nmo_kpts', nmo_kpts)

  # Wrote all the parameter need to creat a dummy EZFIO folder who will containt the integral after.
  # More an implentation detail than a real thing
  with open('param','w') as f:
    f.write(' '.join(map(str,(nelectron, nmo, natom))))

  def trans(v):
    try:
      return v.real
    except:
      return v

  #                             _                             
  # |\ |      _ |  _   _. ._   |_)  _  ._      |  _ o  _  ._  
  # | \| |_| (_ | (/_ (_| |    | \ (/_ |_) |_| | _> | (_) | | 
  #                                    |                      

  print('mf, cell', mf.energy_nuc(), cell.energy_nuc())
  shift = tools.pbc.madelung(cell, numpy.zeros(3))*cell.nelectron * -.5 if PBC else 0
  e_nuc = cell.energy_nuc() + shift

  print('nucl_repul', e_nuc)
  with open('e_nuc','w') as f:
    f.write(str(e_nuc))


  from itertools import product

  # ___                                              
  #  |  ._ _|_  _   _  ._ _. |  _   |\/|  _  ._   _  
  # _|_ | | |_ (/_ (_| | (_| | _>   |  | (_) | | (_) 
  #                 _|                              



  if PBC:
    if Gamma:
      h_ao = ('kinetic', mf.get_hcore() )
      dummy_ao = ('nuclear', numpy.zeros( (nmo,nmo), dtype=numpy.float ))
    else:
      print('Computing mono electronics integrals for solid with KPTs (not special case Gamma)')
      h_ao = ('kinetic', mf.get_hcore(kpts=mf.kpts).reshape(nkpts, nmo_kpts, nmo_kpts)) # Give only one k point ?
      dummy_ao = ('nuclear', numpy.zeros( (nkpts,nmo_kpts,nmo_kpts), dtype=numpy.float ))
  else:
    h_ao = ('kinetic', mf.get_hcore() )
    dummy_ao = ('nuclear', numpy.zeros( (nmo,nmo), dtype=numpy.float ))


  def gen_mono_MO(mo_coeff,nmo,l_int,shift=0):
    # 2Id transfortion Transformation. For now we handle only one or zero K point.
    #print 'l_int.shape=',l_int.shape

    #l_int_mo = reduce(numpy.dot, (mo_coeff.T, l_int, mo_coeff)) #This formula is only right for one kpt.
    l_int_mo = mo_coeff.T @ l_int @ mo_coeff  #This formula is only right for one kpt. requires python 3.5 or later

    for i,j in product(list(range(nmo)), repeat=2):
      int_ = l_int_mo[i,j]
      #yield (i+1+shift,j+1+shift, trans(int_))
      yield (i+1+shift,j+1+shift, int_.real, int_.imag)

  # Print 
  for name, ao in (h_ao,dummy_ao):
    with open('%s_mo' % name,'w') as f:
      print('%s_mo' % name)
      if not PBC:
        for mono in gen_mono_MO(mo_coeff,nmo,ao):
          f.write('%s %s %s %s\n'% mono)
      else:
        if not Gamma:
          for i,(m,a) in enumerate(zip(l_mo_coeff,ao)):
            for mono in gen_mono_MO(m,nmo_kpts,a,shift=nmo_kpts*i):
              f.write('%s %s %s %s\n'% mono)
        else:
          for mono in gen_mono_MO(mo_coeff,nmo,ao):
            f.write('%s %s %s %s\n'% mono)

  # ___                              _    
  #  |  ._ _|_  _   _  ._ _. |  _   |_) o 
  # _|_ | | |_ (/_ (_| | (_| | _>   |_) | 
  #                 _|                    
  #

  def write_amazing(eri_4d, shift=0):
    print('eri', nmo, '*4')
    # HANDLE 8 FOLD by Scemama way. Maybe we can use compact=True
    for l in range(nmo):
      for k in range(nmo):
        for j in range(l,nmo):
          for i in range(max(j,k),nmo):
            if i==j and k>l: continue     
            v = eri_4d[i,k,j,l]
            if abs(v) > int_threshold:
              #f.write('%s %s %s %s %s\n' % (i+1+shift,j+1+shift,k+1+shift,l+1+shift,trans(v)))
              f.write('%s %s %s %s %s %s\n' % (i+1+shift,j+1+shift,k+1+shift,l+1+shift,v.real,v.imag))


  import pickle
  from itertools import chain
  pickle.dump(list(chain.from_iterable(l_mo_coeff)), open("mo_coeff","wb"))


  def gen_eri_full():

    kenec = dict()
    kpts = mf.kpts
    kconserv = tools.get_kconserv(cell, kpts)
    for a, ka in enumerate(kpts):
      for b, kb in enumerate(kpts):
        for c, kc in enumerate(kpts):
          d = kconserv[a,b,c]
          kd=kpts[d]
          eri_4d = mf.with_df.get_eri((ka,kb,kc,kd), compact=False).reshape((nmo_kpts,)*4)

          for i in range(nmo_kpts):
            for j in range(nmo_kpts):
              for k in range(nmo_kpts):
                for l in range(nmo_kpts):
                  v = eri_4d[i,k,j,l]
                  #f.write('%s %s %s %s %s\n' % (a*nmo+i+1,b*nmo+j+1, c*nmo+k+1, d*nmo+l+1, trans(v))) 
                  print((a*nmo_kpts+i,b*nmo_kpts+j, c*nmo_kpts+k, d*nmo_kpts+l))
                  #kenec[(a*nmo_kpts+i,b*nmo_kpts+j, c*nmo_kpts+k, d*nmo_kpts+l) ] =  trans(v) 
                  kenec[(a*nmo_kpts+i,b*nmo_kpts+j, c*nmo_kpts+k, d*nmo_kpts+l) ] =  v

    # 8 folds symmetrie
    f = open('bielec_ao','w')
    for l in range(nmo):
      for k in range(nmo):
        for j in range(l,nmo):
          for i in range(max(j,k),nmo):
            try:
              v = kenec[(i,j,k,l)]
            except KeyError:
              pass
            else:
              if abs(v) > int_threshold:
                #f.write('%s %s %s %s %s\n' % (i+1,j+1,k+1,l+1,trans(v)))    
                f.write('%s %s %s %s %s %s\n' % (i+1,j+1,k+1,l+1,v.real,v.imag))    

  if not PBC:
    f = open('bielec_ao','w')
    eri_4d = mf.with_df.get_eri(compact=False).reshape((nmo,)*4)
    write_amazing(eri_4d, 0)

  else :
    if Gamma:
      f = open('bielec_ao','w')
      eri_4d = mf.with_df.get_eri(compact=False).reshape((nmo,)*4)
      write_amazing(eri_4d, 0)
    else:
      #KPTS 
      gen_eri_full()


  f = open('bielec_mo','w')
  eri_4d= mf.with_df.ao2mo(mo_coeff,compact=False).reshape((nmo,)*4) 
  write_amazing(eri_4d, 0)







