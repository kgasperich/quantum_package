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


nxcell=3
nycell=2
nzcell=1
nmp = [nxcell, nycell, nzcell]

cell = gto.Cell()

cell.a = '''
        4 0 0 
        0 4 0
        0 0 4'''
cell.atom = '''  
   H        0.00000000       0.00000000       0.00000000
   H        2.5 2.2 0.5
   H        0.7 0.9 3.3
   H        1.1 3.9 1.3
            '''
#cell.basis='bfd-vdz'
#cell.ecp = 'bfd'
cell.basis = 'sto-3g'
cell.pseudo = 'gth-pade'


cell.unit='B'
cell.drop_exponent=0.1

cell.verbose = 5


cell.build()


dokpts=True

if dokpts:
  kpts = cell.make_kpts(nmp)
  kpts -= kpts[0]

  supcell=cell
  mydf = df.GDF(supcell,kpts)
  mydf.auxbasis = 'weigend'
#  mydf._cderi_to_save = 'df_ints.h5'   # new
#  mydf.build()                         # new
  mf = scf.KRHF(supcell,kpts).density_fit()


else:

  supcell = super_cell(cell, nmp)
  mydf = df.GDF(supcell)
#  mydf._cderi_to_save = 'df_ints.h5'   # new
#  mydf.build()                         # new
  mf = scf.RHF(supcell).density_fit()
  kpts=[]


#mf.exxdiv = 'ewald'
mf.exxdiv = None
mf.with_df = mydf
mf.chkfile ='H-scf.chk'
dm = mf.from_chk('H-scf.chk')           # restart
mf.with_df._cderi = 'df_ints.h5'
e_scf=mf.kernel(dm)                     # restart
#e_scf=mf.kernel()                      # new


ener = open('e_scf','w')
ener.write('%s\n' % (e_scf))
print('e_scf',e_scf)
ener.close()

title="H-{:}{:}{:}".format(nxcell,nycell,nzcell)

#from PyscfToQmcpack import savetoqmcpack
#savetoqmcpack(supcell,mf,title=title,kpts=kpts)

if dokpts:
  from MolPyscfToQPkpts import pyscf2QP
  pyscf2QP(mf.cell,mf,kpts=mf.kpts,int_threshold = 1E-15)
#else:
#  from MolPyscfToQP import pyscf2QP
#  pyscf2QP(supcell,mf,kpts=kpts,int_threshold = 1E-15)




