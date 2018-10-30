#!/bin/bash

ezfio=$1
# Create the integral
echo 'Create Integral'

echo 'Create EZFIO'
read nel nmo natom <<< $(cat param) 
read e_nucl <<< $(cat e_nuc)
read nao <<< $(cat num_ao)
read nkpts <<< $(cat num_kpts)
read ndf <<< $(cat num_df)
./create_ezfio_complex.py $ezfio $nel $natom $nmo $e_nucl $nao $nkpts $ndf
#Handle the orbital consitensy check
qp_edit -c $ezfio &> /dev/null
cp $ezfio/{ao,mo}_basis/ao_md5 

#Read the integral
echo 'Read Integral'




################################################
##  using AO mono, 3-idx, mo coef from pyscf  ##
################################################

qp_run read_ao_mono_complex $ezfio 
qp_run read_kconserv $ezfio
qp_run read_ao_df_complex $ezfio
qp_run read_mo_coef_complex $ezfio    #start from converged pyscf MOs
#qp_run mo_from_ao_orth $ezfio        #use canonical orthonormalized AOs as initial MO guess


###############################################################
##  using AO mono, full 4-idx AO bielec, mo coef from pyscf  ##
###############################################################

#qp_run read_ao_mono_complex $ezfio 
#qp_run read_kconserv $ezfio
#qp_run read_ao_eri_chunk_complex $ezfio 
#qp_run read_mo_coef_complex $ezfio    #start from converged pyscf MOs
##qp_run mo_from_ao_orth $ezfio        #use canonical orthonormalized AOs as initial MO guess


######################################################
##  using MO mono, full 4-idx MO bielec from pyscf  ##
######################################################

#qp_run read_mo_mono_complex $ezfio 
#qp_run read_kconserv $ezfio
#qp_run read_mo_eri_chunk_complex $ezfio 

