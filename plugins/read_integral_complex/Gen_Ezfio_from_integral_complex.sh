#!/bin/bash

ezfio=$1
# Create the integral
echo 'Create Integral'

echo 'Create EZFIO'
read nel nmo natom <<< $(cat param) 
read e_nucl <<< $(cat e_nuc)
read nao <<< $(cat num_ao)
read nkpts <<< $(cat num_kpts)
./create_ezfio_complex.py $ezfio $nel $natom $nmo $e_nucl $nao $nkpts
#Handle the orbital consitensy check
qp_edit -c $ezfio &> /dev/null
cp $ezfio/{ao,mo}_basis/ao_md5 

#Read the integral
echo 'Read Integral'

qp_run read_mo_mono_complex $ezfio 
qp_run read_mo_eri_chunk_complex $ezfio 
#qp_run read_ao_mono_complex $ezfio 
#qp_run read_ao_eri_chunk_complex $ezfio 
#
##comment out one of the following two lines
#qp_run read_mo_coef_complex $ezfio #start from converged pyscf MOs 
#qp_run mo_from_ao_orth $ezfio #use canonical orthonormalized AOs as initial MO guess
