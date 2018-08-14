BEGIN_PROVIDER [complex*16, ao_kinetic_integral_kpts, (ao_num,ao_num,num_kpts)]
  implicit none
  BEGIN_DOC
  ! array of the priminitve basis kinetic integrals
  !  \langle \chi_i |\hat{T}| \chi_j \rangle
  END_DOC
  
  if (read_ao_one_integrals_kpts) then
    call read_one_e_integrals_complex_kpts('ao_kinetic_integral_kpts', ao_kinetic_integral_kpts,&
        size(ao_kinetic_integral_kpts,1), size(ao_kinetic_integral_kpts,2), size(ao_kinetic_integral_kpts,3))
    print *,  'AO kinetic integrals read from disk'
  else
    print *, 'complex AO kinetic integrals must be provided'
    stop
  endif
END_PROVIDER

BEGIN_PROVIDER [ complex*16, ao_nucl_elec_integral_kpts, (ao_num,ao_num,num_kpts)]
   implicit none
   BEGIN_DOC
   ! interaction nuclear electron
   END_DOC
   
   if (read_ao_one_integrals_kpts) then
    call read_one_e_integrals_complex_kpts('ao_ne_integral_kpts', ao_nucl_elec_integral_kpts,      &
            size(ao_nucl_elec_integral_kpts,1), size(ao_nucl_elec_integral_kpts,2), size(ao_nucl_elec_integral_kpts,3))
     print *,  'AO N-e integrals read from disk'
   else
     print *, 'complex AO N-e integrals must be provided'
     stop
   endif
   
END_PROVIDER

 BEGIN_PROVIDER [ complex*16, ao_nucl_elec_integral_per_atom_kpts, (ao_num,ao_num,nucl_num,num_kpts)]
 BEGIN_DOC
! ao_nucl_elec_integral_per_atom(i,j,k) = -<AO(i)|1/|r-Rk|AO(j)> 
! where Rk is the geometry of the kth atom
 END_DOC
  print *, 'ao_nucl_elec_integral_per_atom not implemented for k-points'
  stop
END_PROVIDER


 BEGIN_PROVIDER [ complex*16, ao_mono_elec_integral_kpts,(ao_num,ao_num,num_kpts)]
&BEGIN_PROVIDER [ double precision, ao_mono_elec_integral_diag_kpts,(ao_num,num_kpts)]
  implicit none
  integer :: i,j
  BEGIN_DOC
 ! array of the mono electronic hamiltonian on the AOs basis
 ! : sum of the kinetic and nuclear electronic potential 
  END_DOC
  do k = 1, num_kpts
    do j = 1, ao_num
      do i = 1, ao_num
        ao_mono_elec_integral_kpts(i,j,k) = ao_nucl_elec_integral_kpts(i,j,k) + ao_kinetic_integral_kpts(i,j,k)
      enddo
      ao_mono_elec_integral_diag_kpts(j,k) = real(ao_mono_elec_integral_kpts(j,j,k))
    enddo
  enddo
END_PROVIDER

