BEGIN_PROVIDER [complex*16, ao_kinetic_integral, (ao_num,ao_num)]
  implicit none
  BEGIN_DOC
  ! array of the priminitve basis kinetic integrals
  !  \langle \chi_i |\hat{T}| \chi_j \rangle
  END_DOC
  
  if (read_ao_one_integrals) then
    call read_one_e_integrals_complex('ao_kinetic_integral', ao_kinetic_integral,&
        size(ao_kinetic_integral,1), size(ao_kinetic_integral,2))
    print *,  'AO kinetic integrals read from disk'
  else
    print *, 'complex AO kinetic integrals must be provided'
    stop
  endif
END_PROVIDER

BEGIN_PROVIDER [ complex*16, ao_nucl_elec_integral, (ao_num,ao_num)]
   implicit none
   BEGIN_DOC
   ! interaction nuclear electron
   END_DOC
   
   if (read_ao_one_integrals) then
    call read_one_e_integrals_complex('ao_ne_integral', ao_nucl_elec_integral,      &
            size(ao_nucl_elec_integral,1), size(ao_nucl_elec_integral,2))
     print *,  'AO N-e integrals read from disk'
   else
     print *, 'complex AO N-e integrals must be provided'
     stop
   endif
   
END_PROVIDER

 BEGIN_PROVIDER [ complex*16, ao_nucl_elec_integral_per_atom, (ao_num,ao_num,nucl_num)]
 BEGIN_DOC
! ao_nucl_elec_integral_per_atom(i,j,k) = -<AO(i)|1/|r-Rk|AO(j)> 
! where Rk is the geometry of the kth atom
 END_DOC
  print *, 'ao_nucl_elec_integral_per_atom not implemented for k-points'
  stop
END_PROVIDER


 BEGIN_PROVIDER [ complex*16, ao_mono_elec_integral,(ao_num,ao_num)]
&BEGIN_PROVIDER [ double precision, ao_mono_elec_integral_diag,(ao_num)]
  implicit none
  integer :: i,j
  BEGIN_DOC
 ! array of the mono electronic hamiltonian on the AOs basis
 ! : sum of the kinetic and nuclear electronic potential 
  END_DOC
  do j = 1, ao_num
   do i = 1, ao_num
    ao_mono_elec_integral(i,j) = ao_nucl_elec_integral(i,j) + ao_kinetic_integral(i,j)
   enddo
   ao_mono_elec_integral_diag(j) = real(ao_mono_elec_integral(j,j))
  enddo
END_PROVIDER

