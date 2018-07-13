BEGIN_PROVIDER [complex*16, ao_kinetic_integral, (ao_num,ao_num)]
  implicit none
  BEGIN_DOC
  ! array of the priminitve basis kinetic integrals
  !  \langle \chi_i |\hat{T}| \chi_j \rangle
  END_DOC
  integer                        :: i,j,k,l
  
  if (read_ao_one_integrals) then
    call read_one_e_integrals_complex('ao_kinetic_integral', ao_kinetic_integral,&
        size(ao_kinetic_integral,1), size(ao_kinetic_integral,2))
    print *,  'AO kinetic integrals read from disk'
  else
    print *, 'complex AO kinetic integrals must be provided'
  endif
END_PROVIDER



