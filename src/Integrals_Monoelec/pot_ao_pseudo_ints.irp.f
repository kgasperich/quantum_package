BEGIN_PROVIDER [ complex*16, ao_pseudo_integral, (ao_num,ao_num)]
  implicit none
  BEGIN_DOC
  ! Pseudo-potential integrals
  END_DOC
  
  if (read_ao_one_integrals) then
    call read_one_e_integrals_complex('ao_pseudo_integral', ao_pseudo_integral,&
        size(ao_pseudo_integral,1), size(ao_pseudo_integral,2))
    print *,  'AO pseudopotential integrals read from disk'
  else
    print *, 'complex AO pseudopotential integrals must be provided'
    stop
  endif
  
END_PROVIDER

