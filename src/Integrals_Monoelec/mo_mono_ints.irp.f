BEGIN_PROVIDER [complex*16, mo_kinetic_integral, (mo_tot_num,mo_tot_num)]
  implicit none
  BEGIN_DOC
  !  Kinetic energy integrals in the MO basis
  END_DOC

  if (read_mo_one_integrals) then
    call read_one_e_integrals_complex('mo_kinetic_integral', mo_kinetic_integral,&
        size(mo_kinetic_integral,1), size(mo_kinetic_integral,2))
    print *,  'MO kinetic integrals read from disk'
  else
    call complex_ao_to_mo(                                           &
        ao_kinetic_integral,                                         &
        size(ao_kinetic_integral,1),                                 &
        mo_kinetic_integral,                                         &
        size(mo_kinetic_integral,1)                                  &
        )
  endif
  if (write_mo_one_integrals) then
    call write_one_e_integrals_complex('mo_kinetic_integral', mo_kinetic_integral,&
        size(mo_kinetic_integral,1), size(mo_kinetic_integral,2))
    print *,  'MO kinetic integrals written to disk'
  endif

END_PROVIDER

BEGIN_PROVIDER [complex*16, mo_nucl_elec_integral, (mo_tot_num,mo_tot_num)]
 implicit none
 BEGIN_DOC
! interaction nuclear electron on the MO basis
 END_DOC

  if (read_mo_one_integrals) then
    call read_one_e_integrals_complex('mo_ne_integral', mo_nucl_elec_integral, &
        size(mo_nucl_elec_integral,1), size(mo_nucl_elec_integral,2))
    print *,  'MO N-e integrals read from disk'
  else
    call complex_ao_to_mo(                                            &
        ao_nucl_elec_integral,                               &
        size(ao_nucl_elec_integral,1),                       &
        mo_nucl_elec_integral,                                       &
        size(mo_nucl_elec_integral,1)                                &
        )
  endif
  if (write_mo_one_integrals) then
    call write_one_e_integrals_complex('mo_ne_integral', mo_nucl_elec_integral, &
        size(mo_nucl_elec_integral,1), size(mo_nucl_elec_integral,2))
    print *,  'MO N-e integrals written to disk'
  endif

END_PROVIDER


BEGIN_PROVIDER [complex*16, mo_nucl_elec_integral_per_atom, (mo_tot_num,mo_tot_num,nucl_num)]
 implicit none
 BEGIN_DOC
! mo_nucl_elec_integral_per_atom(i,j,k) = -<MO(i)|1/|r-Rk|MO(j)> 
! where Rk is the geometry of the kth atom
 END_DOC

 integer :: k
 mo_nucl_elec_integral_per_atom = 0.d0
 do k = 1, nucl_num 
   call complex_ao_to_mo(                                               &
       ao_nucl_elec_integral_per_atom(1,1,k),                        &
       size(ao_nucl_elec_integral_per_atom,1),                       &
       mo_nucl_elec_integral_per_atom(1,1,k),                        &
       size(mo_nucl_elec_integral_per_atom,1)                        &
       )
 enddo

END_PROVIDER

BEGIN_PROVIDER [ complex*16, mo_mono_elec_integral,(mo_tot_num,mo_tot_num)]
  implicit none
  integer                        :: i,j
  BEGIN_DOC
  ! array of the mono electronic hamiltonian on the MOs basis :
  ! sum of the kinetic and nuclear electronic potential (and pseudo potential if needed)
  END_DOC
  print*,'Providing the mono electronic integrals'

  do j = 1, mo_tot_num
    do i = 1, mo_tot_num
      mo_mono_elec_integral(i,j) = mo_nucl_elec_integral(i,j) + mo_kinetic_integral(i,j)
    enddo
  enddo

END_PROVIDER
