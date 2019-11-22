program green
  implicit none
  BEGIN_DOC
! TODO
  END_DOC
  read_wf = .True.
  touch read_wf
  print*,'ref_bitmask_energy = ',ref_bitmask_energy
!  call psicoefprinttest
  call print_lanczos_eigvals
  call print_spec
end

subroutine psicoefprinttest
  implicit none
  integer :: i
  TOUCH psi_coef
  print *, 'printing ndet', N_det
end
subroutine print_lanczos_eigvals
  implicit none
  integer :: i
  print *, 'printing lanczos eigenvalues'
  do i=1,n_lanczos_iter
    print *, i, lanczos_eigvals(i)
  enddo
end
subroutine print_spec
  implicit none
  integer :: i
  print *, 'printing spectral density'
  do i=1,n_omega
    print *, i, omega_list(i), spectral_lanczos(i)
  enddo
end
