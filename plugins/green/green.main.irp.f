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
  integer :: i, iunit
  integer :: getunitandopen
  iunit = getunitandopen('lanczos_eigval_alpha_beta.out','w')
  print *, 'printing lanczos eigenvalues, alpha, beta to "lanczos_eigval_alpha_beta.out"'
  do i=1,n_lanczos_iter
    write(iunit,'(I6,3(E25.15))') i, lanczos_eigvals(i), alpha_lanczos(i), beta_lanczos(i)
  enddo
  close(iunit)
end
subroutine print_spec
  implicit none
  integer :: i, iunit
  integer :: getunitandopen
  iunit = getunitandopen('omega_A.out','w')
  print *, 'printing frequency, spectral density to "omega_A.out"'
  do i=1,n_omega
    write(iunit,'(2(E25.15))') omega_list(i), spectral_lanczos(i)
  enddo
  close(iunit)
end
