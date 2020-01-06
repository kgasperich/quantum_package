program green
  implicit none
  BEGIN_DOC
! TODO
  END_DOC
  read_wf = .True.
  touch read_wf
!  print*,'ref_bitmask_energy = ',ref_bitmask_energy
!  call psicoefprinttest
  call print_lanczos_eigvals_homo
  call print_spec_homo
end

subroutine psicoefprinttest
  implicit none
  integer :: i
  TOUCH psi_coef
  print *, 'printing ndet', N_det
end
subroutine print_lanczos_eigvals_homo
  implicit none
  integer :: i, iunit, j
  integer :: getunitandopen
  character(5) :: jstr

  do j=0,0
    write(jstr,'(I0.3)') j
    iunit = getunitandopen('lanczos_eigval_alpha_beta.out.'//trim(jstr),'w')
    print *, 'printing lanczos eigenvalues, alpha, beta to "lanczos_eigval_alpha_beta.out.'//trim(jstr)//'"'
    do i=1,n_lanczos_iter
      write(iunit,'(I6,3(E25.15))') i, lanczos_eigvals_homo(i), alpha_lanczos_homo(i), beta_lanczos_homo(i)
    enddo
    close(iunit)
  enddo
end
subroutine print_spec_homo
  implicit none
  integer :: i, iunit, j
  integer :: getunitandopen
  character(5) :: jstr
  do j=0,0
    write(jstr,'(I0.3)') j
    iunit = getunitandopen('omega_A.out.'//trim(jstr),'w')
    print *, 'printing frequency, spectral density to "omega_A.out.'//trim(jstr)//'"'
    do i=1,n_omega
      write(iunit,'(3(E25.15))') omega_list(i), spectral_lanczos_homo(i)
    enddo
    close(iunit)
  enddo
end
