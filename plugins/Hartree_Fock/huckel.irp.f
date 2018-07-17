subroutine huckel_guess
  implicit none
  BEGIN_DOC
! Build the MOs using the extended Huckel model
  END_DOC
  integer                        :: i,j
  double precision               :: accu
  double precision               :: c
  character*(64)                 :: label
  complex*16, allocatable  :: A(:,:)
  label = "Guess"
  c = 0.5d0 * 1.75d0

  allocate (A(ao_num, ao_num))
  A = 0.d0
  do j=1,ao_num
    do i=1,ao_num
      A(i,j) = c * ao_overlap(i,j) * (ao_mono_elec_integral_diag(i) + ao_mono_elec_integral_diag(j))
    enddo
    A(j,j) = ao_mono_elec_integral_diag(j) + real(ao_bi_elec_integral_alpha(j,j))
    if (abs(imag(ao_bi_elec_integral_alpha(j,j))) .gt. 1.0d-10) then
      stop 'diagonal elements of ao_bi_elec_integral_alpha should be real'
    endif
  enddo

!  Fock_matrix_ao_alpha(1:ao_num,1:ao_num) = A(1:ao_num,1:ao_num)
!  Fock_matrix_ao_beta (1:ao_num,1:ao_num) = A(1:ao_num,1:ao_num)
  call zlacp2('X', ao_num, ao_num, A, size(A,1), &
         Fock_matrix_ao_alpha, size(Fock_matrix_ao_alpha,1))
  call zlacp2('X', ao_num, ao_num, A, size(A,1), &
         Fock_matrix_ao_beta,  size(Fock_matrix_ao_beta, 1))
  

!  TOUCH mo_coef

  TOUCH Fock_matrix_ao_alpha Fock_matrix_ao_beta
  mo_coef = eigenvectors_fock_matrix_mo
  SOFT_TOUCH mo_coef
  call save_mos
  deallocate(A)

end
