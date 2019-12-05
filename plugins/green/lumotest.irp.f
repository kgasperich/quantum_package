program lumotest
  implicit none
  BEGIN_DOC
! TODO
  END_DOC
  read_wf = .True.
  touch read_wf
  print*,'ref_bitmask_energy = ',ref_bitmask_energy
!  call psicoefprinttest
  call lumo_add
!  call print_spec
end

subroutine lumo_add
  implicit none
  integer, allocatable :: idx_tmp(:)
  integer :: idx_lumo
  integer :: i,j,k
  integer :: n_det_new
  logical :: is_filled
  integer :: spin_lumo
  PROVIDE psi_det
  allocate(idx_tmp(N_det))

  call get_idx_spin_lumo(idx_lumo,spin_lumo)

  k=0
  do i=1,N_det
    call orb_is_filled(psi_det(i),idx_lumo,spin_lumo,N_int,is_filled)
    if (.not.is_filled) then
      k=k+1
      idx_tmp(k)=i
    endif
  enddo

  n_det_new = k
  integer(bit_kind), allocatable :: psi_det_tmp(:,:,:)
  complex*16, allocatable :: psi_coef_tmp(:)
  double precision :: phase_tmp
  
  allocate(psi_det_tmp(N_int,2,n_det_new),psi_coef_tmp(n_det_new))

  do i=1,n_det_new
    j=idx_tmp(i)
    call ac_operator_phase(psi_det_tmp(i),psi_det(j),idx_lumo,spin_lumo,N_int,phase_tmp)
    psi_coef_tmp(i) = phase_tmp*psi_coef(j)
  enddo
  
  deallocate(psi_coef,psi_det)
  N_det = n_det_new
  psi_det_size = n_det_new
  N_states = 1
  if (spin_lumo==1) then
    elec_alpha_num=elec_alpha_num+1
  else
    elec_beta_num=elec_beta_num+1
  endif

  SOFT_TOUCH psi_det_size N_states N_det elec_alpha_num elec_beta_num
  allocate(psi_det(N_int,2,psi_det_size), psi_coef(psi_det_size,N_states))
  psi_coef=psi_coef_tmp
  psi_det=psi_det_tmp


end

BEGIN_PROVIDER [ double precision, mo_energies_ref, (mo_tot_num, 2) ]
  implicit none
  use bitmasks
  implicit none

  integer :: i_int, i_orb
  integer, allocatable :: list_alpha(:), list_beta(:)
  integer :: n_elem_alpha, n_elem_beta
  allocate(list_alpha(N_int*bit_kind_size),list_beta(N_int*bit_kind_size))

  call bitstring_to_list(HF_bitmask(1:N_int,1),list_alpha,n_elem_alpha,N_int)
  call bitstring_to_list(HF_bitmask(1:N_int,2),list_beta,n_elem_beta,N_int)

  integer :: i
  do i=1,mo_tot_num
    

end

subroutine get_idx_spin_lumo(idx_lumo,spin_lumo)
  use bitmasks
  implicit none

  integer :: i_int, i_orb
  integer, allocatable :: list_alpha(:), list_beta(:)
  integer :: n_elem_alpha, n_elem_beta
  allocate(list_alpha(N_int*bit_kind_size),list_beta(N_int*bit_kind_size))

  call bitstring_to_list(HF_bitmask(1:N_int,1),list_alpha,n_elem_alpha,N_int)
  call bitstring_to_list(HF_bitmask(1:N_int,2),list_beta,n_elem_beta,N_int)



end


 BEGIN_PROVIDER [ integer(bit_kind), psi_det_lumo, (N_int,2,psi_det_lumo_size) ]
&BEGIN_PROVIDER [ complex*16, psi_coef_plus, (N_det_new) ]
  implicit none
  BEGIN_DOC
  ! initial lanczos vector
  ! must be normalized
  END_DOC

  integer :: i
  
  do i=1,N_det
    u1_lanczos(i)=1.d0/(dble(i))**0.1d0
  enddo
  call normalize_complex(u1_lanczos,N_det)

END_PROVIDER

subroutine get_phase_lumo(det_in,lumo_idx,phase)
  implicit none
  integer, intent(in) :: lumo_idx
  integer, intent(out) :: phase
end

