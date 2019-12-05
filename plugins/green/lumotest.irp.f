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


subroutine get_list_hp_banned_ab(tmp_det,N_hp,exc_is_banned,spin_hp,sign_hp,idx_hp,nint,all_banned)
  implicit none
  integer, intent(in) :: N_hp
  integer, intent(in) :: spin_hp(N_hp), sign_hp(N_hp), idx_hp(N_hp)
  integer(bit_kind), intent(in) :: tmp_det(nint,2)
  logical, intent(out) :: exc_is_banned(N_hp)
  logical, intent(out) :: all_banned
  
  integer :: i
  logical :: is_filled

  all_banned = .True.
  do i=1,N_hp
    call orb_is_filled(tmp_det,idx_hp(i),spin_hp(i),nint,is_filled)
    if (sign_hp(i).gt.0) then ! particle creation, banned if already filled
      exc_is_banned(i) = is_filled
    else ! hole creation, banned if already empty
      exc_is_banned(i) = (.not.is_filled)
    endif
    all_banned = (all_banned.and.exc_is_banned(i))
  enddo
end

subroutine get_list_hp_banned_single_spin(tmp_spindet,N_hp,exc_is_banned,spin_hp,sign_hp,idx_hp,ispin,nint,all_banned)
  implicit none
  integer, intent(in) :: N_hp, ispin
  integer, intent(in) :: spin_hp(N_hp), sign_hp(N_hp), idx_hp(N_hp)
  integer(bit_kind), intent(in) :: tmp_spindet(nint)
  logical, intent(out) :: exc_is_banned(N_hp)
  logical, intent(out) :: all_banned
  
  integer :: i
  logical :: is_filled

  all_banned = .True.
  do i=1,N_hp
    if (spin_hp(i).eq.ispin) then
      call orb_is_filled_single_spin(tmp_spindet,idx_hp(i),nint,is_filled)
      if (sign_hp(i).gt.0) then ! particle creation, banned if already filled
        exc_is_banned(i) = is_filled
      else ! hole creation, banned if already empty
        exc_is_banned(i) = (.not.is_filled)
      endif
    else
      exc_is_banned(i) = .False.
    endif
    all_banned = (all_banned.and.exc_is_banned(i))
  enddo
end

subroutine get_list_hp_banned_spin(tmp_det,N_hp,exc_is_banned,spin_hp,sign_hp,idx_hp,ispin,nint,all_banned)
  implicit none
  integer, intent(in) :: N_hp, ispin
  integer, intent(in) :: spin_hp(N_hp), sign_hp(N_hp), idx_hp(N_hp)
  integer(bit_kind), intent(in) :: tmp_det(nint,2)
  logical, intent(out) :: exc_is_banned(N_hp)
  logical, intent(out) :: all_banned
  
  integer :: i
  logical :: is_filled
  spindet(1:nint) = tmp_det(1:nint,ispin)

  all_banned = .True.
  do i=1,N_hp
    if (spin_hp(i).eq.ispin) then
      call orb_is_filled(tmp_det,idx_hp(i),ispin,nint,is_filled)
      if (sign_hp(i).gt.0) then ! particle creation, banned if already filled
        exc_is_banned(i) = is_filled
      else ! hole creation, banned if already empty
        exc_is_banned(i) = (.not.is_filled)
      endif
    else
      exc_is_banned(i) = .False.
    endif
    all_banned = (all_banned.and.exc_is_banned(i))
  enddo
end

  
subroutine orb_is_filled_bit_int(key_ref,iorb_int,iorb_bit,ispin,Nint,is_filled)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! apply creation operator to key_ref
  ! add electron with spin ispin to orbital with index iorb
  ! output resulting det and phase in key_new and phase
  END_DOC
  integer, intent(in)            :: iorb_int, iorb_bit, ispin, Nint
  integer(bit_kind), intent(in) :: key_ref(Nint,2)
  logical, intent(out) :: is_filled
  
  integer                        :: k,l
  
  ASSERT (iorb > 0)
  ASSERT (ispin > 0)
  ASSERT (ispin < 3)
  ASSERT (Nint > 0)
  
  ! alpha det is list of Nint 64-bit ints
  ! k is index of the int where iorb is found
  ! l is index of the bit where iorb is found
  k = ishft(iorb-1,-bit_kind_shift)+1
  ASSERT (k >0)
  l = iorb - ishft(k-1,bit_kind_shift)-1
  ASSERT (l >= 0)
  is_filled = btest(key_ref(iorb_int,ispin),iorb_bit)  
end

subroutine get_orb_bit_int(iorb,Nint)
  use bitmasks
  implicit none
  iorb_int = ishft(iorb-1,-bit_kind_shift)+1
  iorb_bit = iorb - ishft(iorb_int-1,bit_kind_shift)-1

end
subroutine orb_is_filled_single_spin(key_ref,iorb,Nint,is_filled)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! apply creation operator to key_ref
  ! add electron with spin ispin to orbital with index iorb
  ! output resulting det and phase in key_new and phase
  END_DOC
  integer, intent(in)            :: iorb, Nint
  integer(bit_kind), intent(in) :: key_ref(Nint)
  logical, intent(out) :: is_filled
  
  integer                        :: k,l
  
  ASSERT (iorb > 0)
  ASSERT (Nint > 0)
  
  ! alpha det is list of Nint 64-bit ints
  ! k is index of the int where iorb is found
  ! l is index of the bit where iorb is found
  k = ishft(iorb-1,-bit_kind_shift)+1
  ASSERT (k >0)
  l = iorb - ishft(k-1,bit_kind_shift)-1
  ASSERT (l >= 0)
  is_filled = btest(key_ref(k),l)  
end

subroutine orb_is_filled(key_ref,iorb,ispin,Nint,is_filled)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! apply creation operator to key_ref
  ! add electron with spin ispin to orbital with index iorb
  ! output resulting det and phase in key_new and phase
  END_DOC
  integer, intent(in)            :: iorb, ispin, Nint
  integer(bit_kind), intent(in) :: key_ref(Nint,2)
  logical, intent(out) :: is_filled
  
  integer                        :: k,l
  
  ASSERT (iorb > 0)
  ASSERT (ispin > 0)
  ASSERT (ispin < 3)
  ASSERT (Nint > 0)
  
  ! alpha det is list of Nint 64-bit ints
  ! k is index of the int where iorb is found
  ! l is index of the bit where iorb is found
  k = ishft(iorb-1,-bit_kind_shift)+1
  ASSERT (k >0)
  l = iorb - ishft(k-1,bit_kind_shift)-1
  ASSERT (l >= 0)
  is_filled = btest(key_ref(k,ispin),l)  
end

subroutine ac_operator_phase(key_new,key_ref,iorb,ispin,Nint,phase)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! apply creation operator to key_ref
  ! add electron with spin ispin to orbital with index iorb
  ! output resulting det and phase in key_new and phase
  END_DOC
  integer, intent(in)            :: iorb, ispin, Nint
  integer(bit_kind), intent(in) :: key_ref(Nint,2)
  integer(bit_kind), intent(out) :: key_new(Nint,2)
  double precision, intent(out) :: phase
  
  integer                        :: k,l,i

  double precision, parameter :: p(0:1) = (/ 1.d0, -1.d0 /)
  
  ASSERT (iorb > 0)
  ASSERT (ispin > 0)
  ASSERT (ispin < 3)
  ASSERT (Nint > 0)
  
  key_new=key_ref

  ! alpha det is list of Nint 64-bit ints
  ! k is index of the int where iorb is found
  ! l is index of the bit where iorb is found
  k = ishft(iorb-1,-bit_kind_shift)+1
  ASSERT (k >0)
  l = iorb - ishft(k-1,bit_kind_shift)-1
  ASSERT (l >= 0)
  key_new(k,ispin) = ibset(key_new(k,ispin),l)
  
  integer(bit_kind) :: parity_filled

  ! I assume here that the ordering is all alpha spinorbs and then all beta spinorbs
  ! If we add an alpha electron, parity is not affected by beta part of determinant
  !   (only need number of alpha occupied orbs below iorb)
 
  ! If we add a beta electron, the parity is affected by alpha part
  !    (need total number of occupied alpha orbs (all of which come before beta)
  !      and total number of beta occupied orbs below iorb)

  if (ispin==1) then
    parity_filled=0_bit_kind
  else
    parity_filled=iand(elec_alpha_num,1_bit_kind)
  endif

  ! get parity due to orbs in other ints (with lower indices)
  do i=1,k-1
    parity_filled = iand(popcnt(key_ref(i,ispin)),parity_filled)
  enddo
  
  ! get parity due to orbs in same int as iorb
  ! ishft(1_bit_kind,l)-1 has its l rightmost bits set to 1, other bits set to 0
  parity_filled = iand(popcnt(iand(ishft(1_bit_kind,l)-1,key_ref(k,ispin))),parity_filled)
  phase = p(iand(1_bit_kind,parity_filled))

end

subroutine a_operator_phase(key_new,key_ref,iorb,ispin,Nint,phase)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! apply annihilation operator to key_ref
  ! remove electron with spin ispin to orbital with index iorb
  ! output resulting det and phase in key_new and phase
  END_DOC
  integer, intent(in)            :: iorb, ispin, Nint
  integer(bit_kind), intent(in) :: key_ref(Nint,2)
  integer(bit_kind), intent(out) :: key_new(Nint,2)
  double precision, intent(out) :: phase
  
  integer                        :: k,l,i

  double precision, parameter :: p(0:1) = (/ 1.d0, -1.d0 /)
  
  ASSERT (iorb > 0)
  ASSERT (ispin > 0)
  ASSERT (ispin < 3)
  ASSERT (Nint > 0)
  
  key_new=key_ref

  ! alpha det is list of Nint 64-bit ints
  ! k is index of the int where iorb is found
  ! l is index of the bit where iorb is found
  k = ishft(iorb-1,-bit_kind_shift)+1
  ASSERT (k >0)
  l = iorb - ishft(k-1,bit_kind_shift)-1
  ASSERT (l >= 0)
  key_new(k,ispin) = ibclr(key_new(k,ispin),l)
  
  integer(bit_kind) :: parity_filled

  ! I assume here that the ordering is all alpha spinorbs and then all beta spinorbs
  ! If we add an alpha electron, parity is not affected by beta part of determinant
  !   (only need number of alpha occupied orbs below iorb)
 
  ! If we add a beta electron, the parity is affected by alpha part
  !    (need total number of occupied alpha orbs (all of which come before beta)
  !      and total number of beta occupied orbs below iorb)

  if (ispin==1) then
    parity_filled=0_bit_kind
  else
    parity_filled=iand(elec_alpha_num,1_bit_kind)
  endif

  ! get parity due to orbs in other ints (with lower indices)
  do i=1,k-1
    parity_filled = iand(popcnt(key_ref(i,ispin)),parity_filled)
  enddo
  
  ! get parity due to orbs in same int as iorb
  ! ishft(1_bit_kind,l)-1 has its l rightmost bits set to 1, other bits set to 0
  parity_filled = iand(popcnt(iand(ishft(1_bit_kind,l)-1,key_ref(k,ispin))),parity_filled)
  phase = p(iand(1_bit_kind,parity_filled))

end
