
subroutine get_list_hp_banned_ab(tmp_det,N_hp,exc_is_banned,spin_hp,sign_hp,idx_hp,nint,all_banned)
  implicit none
  BEGIN_DOC
  ! input determinant tmp_det and list of single holes/particles
  ! for each hole/particle, determine whether it is filled/empty in tmp_det
  ! return which are disallowed in exc_is_banned
  ! if all are banned, set all_banned to true
  END_DOC
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
  BEGIN_DOC
  ! input spindeterminant tmp_spindet and list of single holes/particles
  ! tmp_spindet is only one spin part of a full det, with spin==ispin
  ! for each hole/particle, determine whether it is filled/empty in tmp_det
  ! return which are disallowed in exc_is_banned
  ! if all are banned, set all_banned to true
  END_DOC
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
  BEGIN_DOC
  ! input determinant tmp_det and list of single holes/particles
  ! for each hole/particle, determine whether it is filled/empty in tmp_det
  ! return which are disallowed in exc_is_banned
  ! if all are banned, set all_banned to true
  ! only consider tmp_det(1:N_int, ispin)
  END_DOC
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
  !  todo: not finished
  ! determine whether iorb is filled in key_ref
  ! iorb is specified by int and bit locations within the determinant
  END_DOC
  integer, intent(in)            :: iorb_int, iorb_bit, ispin, Nint
  integer(bit_kind), intent(in) :: key_ref(Nint,2)
  logical, intent(out) :: is_filled
  
  integer                        :: k,l
  
  ASSERT (iorb > 0)
  ASSERT (ispin > 0)
  ASSERT (ispin < 3)
  ASSERT (Nint > 0)
  
  ! k is index of the int where iorb is found
  ! l is index of the bit where iorb is found
  k = ishft(iorb-1,-bit_kind_shift)+1
  ASSERT (k >0)
  l = iorb - ishft(k-1,bit_kind_shift)-1
  ASSERT (l >= 0)
  is_filled = btest(key_ref(iorb_int,ispin),iorb_bit)  
end

subroutine get_orb_bit_int(iorb,Nint,iorb_bit,iorb_int)
  BEGIN_DOC
  !  todo: not finished
  ! get int and bit corresponding to orbital index iorb 
  END_DOC
  use bitmasks
  implicit none
  iorb_int = ishft(iorb-1,-bit_kind_shift)+1
  iorb_bit = iorb - ishft(iorb_int-1,bit_kind_shift)-1
end

subroutine orb_is_filled_single_spin(key_ref,iorb,Nint,is_filled)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! determine whether iorb is filled in key_ref
  ! key_ref is single alpha or beta determinant
  END_DOC
  integer, intent(in)            :: iorb, Nint
  integer(bit_kind), intent(in) :: key_ref(Nint)
  logical, intent(out) :: is_filled
  
  integer                        :: k,l
  
  ASSERT (iorb > 0)
  ASSERT (Nint > 0)
  
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
  ! determine whether iorb, ispin is filled in key_ref
  ! key_ref has alpha and beta parts
  END_DOC
  integer, intent(in)            :: iorb, ispin, Nint
  integer(bit_kind), intent(in) :: key_ref(Nint,2)
  logical, intent(out) :: is_filled
  
  integer                        :: k,l
  
  ASSERT (iorb > 0)
  ASSERT (ispin > 0)
  ASSERT (ispin < 3)
  ASSERT (Nint > 0)
  
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
