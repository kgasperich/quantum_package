
subroutine get_mask_phase_new(det1, pm, Nint)
  use bitmasks
  BEGIN_DOC
  ! phasemask copied from qp2
  ! return phasemask of det1 in pm
  END_DOC
  implicit none
  integer, intent(in) :: Nint
  integer(bit_kind), intent(in) :: det1(Nint,2)
  integer(bit_kind), intent(out) :: pm(Nint,2)
  integer(bit_kind) :: tmp1, tmp2
  integer :: i
  pm(1:Nint,1:2) = det1(1:Nint,1:2)
  tmp1 = 0_8
  tmp2 = 0_8
  do i=1,Nint
    pm(i,1) = ieor(pm(i,1), shiftl(pm(i,1), 1))
    pm(i,2) = ieor(pm(i,2), shiftl(pm(i,2), 1))
    pm(i,1) = ieor(pm(i,1), shiftl(pm(i,1), 2))
    pm(i,2) = ieor(pm(i,2), shiftl(pm(i,2), 2))
    pm(i,1) = ieor(pm(i,1), shiftl(pm(i,1), 4))
    pm(i,2) = ieor(pm(i,2), shiftl(pm(i,2), 4))
    pm(i,1) = ieor(pm(i,1), shiftl(pm(i,1), 8))
    pm(i,2) = ieor(pm(i,2), shiftl(pm(i,2), 8))
    pm(i,1) = ieor(pm(i,1), shiftl(pm(i,1), 16))
    pm(i,2) = ieor(pm(i,2), shiftl(pm(i,2), 16))
    pm(i,1) = ieor(pm(i,1), shiftl(pm(i,1), 32))
    pm(i,2) = ieor(pm(i,2), shiftl(pm(i,2), 32))
    pm(i,1) = ieor(pm(i,1), tmp1)
    pm(i,2) = ieor(pm(i,2), tmp2)
    if(iand(popcnt(det1(i,1)), 1) == 1) tmp1 = not(tmp1)
    if(iand(popcnt(det1(i,2)), 1) == 1) tmp2 = not(tmp2)
  end do
end subroutine

subroutine get_phase_hp(g_idx_int,g_idx_bit,g_spin,g_sign,det_in,g_det_phase,nint,n_g)
  implicit none
  integer, intent(in) :: nint,n_g
  integer, intent(in) :: g_idx_int(n_g), g_idx_bit(n_g),g_spin(n_g)
  double precision, intent(in) :: g_sign(n_g)
  integer(bit_kind), intent(in) :: det_in(nint,2)
  double precision, intent(out) :: g_det_phase(n_g)

  integer(bit_kind) :: tmp_spindet(nint), pm(nint,2)
  double precision, parameter :: phase_dble(0:1) = (/1.d0,-1.d0/)

  integer :: i
  logical :: is_allowed(n_g), all_banned, is_filled

  all_banned=.True.
  do i=1,n_g
    tmp_spindet(1:nint) = det_in(1:nint,g_spin(i))
    call spinorb_is_filled_int_bit(tmp_spindet,g_idx_int(i),g_idx_bit(i),nint,is_filled)
    is_allowed(i) = (.not.(((g_sign(i)<0).and.(.not.is_filled)).or.((g_sign(i)>0).and.(is_filled))))
    all_banned=(all_banned.and.(.not.is_allowed(i)))
  enddo

  if (all_banned)
    g_det_phase(:)=0.d0
  else
    call get_phase_mask_new(det_in,pm,nint)
    do i=1,n_g
      if (is_allowed(i)) then
        g_det_phase(i) = phase_dble(popcnt(iand(ibset(0_bit_kind,g_idx_bit(i)),pm(g_idx_int(i),g_spin(i)))))
      else
        g_det_phase(i)=0.d0
      endif
    enddo
  endif
end

subroutine get_homo_lumo(idx_homo,spin_homo,idx_lumo,spin_lumo)
  implicit none
  integer, intent(out) :: idx_homo, spin_homo, idx_lumo, spin_lumo

  ! dummy vals for testing
  idx_homo=1
  spin_homo=1
  idx_lumo=1
  spin_lumo=1
  
end

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

  
subroutine spinorb_is_filled_int_bit(key_ref,iorb_int,iorb_bit,Nint,is_filled)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! determine whether iorb is filled in key_ref
  ! iorb is specified by int and bit locations within the determinant
  END_DOC
  integer, intent(in)            :: iorb_int, iorb_bit, Nint
  integer(bit_kind), intent(in) :: key_ref(Nint)
  logical, intent(out) :: is_filled
  
  ASSERT (iorb_int > 0)
  ASSERT (iorb_bit >= 0)
  ASSERT (Nint > 0)
  is_filled = btest(key_ref(iorb_int),iorb_bit)
end

subroutine orb_is_filled_int_bit(key_ref,iorb_int,iorb_bit,ispin,Nint,is_filled)
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
  
  ASSERT (iorb_int > 0)
  ASSERT (iorb_bit >= 0)
  ASSERT (ispin > 0)
  ASSERT (ispin < 3)
  ASSERT (Nint > 0)
  is_filled = btest(key_ref(iorb_int,ispin),iorb_bit)
!  call spinorb_is_filled_int_bit(key_ref(1,ispin),iorb_int,iorb_bit,Nint,is_filled)
end

subroutine get_orb_int_bit(iorb,iorb_int,iorb_bit)
  BEGIN_DOC
  ! get int and bit corresponding to orbital index iorb 
  END_DOC
  use bitmasks
  implicit none
  integer, intent(in)            :: iorb
  integer, intent(out)            :: iorb_int, iorb_bit
  ASSERT (iorb > 0)
  iorb_int = ishft(iorb-1,-bit_kind_shift)+1
  ASSERT (iorb_int > 0)
  iorb_bit = iorb - ishft(iorb_int-1,bit_kind_shift)-1
  ASSERT (iorb_bit >= 0)
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
