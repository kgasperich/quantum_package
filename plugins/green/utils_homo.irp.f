
subroutine get_homo_banned_ab(tmp_det,is_banned,spin_hp,sign_hp,idx_hp,nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! input determinant tmp_det and list of single holes/particles
  ! for each hole/particle, determine whether it is filled/empty in tmp_det
  ! return which are disallowed in exc_is_banned
  ! if all are banned, set all_banned to true
  END_DOC
  integer, intent(in) :: nint
  integer, intent(in) :: spin_hp, idx_hp
  double precision, intent(in) :: sign_hp
  integer(bit_kind), intent(in) :: tmp_det(nint,2)
  logical, intent(out) :: is_banned
  
  integer :: i
  logical :: is_filled

    call orb_is_filled(tmp_det,idx_hp,spin_hp,nint,is_filled)
    is_banned = ((sign_hp.gt.0).eqv.is_filled)
end

subroutine get_homo_banned_single_spin(tmp_spindet,is_banned,spin_hp,sign_hp,idx_hp,ispin,nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! input spindeterminant tmp_spindet and list of single holes/particles
  ! tmp_spindet is only one spin part of a full det, with spin==ispin
  ! for each hole/particle, determine whether it is filled/empty in tmp_det
  ! return which are disallowed in exc_is_banned
  ! if all are banned, set all_banned to true
  END_DOC
  integer, intent(in) :: ispin, nint
  integer, intent(in) :: spin_hp, idx_hp
  double precision, intent(in) :: sign_hp
  integer(bit_kind), intent(in) :: tmp_spindet(nint)
  logical, intent(out) :: is_banned
  
  integer :: i
  logical :: is_filled

    if (spin_hp.eq.ispin) then
      call orb_is_filled_single_spin(tmp_spindet,idx_hp,nint,is_filled)
      is_banned = ((sign_hp.gt.0).eqv.is_filled)
    else
      is_banned = .False.
    endif
end

subroutine get_homo_banned_spin(tmp_det,is_banned,spin_hp,sign_hp,idx_hp,ispin,nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! input determinant tmp_det and list of single holes/particles
  ! for each hole/particle, determine whether it is filled/empty in tmp_det
  ! return which are disallowed in exc_is_banned
  ! if all are banned, set all_banned to true
  ! only consider tmp_det(1:N_int, ispin)
  END_DOC
  integer, intent(in) :: ispin, nint
  integer, intent(in) :: spin_hp, idx_hp
  double precision, intent(in) :: sign_hp
  integer(bit_kind), intent(in) :: tmp_det(nint,2)
  logical, intent(out) :: is_banned

  integer(bit_kind) :: spindet(nint)
  
  integer :: i
  logical :: is_filled
  spindet(1:nint) = tmp_det(1:nint,ispin)

    if (spin_hp.eq.ispin) then
      call orb_is_filled(tmp_det,idx_hp,ispin,nint,is_filled)
      is_banned = ((sign_hp.gt.0).eqv.is_filled)
    else
      is_banned = .False.
    endif
end

  
