double precision function diag_S_mat_elem(key_i,Nint)
  implicit none
  use bitmasks
  include 'Utils/constants.include.F'

  integer                        :: Nint
  integer(bit_kind), intent(in)  :: key_i(Nint,2)
  BEGIN_DOC
! Returns <i|S^2|i>
! (returns number of unpaired alpha electrons in |i>; is this correct?)
  END_DOC
  integer                        :: nup, i
  integer(bit_kind)              :: xorvec(N_int_max)

  do i=1,Nint
    xorvec(i) = xor(key_i(i,1),key_i(i,2))
  enddo

  do i=1,Nint
    xorvec(i) = iand(xorvec(i),key_i(i,1))
  enddo

  nup = 0
  do i=1,Nint
    if (xorvec(i) /= 0_bit_kind) then
      nup += popcnt(xorvec(i))
    endif
  enddo
  diag_S_mat_elem = dble(nup)

end

subroutine get_s2(key_i,key_j,Nint,s2)
  implicit none
  use bitmasks
  BEGIN_DOC
  ! Returns <S^2>
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key_i(Nint,2)
  integer(bit_kind), intent(in)  :: key_j(Nint,2)
  double precision, intent(out)  :: s2
  integer                        :: exc(0:2,2,2)
  integer                        :: degree
  double precision               :: phase_spsm
  integer                        :: nup, i
  
  s2 = 0.d0
  !$FORCEINLINE
  call get_excitation_degree(key_i,key_j,degree,Nint)
  select case (degree)
    case(2)
      call get_double_excitation(key_j,key_i,exc,phase_spsm,Nint)
      if (exc(0,1,1) == 1) then   ! Mono alpha + mono-beta
        if ( (exc(1,1,1) == exc(1,2,2)).and.(exc(1,1,2) == exc(1,2,1)) ) then
          s2 =  -phase_spsm
        endif
      endif
    case(0)
      double precision, external :: diag_S_mat_elem
      !DIR$ FORCEINLINE
      s2 = diag_S_mat_elem(key_i,Nint)
  end select
end

BEGIN_PROVIDER [ double precision, S_z ]
&BEGIN_PROVIDER [ double precision, S_z2_Sz ]
 implicit none
 BEGIN_DOC
! z component of the Spin
 END_DOC

 S_z = 0.5d0*dble(elec_alpha_num-elec_beta_num)
 S_z2_Sz = S_z*(S_z-1.d0)

END_PROVIDER

BEGIN_PROVIDER [ double precision, expected_s2]
 implicit none
 BEGIN_DOC
! Expected value of S2 : S*(S+1)
 END_DOC
   logical :: has_expected_s2

   call ezfio_has_determinants_expected_s2(has_expected_s2)
   if (has_expected_s2) then
     call ezfio_get_determinants_expected_s2(expected_s2)
   else
     double precision :: S
     S = (elec_alpha_num-elec_beta_num)*0.5d0 
     expected_s2 = S * (S+1.d0)
   endif

END_PROVIDER 

BEGIN_PROVIDER [ double precision, s2_values, (N_states) ]
 implicit none
 BEGIN_DOC
! array of the averaged values of the S^2 operator on the various states
 END_DOC
 integer :: i
 call u_0_S2_u_0(s2_values,psi_coef,n_det,psi_det,N_int,N_states,psi_det_size)

END_PROVIDER



subroutine u_0_S2_u_0(e_0,u_0,n,keys_tmp,Nint,N_st,sze_8)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes e_0 = <u_0|S2|u_0>/<u_0|u_0>
  !
  ! n : number of determinants
  !
  END_DOC
  integer, intent(in)            :: n,Nint, N_st, sze_8
  double precision, intent(out)  :: e_0(N_st)
  complex*16, intent(in)   :: u_0(sze_8,N_st)
  integer(bit_kind),intent(in)   :: keys_tmp(Nint,2,n)
  
  complex*16, allocatable  :: v_0(:,:)
  double precision               :: u_dot_u_complex
  complex*16                     :: u_dot_v_complex
  integer :: i,j
  allocate (v_0(sze_8,N_st))
  
  call S2_u_0_nstates(v_0,u_0,n,keys_tmp,Nint,N_st,sze_8)
  do i=1,N_st
    e_0(i) = dreal(u_dot_v_complex(u_0(1,i),v_0(1,i),n))/u_dot_u_complex(u_0(1,i),n) + S_z2_Sz
  enddo
end



subroutine S2_u_0(v_0,u_0,n,keys_tmp,Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes v_0 = S^2|u_0>
  !
  ! n : number of determinants
  !
  END_DOC
  integer, intent(in)            :: n,Nint
  complex*16, intent(out)  :: v_0(n)
  complex*16, intent(in)   :: u_0(n)
  integer(bit_kind),intent(in)   :: keys_tmp(Nint,2,n)
  call S2_u_0_nstates(v_0,u_0,n,keys_tmp,Nint,1,n)
end

subroutine S2_u_0_nstates(v_0,u_0,n,keys_tmp,Nint,N_st,sze_8)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes v_0  = S^2|u_0>
  !
  ! n : number of determinants
  !
  END_DOC
  integer, intent(in)            :: N_st,n,Nint, sze_8
  complex*16, intent(out)  :: v_0(sze_8,N_st)
  complex*16, intent(in)   :: u_0(sze_8,N_st)
  integer(bit_kind),intent(in)   :: keys_tmp(Nint,2,n)
  double precision               :: s2_tmp
  complex*16, allocatable  :: vt(:,:)
  integer                        :: i,j,k,l, jj,ii
  integer                        :: i0, j0
  
  integer, allocatable           :: shortcut(:,:), sort_idx(:,:)
  integer(bit_kind), allocatable :: sorted(:,:,:), version(:,:,:)
  integer(bit_kind)              :: sorted_i(Nint)
  
  integer                        :: sh, sh2, ni, exa, ext, org_i, org_j, endi, istate
  

  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  ASSERT (n>0)
  PROVIDE ref_bitmask_energy 

  allocate (shortcut(0:n+1,2), sort_idx(n,2), sorted(Nint,n,2), version(Nint,n,2))
  v_0 = 0.d0

  call sort_dets_ab_v(keys_tmp, sorted(1,1,1), sort_idx(1,1), shortcut(0,1), version(1,1,1), n, Nint)
  call sort_dets_ba_v(keys_tmp, sorted(1,1,2), sort_idx(1,2), shortcut(0,2), version(1,1,2), n, Nint)
  
  !$OMP PARALLEL DEFAULT(NONE)                                       &
      !$OMP PRIVATE(i,s2_tmp,j,k,jj,vt,ii,sh,sh2,ni,exa,ext,org_i,org_j,endi,sorted_i,istate)&
      !$OMP SHARED(n,u_0,keys_tmp,Nint,v_0,sorted,shortcut,sort_idx,version,N_st,sze_8)
  allocate(vt(sze_8,N_st))
  vt = 0.d0
  
  !$OMP DO SCHEDULE(dynamic)
  do sh=1,shortcut(0,1)
    do sh2=sh,shortcut(0,1)
      exa = 0
      do ni=1,Nint
        exa = exa + popcnt(xor(version(ni,sh,1), version(ni,sh2,1)))
      end do
      if(exa > 2) then
        cycle
      end if
      
      do i=shortcut(sh,1),shortcut(sh+1,1)-1
        org_i = sort_idx(i,1)
        if(sh==sh2) then
          endi = i-1
        else
          endi = shortcut(sh2+1,1)-1
        end if
        do ni=1,Nint
          sorted_i(ni) = sorted(ni,i,1)
        enddo
        
        do j=shortcut(sh2,1),endi
          org_j = sort_idx(j,1)
          ext = exa
          do ni=1,Nint
            ext = ext + popcnt(xor(sorted_i(ni), sorted(ni,j,1)))
          end do
          if(ext <= 4) then
            call get_s2(keys_tmp(1,1,org_i),keys_tmp(1,1,org_j),Nint,s2_tmp)
            do istate=1,N_st
              vt (org_i,istate) = vt (org_i,istate) + s2_tmp*u_0(org_j,istate)
              vt (org_j,istate) = vt (org_j,istate) + s2_tmp*u_0(org_i,istate)
            enddo
          endif
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO NOWAIT
  
  !$OMP DO SCHEDULE(dynamic)
  do sh=1,shortcut(0,2)
    do i=shortcut(sh,2),shortcut(sh+1,2)-1
      org_i = sort_idx(i,2)
      do j=shortcut(sh,2),i-1
        org_j = sort_idx(j,2)
        ext = 0
        do ni=1,Nint
          ext = ext + popcnt(xor(sorted(ni,i,2), sorted(ni,j,2)))
        end do
        if(ext == 4) then
          call get_s2(keys_tmp(1,1,org_i),keys_tmp(1,1,org_j),Nint,s2_tmp)
          do istate=1,N_st
            vt (org_i,istate) = vt (org_i,istate) + s2_tmp*u_0(org_j,istate)
            vt (org_j,istate) = vt (org_j,istate) + s2_tmp*u_0(org_i,istate)
          enddo
        end if
      end do
    end do
  enddo
  !$OMP END DO NOWAIT
  
  do istate=1,N_st
    do i=n,1,-1
      !$OMP ATOMIC
      v_0(i,istate) = v_0(i,istate) + vt(i,istate)
    enddo
  enddo

  deallocate(vt)
  !$OMP END PARALLEL
  
  do i=1,n
    call get_s2(keys_tmp(1,1,i),keys_tmp(1,1,i),Nint,s2_tmp)
    do istate=1,N_st
      v_0(i,istate) += s2_tmp * u_0(i,istate)
    enddo
  enddo

  deallocate (shortcut, sort_idx, sorted, version)
end

