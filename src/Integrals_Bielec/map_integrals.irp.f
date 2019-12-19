use map_module

!! AO Map
!! ======

subroutine idx2_sq(i,j,ij)
  implicit none
  integer, intent(in)  :: i,j
  integer, intent(out) :: ij
  if (i<j) then
    ij=(j-1)*(j-1)+2*i-mod(j+1,2)
  else if (i>j) then
    ij=(i-1)*(i-1)+2*j-mod(i,2)
  else
    ij=i*i
  endif
end

subroutine idx2_tri_int(i,j,ij)
  implicit none
  integer, intent(in)  :: i,j
  integer, intent(out) :: ij
  integer :: p,q
  p = max(i,j)
  q = min(i,j)
  ij = q+ishft(p*p-p,-1)
end

subroutine idx2_tri_key(i,j,ij)
  use map_module
  implicit none
  integer, intent(in)  :: i,j
  integer(key_kind), intent(out) :: ij
  integer(key_kind) :: p,q
  p = max(i,j)
  q = min(i,j)
  ij = q+ishft(p*p-p,-1)
end


BEGIN_PROVIDER [ type(map_type), ao_integrals_map ]
  implicit none
  BEGIN_DOC
  ! AO integrals
  END_DOC
  integer(key_kind)              :: key_max
  integer(map_size_kind)         :: sze
  call bielec_integrals_index_2fold(ao_num,ao_num,ao_num,ao_num,key_max)
  sze = key_max
  call map_init(ao_integrals_map,sze)
  print*,  'AO map initialized : ', sze
END_PROVIDER

subroutine bielec_integrals_index(i,j,k,l,i1)
  use map_module
  implicit none
  integer, intent(in)            :: i,j,k,l
  integer(key_kind), intent(out) :: i1
  integer                        :: p,q
  call idx2_tri_int(i,k,p)
  call idx2_tri_int(j,l,q)
  call idx2_tri_key(p,q,i1)
end

subroutine bielec_integrals_index_2fold(i,j,k,l,i1)
  use map_module
  implicit none
  integer, intent(in)            :: i,j,k,l
  integer(key_kind), intent(out) :: i1
  integer                        :: ik,jl

  call idx2_sq(i,k,ik)
  call idx2_sq(j,l,jl)
  call idx2_tri_key(ik,jl,i1)
end

subroutine idx2_tri_rev_key(i,k,ik)
  use map_module
  BEGIN_DOC
  !return i<=k
  END_DOC
  integer(key_kind), intent(in) :: ik
  integer, intent(out)          :: i,k
  integer(key_kind) :: tmp_k
  k = ceiling(0.5d0*(dsqrt(8.d0*dble(ik)+1.d0)-1.d0))
  tmp_k = k
  i = int(ik - ishft(tmp_k*tmp_k-tmp_k,-1))
end

subroutine idx2_tri_rev_int(i,k,ik)
  BEGIN_DOC
  !return i<=k
  END_DOC
  integer, intent(in)           :: ik
  integer, intent(out)          :: i,k
  k = ceiling(0.5d0*(dsqrt(8.d0*dble(ik)+1.d0)-1.d0))
  i = int(ik - ishft(k*k-k,-1))
end

subroutine bielec_integrals_index_reverse(i,j,k,l,i1)
  use map_module
  implicit none
  integer, intent(out)           :: i(8),j(8),k(8),l(8)
  integer(key_kind), intent(in)  :: i1
  integer(key_kind)              :: i0
  integer                        :: i2,i3
  i = 0
  call idx2_tri_rev_key(i3,i2,i1)
  call idx2_tri_rev_int(j(1),l(1),i2)
  call idx2_tri_rev_int(i(1),k(1),i3)

              !ijkl
  i(2) = i(1) !ilkj
  j(2) = l(1)
  k(2) = k(1)
  l(2) = j(1)

  i(3) = k(1) !kjil
  j(3) = j(1)
  k(3) = i(1)
  l(3) = l(1)

  i(4) = k(1) !klij
  j(4) = l(1)
  k(4) = i(1)
  l(4) = j(1)

  i(5) = j(1) !jilk
  j(5) = i(1)
  k(5) = l(1)
  l(5) = k(1)

  i(6) = j(1) !jkli
  j(6) = k(1)
  k(6) = l(1)
  l(6) = i(1)

  i(7) = l(1) !lijk
  j(7) = i(1)
  k(7) = j(1)
  l(7) = k(1)

  i(8) = l(1) !lkji
  j(8) = k(1)
  k(8) = j(1)
  l(8) = i(1)

  integer :: ii, jj
  do ii=2,8
    do jj=1,ii-1
      if ( (i(ii) == i(jj)).and. &
           (j(ii) == j(jj)).and. &
           (k(ii) == k(jj)).and. &
           (l(ii) == l(jj)) ) then
         i(ii) = 0
         exit
      endif
    enddo
  enddo
  do ii=1,8
    if (i(ii) /= 0) then
      call bielec_integrals_index(i(ii),j(ii),k(ii),l(ii),i0)
      if (i1 /= i0) then
        print *,  i1, i0
        print *,  i(ii), j(ii), k(ii), l(ii)
        stop 'bielec_integrals_index_reverse failed'
      endif
    endif
  enddo
end

subroutine idx2_sq_rev(i,k,ik)
  BEGIN_DOC
  ! reverse square compound index
  END_DOC
!  p = ceiling(dsqrt(dble(ik)))
!  q = ceiling(0.5d0*(dble(ik)-dble((p-1)*(p-1))))
!  if (mod(ik,2)==0) then
!    k=p
!    i=q
!  else
!    i=p
!    k=q
!  endif
  integer, intent(in)           :: ik
  integer, intent(out)          :: i,k
  integer                       :: pq(0:1),i1,i2
  pq(0) = ceiling(dsqrt(dble(ik)))
  pq(1) = ceiling(0.5d0*(dble(ik)-dble((pq(0)-1)*(pq(0)-1))))
  i1=mod(ik,2)
  i2=mod(ik+1,2)

  k=pq(i1)
  i=pq(i2)
end

subroutine bielec_integrals_index_reverse_2fold(i,j,k,l,i1)
  use map_module
  implicit none
  integer, intent(out)           :: i(2),j(2),k(2),l(2)
  integer(key_kind), intent(in)  :: i1
  integer(key_kind)              :: i0
  integer                        :: i2,i3
  i = 0
  call idx2_tri_rev_key(i3,i2,i1)

  call idx2_sq_rev(j(1),l(1),i2)
  call idx2_sq_rev(i(1),k(1),i3)

              !ijkl
  i(2) = j(1) !jilk
  j(2) = i(1)
  k(2) = l(1)
  l(2) = k(1)

!  i(3) = k(1) !klij   complex conjugate
!  j(3) = l(1)
!  k(3) = i(1)
!  l(3) = j(1)
!
!  i(4) = l(1) !lkji   complex conjugate
!  j(4) = k(1)
!  k(4) = j(1)
!  l(4) = i(1)

  integer :: ii
  if ( (i(1)==i(2)).and. &
       (j(1)==j(2)).and. &
       (k(1)==k(2)).and. &
       (l(1)==l(2)) ) then
    i(2) = 0
  endif
  do ii=1,2
    if (i(ii) /= 0) then
      call bielec_integrals_index_2fold(i(ii),j(ii),k(ii),l(ii),i0)
      if (i1 /= i0) then
        print *,  i1, i0
        print *,  i(ii), j(ii), k(ii), l(ii)
        stop 'bielec_integrals_index_reverse_2fold failed'
      endif
    endif
  enddo
end

subroutine bielec_integrals_index_reverse_4fold(i,j,k,l,i1,key2)
  use map_module
  implicit none
  integer, intent(out)           :: i(4),j(4),k(4),l(4)
  integer(key_kind), intent(in)  :: i1
  integer(key_kind), intent(out) :: key2
  integer(key_kind)              :: i0
  integer                        :: i2,i3
  i = 0
  call idx2_tri_rev_key(i3,i2,i1)

  call idx2_sq_rev(j(1),l(1),i2)
  call idx2_sq_rev(i(1),k(1),i3)

              !ijkl
  i(2) = j(1) !jilk
  j(2) = i(1)
  k(2) = l(1)
  l(2) = k(1)

  i(3) = k(1) !klij   complex conjugate
  j(3) = l(1)
  k(3) = i(1)
  l(3) = j(1)

  i(4) = l(1) !lkji   complex conjugate
  j(4) = k(1)
  k(4) = j(1)
  l(4) = i(1)

  integer :: ii, jj
  do ii=2,4
    do jj=1,ii-1
      if ( (i(ii) == i(jj)).and. &
           (j(ii) == j(jj)).and. &
           (k(ii) == k(jj)).and. &
           (l(ii) == l(jj)) ) then
         i(ii) = 0
         exit
      endif
    enddo
  enddo
  do ii=1,2
    if (i(ii) /= 0) then
      call bielec_integrals_index_2fold(i(ii),j(ii),k(ii),l(ii),i0)
      if (i1 /= i0) then
        print *,  i1, i0
        print *,  i(ii), j(ii), k(ii), l(ii)
        stop 'bielec_integrals_index_reverse_4fold failed'
      endif
    endif
  enddo
  call bielec_integrals_index_2fold(k(1),l(1),i(1),j(1),key2)
end

 BEGIN_PROVIDER [ integer, ao_integrals_cache_min ]
&BEGIN_PROVIDER [ integer, ao_integrals_cache_max ]
 implicit none
 BEGIN_DOC
 ! Min and max values of the AOs for which the integrals are in the cache
 END_DOC
 ao_integrals_cache_min = max(1,ao_num - 63)
 ao_integrals_cache_max = ao_num

END_PROVIDER

BEGIN_PROVIDER [ complex*16, ao_integrals_cache, (0:64*64*64*64) ]
 implicit none
 BEGIN_DOC
 ! Cache of AO integrals for fast access
 END_DOC
 PROVIDE ao_bielec_integrals_in_map
 integer                        :: i,j,k,l,ii
 integer(key_kind)              :: idx1,idx2
 integer(key_kind)              :: idx_re,idx_im
 real(integral_kind)            :: tmp_re, tmp_im
 complex(integral_kind)         :: integral
 !$OMP PARALLEL DO PRIVATE (i,j,k,l,idx1,idx2,ii,tmp_re,tmp_im,integral,idx_re,idx_im)
 do l=ao_integrals_cache_min,ao_integrals_cache_max
   do k=ao_integrals_cache_min,ao_integrals_cache_max
     do j=ao_integrals_cache_min,ao_integrals_cache_max
       do i=ao_integrals_cache_min,ao_integrals_cache_max
         !DIR$ FORCEINLINE
         call bielec_integrals_index_2fold(i,j,k,l,idx1)
         call bielec_integrals_index_2fold(k,l,i,j,idx2)
         !DIR$ FORCEINLINE
         idx_re = min(idx1,idx2)
         idx_im = max(idx1,idx2)
         call map_get(ao_integrals_map,idx_re,tmp_re)
         if (idx_re /= idx_im) then
           call map_get(ao_integrals_map,idx_im,tmp_im)
           if (idx1 < idx2) then
             integral = dcmplx(tmp_re,tmp_im)
           else
             integral = dcmplx(tmp_re,-tmp_im)
           endif
         else
           tmp_im = 0.d0
           integral = dcmplx(tmp_re,tmp_im)
         endif

         ii = l-ao_integrals_cache_min
         ii = ior( ishft(ii,6), k-ao_integrals_cache_min)
         ii = ior( ishft(ii,6), j-ao_integrals_cache_min)
         ii = ior( ishft(ii,6), i-ao_integrals_cache_min)
         ao_integrals_cache(ii) = integral
       enddo
     enddo
   enddo
 enddo
 !$OMP END PARALLEL DO

END_PROVIDER


complex*16 function get_ao_bielec_integral(i,j,k,l,map) result(result)
  use map_module
  implicit none
  BEGIN_DOC
  ! Gets one AO bi-electronic integral from the AO map <ij|kl>
  END_DOC
  integer, intent(in)            :: i,j,k,l
  integer(key_kind)              :: idx1,idx2
  integer(key_kind)              :: idx_re,idx_im
  type(map_type), intent(inout)  :: map
  integer                        :: ii
  real(integral_kind)            :: tmp_re,tmp_im
  complex(integral_kind)         :: tmp_int
  PROVIDE ao_bielec_integrals_in_map ao_integrals_cache ao_integrals_cache_min
  !DIR$ FORCEINLINE
  if (ao_overlap_abs(i,k)*ao_overlap_abs(j,l) < ao_integrals_threshold ) then
    tmp_int = (0.d0,0.d0)
!  else if (ao_bielec_integral_schwartz(i,k)*ao_bielec_integral_schwartz(j,l) < ao_integrals_threshold) then
!    tmp_int = (0.d0,0.d0)
  else
    ii = l-ao_integrals_cache_min
    ii = ior(ii, k-ao_integrals_cache_min)
    ii = ior(ii, j-ao_integrals_cache_min)
    ii = ior(ii, i-ao_integrals_cache_min)
    if (iand(ii, -64) /= 0) then
      !DIR$ FORCEINLINE
      call bielec_integrals_index_2fold(i,j,k,l,idx1)
      call bielec_integrals_index_2fold(k,l,i,j,idx2)
      !DIR$ FORCEINLINE
      idx_re = min(idx1,idx2)
      idx_im = max(idx1,idx2)
      call map_get(map,idx_re,tmp_re)
      if (idx_re /= idx_im) then
        call map_get(map,idx_im,tmp_im)
        if (idx1 < idx2) then
          tmp_int = dcmplx(tmp_re,tmp_im)
        else
          tmp_int = dcmplx(tmp_re,-tmp_im)
        endif
      else
        tmp_im = 0.d0
        tmp_int = dcmplx(tmp_re,tmp_im)
      endif
    else
      ii = l-ao_integrals_cache_min
      ii = ior( ishft(ii,6), k-ao_integrals_cache_min)
      ii = ior( ishft(ii,6), j-ao_integrals_cache_min)
      ii = ior( ishft(ii,6), i-ao_integrals_cache_min)
      tmp_int = ao_integrals_cache(ii)
    endif
  endif
  result = tmp_int
end


subroutine get_ao_bielec_integrals(j,k,l,sze,out_val)
  use map_module
  BEGIN_DOC
  ! Gets multiple AO bi-electronic integral from the AO map .
  ! All i are retrieved for j,k,l fixed.
  END_DOC
  implicit none
  integer, intent(in)            :: j,k,l, sze
  complex(integral_kind), intent(out) :: out_val(sze)
  
  integer                        :: i
  integer(key_kind)              :: hash
  double precision               :: thresh
  PROVIDE ao_bielec_integrals_in_map ao_integrals_map
  thresh = ao_integrals_threshold
  
  if (ao_overlap_abs(j,l) < thresh) then
    out_val = (0.d0,0.d0)
    return
  endif
  
  complex*16 :: get_ao_bielec_integral
  do i=1,sze
    out_val(i) = get_ao_bielec_integral(i,j,k,l,ao_integrals_map)
  enddo
  
end

subroutine get_ao_bielec_integrals_non_zero(j,k,l,sze,out_val,out_val_index,non_zero_int)
  use map_module
  implicit none
  BEGIN_DOC
  ! Gets multiple AO bi-electronic integral from the AO map .
  ! All non-zero i are retrieved for j,k,l fixed.
  END_DOC
  integer, intent(in)            :: j,k,l, sze
  complex(integral_kind), intent(out) :: out_val(sze)
  integer, intent(out)           :: out_val_index(sze),non_zero_int
  
  integer                        :: i
  double precision               :: thresh, schwarz_jl
  complex*16               :: tmp
  PROVIDE ao_bielec_integrals_in_map
  thresh = ao_integrals_threshold
  
  non_zero_int = 0
  if (ao_overlap_abs(j,l) < thresh) then
    out_val = (0.d0,0.d0)
    return
  endif
 
  non_zero_int = 0
  schwarz_jl = ao_bielec_integral_schwartz(j,l)
  do i=1,sze
    !DIR$ FORCEINLINE
    if (ao_bielec_integral_schwartz(i,k)*schwarz_jl < thresh) then
      cycle
    endif
    complex*16 :: get_ao_bielec_integral
    tmp = get_ao_bielec_integral(i,j,k,l,ao_integrals_map)
    if (cdabs(tmp) < thresh ) cycle
    non_zero_int = non_zero_int+1
    out_val_index(non_zero_int) = i
    out_val(non_zero_int) = tmp
  enddo
  
end


function get_ao_map_size()
  implicit none
  integer (map_size_kind) :: get_ao_map_size
  BEGIN_DOC
  ! Returns the number of elements in the AO map
  END_DOC
  get_ao_map_size = ao_integrals_map % n_elements
end

subroutine clear_ao_map
  implicit none
  BEGIN_DOC
  ! Frees the memory of the AO map
  END_DOC
  call map_deinit(ao_integrals_map)
  FREE ao_integrals_map
end


!! MO Map
!! ======

BEGIN_PROVIDER [ type(map_type), mo_integrals_map ]
  implicit none
  BEGIN_DOC
  ! MO integrals
  END_DOC
  integer(key_kind)              :: key_max
  integer(map_size_kind)         :: sze
  call bielec_integrals_index_2fold(mo_tot_num,mo_tot_num,mo_tot_num,mo_tot_num,key_max)
  sze = key_max
  call map_init(mo_integrals_map,sze)
  print*, 'complex MO map initialized: ', sze
END_PROVIDER

subroutine insert_into_ao_integrals_map(n_integrals,buffer_i, buffer_values)
  use map_module
  implicit none
  BEGIN_DOC
  ! Create new entry into AO map
  ! TODO: if using 3-index integrals (density fitting) to construct AO ints, change this to map_update
  END_DOC
  
  integer, intent(in)                :: n_integrals
  integer(key_kind), intent(inout)   :: buffer_i(n_integrals)
  real(integral_kind), intent(inout) :: buffer_values(n_integrals)
  
  call map_append(ao_integrals_map, buffer_i, buffer_values, n_integrals)
end

subroutine insert_into_mo_integrals_map(n_integrals,                 &
      buffer_i, buffer_values, thr)
  use map_module
  implicit none
  
  BEGIN_DOC
  ! Create new entry into MO map, or accumulate in an existing entry
  END_DOC
  
  integer, intent(in)                :: n_integrals
  integer(key_kind), intent(inout)   :: buffer_i(n_integrals)
  real(integral_kind), intent(inout) :: buffer_values(n_integrals)
  real(integral_kind), intent(in)    :: thr
  call map_update(mo_integrals_map, buffer_i, buffer_values, n_integrals, thr)
end

 BEGIN_PROVIDER [ integer*4, mo_integrals_cache_min ]
&BEGIN_PROVIDER [ integer*4, mo_integrals_cache_max ]
&BEGIN_PROVIDER [ integer*8, mo_integrals_cache_min_8 ]
&BEGIN_PROVIDER [ integer*8, mo_integrals_cache_max_8 ]
 implicit none
 BEGIN_DOC
 ! Min and max values of the MOs for which the integrals are in the cache
 END_DOC
 mo_integrals_cache_min_8 = max(1_8,elec_alpha_num - 63_8)
 mo_integrals_cache_max_8 = min(int(mo_tot_num,8),mo_integrals_cache_min_8+127_8)
 mo_integrals_cache_min   = max(1,elec_alpha_num - 63)
 mo_integrals_cache_max   = min(mo_tot_num,mo_integrals_cache_min+127)

END_PROVIDER

BEGIN_PROVIDER [ complex*16, mo_integrals_cache, (0_8:128_8*128_8*128_8*128_8) ]
 implicit none
 BEGIN_DOC
 ! Cache of MO integrals for fast access
 END_DOC
 !TODO: make this more efficient (can save at least a factor of 2)
 PROVIDE mo_bielec_integrals_in_map
 integer*8                      :: i,j,k,l
 integer*4                      :: i4,j4,k4,l4
 integer*8                      :: ii
 integer(key_kind)              :: idx1,idx2
 real(integral_kind)            :: tmp1,tmp2
 complex(integral_kind)            :: integral
 FREE ao_integrals_cache
 do l=mo_integrals_cache_min_8,mo_integrals_cache_max_8
   l4 = int(l,4)
   do k=mo_integrals_cache_min_8,mo_integrals_cache_max_8
     k4 = int(k,4)
     do j=mo_integrals_cache_min_8,mo_integrals_cache_max_8
       j4 = int(j,4)
       do i=mo_integrals_cache_min_8,mo_integrals_cache_max_8
         i4 = int(i,4)
         !DIR$ FORCEINLINE
         call bielec_integrals_index_2fold(i4,j4,k4,l4,idx1)
         call bielec_integrals_index_2fold(k4,l4,i4,j4,idx2)
         !DIR$ FORCEINLINE
         call map_get(mo_integrals_map,idx1,tmp1)
         if (idx1==idx2) then
           integral=dcmplx(tmp1,0.d0)
         else if (idx1.lt.idx2) then
           call map_get(mo_integrals_map,idx2,tmp2)
           integral=dcmplx(tmp1,tmp2)
         else 
           call map_get(mo_integrals_map,idx2,tmp2)
           integral=dcmplx(tmp2,-tmp1)
         endif

         ii = l-mo_integrals_cache_min_8
         ii = ior( ishft(ii,7), k-mo_integrals_cache_min_8)
         ii = ior( ishft(ii,7), j-mo_integrals_cache_min_8)
         ii = ior( ishft(ii,7), i-mo_integrals_cache_min_8)
         mo_integrals_cache(ii) = integral
       enddo
     enddo
   enddo
 enddo
END_PROVIDER

complex*16 function get_mo_bielec_integral(i,j,k,l,map)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns one complex integral <ij|kl> in the MO basis
  END_DOC
  integer, intent(in)            :: i,j,k,l
  integer(key_kind)              :: idx1,idx2
  integer                        :: ii
  integer*8                      :: ii_8
  type(map_type), intent(inout)  :: map
  real(integral_kind)            :: tmp1,tmp2

  PROVIDE mo_bielec_integrals_in_map mo_integrals_cache
  ii = l-mo_integrals_cache_min
  ii = ior(ii, k-mo_integrals_cache_min)
  ii = ior(ii, j-mo_integrals_cache_min)
  ii = ior(ii, i-mo_integrals_cache_min)
  if (iand(ii, -128) /= 0) then
    !DIR$ FORCEINLINE
    call bielec_integrals_index_2fold(i,j,k,l,idx1)
    call bielec_integrals_index_2fold(k,l,i,j,idx2)

    !DIR$ FORCEINLINE
    call map_get(map,idx1,tmp1)
    call map_get(map,idx2,tmp2) !todo:move this inside conditional below
    if (idx1.eq.idx2) then
      get_mo_bielec_integral = dcmplx(tmp1,0.d0)
    else if (idx1.lt.idx2) then
      get_mo_bielec_integral = dcmplx(tmp1,tmp2)
    else
      get_mo_bielec_integral = dcmplx(tmp2,-tmp1)
    endif
  else
    ii_8 = int(l,8)-mo_integrals_cache_min_8
    ii_8 = ior( ishft(ii_8,7), int(k,8)-mo_integrals_cache_min_8)
    ii_8 = ior( ishft(ii_8,7), int(j,8)-mo_integrals_cache_min_8)
    ii_8 = ior( ishft(ii_8,7), int(i,8)-mo_integrals_cache_min_8)
    get_mo_bielec_integral = mo_integrals_cache(ii_8)
  endif
  return
end


complex*16 function mo_bielec_integral(i,j,k,l)
  implicit none
  BEGIN_DOC
  ! Returns one integral <ij|kl> in the MO basis
  END_DOC
  integer, intent(in)            :: i,j,k,l
  complex*16                 :: get_mo_bielec_integral
  PROVIDE mo_bielec_integrals_in_map mo_integrals_cache
  !DIR$ FORCEINLINE
  PROVIDE mo_bielec_integrals_in_map
  mo_bielec_integral = get_mo_bielec_integral(i,j,k,l,mo_integrals_map)
  return
end

subroutine get_mo_bielec_integrals(j,k,l,sze,out_val,map)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns multiple integrals <ij|kl> in the MO basis, all
  ! i for j,k,l fixed.
  END_DOC
  integer, intent(in)            :: j,k,l, sze
  complex(integral_kind), intent(out)    :: out_val(sze)
  real(integral_kind)            :: out_val1(sze),out_val2(sze)
  type(map_type), intent(inout)  :: map
  integer                        :: i
  integer(key_kind)              :: hash(2,sze)
  real(integral_kind)            :: tmp1,tmp2
  real(integral_kind)         :: tmp_val1(sze),tmp_val2(sze)
  PROVIDE mo_bielec_integrals_in_map
  
  do i=1,sze
    !DIR$ FORCEINLINE
    call bielec_integrals_index_2fold(i,j,k,l,hash(1,i))
    call bielec_integrals_index_2fold(k,l,i,j,hash(2,i))
  enddo
  
  if (key_kind == 8) then
    call map_get_many(map, hash(1,:), out_val1, sze)
    call map_get_many(map, hash(2,:), out_val2, sze)
    do i=1,sze
      if (hash(1,i)==hash(2,i)) then
        out_val(i) = dcmplx(out_val1(i),0.d0)
      else if (hash(1,i).lt.hash(2,i)) then
        out_val(i) = dcmplx(out_val1(i),out_val2(i))
      else
        out_val(i) = dcmplx(out_val2(i),-out_val1(i))
      endif
    enddo
  else
    call map_get_many(map, hash(1,:), tmp_val1, sze)
    call map_get_many(map, hash(2,:), tmp_val2, sze)
    ! Conversion to double precision 
    do i=1,sze
      if (hash(1,i)==hash(2,i)) then
        out_val(i) = dcmplx(tmp_val1(i))
      else if (hash(1,i).lt.hash(2,i)) then
        out_val(i) = dcmplx(tmp_val1(i),tmp_val2(i))
      else
        out_val(i) = dcmplx(tmp_val2(i),-tmp_val1(i))
      endif
    enddo
  endif
end

!TODO: make complex implementation more efficient?
subroutine get_mo_bielec_integrals_ij(k,l,sze,out_array,map)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns multiple integrals <ij|kl> in the MO basis, all
  ! i(1)j(2) 1/r12 k(1)l(2)
  ! i, j for k,l fixed.
  END_DOC
  integer, intent(in)            :: k,l, sze
  complex*16, intent(out)  :: out_array(sze,sze)
  type(map_type), intent(inout)  :: map
  integer                        :: i,j,kk,ll,m1,m2
  integer(key_kind)              :: idx1,idx2
  real(integral_kind)            :: p,q
  integer(key_kind),allocatable  :: hash1(:),hash2(:)
  integer  ,allocatable          :: pairs(:,:), iorder1(:),iorder2(:)
  integer  ,allocatable          :: jorder1(:),jorder2(:)
  real(integral_kind), allocatable :: tmp_val1(:),tmp_val2(:)

  PROVIDE mo_bielec_integrals_in_map
  allocate (hash1(sze*sze), hash2(sze*sze), pairs(2,sze*sze), &
  iorder1(sze*sze), iorder2(sze*sze), jorder1(sze*sze), jorder2(sze*sze), &
  tmp_val1(sze*sze), tmp_val2(sze*sze))
  
  kk=0
  out_array = 0.d0
  do j=1,sze
   do i=1,sze
    kk += 1
    !DIR$ FORCEINLINE
    call bielec_integrals_index_2fold(i,j,k,l,hash1(kk))
    call bielec_integrals_index_2fold(k,l,i,j,hash2(kk))
    pairs(1,kk) = i
    pairs(2,kk) = j
    iorder1(kk) = kk
    iorder2(kk) = kk
    jorder1(kk) = kk
    jorder2(kk) = kk
   enddo
  enddo

  logical :: integral_is_in_map
!  if (key_kind == 8) then
    call i8radix_sort(hash1,iorder1,kk,-1)
    call i8radix_sort(hash2,iorder2,kk,-1)
!  else if (key_kind == 4) then
!    call iradix_sort(hash1,iorder1,kk,-1)
!    call iradix_sort(hash2,iorder2,kk,-1)
!  else if (key_kind == 2) then
!    call i2radix_sort(hash1,iorder1,kk,-1)
!    call i2radix_sort(hash2,iorder2,kk,-1)
!  endif
  call iradix_sort(iorder1,jorder1,kk,-1)
  call iradix_sort(iorder2,jorder2,kk,-1)

  call map_get_many(mo_integrals_map, hash1, tmp_val1, kk)
  call map_get_many(mo_integrals_map, hash2, tmp_val2, kk)

  do ll=1,kk
    m1=jorder1(ll)
    m2=jorder2(ll)
    idx1=hash1(m1)
    idx2=hash2(m2)
    i=pairs(1,ll)
    j=pairs(2,ll)
    p=tmp_val1(m1)
    q=tmp_val2(m2)
    if (idx1==idx2) then
      out_array(i,j) = dcmplx(p,0.d0)
    else if (idx1.lt.idx2) then
      out_array(i,j) = dcmplx(p,q)
    else
      out_array(i,j) = dcmplx(q,-p)
    endif
  enddo  

  deallocate(pairs,hash1,hash2,iorder1,iorder2,jorder1,jorder2,tmp_val1,tmp_val2)
end

!TODO: make complex implementation more efficient?
!      
subroutine get_mo_bielec_integrals_coulomb_ii(k,l,sze,out_val,map)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns multiple integrals <ki|li> 
  ! k(1)i(2) 1/r12 l(1)i(2) :: out_val(i1)
  ! for k,l fixed.
  END_DOC
  integer, intent(in)            :: k,l, sze
  complex*16, intent(out)        :: out_val(sze)
  type(map_type), intent(inout)  :: map
  integer                        :: i
  integer(key_kind)              :: idx1,idx2
  integer(key_kind)              :: hash1(sze),hash2(sze)
  real(integral_kind)            :: tmp_val1(sze),tmp_val2(sze)
  PROVIDE mo_bielec_integrals_in_map
  
  do i=1,sze
    !DIR$ FORCEINLINE
    call bielec_integrals_index_2fold(k,i,l,i,hash1(i))
    call bielec_integrals_index_2fold(l,i,k,i,hash2(i))
  enddo
  
  call map_get_many(map, hash1, tmp_val1, sze)
  call map_get_many(map, hash2, tmp_val2, sze)
  do i=1,sze
    idx1=hash1(i)
    idx2=hash2(i)
    !TODO: this set of conditionals only needs to be evaluated once
    if (idx1==idx2) then
      out_val(i) = dcmplx(tmp_val1(i),0.d0)
    else if (idx1.lt.idx2) then
      out_val(i) = dcmplx(tmp_val1(i),tmp_val2(i))
    else
      out_val(i) = dcmplx(tmp_val2(i),-tmp_val1(i))
    endif
  enddo
end

!TODO: make complex implementation more efficient 
subroutine get_mo_bielec_integrals_exch_ii(k,l,sze,out_val,map)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns multiple integrals <ki|il> 
  ! k(1)i(2) 1/r12 i(1)l(2) :: out_val(i1)
  ! for k,l fixed.
  END_DOC
  integer, intent(in)            :: k,l, sze
  complex*16, intent(out)        :: out_val(sze)
  type(map_type), intent(inout)  :: map
  integer                        :: i
  integer(key_kind)              :: idx1,idx2
  integer(key_kind)              :: hash1(sze),hash2(sze)
  real(integral_kind)            :: tmp_val1(sze),tmp_val2(sze)
  PROVIDE mo_bielec_integrals_in_map
  
  do i=1,sze
    !DIR$ FORCEINLINE
    call bielec_integrals_index_2fold(k,i,i,l,hash1(i))
    call bielec_integrals_index_2fold(i,l,k,i,hash2(i))
  enddo
  
  call map_get_many(map, hash1, tmp_val1, sze)
  call map_get_many(map, hash2, tmp_val2, sze)
  do i=1,sze
    idx1=hash1(i)
    idx2=hash2(i)
    !TODO: this set of conditionals may only need to be evaluated a few times
    ! (probably better to replace with a different set of conditionals)
    ! which branch to take depends on whether i is less than both l and k,
    ! between l and k, or greater than l and k
    if (idx1==idx2) then
      out_val(i) = dcmplx(tmp_val1(i),0.d0)
    else if (idx1.lt.idx2) then
      out_val(i) = dcmplx(tmp_val1(i),tmp_val2(i))
    else
      out_val(i) = dcmplx(tmp_val2(i),-tmp_val1(i))
    endif
  enddo
end


integer*8 function get_mo_map_size()
  implicit none
  BEGIN_DOC
  ! Return the number of elements in the MO map
  END_DOC
  get_mo_map_size = mo_integrals_map % n_elements
end

BEGIN_TEMPLATE

subroutine dump_$ao_integrals(filename)
  use map_module
  implicit none
  BEGIN_DOC
  ! Save to disk the $ao integrals
  END_DOC
  character*(*), intent(in)      :: filename
  integer(cache_key_kind), pointer :: key(:)
  real(integral_kind), pointer   :: val(:)
  integer*8                      :: i,j, n
  if (.not.mpi_master) then
    return
  endif
  call ezfio_set_work_empty(.False.)
  open(unit=66,file=filename,FORM='unformatted')
  write(66) integral_kind, key_kind
  write(66) $ao_integrals_map%sorted, $ao_integrals_map%map_size,    &
      $ao_integrals_map%n_elements
  do i=0_8,$ao_integrals_map%map_size
    write(66) $ao_integrals_map%map(i)%sorted, $ao_integrals_map%map(i)%map_size,&
        $ao_integrals_map%map(i)%n_elements
  enddo
  do i=0_8,$ao_integrals_map%map_size
    key => $ao_integrals_map%map(i)%key
    val => $ao_integrals_map%map(i)%value
    n = $ao_integrals_map%map(i)%n_elements
    write(66) (key(j), j=1,n), (val(j), j=1,n)
  enddo
  close(66)
  
end

IRP_IF COARRAY
subroutine communicate_$ao_integrals()
  use map_module
  implicit none
  BEGIN_DOC
  ! Communicate the $ao integrals with co-array
  END_DOC
  integer(cache_key_kind), pointer :: key(:)
  real(integral_kind), pointer   :: val(:)
  integer*8                      :: i,j, k, nmax
  integer*8, save                :: n[*]
  integer                        :: copy_n

  real(integral_kind), allocatable            :: buffer_val(:)[:]
  integer(cache_key_kind), allocatable        :: buffer_key(:)[:]
  real(integral_kind), allocatable            :: copy_val(:)
  integer(key_kind), allocatable              :: copy_key(:)

  n = 0_8
  do i=0_8,$ao_integrals_map%map_size
    n = max(n,$ao_integrals_map%map(i)%n_elements)
  enddo
  sync all
  nmax = 0_8
  do j=1,num_images()
    nmax = max(nmax,n[j])
  enddo
  allocate( buffer_key(nmax)[*], buffer_val(nmax)[*])
  allocate( copy_key(nmax), copy_val(nmax))
  do i=0_8,$ao_integrals_map%map_size
    key => $ao_integrals_map%map(i)%key
    val => $ao_integrals_map%map(i)%value
    n = $ao_integrals_map%map(i)%n_elements
    do j=1,n
      buffer_key(j) = key(j)
      buffer_val(j) = val(j)
    enddo
    sync all
    do j=1,num_images()
      if (j /= this_image()) then
        copy_n = n[j]
        do k=1,copy_n
          copy_val(k) = buffer_val(k)[j]
          copy_key(k) = buffer_key(k)[j]
          copy_key(k) = copy_key(k)+ishft(i,-map_shift)
        enddo
        call map_append($ao_integrals_map, copy_key, copy_val, copy_n )
      endif
    enddo
    sync all
  enddo
  deallocate( buffer_key, buffer_val, copy_val, copy_key)
  
end
IRP_ENDIF 


integer function load_$ao_integrals(filename)
  implicit none
  BEGIN_DOC
  ! Read from disk the $ao integrals
  END_DOC
  character*(*), intent(in)      :: filename
  integer*8                      :: i
  integer(cache_key_kind), pointer :: key(:)
  real(integral_kind), pointer   :: val(:)
  integer                        :: iknd, kknd
  integer*8                      :: n, j
  load_$ao_integrals = 1
  open(unit=66,file=filename,FORM='unformatted',STATUS='UNKNOWN')
  read(66,err=98,end=98) iknd, kknd
  if (iknd /= integral_kind) then
    print *,  'Wrong integrals kind in file :', iknd
    stop 1
  endif
  if (kknd /= key_kind) then
    print *,  'Wrong key kind in file :', kknd
    stop 1
  endif
  read(66,err=98,end=98) $ao_integrals_map%sorted, $ao_integrals_map%map_size,&
      $ao_integrals_map%n_elements
  do i=0_8, $ao_integrals_map%map_size
    read(66,err=99,end=99) $ao_integrals_map%map(i)%sorted,          &
        $ao_integrals_map%map(i)%map_size, $ao_integrals_map%map(i)%n_elements
    call cache_map_reallocate($ao_integrals_map%map(i),$ao_integrals_map%map(i)%map_size)
  enddo
  do i=0_8, $ao_integrals_map%map_size
    key => $ao_integrals_map%map(i)%key
    val => $ao_integrals_map%map(i)%value
    n = $ao_integrals_map%map(i)%n_elements
    read(66,err=99,end=99) (key(j), j=1,n), (val(j), j=1,n)
  enddo
  call map_sort($ao_integrals_map)
  load_$ao_integrals = 0
  return
  99 continue
  call map_deinit($ao_integrals_map)
  98 continue
  stop 'Problem reading $ao_integrals_map file in work/'
  
end

SUBST [ ao_integrals_map, ao_integrals, ao_num ]
ao_integrals_map ; ao_integrals ; ao_num ;;
mo_integrals_map ; mo_integrals ; mo_tot_num ;;
END_TEMPLATE
