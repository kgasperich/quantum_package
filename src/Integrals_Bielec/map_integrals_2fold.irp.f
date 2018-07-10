use map_module

subroutine bielec_integrals_index_2fold(i,j,k,l,i1)
  use map_module
  implicit none
  BEGIN_DOC
  ! index of <ij|kl> with 2-fold symmetry
  END_DOC
  integer, intent(in)            :: i,j,k,l
  integer(key_kind), intent(out) :: i1
  integer(key_kind)              :: p,q,i2
  if (i==k) then
    p=i*i
  else if (i.lt.k) then
    p=(k-1)*(k-1)+2*i-mod(k+1,2)
  else
    p=(i-1)*(i-1)+2*k-mod(i,2)
  endif
  if (j==l) then
    q=j*j
  else if (j.lt.l) then
    q=(l-1)*(l-1)+2*j-mod(l+1,2)
  else
    q=(j-1)*(j-1)+2*l-mod(j,2)
  endif
  i1 = min(p,q)
  i2 = max(p,q)
  i1 = i1+ishft(i2*i2-i2,-1)
end

subroutine bielec_integrals_index_reverse_2fold(i,j,k,l,i1)
  use map_module
  implicit none
  integer, intent(out)           :: i(2),j(2),k(2),l(2)
  integer(key_kind), intent(in)  :: i1
  integer(key_kind)              :: i2,i3,p,q
  i = 0
  i2   = ceiling(0.5d0*(dsqrt(8.d0*dble(i1)+1.d0)-1.d0))
  i3   = i1 - ishft(i2*i2-i2,-1)
  p = ceiling(dsqrt(dble(i2)))
  q = ceiling(0.5d0*(dble(i2)-dble((p-1)*(p-1))))
  if (mod(i2,2)==0) then
    l(1)=p
    j(1)=q
  else
    j(1)=p
    l(1)=q
  endif
  p = ceiling(dsqrt(dble(i3)))
  q = ceiling(0.5d0*(dble(i3)-dble((p-1)*(p-1))))
  if (mod(i3,2)==0) then
    k(1)=p
    i(1)=q
  else
    i(1)=p
    k(1)=q
  endif
              !ijkl
  i(2) = j(1) !jilk
  j(2) = i(1)
  k(2) = l(1)
  l(2) = k(1)

  integer :: ii, jj
  do ii=2,2
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
      call bielec_integrals_index_2fold(i(ii),j(ii),k(ii),l(ii),i2)
      if (i1 /= i2) then
        print *,  i1, i2
        print *,  i(ii), j(ii), k(ii), l(ii)
        stop 'bielec_integrals_index_reverse_2fold failed'
      endif
    endif
  enddo
end

subroutine bielec_integrals_index_reverse_4fold(i,j,k,l,i1)
  use map_module
  implicit none
  integer, intent(out)           :: i(4),j(4),k(4),l(4)
  integer(key_kind), intent(in)  :: i1
  integer(key_kind)              :: i2,i3,p,q
  i = 0
  i2   = ceiling(0.5d0*(dsqrt(8.d0*dble(i1)+1.d0)-1.d0))
  i3   = i1 - ishft(i2*i2-i2,-1)
  p = ceiling(dsqrt(dble(i2)))
  q = ceiling(0.5d0*(dble(i2)-dble((p-1)*(p-1))))
  if (mod(i2,2)==0) then
    l(1)=p
    j(1)=q
  else
    j(1)=p
    l(1)=q
  endif
  p = ceiling(dsqrt(dble(i3)))
  q = ceiling(0.5d0*(dble(i3)-dble((p-1)*(p-1))))
  if (mod(i3,2)==0) then
    k(1)=p
    i(1)=q
  else
    i(1)=p
    k(1)=q
  endif
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
  do ii=1,4
    if (i(ii) /= 0) then
      call complex_bielec_integrals_index(i(ii),j(ii),k(ii),l(ii),i2)
      if (i1 /= i2) then
        print *,  i1, i2
        print *,  i(ii), j(ii), k(ii), l(ii)
        stop 'bielec_integrals_index_reverse_4fold failed (this has not been implemented yet)'
      endif
    endif
  enddo
end




!! MO Map
!! ======

BEGIN_PROVIDER [ type(map_type), mo_integrals_map_2fold ]
  implicit none
  BEGIN_DOC
  ! MO integrals
  END_DOC
  integer(key_kind)              :: key_max
  integer(map_size_kind)         :: sze
  call bielec_integrals_index_2fold(mo_tot_num,mo_tot_num,mo_tot_num,mo_tot_num,key_max)
  sze = key_max
  call map_init(mo_integrals_map_2fold,sze)
  print*, 'MO map initialized (2-fold): ', sze
END_PROVIDER


subroutine insert_into_mo_integrals_map_2fold(n_integrals,                 &
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
  call map_update(mo_integrals_map_2fold, buffer_i, buffer_values, n_integrals, thr)
end

 BEGIN_PROVIDER [ integer*4, mo_integrals_cache_min_2fold ]
&BEGIN_PROVIDER [ integer*4, mo_integrals_cache_max_2fold ]
&BEGIN_PROVIDER [ integer*8, mo_integrals_cache_min_8_2fold ]
&BEGIN_PROVIDER [ integer*8, mo_integrals_cache_max_8_2fold ]
 implicit none
 BEGIN_DOC
 ! Min and max values of the MOs for which the integrals are in the cache
 END_DOC
 mo_integrals_cache_min_8_2fold = max(1_8,elec_alpha_num - 63_8)
 mo_integrals_cache_max_8_2fold = min(int(mo_tot_num,8),mo_integrals_cache_min_8_2fold+127_8)
 mo_integrals_cache_min_2fold   = max(1,elec_alpha_num - 63)
 mo_integrals_cache_max_2fold   = min(mo_tot_num,mo_integrals_cache_min_2fold+127)

END_PROVIDER

BEGIN_PROVIDER [ double precision, mo_integrals_cache_2fold, (0_8:128_8*128_8*128_8*128_8) ]
 implicit none
 BEGIN_DOC
 ! Cache of MO integrals for fast access
 END_DOC
 PROVIDE mo_bielec_integrals_in_map_2fold
 integer*8                      :: i,j,k,l
 integer*4                      :: i4,j4,k4,l4
 integer*8                      :: ii
 integer(key_kind)              :: idx,idx1,idx2
 real(integral_kind)            :: integral
 FREE ao_integrals_cache
 !$OMP PARALLEL DO PRIVATE (i,j,k,l,i4,j4,k4,l4,idx,ii,integral)
 do l=mo_integrals_cache_min_8_2fold,mo_integrals_cache_max_8_2fold
   l4 = int(l,4)
   do k=mo_integrals_cache_min_8_2fold,mo_integrals_cache_max_8_2fold
     k4 = int(k,4)
     do j=mo_integrals_cache_min_8_2fold,mo_integrals_cache_max_8_2fold
       j4 = int(j,4)
       do i=mo_integrals_cache_min_8_2fold,mo_integrals_cache_max_8_2fold
         i4 = int(i,4)
         !DIR$ FORCEINLINE
         call bielec_integrals_index_2fold(i4,j4,k4,l4,idx1)
         call bielec_integrals_index_2fold(k4,l4,i4,j4,idx2)
         !DIR$ FORCEINLINE
         idx=min(idx1,idx2)
         call map_get(mo_integrals_map_2fold,idx,integral)
         ii = l-mo_integrals_cache_min_8_2fold
         ii = ior( ishft(ii,7), k-mo_integrals_cache_min_8_2fold)
         ii = ior( ishft(ii,7), j-mo_integrals_cache_min_8_2fold)
         ii = ior( ishft(ii,7), i-mo_integrals_cache_min_8_2fold)
         mo_integrals_cache_2fold(ii) = integral
       enddo
     enddo
   enddo
 enddo
 !$OMP END PARALLEL DO

END_PROVIDER


double precision function get_mo_bielec_integral_2fold(i,j,k,l,map)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns one integral <ij|kl> in the MO basis
  END_DOC
  integer, intent(in)            :: i,j,k,l
  integer(key_kind)              :: idx,idx1,idx2
  integer                        :: ii
  integer*8                      :: ii_8
  type(map_type), intent(inout)  :: map
  real(integral_kind)            :: tmp
  PROVIDE mo_bielec_integrals_in_map_2fold mo_integrals_cache_2fold
  ii = l-mo_integrals_cache_min_2fold
  ii = ior(ii, k-mo_integrals_cache_min_2fold)
  ii = ior(ii, j-mo_integrals_cache_min_2fold)
  ii = ior(ii, i-mo_integrals_cache_min_2fold)
  if (iand(ii, -128) /= 0) then
    !DIR$ FORCEINLINE
    call bielec_integrals_index_2fold(i,j,k,l,idx1)
    call bielec_integrals_index_2fold(k,l,i,j,idx2)
    !DIR$ FORCEINLINE
    idx=min(idx1,idx2)
    call map_get(map,idx,tmp)
    get_mo_bielec_integral_2fold = dble(tmp)
  else
    ii_8 = int(l,8)-mo_integrals_cache_min_8_2fold
    ii_8 = ior( ishft(ii_8,7), int(k,8)-mo_integrals_cache_min_8_2fold)
    ii_8 = ior( ishft(ii_8,7), int(j,8)-mo_integrals_cache_min_8_2fold)
    ii_8 = ior( ishft(ii_8,7), int(i,8)-mo_integrals_cache_min_8_2fold)
    get_mo_bielec_integral_2fold = mo_integrals_cache_2fold(ii_8)
  endif
end


double precision function mo_bielec_integral_2fold(i,j,k,l)
  implicit none
  BEGIN_DOC
  ! Returns one integral <ij|kl> in the MO basis
  END_DOC
  integer, intent(in)            :: i,j,k,l
  double precision               :: get_mo_bielec_integral_2fold
  PROVIDE mo_bielec_integrals_in_map_2fold mo_integrals_cache_2fold
  !DIR$ FORCEINLINE
  PROVIDE mo_bielec_integrals_in_map_2fold
  mo_bielec_integral_2fold = get_mo_bielec_integral_2fold(i,j,k,l,mo_integrals_map_2fold)
  return
end

subroutine get_mo_bielec_integrals_2fold(j,k,l,sze,out_val,map)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns multiple integrals <ij|kl> in the MO basis, all
  ! i for j,k,l fixed.
  END_DOC
  integer, intent(in)            :: j,k,l, sze
  double precision, intent(out)  :: out_val(sze)
  type(map_type), intent(inout)  :: map
  integer                        :: i
  integer(key_kind)              :: tmp_hash1,tmp_hash2
  integer(key_kind)              :: hash(sze)
  PROVIDE mo_bielec_integrals_in_map_2fold
  
  do i=1,sze
    !DIR$ FORCEINLINE
    call bielec_integrals_index_2fold(i,j,k,l,tmp_hash1)
    call bielec_integrals_index_2fold(k,l,i,j,tmp_hash2)
    hash(i) = min(tmp_hash1,tmp_hash2)
  enddo
  
  if (integral_kind == 8) then
    call map_get_many(map, hash(:), out_val, sze)
  else
    call map_get_many(map, hash(:), tmp_val, sze)
    ! Conversion to double precision 
    do i=1,sze
      out_val(i) = dble(tmp_val(i))
    enddo
  endif
end
!TODO: remove map argument.  always uses same map
subroutine get_mo_bielec_integrals_ij_2fold(k,l,sze,out_array,map)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns multiple integrals <ij|kl> in the MO basis, all
  ! i(1)j(2) 1/r12 k(1)l(2)
  ! i, j for k,l fixed.
  END_DOC
  integer, intent(in)            :: k,l, sze
  double precision, intent(out)  :: out_array(sze,sze)
  type(map_type), intent(inout)  :: map
  integer                        :: i,j,kk,ll,m
  integer(key_kind),allocatable  :: hash(:)
  integer(key_kind)              :: tmp_hash1,tmp_hash2
  integer  ,allocatable          :: pairs(:,:), iorder(:)
  real(integral_kind), allocatable :: tmp_val(:)

  PROVIDE mo_bielec_integrals_in_map_2fold
  allocate (hash(sze*sze), pairs(2,sze*sze), &
  iorder(sze*sze), tmp_val(sze*sze))
  
  kk=0
  out_array = 0.d0
  do j=1,sze
   do i=1,sze
    kk += 1
    !DIR$ FORCEINLINE
    call bielec_integrals_index_2fold(i,j,k,l,tmp_hash1)
    call bielec_integrals_index_2fold(k,l,i,j,tmp_hash2)
    hash(kk) = min(tmp_hash1,tmp_hash2)
    pairs(1,kk) = i
    pairs(2,kk) = j
   enddo
  enddo

  logical :: integral_is_in_map
  if (key_kind == 8) then
    call i8radix_sort(hash,iorder,kk,-1)
  else if (key_kind == 4) then
    call iradix_sort(hash,iorder,kk,-1)
  else if (key_kind == 2) then
    call i2radix_sort(hash,iorder,kk,-1)
  endif

  call map_get_many(mo_integrals_map_2fold, hash, tmp_val, kk)

  do ll=1,kk
    m = iorder(ll)
    i = pairs(1,m)
    j = pairs(2,m)
    out_array(i,j) = tmp_val(ll)
  enddo  

  deallocate(pairs,hash,iorder,tmp_val)
end

subroutine get_mo_bielec_integrals_coulomb_ii_2fold(k,l,sze,out_val,map)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns multiple integrals <ki|li> 
  ! k(1)i(2) 1/r12 l(1)i(2) :: out_val(i1)
  ! for k,l fixed.
  END_DOC
  integer, intent(in)            :: k,l, sze
  double precision, intent(out)  :: out_val(sze)
  type(map_type), intent(inout)  :: map
  integer                        :: i
  integer(key_kind)              :: tmp_hash1,tmp_hash2
  integer(key_kind)              :: hash(sze)
  real(integral_kind)            :: tmp_val(sze)
  PROVIDE mo_bielec_integrals_in_map_2fold
  
  do i=1,sze
    !DIR$ FORCEINLINE
    call bielec_integrals_index_2fold(k,i,l,i,tmp_hash1)
    call bielec_integrals_index_2fold(l,i,k,i,tmp_hash2)
    hash(i) = min(tmp_hash1,tmp_hash2)
  enddo
  
  if (integral_kind.eq.8) then
    call map_get_many(mo_integrals_map_2fold, hash, out_val, sze)
  else
    call map_get_many(mo_integrals_map_2fold, hash, tmp_val, sze)
    do i=1,sze
      out_val(i) = dble(tmp_val(i))
    enddo
  endif
end

subroutine get_mo_bielec_integrals_exch_ii_2fold(k,l,sze,out_val,map)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns multiple integrals <ki|il> 
  ! k(1)i(2) 1/r12 i(1)l(2) :: out_val(i1)
  ! for k,l fixed.
  END_DOC
  integer, intent(in)            :: k,l, sze
  double precision, intent(out)  :: out_val(sze)
  type(map_type), intent(inout)  :: map
  integer                        :: i
  integer(key_kind)              :: tmp_hash1,tmp_hash2
  integer(key_kind)              :: hash(sze)
  real(integral_kind)            :: tmp_val(sze)
  PROVIDE mo_bielec_integrals_in_map_2fold
  
  do i=1,sze
    !DIR$ FORCEINLINE
    call bielec_integrals_index_2fold(k,i,i,l,tmp_hash1)
    call bielec_integrals_index_2fold(i,l,k,i,tmp_hash2)
    hash(i) = min(tmp_hash1,tmp_hash2)
  enddo
  
  if (integral_kind.eq.8) then
    call map_get_many(mo_integrals_map_2fold, hash, out_val, sze)
  else
    call map_get_many(mo_integrals_map_2fold, hash, tmp_val, sze)
    do i=1,sze
      out_val(i) = dble(tmp_val(i))
    enddo
  endif
end


integer*8 function get_mo_map_size_2fold()
  implicit none
  BEGIN_DOC
  ! Return the number of elements in the MO map
  END_DOC
  get_mo_map_size_2fold = mo_integrals_map_2fold % n_elements
end


