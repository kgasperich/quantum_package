complex*16 function ao_bielec_integral(i,j,k,l)
  use map_module
  implicit none
  BEGIN_DOC
  !  integral of the AO basis <ik|jl> or (ij|kl)
  !     i(r1) j(r1) 1/r12 k(r2) l(r2)
  END_DOC

  integer,intent(in)             :: i,j,k,l
  complex*16                     :: get_ao_bielec_integral

  PROVIDE ao_bielec_integrals_in_map

  ao_bielec_integral = get_ao_bielec_integral(i,k,j,l,ao_integrals_map)
  return
end


BEGIN_PROVIDER [ logical, ao_bielec_integrals_in_map ]
  implicit none
  use f77_zmq
  use map_module
  BEGIN_DOC
  !  Map of Atomic integrals
  !     i(r1) j(r2) 1/r12 k(r1) l(r2)
  END_DOC
  
  PROVIDE read_ao_integrals disk_access_ao_integrals
  if (read_ao_integrals) then
    print*,'Reading the AO integrals'
    call map_load_from_disk(trim(ezfio_filename)//'/work/ao_ints',ao_integrals_map)
    print*, 'AO integrals provided'
    ao_bielec_integrals_in_map = .True.
    return
  else if (read_df_integral_array) then
    call ao_map_fill_from_df
    ao_bielec_integrals_in_map = .True.
    return
  else
    print*,'ERROR: complex AO integrals must be provided'
    stop
  endif
  
END_PROVIDER
 
subroutine ao_map_fill_from_df
  use map_module
  implicit none
  BEGIN_DOC
  ! fill ao bielec integral map using 3-index df integrals
  ! todo: add OMP directives, do smarter things with pointers/reshaping/transposing
  END_DOC

  integer :: i,k,j,l
  integer :: ki,kk,kj,kl
  integer :: ii,ik,ij,il
  integer :: kikk2,kjkl2,jl2,ik2

  complex*16,allocatable :: ints_ik(:,:), ints_jl(:,:), ints_ikjl_tmp(:,:), ints_ikjl(:,:,:,:)

  complex*16 :: integral
  integer                        :: n_integrals
  integer                        :: size_buffer
  integer(key_kind),allocatable  :: buffer_i1(:), buffer_i2(:)
  real(integral_kind),allocatable :: buffer_value1(:),buffer_value2(:)
  integer(key_kind)              :: tmp_idx1,tmp_idx2
  double precision               :: tmp_re,tmp_im

  size_buffer = min(ao_num*ao_num*ao_num,16000000)
  print*, 'Providing the ao_bielec integrals from 3-index df integrals'
  call ezfio_set_integrals_bielec_disk_access_ao_integrals('Write')

  !$OMP PARALLEL PRIVATE(i,k,j,l,ki,kk,kj,kl,ii,ik,ij,il,kikk2,kjkl2,jl2,ik2, &
      !$OMP  ints_ik, ints_jl, ints_ikjl_tmp, ints_ikjl, &
      !$OMP  n_integrals, buffer_i1, buffer_i2, buffer_value1, buffer_value2, &
      !$OMP  tmp_idx1, tmp_idx2, tmp_re, tmp_im, integral) &
      !$OMP  DEFAULT(NONE)  &
      !$OMP  SHARED(size_buffer, ao_num, num_kpts, df_num, ao_num_per_kpt, &
      !$OMP  kconserv, df_integral_array, ao_integrals_threshold)
  
  allocate( &
    ints_ik(ao_num_per_kpt**2,df_num),&
    ints_jl(ao_num_per_kpt**2,df_num),&
    ints_ikjl_tmp(ao_num_per_kpt**2,ao_num_per_kpt**2),&
    ints_ikjl(ao_num_per_kpt,ao_num_per_kpt,ao_num_per_kpt,ao_num_per_kpt),&
    buffer_i1(size_buffer),                                         &
    buffer_i2(size_buffer),                                         &
    buffer_value1(size_buffer),                                     &
    buffer_value2(size_buffer) & 
  )

  n_integrals=0

  !$OMP DO SCHEDULE(guided)
  do kl=1, num_kpts
    do kj=1, kl
      call idx2_tri_int(kj,kl,kjkl2)
      do kk=1,kl
        ki=kconserv(kl,kk,kj)
        call idx2_tri_int(ki,kk,kikk2)
        if (kikk2 > kjkl2) cycle
        if ((kl == kj) .and. (ki > kk)) cycle
        ! maybe use pointers instead of reshaping?
        if (ki >= kk) then
          ints_ik = reshape(reshape(df_integral_array(:,:,:,kikk2),(/ao_num_per_kpt,ao_num_per_kpt,df_num/),order=(/1,2,3/)),&
                    (/ao_num_per_kpt**2,df_num/))
        else
          ints_ik = reshape( &
                   conjg(reshape(df_integral_array(:,:,:,kikk2),(/ao_num_per_kpt,ao_num_per_kpt,df_num/),order=(/2,1,3/))),&
                   (/ao_num_per_kpt**2,df_num/))
        endif

        ints_jl = reshape( &
                 conjg(reshape(df_integral_array(:,:,:,kjkl2),(/ao_num_per_kpt,ao_num_per_kpt,df_num/),order=(/2,1,3/))),&
                 (/ao_num_per_kpt**2,df_num/))


        ! todo: option 1: change bounds so that kl <= kj
        !       option 2: change df_integral array so that it is the conjugate transpose of what it is now (transpose first 2 dimensions)

        ! todo: maybe just use 'C' instead of 'T' rather than doing conjugation above? (some cases will still require conjg on ikjl array)
        !       figure this out in conjunction with deciding how to structure df integrals when first constructed
        call zgemm('N','T', ao_num_per_kpt**2, ao_num_per_kpt**2, df_num, &
               (1.d0,0.d0), ints_ik, size(ints_ik,1), &
               ints_jl, size(ints_jl,1), &
               (0.d0,0.d0), ints_ikjl_tmp, size(ints_ikjl_tmp,1))

        ! this is bad
        ! use a pointer instead?
        ints_ikjl = reshape(ints_ikjl_tmp,(/ao_num_per_kpt,ao_num_per_kpt,ao_num_per_kpt,ao_num_per_kpt/))

        do il=1,ao_num_per_kpt
          l=il+(kl-1)*ao_num_per_kpt
          do ij=1,ao_num_per_kpt
            j=ij+(kj-1)*ao_num_per_kpt
            if (j>l) exit
            call idx2_tri_int(j,l,jl2)
            do ik=1,ao_num_per_kpt
              k=ik+(kk-1)*ao_num_per_kpt
              if (k>l) exit
              do ii=1,ao_num_per_kpt
                i=ii+(ki-1)*ao_num_per_kpt
                call idx2_tri_int(i,k,ik2)
                if (ik2 > jl2) exit
                if ((j==l) .and. (i>k)) exit
                integral = ints_ikjl(ii,ik,ij,il)
                if (cdabs(integral) < ao_integrals_threshold) then
                  cycle
                endif
                n_integrals += 1
                tmp_re = real(integral)
                tmp_im = imag(integral)
                call mo_bielec_integrals_index(i,j,k,l,tmp_idx1)
                call mo_bielec_integrals_index(k,l,i,j,tmp_idx2)
                if (tmp_idx1.eq.tmp_idx2) then
                  buffer_i1(n_integrals) = tmp_idx1
                  buffer_i2(n_integrals) = tmp_idx1
                  buffer_value1(n_integrals) = tmp_re
                  buffer_value2(n_integrals) = 0.d0
                else if (tmp_idx1 .lt. tmp_idx2) then
                  buffer_i1(n_integrals) = tmp_idx1
                  buffer_i2(n_integrals) = tmp_idx2
                  buffer_value1(n_integrals) = tmp_re
                  buffer_value2(n_integrals) = tmp_im
                else
                  buffer_i1(n_integrals) = tmp_idx2
                  buffer_i2(n_integrals) = tmp_idx1
                  buffer_value1(n_integrals) = tmp_re
                  buffer_value2(n_integrals) = -tmp_im
                endif


                if (n_integrals == size_buffer) then
                  call insert_into_ao_integrals_map(n_integrals,buffer_i1,buffer_value1)
                  call insert_into_ao_integrals_map(n_integrals,buffer_i2,buffer_value2)
                  n_integrals = 0
                endif
              enddo !ii
            enddo !ik
          enddo !ij
        enddo !il
      enddo !kk
    enddo !kj
  enddo !kl
 
  !$OMP END DO NOWAIT

  deallocate( &
    ints_ik,&
    ints_jl,&
    ints_ikjl_tmp,&
    ints_ikjl&
    )
  if (n_integrals /= 0) then
    call insert_into_ao_integrals_map(n_integrals,buffer_i1,buffer_value1)
    call insert_into_ao_integrals_map(n_integrals,buffer_i2,buffer_value2)
    n_integrals=0
  endif
  
  deallocate( &
    buffer_i1,&
    buffer_i2,&
    buffer_value1,&
    buffer_value2&
    )
  !$OMP END PARALLEL

  call map_sort(ao_integrals_map)
  call map_unique(ao_integrals_map)

  call map_save_to_disk(trim(ezfio_filename)//'/work/ao_ints',ao_integrals_map)
  call ezfio_set_integrals_bielec_disk_access_ao_integrals('Read')
  

end subroutine ao_map_fill_from_df


BEGIN_PROVIDER [ double precision, ao_bielec_integral_schwartz,(ao_num,ao_num)  ]
  implicit none
  BEGIN_DOC
  !  Needed to compute Schwartz inequalities
  END_DOC
  
  integer                        :: i,k
  double precision               :: cpu_1,cpu_2, wall_1, wall_2
  complex*16                     :: ao_bielec_integral
  
!  ao_bielec_integral_schwartz(1,1) = real(ao_bielec_integral(1,1,1,1))
  !$OMP PARALLEL DO PRIVATE(i,k)                                     &
      !$OMP DEFAULT(NONE)                                            &
      !$OMP SHARED (ao_num,ao_bielec_integral_schwartz)              &
      !$OMP SCHEDULE(dynamic)
  do i=1,ao_num
    do k=1,i
      ao_bielec_integral_schwartz(i,k) = dsqrt(real(ao_bielec_integral(i,k,k,i)))
      ao_bielec_integral_schwartz(k,i) = ao_bielec_integral_schwartz(i,k)
    enddo
  enddo
  !$OMP END PARALLEL DO

END_PROVIDER



