BEGIN_PROVIDER [complex*16, df_ao_integral_array, (ao_num_per_kpt, ao_num_per_kpt, df_num, num_kpt_pairs)]
  implicit none
  BEGIN_DOC
  !
  END_DOC
  
  if (read_df_ao_integral_array) then
    print *, 'reading df_ao_integral_array from disk'
    call read_df_integral_array_file('df_ao_integral_array',df_ao_integral_array,ao_num_per_kpt,df_num,num_kpt_pairs)
    print *, 'read df_ao_integral_array from disk'
  else
    print *, 'df_ao_integral_array must be provided'
  endif

END_PROVIDER

BEGIN_PROVIDER [complex*16, df_mo_integral_array, (mo_num_per_kpt, mo_num_per_kpt, df_num, num_kpt_pairs)]
  implicit none
  BEGIN_DOC
  !
  END_DOC
 
  if (read_df_mo_integral_array) then
    print *, 'reading df_mo_integral_array from disk'
    call read_df_integral_array_file('df_mo_integral_array',df_mo_integral_array,mo_num_per_kpt,df_num,num_kpt_pairs)
    print *, 'read df_mo_integral_array from disk'
  else
    call df_mo_from_df_ao(df_mo_integral_array,mo_num_per_kpt,df_num,num_kpt_pairs)
    if (write_df_mo_integral_array) then
      print*, 'saving df mo integrals to disk'
      call write_df_integral_array_file('df_mo_integral_array',df_mo_integral_array,mo_num_per_kpt,df_num,num_kpt_pairs)
      print*, 'setting integrals_bielec/disk_access_df_mo_integral_array = Read'
      call ezfio_set_integrals_bielec_disk_access_df_mo_integral_array("Read")
      disk_access_df_mo_integral_array = 'Read'
      SOFT_TOUCH disk_access_df_mo_integral_array
!      print*, 'read/write = ',read_df_mo_integral_array,write_df_mo_integral_array
!      print*, 'touching read/write_df_mo_integral_array'
!      TOUCH read_df_mo_integral_array write_df_mo_integral_array
!      print*, 'read/write = ',read_df_mo_integral_array,write_df_mo_integral_array
!      print*, 'df mo integrals saved to disk'
    endif
  endif

END_PROVIDER


 BEGIN_PROVIDER [ logical, read_df_ao_integral_array ]
&BEGIN_PROVIDER [ logical, write_df_ao_integral_array ]
   
   BEGIN_DOC
   ! One level of abstraction for disk_access_df_ao_integral_array
   END_DOC
   implicit none
   
   if (disk_access_df_ao_integral_array.EQ.'Read') then
     read_df_ao_integral_array =  .True.
     write_df_ao_integral_array = .False.
     
   else if  (disk_access_df_ao_integral_array.EQ.'Write') then
     read_df_ao_integral_array = .False.
     write_df_ao_integral_array =  .True.
     
   else if (disk_access_df_ao_integral_array.EQ.'None') then
     read_df_ao_integral_array = .False.
     write_df_ao_integral_array = .False.
     
   else
     print *, 'integrals_bielec/disk_access_df_ao_integral_array has a wrong type'
     stop 1
     
   endif
   
END_PROVIDER

 BEGIN_PROVIDER [ logical, read_df_mo_integral_array ]
&BEGIN_PROVIDER [ logical, write_df_mo_integral_array ]
   
   BEGIN_DOC
   ! One level of abstraction for disk_access_df_mo_integral_array
   END_DOC
   implicit none
   
   if (disk_access_df_mo_integral_array.EQ.'Read') then
     read_df_mo_integral_array =  .True.
     write_df_mo_integral_array = .False.
     
   else if  (disk_access_df_mo_integral_array.EQ.'Write') then
     read_df_mo_integral_array = .False.
     write_df_mo_integral_array =  .True.
     
   else if (disk_access_df_mo_integral_array.EQ.'None') then
     read_df_mo_integral_array = .False.
     write_df_mo_integral_array = .False.
     
   else
     print *, 'integrals_bielec/disk_access_df_mo_integral_array has a wrong type'
     stop 1
     
   endif
   
END_PROVIDER

subroutine write_df_integral_array_file(filename, A, n_ao, n_df, n_kpt_pairs)
  implicit none
  BEGIN_DOC
! Write the 1-electron integrals stored in A(m,n) into file 'filename'
  END_DOC
  character(len=*), intent(in)   :: filename
  integer, intent(in)            :: n_df, n_ao, n_kpt_pairs
  complex*16, intent(in)         :: A(n_ao,n_ao,n_df,n_kpt_pairs)
  double precision, allocatable  :: A_re(:,:,:,:), A_im(:,:,:,:)

  integer                        :: iunit
  integer, external              :: getUnitAndOpen
  character*(256)                :: f
  
  integer                        :: i,j,k,l

  allocate( A_re(n_ao,n_ao,n_df,n_kpt_pairs), A_im(n_ao,n_ao,n_df,n_kpt_pairs))

  do i=1,n_kpt_pairs
    do j=1,n_df
      do k=1,n_ao
        do l=1,n_ao
          A_re(l,k,j,i) = real(A(l,k,j,i))
          A_im(l,k,j,i) = imag(A(l,k,j,i))
        enddo
      enddo
    enddo
  enddo

  iunit = getUnitAndOpen( trim(ezfio_work_dir)//trim(filename)//'_real', 'W' )
  write(iunit) A_re
  close(iunit)

  iunit = getUnitAndOpen( trim(ezfio_work_dir)//trim(filename)//'_imag', 'W' )
  write(iunit) A_im
  close(iunit)

  deallocate(A_re,A_im)
end

subroutine read_df_integral_array_file(filename, A, n_ao, n_df, n_kpt_pairs)
  implicit none
  BEGIN_DOC
! Read the 1-electron integrals into in A(m,n) from file 'filename'
  END_DOC
  character(len=*), intent(in)   :: filename
  integer, intent(in)            :: n_df, n_ao, n_kpt_pairs
  complex*16, intent(out)        :: A(n_ao,n_ao,n_df,n_kpt_pairs)
  double precision, allocatable  :: A_re(:,:,:,:), A_im(:,:,:,:)

  integer                        :: iunit
  integer, external              :: getUnitAndOpen
  character*(256)                :: f
  
  integer                        :: i,j,k,l

  allocate( A_re(n_ao,n_ao,n_df,n_kpt_pairs), A_im(n_ao,n_ao,n_df,n_kpt_pairs))

  iunit = getUnitAndOpen( trim(ezfio_work_dir)//trim(filename)//'_real', 'R' )
  read(iunit) A_re
  close(iunit)

  iunit = getUnitAndOpen( trim(ezfio_work_dir)//trim(filename)//'_imag', 'R' )
  read(iunit) A_im
  close(iunit)

  do i=1,n_kpt_pairs
    do j=1,n_df
      do k=1,n_ao
        do l=1,n_ao
          A(l,k,j,i) = cmplx(A_re(l,k,j,i),A_im(l,k,j,i))
        enddo
      enddo
    enddo
  enddo
  deallocate(A_re,A_im)

end


subroutine df_mo_from_df_ao(df_mo,n_mo,n_df,n_k_pairs)
  use map_module
  implicit none
  BEGIN_DOC
  ! create 3-idx mo ints from 3-idx ao ints
  END_DOC
  integer,intent(in) :: n_mo,n_df,n_k_pairs
  complex*16,intent(out) :: df_mo(n_mo,n_mo,n_df,n_k_pairs)
  integer :: kl,kj,kjkl2,mu
  complex*16,allocatable :: coef_l(:,:), coef_j(:,:), ints_jl(:,:), ints_tmp(:,:)
  double precision :: wall_1,wall_2,cpu_1,cpu_2

  print*,'providing 3-index MO integrals from 3-index AO integrals'

  call wall_time(wall_1)
  call cpu_time(cpu_1)
  allocate( &
            coef_l(ao_num_per_kpt,mo_num_per_kpt),&
            coef_j(ao_num_per_kpt,mo_num_per_kpt),&
            ints_jl(ao_num_per_kpt,ao_num_per_kpt),&
            ints_tmp(mo_num_per_kpt,ao_num_per_kpt)&
          )

  do kl=1, num_kpts
    coef_l = mo_coef_kpts_trunc(:,:,kl)
    do kj=1, kl
      coef_j = mo_coef_kpts_trunc(:,:,kj)
      call idx2_tri_int(kj,kl,kjkl2)
      do mu=1, df_num
        ints_jl = df_ao_integral_array(:,:,mu,kjkl2)
        call zgemm('C','N',mo_num_per_kpt,ao_num_per_kpt,ao_num_per_kpt, &
              (1.d0,0.d0), coef_j, ao_num_per_kpt, &
              ints_jl, ao_num_per_kpt, &
              (0.d0,0.d0), ints_tmp, mo_num_per_kpt)

        call zgemm('N','N',mo_num_per_kpt,mo_num_per_kpt,ao_num_per_kpt, &
              (1.d0,0.d0), ints_tmp, mo_num_per_kpt, &
              coef_l, ao_num_per_kpt, &
              (0.d0,0.d0), df_mo(:,:,mu,kjkl2), mo_num_per_kpt)
      enddo
    enddo
    call wall_time(wall_2)
    print*,100.*float(kl*(kl+1))/(2.*n_k_pairs), '% in ', &
                wall_2-wall_1, 's'
  enddo

  deallocate( &
            coef_l, &
            coef_j, &
            ints_jl, &
            ints_tmp &
          )
  call wall_time(wall_2)
  call cpu_time(cpu_2)
  print*,' 3-idx MO provided'
  print*,'  cpu  time:',cpu_2-cpu_1,'s'
  print*,'  wall time:',wall_2-wall_1,'s  ( x ',(cpu_2-cpu_1)/(wall_2-wall_1),')'

end subroutine df_mo_from_df_ao
