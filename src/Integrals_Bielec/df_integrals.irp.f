
BEGIN_PROVIDER [integer, ao_num_per_kpt ]
  implicit none
  ao_num_per_kpt = ao_num / num_kpts
END_PROVIDER

BEGIN_PROVIDER [integer, num_kpt_pairs ]
  implicit none
  num_kpt_pairs = num_kpts * (num_kpts + 1) / 2
END_PROVIDER

BEGIN_PROVIDER [complex*16, df_integral_array, (num_kpt_pairs, df_tot_num, ao_num_per_kpt, ao_num_per_kpt)]
  implicit none
  BEGIN_DOC
  !  for k-points i,j,k, give k-point l for which the total momentum is conserved over k_i* k_j* k_k k_l
  END_DOC
  
  if (read_df_integral_array) then
    print *, 'reading df_integral_array from disk'
    call read_df_integral_array_file('df_integral_array',df_integral_array,df_tot_num,ao_num_per_kpt,num_kpt_pairs)
    print *, 'read df_integral_array from disk'
  else
    print *, 'df_integral_array must be provided'
  endif

END_PROVIDER


 BEGIN_PROVIDER [ logical, read_df_integral_array ]
&BEGIN_PROVIDER [ logical, write_df_integral_array ]
   
   BEGIN_DOC
   ! One level of abstraction for disk_access_df_integral_array
   END_DOC
   implicit none
   
   if (disk_access_df_integral_array.EQ.'Read') then
     read_df_integral_array =  .True.
     write_df_integral_array = .False.
     
   else if  (disk_access_df_integral_array.EQ.'Write') then
     read_df_integral_array = .False.
     write_df_integral_array =  .True.
     
   else if (disk_access_df_integral_array.EQ.'None') then
     read_df_integral_array = .False.
     write_df_integral_array = .False.
     
   else
     print *, 'utils/disk_access_df_integral_array has a wrong type'
     stop 1
     
   endif
   
END_PROVIDER

subroutine write_df_integral_array_file(filename, A, n_df, n_ao, n_kpt_pairs)
  implicit none
  BEGIN_DOC
! Write the 1-electron integrals stored in A(m,n) into file 'filename'
  END_DOC
  character(len=*), intent(in)   :: filename
  integer, intent(in)            :: n_df, n_ao, n_kpt_pairs
  complex*16, intent(in)         :: A(n_kpt_pairs,n_df,n_ao,n_ao)
  double precision, allocatable  :: A_re(:,:,:), A_im(:,:,:)

  integer                        :: iunit
  integer, external              :: getUnitAndOpen
  character*(256)                :: f
  
  integer                        :: i,j,k,l

  allocate( A_re(n_kpt_pairs,n_df,n_ao,n_ao), A_im(n_kpt_pairs,n_df,n_ao,n_ao))

  do i=1,n_kpt_pairs
    do j=1,n_df
      do k=1,n_ao
        do l=1,n_ao
          A_re(i,j,k,l) = real(A(i,j,k,l))
          A_im(i,j,k,l) = imag(A(i,j,k,l))
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

subroutine read_df_integral_array_file(filename, A, n_df, n_ao, n_kpt_pairs)
  implicit none
  BEGIN_DOC
! Read the 1-electron integrals into in A(m,n) from file 'filename'
  END_DOC
  character(len=*), intent(in)   :: filename
  integer, intent(in)            :: n_df, n_ao, n_kpt_pairs
  complex*16, intent(out)        :: A(n_kpt_pairs,n_df,n_ao,n_ao)
  double precision, allocatable  :: A_re(:,:,:,:), A_im(:,:,:,:)

  integer                        :: iunit
  integer, external              :: getUnitAndOpen
  character*(256)                :: f
  
  integer                        :: i,j,k,l

  allocate( A_re(n_kpt_pairs,n_df,n_ao,n_ao), A_im(n_kpt_pairs,n_df,n_ao,n_ao))

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
          A(i,j,k,l) = cmplx(A_re(i,j,k,l),A_im(i,j,k,l))
        enddo
      enddo
    enddo
  enddo
  deallocate(A_re,A_im)

end


