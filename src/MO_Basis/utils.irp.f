subroutine save_mos
  implicit none
  double precision, allocatable  :: buffer_re(:,:), buffer_im(:,:)
  integer                        :: i,j
  
  call system('$QP_ROOT/scripts/save_current_mos.sh '//trim(ezfio_filename))
  
  call ezfio_set_mo_basis_mo_tot_num(mo_tot_num)
  call ezfio_set_mo_basis_mo_label(mo_label)
  call ezfio_set_mo_basis_ao_md5(ao_md5)
  allocate ( buffer_re(ao_num,mo_tot_num), buffer_im(ao_num,mo_tot_num) )
  buffer_re = 0.d0
  buffer_im = 0.d0
  do j = 1, mo_tot_num
    do i = 1, ao_num
      buffer_re(i,j) = real(mo_coef(i,j))
      buffer_im(i,j) = imag(mo_coef(i,j))
    enddo
  enddo
  call ezfio_set_mo_basis_mo_coef_real(buffer_re)
  call ezfio_set_mo_basis_mo_coef_imag(buffer_im)
  call ezfio_set_mo_basis_mo_occ(mo_occ)
  deallocate (buffer_re, buffer_im)
  
end

subroutine save_mos_truncated(n)
  implicit none
  double precision, allocatable  :: buffer_re(:,:), buffer_im(:,:)
  integer                        :: i,j,n
  
  call system('$QP_ROOT/scripts/save_current_mos.sh '//trim(ezfio_filename))
  
  call ezfio_set_mo_basis_mo_tot_num(n)
  call ezfio_set_mo_basis_mo_label(mo_label)
  call ezfio_set_mo_basis_ao_md5(ao_md5)
  allocate ( buffer_re(ao_num,n), buffer_im(ao_num,n) )
  buffer = 0.d0
  do j = 1, n
    do i = 1, ao_num
      buffer_re(i,j) = real(mo_coef(i,j))
      buffer_im(i,j) = imag(mo_coef(i,j))
    enddo
  enddo
  call ezfio_set_mo_basis_mo_coef_real(buffer_re)
  call ezfio_set_mo_basis_mo_coef_imag(buffer_im)
  call ezfio_set_mo_basis_mo_occ(mo_occ)
  deallocate (buffer_re, buffer_im)
  
end

subroutine mo_as_eigvectors_of_mo_matrix(matrix,n,m,label,sign,output)
  implicit none
  integer,intent(in)             :: n,m, sign
  character*(64), intent(in)     :: label
  complex*16, intent(in)   :: matrix(n,m)
  logical, intent(in)            :: output
  
  integer :: i,j
  double precision, allocatable  :: eigvalues(:)
  complex*16, allocatable  :: mo_coef_new(:,:), R(:,:), A(:,:)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: mo_coef_new, R
  
  call write_time(6)
  if (m /= mo_tot_num) then
    print *, irp_here, ': Error : m/= mo_tot_num'
    stop 1
  endif
  allocate(A(n,m),R(n,m),mo_coef_new(ao_num,m),eigvalues(m))
  if (sign == -1) then
    do j=1,m
      do i=1,n
        A(i,j) = -matrix(i,j)
      enddo
    enddo
  else
    do j=1,m
      do i=1,n
        A(i,j) = matrix(i,j)
      enddo
    enddo
  endif
  mo_coef_new = mo_coef
  
  call lapack_diag_z(eigvalues,R,A,n,m)
  if (output) then
    write (6,'(A)')  'MOs are now **'//trim(label)//'**'
    write (6,'(A)') ''
    write (6,'(A)')  'Eigenvalues'
    write (6,'(A)') '-----------'
    write (6,'(A)')  ''
    write (6,'(A)') '======== ================'
  endif
  if (sign == -1) then
    do i=1,m
      eigvalues(i) = -eigvalues(i)
    enddo
  endif
  if (output) then
    do i=1,m
      write (6,'(I8,1X,F16.10)')  i,eigvalues(i)
    enddo
    write (6,'(A)') '======== ================'
    write (6,'(A)')  ''
  endif
  
  call zgemm('N','N',ao_num,m,m,(1.d0,0.d0),mo_coef_new,size(mo_coef_new,1),&
             R,size(R,1),(0.d0,0.d0),mo_coef,size(mo_coef,1))
  deallocate(A,mo_coef_new,R,eigvalues)
  call write_time(6)
  
  mo_label = label
end

subroutine mo_as_svd_vectors_of_mo_matrix(matrix,lda,m,n,label)
  implicit none
  integer,intent(in)             :: lda,m,n
  character*(64), intent(in)     :: label
  complex*16, intent(in)   :: matrix(lda,n)
  
  integer :: i,j
  double precision, allocatable  :: D(:)
  complex*16, allocatable  :: mo_coef_new(:,:), U(:,:), A(:,:), Vt(:,:)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: mo_coef_new, U, Vt, A
  
  call write_time(6)
  if (m /= mo_tot_num) then
    print *, irp_here, ': Error : m/= mo_tot_num'
    stop 1
  endif

  allocate(A(lda,n),U(lda,n),mo_coef_new(ao_num,m),D(m),Vt(lda,n))

  do j=1,n
    do i=1,m
      A(i,j) = matrix(i,j)
    enddo
  enddo
  mo_coef_new = mo_coef
  
  call svd_z(A,lda,U,lda,D,Vt,lda,m,n)

  write (6,'(A)') 'MOs are now **'//trim(label)//'**'
  write (6,'(A)')  ''
  write (6,'(A)') 'Eigenvalues'
  write (6,'(A)')  '-----------'
  write (6,'(A)') ''
  write (6,'(A)')  '======== ================'

  do i=1,m
    write (6,'(I8,1X,F16.10)')  i,D(i)
  enddo
  write (6,'(A)')  '======== ================'
  write (6,'(A)')  ''
  
  call zgemm('N','N',ao_num,m,m,(1.d0,0.d0),mo_coef_new,size(mo_coef_new,1),&
             U,size(U,1),(0.d0,0.d0),mo_coef,size(mo_coef,1))
  deallocate(A,mo_coef_new,U,Vt,D)
  call write_time(6)
  
  mo_label = label
end



subroutine give_all_mos_at_r(r,mos_array)
 implicit none
 double precision, intent(in) :: r(3)
 complex*16, intent(out) :: mos_array(mo_tot_num)
 call give_specific_mos_at_r(r,mos_array, mo_coef)
end

subroutine give_specific_mos_at_r(r,mos_array, mo_coef_specific)
 implicit none
 double precision, intent(in)  :: r(3)
 complex*16, intent(in)        :: mo_coef_specific(ao_num, mo_tot_num)
 complex*16, intent(out)       :: mos_array(mo_tot_num)
 double precision              :: aos_array(ao_num)
 complex*16                    :: accu
 integer :: i,j
 call give_all_aos_at_r(r,aos_array)
 do i = 1, mo_tot_num
  accu = (0.d0,0.d0)
  do j = 1, ao_num
   accu += mo_coef_specific(j,i) * aos_array(j) 
  enddo
  mos_array(i) = accu
 enddo
end
