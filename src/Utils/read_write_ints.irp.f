 BEGIN_PROVIDER [ logical, read_ao_overlap_integrals ]
&BEGIN_PROVIDER [ logical, write_ao_overlap_integrals ]
   
   BEGIN_DOC
   ! One level of abstraction for disk_access_ao_integrals and disk_access_mo_integrals
   END_DOC
   implicit none
   
   if (disk_access_ao_overlap_integrals.EQ.'Read') then
     read_ao_overlap_integrals =  .True.
     write_ao_overlap_integrals = .False.
     
   else if  (disk_access_ao_overlap_integrals.EQ.'Write') then
     read_ao_overlap_integrals = .False.
     write_ao_overlap_integrals =  .True.
     
   else if (disk_access_ao_overlap_integrals.EQ.'None') then
     read_ao_overlap_integrals = .False.
     write_ao_overlap_integrals = .False.
     
   else
     print *, 'utils/disk_access_ao_overlap_integrals has a wrong type'
     stop 1
     
   endif
   
END_PROVIDER

subroutine write_one_e_integrals(filename, A, m, n)
  implicit none
  BEGIN_DOC
! Write the 1-electron integrals stored in A(m,n) into file 'filename'
  END_DOC
  character(len=*), intent(in)   :: filename
  integer, intent(in)            :: m,n
  double precision, intent(in)   :: A(m,n)

  integer                        :: iunit
  integer, external              :: getUnitAndOpen
  character*(256)                :: f

  iunit = getUnitAndOpen( trim(ezfio_work_dir)//trim(filename), 'W' )
  write(iunit) A
  close(iunit)
end

subroutine read_one_e_integrals(filename, A, m, n)
  implicit none
  BEGIN_DOC
! Read the 1-electron integrals into in A(m,n) from file 'filename'
  END_DOC
  character(len=*), intent(in)   :: filename
  integer, intent(in)            :: m,n
  double precision, intent(out)  :: A(m,n)

  integer                        :: iunit
  integer, external              :: getUnitAndOpen
  character*(256)                :: f

  iunit = getUnitAndOpen( trim(ezfio_work_dir)//trim(filename), 'R' )
  read(iunit) A
  close(iunit)
end

subroutine write_one_e_integrals_complex(filename, A, m, n)
  implicit none
  BEGIN_DOC
! Write the 1-electron integrals stored in A(m,n) into file 'filename'
  END_DOC
  character(len=*), intent(in)   :: filename
  integer, intent(in)            :: m,n
  complex*16, intent(in)   :: A(m,n)
  double precision, allocatable :: A_re(:,:), A_im(:,:)

  integer                        :: i,j
  integer                        :: iunit
  integer, external              :: getUnitAndOpen
  character*(256)                :: f

  allocate( A_re(m,n), A_im(m,n) )

  do i=1,m
    do j=1,n
      A_re(i,j) = real(A(i,j))
      A_im(i,j) = imag(A(i,j))
    enddo
  enddo

  iunit = getUnitAndOpen( trim(ezfio_work_dir)//trim(filename)//'_real', 'W' )
  write(iunit) A_re
  close(iunit)

  iunit = getUnitAndOpen( trim(ezfio_work_dir)//trim(filename)//'_imag', 'W' )
  write(iunit) A_im
  close(iunit)
  deallocate( A_re, A_im )
end

subroutine read_one_e_integrals_complex(filename, A, m, n)
  implicit none
  BEGIN_DOC
! Read the 1-electron integrals into in A(m,n) from file 'filename'
  END_DOC
  character(len=*), intent(in)   :: filename
  integer, intent(in)            :: m,n
  complex*16, intent(out)  :: A(m,n)
  double precision, allocatable  :: A_re(:,:), A_im(:,:)

  integer                        :: i,j
  integer                        :: iunit
  integer, external              :: getUnitAndOpen
  character*(256)                :: f

  allocate( A_re(m,n), A_im(m,n) )

  iunit = getUnitAndOpen( trim(ezfio_work_dir)//trim(filename)//'_real', 'R' )
  read(iunit) A_re
  close(iunit)
  
  iunit = getUnitAndOpen( trim(ezfio_work_dir)//trim(filename)//'_imag', 'R' )
  read(iunit) A_im
  close(iunit)

  do i=1,m
    do j=1,n
      A(i,j) = cmplx(A_re(i,j),A_im(i,j))
    enddo
  enddo
  deallocate( A_re, A_im )

end

subroutine write_one_e_integrals_complex_kpts(filename, A, m, n, p)
  implicit none
  BEGIN_DOC
! Write the 1-electron integrals stored in A(m,n) into file 'filename'
  END_DOC
  character(len=*), intent(in)   :: filename
  integer, intent(in)            :: m,n,p
  complex*16, intent(in)   :: A(m,n,p)
  double precision, allocatable :: A_re(:,:,:), A_im(:,:,:)

  integer                        :: i,j,k
  integer                        :: iunit
  integer, external              :: getUnitAndOpen
  character*(256)                :: f

  allocate( A_re(m,n,p), A_im(m,n,p) )

  do k=1,p
    do j=1,n
      do i=1,m
        A_re(i,j,k) = real(A(i,j,k))
        A_im(i,j,k) = imag(A(i,j,k))
      enddo
    enddo
  enddo

  iunit = getUnitAndOpen( trim(ezfio_work_dir)//trim(filename)//'_real', 'W' )
  write(iunit) A_re
  close(iunit)

  iunit = getUnitAndOpen( trim(ezfio_work_dir)//trim(filename)//'_imag', 'W' )
  write(iunit) A_im
  close(iunit)
  deallocate( A_re, A_im )
end

subroutine read_one_e_integrals_complex_kpts(filename, A, m, n, p)
  implicit none
  BEGIN_DOC
! Read the 1-electron integrals into in A(m,n) from file 'filename'
  END_DOC
  character(len=*), intent(in)   :: filename
  integer, intent(in)            :: m,n,p
  complex*16, intent(out)  :: A(m,n,p)
  double precision, allocatable  :: A_re(:,:,:), A_im(:,:,:)

  integer                        :: i,j,k
  integer                        :: iunit
  integer, external              :: getUnitAndOpen
  character*(256)                :: f

  allocate( A_re(m,n,p), A_im(m,n,p) )

  iunit = getUnitAndOpen( trim(ezfio_work_dir)//trim(filename)//'_real', 'R' )
  read(iunit) A_re
  close(iunit)
  
  iunit = getUnitAndOpen( trim(ezfio_work_dir)//trim(filename)//'_imag', 'R' )
  read(iunit) A_im
  close(iunit)

  do k=1,p
    do j=1,n
      do i=1,n
        A(i,j,k) = cmplx(A_re(i,j,k),A_im(i,j,k))
      enddo
    enddo
  enddo
  deallocate( A_re, A_im )

end

