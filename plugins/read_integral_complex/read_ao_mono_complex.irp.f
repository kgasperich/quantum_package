program read_integrals

  PROVIDE ezfio_filename
  call ezfio_set_integrals_monoelec_disk_access_ao_one_integrals("None")
  call ezfio_set_utils_disk_access_ao_overlap_integrals("None")
  call run
end

subroutine run
  use map_module
  implicit none
  
  integer :: iunit
  integer :: getunitandopen

  integer ::i,j
  double precision :: int_re, int_im
  complex*16, allocatable :: A(:,:)

  allocate (A(ao_num,ao_num))


  A = 0.d0
  iunit = getunitandopen('kinetic_ao_complex','r')
  do 
    read (iunit,*,end=10) i,j, int_re, int_im
    if (i.eq.j) then
      int_im = 0.d0
      A(i,i) = dcmplx(int_re,int_im)
    else
      A(i,j) = dcmplx(int_re,int_im)
      A(j,i) = dcmplx(int_re,-int_im)
    endif
  enddo
  10 continue
  close(iunit)
  call write_one_e_integrals_complex('ao_kinetic_integral', A, size(A,1), size(A,2))


!  A = 0.d0
!  print*,'opening ne_ao_complex'
!  iunit = getunitandopen('ne_ao_complex','r')
!  print*,'getting integrals from ne_ao_complex'
!  do
!    read (iunit,*,end=11) i,j, int_re, int_im
!    if (i.eq.j) then
!      int_im = 0.d0
!      A(i,i) = dcmplx(int_re,int_im)
!    else
!      A(i,j) = dcmplx(int_re,int_im)
!      A(j,i) = dcmplx(int_re,-int_im)
!    endif
!  enddo
!  11 continue
!  print*,'closing ne_ao_complex'
!  close(iunit)
!  print*,'saving ne_ao_integrals to work dir'
!  call write_one_e_integrals_complex('ao_ne_integral', A, size(A,1), size(A,2))
!  print*,'ne_ao integrals saved'


  A = 0.d0
  iunit = getunitandopen('overlap_ao_complex','r')
  do 
    read (iunit,*,end=12) i,j, int_re, int_im
    if (i.eq.j) then
      int_im = 0.d0
      A(i,i) = dcmplx(int_re,int_im)
    else
      A(i,j) = dcmplx(int_re,int_im)
      A(j,i) = dcmplx(int_re,-int_im)
    endif
  enddo
  12 continue
  close(iunit)
  call write_one_e_integrals_complex('ao_overlap', A, size(A,1), size(A,2))


  A = 0.d0
  print*,'opening ne_ao_complex'
  iunit = getunitandopen('ne_ao_complex','r')
  print*,'getting integrals from ne_ao_complex'
  do
    read (iunit,*,end=11) i,j, int_re, int_im
    if (i.eq.j) then
      int_im = 0.d0
      A(i,i) = dcmplx(int_re,int_im)
    else
      A(i,j) = dcmplx(int_re,int_im)
      A(j,i) = dcmplx(int_re,-int_im)
    endif
  enddo
  11 continue
  print*,'closing ne_ao_complex'
  close(iunit)
  print*,'saving ne_ao_integrals to work dir'
  call write_one_e_integrals_complex('ao_ne_integral', A, size(A,1), size(A,2))
  print*,'ne_ao integrals saved'

  deallocate(A)

  call ezfio_set_integrals_monoelec_disk_access_ao_one_integrals("Read")
  call ezfio_set_utils_disk_access_ao_overlap_integrals("Read")

end
