program read_integrals

  PROVIDE ezfio_filename
  call ezfio_set_integrals_monoelec_disk_access_ao_one_integrals("None")
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
    A(i,j) = dcmplx(int_re,int_im)
    A(j,i) = dcmplx(int_re,-int_im)
  enddo
  10 continue
  close(iunit)
  call write_one_e_integrals_complex('ao_kinetic_integral', A, size(A,1), size(A,2))


  A = 0.d0
  iunit = getunitandopen('nuclear_ao_complex','r')
  do 
    read (iunit,*,end=12) i,j, int_re, int_im
    A(i,j) = dcmplx(int_re,int_im)
    A(j,i) = dcmplx(int_re,-int_im)
  enddo
  12 continue
  close(iunit)
  call write_one_e_integrals_complex('ao_ne_integral', A, size(A,1), size(A,2))

  deallocate(A)

  call write_one_e_integrals_complex('ao_pseudo_integral', ao_pseudo_integral,&
        size(ao_pseudo_integral,1), size(ao_pseudo_integral,2))

  call ezfio_set_integrals_monoelec_disk_access_ao_one_integrals("Read")
end
