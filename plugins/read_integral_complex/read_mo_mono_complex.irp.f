program read_integrals

  PROVIDE ezfio_filename
  call ezfio_set_integrals_monoelec_disk_access_mo_one_integrals("None")
!  call ezfio_set_utils_disk_access_ao_overlap_integrals("None")
  call run_mo_mono_complex
end

subroutine run_mo_mono_complex
  use map_module
  implicit none
  
  integer :: iunit
  integer :: getunitandopen

  integer ::i,j
  double precision :: int_re, int_im
  complex*16, allocatable :: A(:,:)

  call ezfio_get_mo_basis_mo_tot_num(mo_tot_num)

  allocate (A(mo_tot_num,mo_tot_num))


  A = 0.d0
  iunit = getunitandopen('kinetic_mo_complex','r')
  do 
    read (iunit,*,end=10) i,j, int_re, int_im
    A(i,j) = dcmplx(int_re,int_im)
    A(j,i) = dcmplx(int_re,-int_im)
  enddo
  10 continue
  close(iunit)
  call write_one_e_integrals_complex('mo_kinetic_integral', A, size(A,1), size(A,2))


  A = 0.d0
  iunit = getunitandopen('ne_mo_complex','r')
  do
    read (iunit,*,end=11) i,j, int_re, int_im
    A(i,j) = dcmplx(int_re,int_im)
    A(j,i) = dcmplx(int_re,-int_im)
  enddo
  11 continue
  close(iunit)
  call write_one_e_integrals_complex('mo_ne_integral', A, size(A,1), size(A,2))


  A = 0.d0
  iunit = getunitandopen('overlap_mo_complex','r')
  do 
    read (iunit,*,end=12) i,j, int_re, int_im
    A(i,j) = dcmplx(int_re,int_im)
    A(j,i) = dcmplx(int_re,-int_im)
  enddo
  12 continue
  close(iunit)
  call write_one_e_integrals_complex('mo_overlap', A, size(A,1), size(A,2))


  deallocate(A)

  call ezfio_set_integrals_monoelec_disk_access_mo_one_integrals("Read")
!  call ezfio_set_utils_disk_access_ao_overlap_integrals("Read")

end
