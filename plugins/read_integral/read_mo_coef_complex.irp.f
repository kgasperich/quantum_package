program read_integrals

  PROVIDE ezfio_filename
!  call ezfio_set_integrals_monoelec_disk_access_ao_one_integrals("None")
  call run
end

subroutine run
  use map_module
  implicit none
  
  integer :: iunit
  integer :: getunitandopen

  integer ::i,j
  double precision :: int_re, int_im

  
  iunit = getunitandopen('mo_coef_complex','r')
  provide mo_coef
  do 
    read (iunit,*,end=10) i,j, int_re, int_im
    mo_coef(i,j) = dcmplx(int_re,int_im)
  enddo
  10 continue
  close(iunit)
  call save_mos

end
