program read_integrals

  PROVIDE ezfio_filename
  call ezfio_set_integrals_bielec_disk_access_df_integral_array("None")
  call run
end

subroutine run
  use map_module
  implicit none
  
  integer :: iunit
  integer :: getunitandopen

  integer :: kikj,mu,i,j
  double precision :: int_re, int_im
  complex*16, allocatable :: A(:,:,:,:)

  allocate (A(ao_num_per_kpt,ao_num_per_kpt,df_num,num_kpt_pairs))


  A = 0.d0
  iunit = getunitandopen('df_integral_array','r')
  do 
    read (iunit,*,end=10) i,j,mu,kikj, int_re, int_im
    A(i,j,mu,kikj) = dcmplx(int_re,int_im)
  enddo
  10 continue
  close(iunit)
  call write_df_integral_array_file('df_integral_array', A, size(A,1), size(A,3), size(A,4))

  deallocate(A)

  call ezfio_set_integrals_bielec_disk_access_df_integral_array("Read")

end
