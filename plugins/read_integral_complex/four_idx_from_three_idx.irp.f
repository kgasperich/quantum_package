program print_integrals

  PROVIDE ezfio_filename
  call ezfio_set_integrals_bielec_disk_access_ao_integrals('Write')
  call run_34idx
end

subroutine run_34idx
  implicit none
  
  integer :: iunit
  integer :: getunitandopen

  integer ::i,j,k,l
  complex*16 :: integral
!  call ao_map_fill_from_df
  print*,'transforming 3-idx to 4-idx (AOs)'
  PROVIDE ao_bielec_integrals_in_map
  print*,'transformation complete'
end
