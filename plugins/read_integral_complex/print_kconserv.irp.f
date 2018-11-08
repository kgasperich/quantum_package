program print_integrals

  PROVIDE ezfio_filename
  call ezfio_set_utils_disk_access_kconserv('Read')
  call run
end

subroutine run
  implicit none
  
  integer :: iunit
  integer :: getunitandopen

  integer ::i,j,k,l

  iunit = getunitandopen('kconserv','w')
  do i=1,num_kpts
    do j=1,num_kpts
      do k=1,num_kpts
        write(iunit,*) i,j,k,kconserv(i,j,k)
      enddo
    enddo
  enddo
  close(iunit)
  
end
