program read_integrals

  PROVIDE ezfio_filename
  call ezfio_set_utils_disk_access_kconserv("None")
  call run
end

subroutine run
  use map_module
  implicit none
  BEGIN_DOC
  ! read kconserv in physicist notation order <ij|kl>
  ! if kconserv(i,j,k)=l, then <ij|kl> is allowed by symmetry
  ! pyscf stores this internally in the order of chemist notation (ik|jl)
  END_DOC
  integer :: iunit
  integer :: getunitandopen

  integer ::i,j,k,l
  integer, allocatable :: A(:,:,:)

  call ezfio_get_utils_num_kpts(num_kpts)

  allocate (A(num_kpts,num_kpts,num_kpts))


  A = 0
  iunit = getunitandopen('kconserv_complex','r')
  do 
    read (iunit,*,end=10) i,j,k,l
    A(i,j,k) = l
  enddo
  10 continue
  close(iunit)
  call write_kconserv_file('kconserv', A, size(A,1))

  deallocate(A)

  call ezfio_set_utils_disk_access_kconserv("Read")

end
