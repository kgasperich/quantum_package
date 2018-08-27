 BEGIN_PROVIDER [ logical, read_kconserv ]
&BEGIN_PROVIDER [ logical, write_kconserv ]
   
   BEGIN_DOC
   ! One level of abstraction for disk_access_kconserv
   END_DOC
   implicit none
   
   if (disk_access_kconserv.EQ.'Read') then
     read_kconserv =  .True.
     write_kconserv = .False.
     
   else if  (disk_access_kconserv.EQ.'Write') then
     read_kconserv = .False.
     write_kconserv =  .True.
     
   else if (disk_access_kconserv.EQ.'None') then
     read_kconserv = .False.
     write_kconserv = .False.
     
   else
     print *, 'utils/disk_access_kconserv has a wrong type'
     stop 1
     
   endif
   
END_PROVIDER

BEGIN_PROVIDER [integer, kconserv, (num_kpts,num_kpts,num_kpts)]
  implicit none
  BEGIN_DOC
  !  for k-points i,j,k, give k-point l for which the total momentum is conserved over k_i* k_j* k_k k_l
  END_DOC
  
  if (read_kconserv) then
    print *, 'reading kconserv from disk'
    call read_kconserv_file('kconserv',kconserv,num_kpts)
    print *, 'read kconserv from disk'
  else
    print *, 'kconserv must be provided'
  endif

END_PROVIDER

subroutine write_kconserv_file(filename, A, n_k)
  implicit none
  BEGIN_DOC
! Write the 1-electron integrals stored in A(m,n) into file 'filename'
  END_DOC
  character(len=*), intent(in)   :: filename
  integer, intent(in)            :: n_k
  integer, intent(in)            :: A(n_k,n_k,n_k)

  integer                        :: iunit
  integer, external              :: getUnitAndOpen
  character*(256)                :: f

  iunit = getUnitAndOpen( trim(ezfio_work_dir)//trim(filename), 'W' )
  write(iunit) A
  close(iunit)
end

subroutine read_kconserv_file(filename, A, n_k)
  implicit none
  BEGIN_DOC
! Read the 1-electron integrals into in A(m,n) from file 'filename'
  END_DOC
  character(len=*), intent(in)   :: filename
  integer, intent(in)            :: n_k
  integer, intent(out)           :: A(n_k,n_k,n_k)

  integer                        :: iunit
  integer, external              :: getUnitAndOpen
  character*(256)                :: f

  iunit = getUnitAndOpen( trim(ezfio_work_dir)//trim(filename), 'R' )
  read(iunit) A
  close(iunit)
end


