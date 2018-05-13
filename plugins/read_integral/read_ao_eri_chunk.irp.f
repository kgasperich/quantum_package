program read_integrals

  PROVIDE ezfio_filename
  call run
end

subroutine run
  use map_module
  implicit none
  
  integer :: iunit
  integer :: getunitandopen

  integer ::i,j,k,l
  double precision :: integral

  integer   :: n_integrals 
  integer   :: n_chunk_integrals

  integer(key_kind), allocatable   :: buffer_i(:) 
  real(integral_kind), allocatable :: buffer_values(:)

  n_chunk_integrals = 1024
  allocate(buffer_i(n_chunk_integrals), buffer_values(n_chunk_integrals))

  iunit = getunitandopen('bielec_ao','r')
  n_integrals=0
  do 
    read (iunit,*,end=13) i,j,k,l, integral
    n_integrals += 1
    call bielec_integrals_index(i, j, k, l, buffer_i(n_integrals) )
    buffer_values(n_integrals) = integral

    IF ( MOD(n_integrals,n_chunk_integrals).EQ.0) THEN
          call insert_into_ao_integrals_map(n_integrals,buffer_i,buffer_values)
          n_integrals = 0; buffer_i(:) = 0; buffer_values(:) = 0.d0
    ENDIF
  enddo
  13 continue
  close(iunit)
  
  !Tail loop
  IF (n_integrals.ne.0) THEN
    call insert_into_ao_integrals_map(n_integrals,buffer_i,buffer_values)
  END IF

  call map_sort(ao_integrals_map)
  call map_unique(ao_integrals_map)

  call map_save_to_disk(trim(ezfio_filename)//'/work/ao_ints',ao_integrals_map)
  call ezfio_set_integrals_bielec_disk_access_ao_integrals('Read')

end
