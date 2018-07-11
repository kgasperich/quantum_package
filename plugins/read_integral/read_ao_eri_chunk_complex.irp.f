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
  double precision :: int_re,int_im

  integer   :: n_integrals 
  integer   :: n_chunk_integrals

  integer(key_kind) :: tmpkey1,tmpkey2

  integer(key_kind), allocatable   :: buffer_re(:), buffer_im(:)
  real(integral_kind), allocatable :: buffer_values_re(:), buffer_values_im(:)

  n_chunk_integrals = 1024
  allocate(buffer_i(n_chunk_integrals), buffer_values(n_chunk_integrals))

  iunit = getunitandopen('bielec_ao_complex','r')
  n_integrals=0
  do 
    read (iunit,*,end=13) i,j,k,l, int_re, int_im
    n_integrals += 1
    call bielec_integrals_index_2fold(i, j, k, l, tmpkey1)
    call bielec_integrals_index_2fold(k, l, i, j, tmpkey2)
    buffer_re(n_integrals) = min(tmpkey1,tmpkey2)
    buffer_im(n_integrals) = max(tmpkey1,tmpkey2)
    if ( tmpkey1 > tmpkey2 ) then
      int_im = -int_im
    endif
    
    buffer_values_re(n_integrals) = int_re
    buffer_values_im(n_integrals) = int_im

    IF ( MOD(n_integrals,n_chunk_integrals).EQ.0) THEN
          call insert_into_ao_integrals_map_2fold(n_integrals,buffer_re,buffer_values_re)
          call insert_into_ao_integrals_map_2fold(n_integrals,buffer_im,buffer_values_im)
          n_integrals = 0; buffer_re(:) = 0; buffer_values_re(:) = 0.d0
          buffer_im(:) = 0; buffer_values_im(:) = 0.d0
    ENDIF
  enddo
  13 continue
  close(iunit)
  
  !Tail loop
  IF (n_integrals.ne.0) THEN
    call insert_into_ao_integrals_map_2fold(n_integrals,buffer_re,buffer_values_re)
    call insert_into_ao_integrals_map_2fold(n_integrals,buffer_im,buffer_values_im)
  END IF

  call map_sort(ao_integrals_map_2fold)
  call map_unique(ao_integrals_map_2fold)

  call map_save_to_disk(trim(ezfio_filename)//'/work/ao_ints_2fold',ao_integrals_map_2fold)
  call ezfio_set_integrals_bielec_disk_access_ao_integrals('Read')

end
