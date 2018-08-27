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
  allocate( &
    buffer_re(n_chunk_integrals), &
    buffer_im(n_chunk_integrals), &
    buffer_values_re(n_chunk_integrals), &
    buffer_values_im(n_chunk_integrals))

  iunit = getunitandopen('bielec_mo_complex','r')
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
      call map_append(mo_integrals_map, buffer_re, buffer_values_re, n_integrals)
      call map_append(mo_integrals_map, buffer_im, buffer_values_im, n_integrals)
      n_integrals = 0; buffer_re(:) = 0; buffer_values_re(:) = 0.d0
      buffer_im(:) = 0; buffer_values_im(:) = 0.d0
    ENDIF
  enddo
  13 continue
  close(iunit)
  
  !Tail loop
  IF (n_integrals.ne.0) THEN
    call map_append(mo_integrals_map, buffer_re, buffer_values_re, n_integrals)
    call map_append(mo_integrals_map, buffer_im, buffer_values_im, n_integrals)
  END IF

  call map_sort(mo_integrals_map)
  call map_unique(mo_integrals_map)

  call map_save_to_disk(trim(ezfio_filename)//'/work/mo_ints',mo_integrals_map)
  call ezfio_set_integrals_bielec_disk_access_mo_integrals('Read')

end
