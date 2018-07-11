program read_integrals

  PROVIDE ezfio_filename
  call ezfio_set_integrals_monoelec_disk_access_ao_one_integrals("None")
  call run
end

subroutine run
  use map_module
  implicit none
  
  integer :: iunit
  integer :: getunitandopen

  integer ::i,j,k,l
  double precision :: int_re, int_im
  complex*16, allocatable :: A(:,:)

  integer             :: n_integrals 

  integer(key_kind)                :: tmpkey1, tmpkey2
  integer(key_kind), allocatable   :: buffer_re(:), buffer_im(:)
  real(integral_kind), allocatable :: buffer_values_re(:), buffer_values_im(:)
  integer(key_kind)  :: key
   
  allocate (A(ao_num,ao_num))
  A = 0.d0
  
  iunit = getunitandopen('kinetic_ao_complex','r')
  do 
    read (iunit,*,end=10) i,j, int_re, int_im
    A(i,j) = dcmplx(int_re,int_im)
    A(j,i) = dcmplx(int_re,-int_im)
  enddo
  10 continue
  close(iunit)
  call write_one_e_integrals('ao_kinetic_integral', A, size(A,1), size(A,2))


  A = 0.d0
  iunit = getunitandopen('nuclear_ao_complex','r')
  do 
    read (iunit,*,end=12) i,j, int_re, int_im
    A(i,j) = dcmplx(int_re,int_im)
    A(j,i) = dcmplx(int_re,-int_im)
  enddo
  12 continue
  close(iunit)
  call write_one_e_integrals('ao_ne_integral', A, size(A,1), size(A,2))

  call write_one_e_integrals('ao_pseudo_integral', ao_pseudo_integral,&
        size(ao_pseudo_integral,1), size(ao_pseudo_integral,2))


  call ezfio_set_integrals_monoelec_disk_access_ao_one_integrals("Read")

  allocate(buffer_i(ao_num**4), buffer_values(ao_num**4))
   
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
  enddo
  13 continue
  close(iunit)
  
  call insert_into_ao_integrals_map_2fold(n_integrals,buffer_re,buffer_values_re)
  call insert_into_ao_integrals_map_2fold(n_integrals,buffer_im,buffer_values_im)

  call map_sort(ao_integrals_map_2fold)
  call map_unique(ao_integrals_map_2fold)

  call map_save_to_disk(trim(ezfio_filename)//'/work/ao_ints',ao_integrals_map_2fold)
  call ezfio_set_integrals_bielec_disk_access_ao_integrals('Read')

end
