program read_integrals

  PROVIDE ezfio_filename

  call ezfio_set_integrals_monoelec_disk_access_ao_one_integrals("None")
  call ezfio_set_utils_disk_access_ao_overlap_integrals("None")
  disk_access_ao_one_integrals='None'
  disk_access_ao_overlap_integrals='None'
  TOUCH disk_access_ao_one_integrals disk_access_ao_overlap_integrals
  call run_read_ao_mono_complex
  call ezfio_set_integrals_monoelec_disk_access_ao_one_integrals("Read")
  call ezfio_set_utils_disk_access_ao_overlap_integrals("Read")
  disk_access_ao_one_integrals='Read'
  disk_access_ao_overlap_integrals='Read'
  TOUCH disk_access_ao_one_integrals disk_access_ao_overlap_integrals

  call ezfio_set_utils_disk_access_kconserv("None")
  disk_access_kconserv='None'
  TOUCH disk_access_kconserv
  call run_read_kconserv
  call ezfio_set_utils_disk_access_kconserv("Read")
  disk_access_kconserv='Read'
  TOUCH disk_access_kconserv

  call ezfio_set_integrals_bielec_disk_access_df_ao_integral_array("None")
  disk_access_df_ao_integral_array='None'
  TOUCH disk_access_df_ao_integral_array
  call run_read_ao_df_complex
  call ezfio_set_integrals_bielec_disk_access_df_ao_integral_array("Read")
  disk_access_df_ao_integral_array='Read'
  TOUCH disk_access_df_ao_integral_array

!  call ezfio_set_integrals_monoelec_disk_access_ao_one_integrals("None")
  call run_read_mo_coef_complex
  call run_make_mo_df_integrals_complex
end

subroutine run_read_ao_mono_complex
  use map_module
  implicit none
  
  integer :: iunit
  integer :: getunitandopen

  integer ::i,j
  double precision :: int_re, int_im
  complex*16, allocatable :: A(:,:)

  allocate (A(ao_num,ao_num))


  A = 0.d0
  iunit = getunitandopen('kinetic_ao_complex','r')
  do 
    read (iunit,*,end=10) i,j, int_re, int_im
    if (i.eq.j) then
      int_im = 0.d0
      A(i,i) = dcmplx(int_re,int_im)
    else
      A(i,j) = dcmplx(int_re,int_im)
      A(j,i) = dcmplx(int_re,-int_im)
    endif
  enddo
  10 continue
  close(iunit)
  call write_one_e_integrals_complex('ao_kinetic_integral', A, size(A,1), size(A,2))


  A = 0.d0
  iunit = getunitandopen('ne_ao_complex','r')
  do
    read (iunit,*,end=11) i,j, int_re, int_im
    if (i.eq.j) then
      int_im = 0.d0
      A(i,i) = dcmplx(int_re,int_im)
    else
      A(i,j) = dcmplx(int_re,int_im)
      A(j,i) = dcmplx(int_re,-int_im)
    endif
  enddo
  11 continue
  close(iunit)
  call write_one_e_integrals_complex('ao_ne_integral', A, size(A,1), size(A,2))


  A = 0.d0
  iunit = getunitandopen('overlap_ao_complex','r')
  do 
    read (iunit,*,end=12) i,j, int_re, int_im
    if (i.eq.j) then
      int_im = 0.d0
      A(i,i) = dcmplx(int_re,int_im)
    else
      A(i,j) = dcmplx(int_re,int_im)
      A(j,i) = dcmplx(int_re,-int_im)
    endif
  enddo
  12 continue
  close(iunit)
  call write_one_e_integrals_complex('ao_overlap', A, size(A,1), size(A,2))


  deallocate(A)

!  call ezfio_set_integrals_monoelec_disk_access_ao_one_integrals("Read")
!  call ezfio_set_utils_disk_access_ao_overlap_integrals("Read")

end

subroutine run_read_kconserv
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

!  call ezfio_set_utils_disk_access_kconserv("Read")

end

subroutine run_read_ao_df_complex
  use map_module
  implicit none
  
  integer :: iunit
  integer :: getunitandopen

  integer :: kikj,mu,i,j
  double precision :: int_re, int_im
  complex*16, allocatable :: A(:,:,:,:)

  allocate (A(ao_num_per_kpt,ao_num_per_kpt,df_num,num_kpt_pairs))


  A = 0.d0
  iunit = getunitandopen('df_ao_integral_array','r')
  do 
  !  read (iunit,*,end=10) i,j,mu,kikj, int_re, int_im
  !  A(i,j,mu,kikj) = dcmplx(int_re,int_im)
    read (iunit,*,end=10) j,i,mu,kikj, int_re, int_im
    A(i,j,mu,kikj) = dcmplx(int_re,-int_im)
  enddo
  10 continue
  close(iunit)
  call write_df_integral_array_file('df_ao_integral_array', A, size(A,1), size(A,3), size(A,4))

  deallocate(A)

!  call ezfio_set_integrals_bielec_disk_access_df_ao_integral_array("Read")

end

subroutine run_read_mo_coef_complex
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
  mo_label = "None"
  call save_mos

end

subroutine run_make_mo_df_integrals_complex
  use map_module
  implicit none
  disk_access_df_mo_integral_array='Write'
  TOUCH disk_access_df_mo_integral_array
  PROVIDE df_mo_integral_array
end
