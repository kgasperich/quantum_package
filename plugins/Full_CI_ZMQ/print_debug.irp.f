subroutine print_debug_wf
  use bitmasks
  implicit none
  integer :: i,j
  integer :: iunit
  integer :: getunitandopen
  integer, allocatable :: H_matrix_degree(:,:)
  double precision, allocatable :: H_matrix_phase(:,:)
  integer :: degree
  integer(bit_kind), allocatable :: keys_tmp(:,:,:)
  character*(2048)                :: output(2)

  allocate(keys_tmp(N_int,2,N_det))
  
  iunit = getunitandopen('wfdump.txt','a')
!  write(iunit,fmtstring) arg1,arg2,...
  write(iunit,*) '# ==================================================='
  write(iunit,*) 'N_det      = ',N_det
  write(iunit,*) 'N_states   = ',N_states

!  print dets
  do i = 1, N_det
!    print*,''
    call bitstring_to_hexa( output(1), psi_det(1,1,i), N_int )
    call bitstring_to_hexa( output(2), psi_det(1,2,i), N_int )
    write(iunit,*) i,' | ', trim(output(1)) , ' | ', trim(output(2))
  !  call debug_det(psi_det(1,1,i),N_int)
    do j = 1, N_int
      keys_tmp(j,1,i) = psi_det(j,1,i)
      keys_tmp(j,2,i) = psi_det(j,2,i)
    enddo
  enddo
  if(N_det.ge.10000)then
    print*,'Warning !!!'
    print*,'Number of determinants is ',N_det
    print*,'It means that the H matrix will be enormous !'
    print*,'stoppping ..'
    stop
  endif
 
   

  allocate(H_matrix_degree(N_det,N_det),H_matrix_phase(N_det,N_det))
  integer         :: exc(0:2,2,2)
  double precision  :: phase
  do i = 1, N_det
   do j = i, N_det 
    call get_excitation_degree(psi_det(1,1,i),psi_det(1,1,j),degree,N_int)
    H_matrix_degree(i,j) = degree
    H_matrix_degree(j,i) = degree
    phase = 0.d0
    if(degree==1.or.degree==2)then
     call get_excitation(psi_det(1,1,i),psi_det(1,1,j),exc,degree,phase,N_int)
    endif
    H_matrix_phase(i,j) = phase
    H_matrix_phase(j,i) = phase
   enddo
  enddo
  write(iunit,*) 'H matrix '
  double precision :: ref_h_matrix,s2
  ref_h_matrix = H_matrix_all_dets(1,1)
  write(iunit,*)'HF like determinant energy = ',ref_bitmask_energy+nuclear_repulsion
  write(iunit,*)'Ref element of H_matrix    = ',ref_h_matrix+nuclear_repulsion
  write(iunit,*)'Printing the H matrix ...'
  write(iunit,*)''
  write(iunit,*)''
 !do i = 1, N_det
 ! H_matrix_all_dets(i,i) -= ref_h_matrix
 !enddo
 
  do i = 1, N_det
   H_matrix_all_dets(i,i) += nuclear_repulsion
  enddo
  
  do i = 1, N_det
   write(iunit,'(I5,X,A3,1000(F16.7))')i,' | ',H_matrix_all_dets(i,:)
  enddo
 
  write(iunit,*)''
  write(iunit,*)''
  write(iunit,*)''
  write(iunit,*)'Printing the degree of excitations within the H matrix'
  write(iunit,*)''
  write(iunit,*)''
  do i = 1, N_det
   write(iunit,'(I5,X,A3,X,1000(I1,X))')i,' | ',H_matrix_degree(i,:)
  enddo
 
 
  write(iunit,*)''
  write(iunit,*)''
  write(iunit,*)'Printing the phase of the Hamiltonian matrix elements '
  write(iunit,*)''
  write(iunit,*)''
  do i = 1, N_det
   write(iunit,'(I5,X,A3,X,1000(F3.0,X))')i,' | ',H_matrix_phase(i,:)
  enddo
  write(iunit,*)''
 
 
  complex*16, allocatable  :: eigenvectors(:,:)
  double precision, allocatable  :: eigenvalues(:)
  double precision, allocatable  :: s2_eigvalues(:)
  allocate (eigenvectors(size(H_matrix_all_dets,1),N_det))
  allocate (eigenvalues(N_det),s2_eigvalues(N_det))
  call lapack_diag_z(eigenvalues,eigenvectors,                       &
      H_matrix_all_dets,size(H_matrix_all_dets,1),N_det)
  write(iunit,*)'Two first eigenvectors '
  call u_0_S2_u_0(s2_eigvalues,eigenvectors,n_det,keys_tmp,N_int,N_det,size(eigenvectors,1))
  write(iunit,*)'s2,e : '
  do j =1, N_det
    write(iunit,*)s2_eigvalues(j),' | ',eigenvalues(j)
  enddo
  write(iunit,*)'coefs : '
  do i = 1, N_det
    write(iunit,'(I5,X,A3,1000(F16.7))')i,' | ',eigenvectors(i,:)
!     write(iunit,*)'i = ',i,eigenvectors(i,j)
  enddo
  close(iunit)
end
