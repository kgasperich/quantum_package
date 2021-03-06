BEGIN_PROVIDER [ integer, mo_tot_num ]
  implicit none
  BEGIN_DOC
  ! Number of MOs
  END_DOC
  
  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    call ezfio_has_mo_basis_mo_tot_num(has)
  endif
  IRP_IF MPI
    include 'mpif.h'
    integer                        :: ierr
    call MPI_BCAST( has, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read mo_tot_num with MPI'
    endif
  IRP_ENDIF
  if (.not.has) then
    mo_tot_num = ao_ortho_canonical_num
  else
    if (mpi_master) then
      call ezfio_get_mo_basis_mo_tot_num(mo_tot_num)
    endif
    IRP_IF MPI
      call MPI_BCAST( mo_tot_num, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
        stop 'Unable to read mo_tot_num with MPI'
      endif
    IRP_ENDIF
  endif
  call write_int(6,mo_tot_num,'mo_tot_num')
  ASSERT (mo_tot_num > 0)
  
END_PROVIDER

BEGIN_PROVIDER [integer, mo_tot_num_per_kpt ]
  implicit none
  mo_tot_num_per_kpt = mo_tot_num / num_kpts
END_PROVIDER

BEGIN_PROVIDER [ integer, mo_num ]
 implicit none
 BEGIN_DOC
 ! mo_tot_num without the highest deleted MOs
 END_DOC
 mo_num = mo_tot_num
 integer :: i
 mo_num = mo_tot_num
 do i=mo_tot_num,1,-1
   if (mo_class(i) == 'Deleted') then
     mo_num -= 1
!   else
!     exit
   endif
 enddo
END_PROVIDER

BEGIN_PROVIDER [integer, mo_num_per_kpt ]
  implicit none
  mo_num_per_kpt = mo_num / num_kpts
END_PROVIDER

 BEGIN_PROVIDER [ double precision, mo_coef_real, (ao_num,mo_tot_num) ]
&BEGIN_PROVIDER [ double precision, mo_coef_imag, (ao_num,mo_tot_num) ]
  implicit none
  BEGIN_DOC
  ! Molecular orbital coefficients on AO basis set
  ! mo_coef(i,j) = coefficient of the ith ao on the jth mo
  ! mo_label : Label characterizing the MOS (local, canonical, natural, etc)
  END_DOC
  integer                        :: i, j
  double precision, allocatable  :: buffer(:,:)
  logical                        :: exists_real,exists_imag
  PROVIDE ezfio_filename 
  

  if (mpi_master) then
    ! Coefs
    call ezfio_has_mo_basis_mo_coef_real(exists_real)
    call ezfio_has_mo_basis_mo_coef_imag(exists_imag)
  endif
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST(exists_real, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read mo_coef_real with MPI'
    endif
    call MPI_BCAST(exists_imag, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read mo_coef_imag with MPI'
    endif
  IRP_ENDIF

  if (exists_real) then
    if (mpi_master) then
      call ezfio_get_mo_basis_mo_coef_real(mo_coef_real)
      write(*,*) 'Read  mo_coef_real'
    endif
    IRP_IF MPI
      call MPI_BCAST( mo_coef_real, mo_tot_num*ao_num, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
        stop 'Unable to read mo_coef_real with MPI'
      endif
    IRP_ENDIF
    if (exists_imag) then
      if (mpi_master) then
        call ezfio_get_mo_basis_mo_coef_imag(mo_coef_imag)
        write(*,*) 'Read  mo_coef_imag'
      endif
      IRP_IF MPI
        call MPI_BCAST( mo_coef_imag, mo_tot_num*ao_num, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        if (ierr /= MPI_SUCCESS) then
          stop 'Unable to read mo_coef_imag with MPI'
        endif
      IRP_ENDIF
    else
      mo_coef_imag = 0.d0
    endif
  else
    ! Orthonormalized AO basis
    do i=1,mo_tot_num
      do j=1,ao_num
        mo_coef_real(j,i) = ao_ortho_canonical_coef(j,i)
      enddo
    enddo
    mo_coef_imag = 0.d0
  endif

END_PROVIDER


BEGIN_PROVIDER [ complex*16, mo_coef, (ao_num,mo_tot_num) ]
  implicit none
  BEGIN_DOC
  ! Molecular orbital coefficients on AO basis set
  ! mo_coef(i,j) = coefficient of the ith ao on the jth mo
  ! mo_label : Label characterizing the MOS (local, canonical, natural, etc)
  END_DOC
  integer                        :: i, j
  PROVIDE mo_coef_real mo_coef_imag 
    ! Orthonormalized AO basis
    do i=1,mo_tot_num
      do j=1,ao_num
        mo_coef(j,i) = cmplx(mo_coef_real(j,i),mo_coef_imag(j,i))
      enddo
    enddo

END_PROVIDER

BEGIN_PROVIDER [ complex*16, mo_coef_kpts_trunc, (ao_num_per_kpt,mo_num_per_kpt,num_kpts) ]
  implicit none
  BEGIN_DOC
  ! Molecular orbital coefficients on AO basis set
  ! mo_coef_kpts(i,j,k) = coefficient of the ith ao on the jth mo in kpt k
  ! mo_label : Label characterizing the MOS (local, canonical, natural, etc)
  END_DOC
  integer                        :: i, j, k
  integer                        :: ii, jj
  PROVIDE mo_coef
    ! Orthonormalized AO basis
    print*,'mo_coef_trunc'
    do k = 1,num_kpts
      do i=1,mo_num_per_kpt
        ii = i + (k-1)*mo_tot_num_per_kpt
        do j=1,ao_num_per_kpt
          jj = j + (k-1)*ao_num_per_kpt
          mo_coef_kpts_trunc(j,i,k) = mo_coef(jj,ii)
!          print*,j,i,k,jj,ii,mo_coef(jj,ii)
        enddo
      enddo
    enddo

END_PROVIDER

BEGIN_PROVIDER [ complex*16, mo_coef_kpts, (ao_num_per_kpt,mo_tot_num_per_kpt,num_kpts) ]
  implicit none
  BEGIN_DOC
  ! Molecular orbital coefficients on AO basis set
  ! mo_coef_kpts(i,j,k) = coefficient of the ith ao on the jth mo in kpt k
  ! mo_label : Label characterizing the MOS (local, canonical, natural, etc)
  END_DOC
  integer                        :: i, j, k
  integer                        :: ii, jj
  PROVIDE mo_coef
    ! Orthonormalized AO basis
    do k = 1,num_kpts
      do i=1,mo_tot_num_per_kpt
        ii = i + (k-1)*mo_tot_num_per_kpt
        do j=1,ao_num_per_kpt
          jj = j + (k-1)*ao_num_per_kpt
          mo_coef_kpts(j,i,k) = mo_coef(jj,ii)
        enddo
      enddo
    enddo 

END_PROVIDER

BEGIN_PROVIDER [ complex*16, mo_coef_in_ao_ortho_basis, (ao_num, mo_tot_num) ]
 implicit none
 BEGIN_DOC
 ! MO coefficients in orthogonalized AO basis
 !
 ! C^(-1).C_mo
 END_DOC

 call zgemm('N','N',ao_num,mo_tot_num,ao_num,(1.d0,0.d0),            &
     ao_ortho_canonical_coef_inv,                                    &
     size(ao_ortho_canonical_coef_inv,1),                            &
     mo_coef, size(mo_coef,1), (0.d0,0.d0),                          &
     mo_coef_in_ao_ortho_basis, size(mo_coef_in_ao_ortho_basis,1))

END_PROVIDER

BEGIN_PROVIDER [ character*(64), mo_label ]
  implicit none
  BEGIN_DOC
  ! Molecular orbital coefficients on AO basis set
  ! mo_coef(i,j) = coefficient of the ith ao on the jth mo
  ! mo_label : Label characterizing the MOS (local, canonical, natural, etc)
  END_DOC

  logical                        :: exists
  PROVIDE ezfio_filename
  if (mpi_master) then
    call ezfio_has_mo_basis_mo_label(exists)
    if (exists) then
      call ezfio_get_mo_basis_mo_label(mo_label)
      mo_label = trim(mo_label)
    else
      mo_label = 'no_label'
    endif
    write(*,*) '* mo_label          ', trim(mo_label)
  endif
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( mo_label, 64, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read mo_label with MPI'
    endif
  IRP_ENDIF

END_PROVIDER

BEGIN_PROVIDER [ complex*16, mo_coef_transp, (mo_tot_num,ao_num) ]
  implicit none
  BEGIN_DOC
  ! Molecular orbital coefficients on AO basis set
  END_DOC
  integer                        :: i, j
  
  do j=1,ao_num
    do i=1,mo_tot_num
      mo_coef_transp(i,j) = mo_coef(j,i)
    enddo
  enddo
  
END_PROVIDER

BEGIN_PROVIDER [ complex*16, mo_coef_conjg_transp, (mo_tot_num,ao_num) ]
  implicit none
  BEGIN_DOC
  ! Molecular orbital coefficients on AO basis set
  END_DOC
  integer                        :: i, j
  
  do j=1,ao_num
    do i=1,mo_tot_num
      mo_coef_conjg_transp(i,j) = conjg(mo_coef(j,i))
    enddo
  enddo
  
END_PROVIDER

 BEGIN_PROVIDER [ complex*16, mo_coef_kpts_transp, (mo_tot_num_per_kpt,ao_num_per_kpt,num_kpts) ]
&BEGIN_PROVIDER [ complex*16, mo_coef_kpts_conjg_transp, (mo_tot_num_per_kpt,ao_num_per_kpt,num_kpts) ]
  implicit none
  BEGIN_DOC
  ! Molecular orbital coefficients on AO basis set
  END_DOC
  
  mo_coef_kpts_transp = reshape(mo_coef_kpts,(/mo_tot_num_per_kpt,ao_num_per_kpt,num_kpts/),order=(/2,1,3/))
  mo_coef_kpts_conjg_transp = conjg(mo_coef_kpts_transp)
  
END_PROVIDER

BEGIN_PROVIDER [ complex*16, S_mo_coef, (ao_num, mo_tot_num) ]
  implicit none
  BEGIN_DOC
  ! Product S.C where S is the overlap matrix in the AO basis and C the mo_coef matrix.
  END_DOC

  call zgemm('N','N',ao_num, mo_tot_num, ao_num, (1.d0,0.d0), &
        ao_overlap, size(ao_overlap,1), &
        mo_coef, size(mo_coef,1), &
        (0.d0,0.d0), &
        S_mo_coef, size(S_mo_coef,1))

END_PROVIDER

BEGIN_PROVIDER [ double precision, mo_occ, (mo_tot_num) ]
  implicit none
  BEGIN_DOC
  ! MO occupation numbers
  END_DOC
  PROVIDE ezfio_filename elec_beta_num elec_alpha_num 
  if (mpi_master) then
    logical :: exists
    call ezfio_has_mo_basis_mo_occ(exists)
    if (exists) then
      call ezfio_get_mo_basis_mo_occ(mo_occ)
    else
      mo_occ = 0.d0
      integer :: i
      do i=1,elec_beta_num
        mo_occ(i) = 2.d0
      enddo
      do i=elec_beta_num+1,elec_alpha_num
        mo_occ(i) = 1.d0
      enddo
    endif
    write(*,*) 'Read mo_occ'
  endif
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( mo_occ, mo_tot_num, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read mo_occ with MPI'
    endif
  IRP_ENDIF

END_PROVIDER


subroutine complex_ao_to_mo(A_ao,LDA_ao,A_mo,LDA_mo)
  implicit none
  BEGIN_DOC
  ! Transform A from the AO basis to the MO basis
  ! where A is complex in the AO basis
  !
  ! Ct.A_ao.C
  END_DOC
  integer, intent(in)            :: LDA_ao,LDA_mo
  complex*16, intent(in)   :: A_ao(LDA_ao,ao_num)
  complex*16, intent(out)  :: A_mo(LDA_mo,mo_tot_num)
  complex*16, allocatable  :: T(:,:)
  
  allocate ( T(ao_num,mo_tot_num) )
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: T
  
  call zgemm('N','N', ao_num, mo_tot_num, ao_num,                    &
      (1.d0,0.d0), A_ao,LDA_ao,                                      &
      mo_coef, size(mo_coef,1),                                      &
      (0.d0,0.d0), T, size(T,1))
  
  call zgemm('C','N', mo_tot_num, mo_tot_num, ao_num,                &
      (1.d0,0.d0), mo_coef,size(mo_coef,1),                          &
      T, ao_num,                                                     &   
      (0.d0,0.d0), A_mo, size(A_mo,1))
  
  deallocate(T)
end

subroutine real_ao_to_mo(A_ao,LDA_ao,A_mo,LDA_mo)
  implicit none
  BEGIN_DOC
  ! Transform A from the AO basis to the MO basis
  ! where A is real in the AO basis
  !
  ! Ct.A_ao.C
  END_DOC
  integer, intent(in)            :: LDA_ao,LDA_mo
  double precision, intent(in)   :: A_ao(LDA_ao,ao_num)
  complex*16, intent(out)  :: A_mo(LDA_mo,mo_tot_num)
  complex*16, allocatable  :: T(:,:)
  double precision, allocatable :: work(:)
  integer :: sze_work

  sze_work = 2 * ao_num * mo_tot_num
  allocate ( T(ao_num,mo_tot_num) , work(sze_work) )
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: T
  
  call zlarcm(ao_num, mo_tot_num, A_ao,LDA_ao,                       &
              mo_coef, size(mo_coef,1),                              &
              T, size(T,1),                                          &
              work)
  
  call zgemm('C','N', mo_tot_num, mo_tot_num, ao_num,                &
      (1.d0,0.d0), mo_coef,size(mo_coef,1),                          &
      T, ao_num,                                                     &   
      (0.d0,0.d0), A_mo, size(A_mo,1))
  
  deallocate(T,work)
end

subroutine mo_to_ao(A_mo,LDA_mo,A_ao,LDA_ao)
  implicit none
  BEGIN_DOC
  ! Transform A from the MO basis to the AO basis
  !
  ! (S.C).A_mo.(S.C)t
  END_DOC
  integer, intent(in)            :: LDA_ao,LDA_mo
  complex*16, intent(in)         :: A_mo(LDA_mo,mo_tot_num)
  complex*16, intent(out)        :: A_ao(LDA_ao,ao_num)
  complex*16, allocatable        :: T(:,:)
  
  allocate ( T(mo_tot_num,ao_num) )
  
  call zgemm('N','C', mo_tot_num, ao_num, mo_tot_num,                &
      (1.d0,0.d0), A_mo,size(A_mo,1),                                &
      S_mo_coef, size(S_mo_coef,1),                                  &
      (0.d0,0.d0), T, size(T,1))
  
  call zgemm('N','N', ao_num, ao_num, mo_tot_num,                    &
      (1.d0,0.d0), S_mo_coef, size(S_mo_coef,1),                     &
      T, size(T,1),                                                  &
      (0.d0,0.d0), A_ao, size(A_ao,1))
  
  deallocate(T)
end

subroutine mo_to_ao_no_overlap(A_mo,LDA_mo,A_ao,LDA_ao)
  implicit none
  BEGIN_DOC
  ! Transform A from the MO basis to the S^-1 AO basis
  ! Useful for density matrix
  END_DOC
  integer, intent(in)            :: LDA_ao,LDA_mo
  complex*16, intent(in)         :: A_mo(LDA_mo,mo_tot_num)
  complex*16, intent(out)        :: A_ao(LDA_ao,ao_num)
  complex*16, allocatable        :: T(:,:)
  
  allocate ( T(mo_tot_num,ao_num) )
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: T
  
  call zgemm('N','C', mo_tot_num, ao_num, mo_tot_num,                &
      (1.d0,0.d0), A_mo,size(A_mo,1),                                &
      mo_coef, size(mo_coef,1),                                      &
      (0.d0,0.d0), T, size(T,1))
  
  call zgemm('N','N', ao_num, ao_num, mo_tot_num,                    &
      (1.d0,0.d0), mo_coef,size(mo_coef,1),                          &
      T, size(T,1),                                                  &
      (0.d0,0.d0), A_ao, size(A_ao,1))
  
  deallocate(T)
end

subroutine mix_mo_jk(j,k)
 implicit none
 integer, intent(in) :: j,k
 integer :: i,i_plus,i_minus
 BEGIN_DOC
! subroutine that rotates the jth MO with the kth MO
! to give two new MO's that are 
!         '+' = 1/sqrt(2) (|j> + |k>) 
!         '-' = 1/sqrt(2) (|j> - |k>)
! by convention, the '+' MO is in the lower index (min(j,k))
! by convention, the '-' MO is in the greater index (max(j,k))
 END_DOC
 double precision :: dsqrt_2
 complex*16       :: array_tmp(ao_num,2)
 if(j==k)then
  print*,'You want to mix two orbitals that are the same !'
  print*,'It does not make sense ... '
  print*,'Stopping ...'
  stop
 endif
 array_tmp = (0.d0,0.d0)
 dsqrt_2 = 1.d0/dsqrt(2.d0)
 do i = 1, ao_num
  array_tmp(i,1) = dsqrt_2 * (mo_coef(i,j) + mo_coef(i,k))
  array_tmp(i,2) = dsqrt_2 * (mo_coef(i,j) - mo_coef(i,k))
 enddo
 i_plus = min(j,k)
 i_minus = max(j,k)
 do i = 1, ao_num
  mo_coef(i,i_plus) = array_tmp(i,1)
  mo_coef(i,i_minus) = array_tmp(i,2)
 enddo

end

subroutine ao_ortho_cano_to_ao(A_ao,LDA_ao,A,LDA)
  implicit none
  BEGIN_DOC
  ! Transform A from the AO basis to the orthogonal AO basis
  ! where A is complex in both bases
  !
  ! C^(-1).A_ao.Ct^(-1)
  END_DOC
  integer, intent(in)            :: LDA_ao,LDA
  complex*16, intent(in)         :: A_ao(LDA_ao,*)
  complex*16, intent(out)        :: A(LDA,*)
  complex*16, allocatable        :: T(:,:)
  
  allocate ( T(ao_num,ao_num) )
  
  call zgemm('C','N', ao_num, ao_num, ao_num,                        &
      (1.d0,0.d0), ao_ortho_canonical_coef_inv,              &
      size(ao_ortho_canonical_coef_inv,1),                   &
      A_ao,size(A_ao,1),                                             &
      (0.d0,0.d0), T, size(T,1))
  
  call zgemm('N','N', ao_num, ao_num, ao_num, (1.d0,0.d0),           &
      T, size(T,1), ao_ortho_canonical_coef_inv,             &
      size(ao_ortho_canonical_coef_inv,1),                   &
      (0.d0,0.d0), A, size(A,1))
  
  deallocate(T)
end

