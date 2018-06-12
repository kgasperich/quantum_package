BEGIN_PROVIDER [ double precision, threshold_DIIS_nonzero ]
 implicit none
 BEGIN_DOC
 ! If threshold_DIIS is zero, choose sqrt(thresh_scf)
 END_DOC
 if (threshold_DIIS == 0.d0) then
   threshold_DIIS_nonzero = dsqrt(thresh_scf)
 else
   threshold_DIIS_nonzero = threshold_DIIS
 endif
 ASSERT (threshold_DIIS_nonzero >= 0.d0)

END_PROVIDER

BEGIN_PROVIDER [complex*16, FPS_SPF_Matrix_AO, (AO_num, AO_num)]
  implicit none
  BEGIN_DOC
  !   Commutator FPS - SPF
  END_DOC
  complex*16, allocatable  :: scratch(:,:)
  double precision, allocatable :: rwork(:)
  allocate(                                                          &
      rwork(2*AO_num*AO_num),                                  &
      scratch(AO_num, AO_num)                                  &
      )
  
  ! Compute FP
  
  call zgemm('N','N',AO_num,AO_num,AO_num,                           &
      (1.d0,0.d0),                                                   &
      Fock_Matrix_AO,Size(Fock_Matrix_AO,1),                         &
      HF_Density_Matrix_AO,Size(HF_Density_Matrix_AO,1),             &
      (0.d0,0.d0),                                                   &
      scratch,Size(scratch,1))
  
  ! Compute FPS
  
  call zlacrm(AO_num,AO_num,                                         &
        scratch,size(scratch,1),                                     &
        AO_Overlap,size(AO_Overlap,1),                               &
        FPS_SPF_Matrix_AO,size(FPS_SPF_Matrix_AO,1),rwork)


  ! Compute SP
  
  call zlarcm(AO_num,AO_num,                                         &
        AO_Overlap,size(AO_Overlap,1),                               &
        HF_Density_Matrix_AO,size(HF_Density_Matrix_AO,1),           &
        scratch,size(scratch,1),rwork)

  
  ! Compute FPS - SPF
  
  call zgemm('N','N',AO_num,AO_num,AO_num,                           &
      (-1.d0,0.d0),                                                  &
      scratch,Size(scratch,1),                                       &
      Fock_Matrix_AO,Size(Fock_Matrix_AO,1),                         &
      (1.d0,0.d0),                                                   &
      FPS_SPF_Matrix_AO,Size(FPS_SPF_Matrix_AO,1))

  deallocate(scratch, rwork)

END_PROVIDER

bEGIN_PROVIDER [complex*16, FPS_SPF_Matrix_MO, (mo_tot_num, mo_tot_num)]
  implicit none
  begin_doc 
!   Commutator FPS - SPF in MO basis
  end_doc
  call complex_ao_to_mo(FPS_SPF_Matrix_AO, size(FPS_SPF_Matrix_AO,1), &
     FPS_SPF_Matrix_MO, size(FPS_SPF_Matrix_MO,1))
END_PROVIDER


 BEGIN_PROVIDER [ double precision, eigenvalues_Fock_matrix_AO, (AO_num) ]
&BEGIN_PROVIDER [ complex*16, eigenvectors_Fock_matrix_AO, (AO_num,AO_num) ]

   BEGIN_DOC
   ! Eigenvalues and eigenvectors of the Fock matrix over the AO basis
   END_DOC

   implicit none
   
   double precision, allocatable  :: work(:),Xt(:,:)
   integer                        :: lwork,info,lrwork
   integer                        :: i,j

   double precision, allocatable :: rwork(:)
   complex*16, allocatable       :: scratch(:,:)
   
   lrwork = 2*ao_num*ao_num
   allocate(                                                         &
       scratch(AO_num,AO_num),                                 &
       rwork(lrwork),                                                  &
       Xt(AO_num,AO_num)                                             &
       )
 
! Calculate Xt

  do i=1,AO_num
    do j=1,AO_num
      Xt(i,j) = S_half_inv(j,i)
    enddo
  enddo

! Calculate Fock matrix in orthogonal basis: F' = Xt.F.X

  call zlacrm(ao_num,ao_num, &
         Fock_matrix_AO,size(Fock_matrix_AO,1), &
         S_half_inv,size(S_half_inv,1),        &
         eigenvectors_Fock_matrix_AO,size(eigenvectors_Fock_matrix_AO,1), &
         rwork)

  call zlarcm(ao_num,ao_num, &
         Xt,size(Xt,1), &
         eigenvectors_Fock_matrix_AO,size(eigenvectors_Fock_matrix_AO,1), &
         scratch,size(scratch,1), &
         rwork)
     
! Diagonalize F' to obtain eigenvectors in orthogonal basis C' and eigenvalues
 
   deallocate(rwork)
   lrwork = 3*ao_num - 2
   allocate( rwork(lrwork), work(1) )
   lwork = -1

   call zheev('V','U',AO_num,       &
        scratch,size(scratch,1),    &
        eigenvalues_Fock_matrix_AO, &
        work,lwork,rwork,info)

   lwork = work(1)
   deallocate(work)
   allocate(work(lwork))

   call zheev('V','U',AO_num,       &
        scratch,size(scratch,1),    &
        eigenvalues_Fock_matrix_AO, &
        work,lwork,rwork,info)

   deallocate(work,rwork)

   if(info /= 0) then
     print *,  irp_here//' failed : ', info
     stop 1
   endif

! Back-transform eigenvectors: C =X.C'

  lrwork = 2*ao_num*ao_num
  allocate(rwork(lrwork))
  call zlarcm(ao_num,ao_num,     &
       S_half_inv,size(S_half_inv,1),        &
       scratch,size(scratch,1),                &
       eigenvectors_Fock_matrix_AO,size(eigenvectors_Fock_matrix_AO,1), &
       rwork)
  deallocate(rwork, scratch, Xt)
   
END_PROVIDER

