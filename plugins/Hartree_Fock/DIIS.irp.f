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
  allocate(                                                          &
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
  
  call zgemm('N','N',AO_num,AO_num,AO_num,                           &
        (1.d0,0.d0), &
        scratch,size(scratch,1),                                     &
        AO_Overlap,size(AO_Overlap,1),                               &
        (0.d0,0.d0), &
        FPS_SPF_Matrix_AO,size(FPS_SPF_Matrix_AO,1))


  ! Compute SP

  call zgemm('N','N',AO_num,AO_num,AO_num,                           &
        (1.d0,0.d0), &
        AO_Overlap,size(AO_Overlap,1),                               &
        HF_Density_Matrix_AO,size(HF_Density_Matrix_AO,1),           &
        (0.d0,0.d0), &
        scratch,size(scratch,1))

  
  ! Compute FPS - SPF
  
  call zgemm('N','N',AO_num,AO_num,AO_num,                           &
      (-1.d0,0.d0),                                                  &
      scratch,Size(scratch,1),                                       &
      Fock_Matrix_AO,Size(Fock_Matrix_AO,1),                         &
      (1.d0,0.d0),                                                   &
      FPS_SPF_Matrix_AO,Size(FPS_SPF_Matrix_AO,1))

  deallocate(scratch)

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
   
   double precision, allocatable  :: work(:)
   integer                        :: lwork,info,lrwork
   integer                        :: i,j

   double precision, allocatable :: rwork(:)
   complex*16, allocatable       :: scratch(:,:)
   
   allocate( scratch(AO_num,AO_num) )
 
! Calculate Fock matrix in orthogonal basis: F' = Xt.F.X
  
  call zgemm('C','N',AO_num,AO_num,AO_num,                           &
         (1.d0,0.d0),                                                   &
         S_half_inv,size(S_half_inv,1),        &
         Fock_Matrix_AO,Size(Fock_Matrix_AO,1),                         &
         (0.d0,0.d0),                                                   &
         eigenvectors_Fock_matrix_AO,size(eigenvectors_Fock_matrix_AO,1))

  call zgemm('N','N',AO_num,AO_num,AO_num,                           &
         (1.d0,0.d0),                                                   &
         eigenvectors_Fock_matrix_AO,size(eigenvectors_Fock_matrix_AO,1), &
         S_half_inv,size(S_half_inv,1),        &
         (0.d0,0.d0),                                                   &
         scratch,size(scratch,1))

! Diagonalize F' to obtain eigenvectors in orthogonal basis C' and eigenvalues
 
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
  
  call zgemm('N','N',AO_num,AO_num,AO_num,                           &
         (1.d0,0.d0),                                                   &
         S_half_inv,size(S_half_inv,1),        &
         scratch,size(scratch,1),                &
         (0.d0,0.d0),                                                   &
         eigenvectors_Fock_matrix_AO,size(eigenvectors_Fock_matrix_AO,1))

  deallocate(rwork, scratch)
   
END_PROVIDER

