! modified from H_S2_u_0_nstates_openmp in Davidson/u0Hu0.irp.f

subroutine h_u_0_homo_openmp(v_0,u_0,sze,spin_hp,sign_hp,idx_hp)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes v_0 = H|u_0> 
  !
  ! Assumes that the determinants are in psi_det
  !
  ! istart, iend, ishift, istep are used in ZMQ parallelization.
  ! 
  ! N_hp is number of holes and particles to be applied
  ! each element of spin_hp is either 1(alpha) or 2(beta)
  ! each element of sign_hp is either 1(particle) or -1(hole)
  ! idx_hp contains orbital indices for holes and particles
  END_DOC
  integer, intent(in)            :: sze
  complex*16, intent(inout)  :: v_0(sze), u_0(sze)
  integer :: k
  integer, intent(in) :: spin_hp, idx_hp
  double precision, intent(in) :: sign_hp

  call cdset_order(u_0(1),psi_bilinear_matrix_order,N_det)
  v_0 = (0.d0,0.d0)
  call h_u_0_homo_openmp_work(v_0,u_0,sze,spin_hp,sign_hp,idx_hp,1,N_det,0,1)

  call cdset_order(v_0(1),psi_bilinear_matrix_order_reverse,N_det)
  call cdset_order(u_0(1),psi_bilinear_matrix_order_reverse,N_det)

end


subroutine h_u_0_homo_openmp_work(v_t,u_t,sze,spin_hp,sign_hp,idx_hp,istart,iend,ishift,istep)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes v_t = H|u_t>
  !
  ! Default should be 1,N_det,0,1
  END_DOC
  integer, intent(in)            :: sze,istart,iend,ishift,istep
  complex*16, intent(in)   :: u_t(N_det)
  complex*16, intent(out)  :: v_t(sze)
  integer, intent(in) :: spin_hp, idx_hp
  double precision, intent(in) :: sign_hp

  
  PROVIDE ref_bitmask_energy N_int

  select case (N_int)
    case (1)
      call H_u_0_homo_openmp_work_1(v_t,u_t,sze,spin_hp,sign_hp,idx_hp,istart,iend,ishift,istep)
    case (2)
      call H_u_0_homo_openmp_work_2(v_t,u_t,sze,spin_hp,sign_hp,idx_hp,istart,iend,ishift,istep)
    case (3)
      call H_u_0_homo_openmp_work_3(v_t,u_t,sze,spin_hp,sign_hp,idx_hp,istart,iend,ishift,istep)
    case (4)
      call H_u_0_homo_openmp_work_4(v_t,u_t,sze,spin_hp,sign_hp,idx_hp,istart,iend,ishift,istep)
    case default
      call H_u_0_homo_openmp_work_N_int(v_t,u_t,sze,spin_hp,sign_hp,idx_hp,istart,iend,ishift,istep)
  end select
end
BEGIN_TEMPLATE

subroutine h_u_0_homo_openmp_work_$N_int(v_t,u_t,sze,spin_hp,sign_hp,idx_hp,istart,iend,ishift,istep)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes v_t = H|u_t> and s_t = S^2 |u_t>
  !
  ! Default should be 1,N_det,0,1
  END_DOC
  integer, intent(in)            :: sze,istart,iend,ishift,istep
  complex*16, intent(in)   :: u_t(N_det)
  complex*16, intent(out)  :: v_t(sze)
  integer, intent(in) :: spin_hp, idx_hp
  double precision, intent(in) :: sign_hp

  complex*16               :: hij
  double precision :: hii 
  integer                        :: i,j,k,l
  integer                        :: k_a, k_b, l_a, l_b, m_a, m_b
  integer                        :: istate
  integer                        :: krow, kcol, krow_b, kcol_b
  integer                        :: lrow, lcol
  integer                        :: mrow, mcol
  integer(bit_kind)              :: spindet($N_int)
  integer(bit_kind)              :: tmp_det($N_int,2)
  integer(bit_kind)              :: tmp_det2($N_int,2)
  integer(bit_kind)              :: tmp_det3($N_int,2)
  integer(bit_kind), allocatable :: buffer(:,:)
  integer                        :: n_doubles
  integer, allocatable           :: doubles(:)
  integer, allocatable           :: singles_a(:)
  integer, allocatable           :: singles_b(:)
  integer, allocatable           :: idx(:), idx0(:)
  integer                        :: maxab, n_singles_a, n_singles_b, kcol_prev
  integer*8                      :: k8

  
  logical :: banned_a1,banned_b1,banned_ab1
  logical :: banned_a2,banned_b2,banned_ab12
  integer                        :: ii,na,nb
  double precision :: hii_hp
  complex*16 :: hij_hp

  maxab = max(N_det_alpha_unique, N_det_beta_unique)+1
  allocate(idx0(maxab))
  
  do i=1,maxab
    idx0(i) = i
  enddo

  ! Prepare the array of all alpha single excitations
  ! -------------------------------------------------

  PROVIDE N_int nthreads_davidson elec_num homo_allowed
  !$OMP PARALLEL DEFAULT(NONE) NUM_THREADS(nthreads_davidson)        &
      !$OMP   SHARED(psi_bilinear_matrix_rows, N_det,                &
      !$OMP          psi_bilinear_matrix_columns,                    &
      !$OMP          psi_det_alpha_unique, psi_det_beta_unique,      &
      !$OMP          n_det_alpha_unique, n_det_beta_unique, N_int,   &
      !$OMP          psi_bilinear_matrix_transp_rows,                &
      !$OMP          psi_bilinear_matrix_transp_columns,             &
      !$OMP          psi_bilinear_matrix_transp_order,         &
      !$OMP          psi_bilinear_matrix_order_transp_reverse,       &
      !$OMP          psi_bilinear_matrix_columns_loc,                &
      !$OMP          psi_bilinear_matrix_transp_rows_loc,            &
      !$OMP          istart, iend, istep, irp_here, v_t,        &
      !$OMP          spin_hp,sign_hp,idx_hp,        &
      !$OMP          elec_num_tab,        &
      !$OMP          ishift, idx0, u_t, maxab)                       &
      !$OMP   PRIVATE(krow, kcol, tmp_det, spindet, k_a, k_b, i,     &
      !$OMP          lcol, lrow, l_a, l_b,                           &
      !$OMP          buffer, doubles, n_doubles,                     &
      !$OMP          tmp_det2, hii, hij, idx, l, kcol_prev,     &
      !$OMP          singles_a, n_singles_a, singles_b,              &
      !$OMP          banned_a1,banned_b1,banned_ab1, &
      !$OMP          banned_a2,banned_b2,banned_ab12, &
      !$OMP          ii, hij_hp, j, hii_hp,na,nb, &
      !$OMP          n_singles_b, k8)
  
  ! Alpha/Beta double excitations
  ! =============================
    
  allocate( buffer($N_int,maxab),                                     &
      singles_a(maxab),                                              &
      singles_b(maxab),                                              &
      doubles(maxab),                                                &
      idx(maxab))

  kcol_prev=-1
  banned_b1=.False.
  ASSERT (iend <= N_det)
  ASSERT (istart > 0)
  ASSERT (istep  > 0)

  !$OMP DO SCHEDULE(dynamic,64)
  do k_a=istart+ishift,iend,istep
  ! iterate over dets in psi

    krow = psi_bilinear_matrix_rows(k_a)
    ASSERT (krow <= N_det_alpha_unique)

    kcol = psi_bilinear_matrix_columns(k_a)
    ASSERT (kcol <= N_det_beta_unique)

    tmp_det(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, krow)
    tmp_det(1:$N_int,2) = psi_det_beta_unique (1:$N_int, kcol)

    if (kcol /= kcol_prev) then !if we've moved to a new unique beta determinant
      call get_homo_banned_spin(tmp_det,banned_b1,spin_hp,sign_hp,idx_hp,2,$N_int)
!      call get_list_hp_banned_spin(tmp_det,N_hp,exc_is_banned_b1,spin_hp,sign_hp,idx_hp,2,$N_int,all_banned_b1)
      if (banned_b1) then
        kcol_prev = kcol
        cycle
      else ! get all unique beta dets connected to this one by a single excitation
        call get_all_spin_singles_$N_int(                              &
            psi_det_beta_unique, idx0,                                 &
            tmp_det(1,2), N_det_beta_unique,                           &
            singles_b, n_singles_b)
        kcol_prev = kcol
      endif
    else
      if (banned_b1) cycle
    endif

    ! at least some beta allowed
    ! check alpha
    call get_homo_banned_spin(tmp_det,banned_a1,spin_hp,sign_hp,idx_hp,1,$N_int)
    if (banned_a1) cycle

!    kcol_prev = kcol ! keep track of old col to see when we've moved to a new one

    ! Loop over singly excited beta columns
    ! -------------------------------------

    do i=1,n_singles_b ! loop over other columns in this row
      lcol = singles_b(i)

      tmp_det2(1:$N_int,2) = psi_det_beta_unique(1:$N_int, lcol)

      call get_homo_banned_spin(tmp_det2,banned_b2,spin_hp,sign_hp,idx_hp,2,$N_int)
      if (banned_b2) cycle

      l_a = psi_bilinear_matrix_columns_loc(lcol) ! location of start of this column within psi_bilinear_mat
      ASSERT (l_a <= N_det)

      do j=1,psi_bilinear_matrix_columns_loc(lcol+1) - l_a ! loop over rows in this column
        lrow = psi_bilinear_matrix_rows(l_a) ! get row (index of unique alpha det)
        ASSERT (lrow <= N_det_alpha_unique)

        buffer(1:$N_int,j) = psi_det_alpha_unique(1:$N_int, lrow) ! get alpha det

        ASSERT (l_a <= N_det)
        idx(j) = l_a ! indices of dets within psi_bilinear_mat
        l_a = l_a+1
      enddo
      j = j-1
      ! get all alpha dets in this column that are connected to ref alpha by a single exc.
      call get_all_spin_singles_$N_int(                              &
          buffer, idx, tmp_det(1,1), j,                              &
          singles_a, n_singles_a )

      ! Loop over alpha singles
      ! -----------------------

      do k = 1,n_singles_a
        l_a = singles_a(k)
        ASSERT (l_a <= N_det)

        lrow = psi_bilinear_matrix_rows(l_a)
        ASSERT (lrow <= N_det_alpha_unique)

        tmp_det2(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, lrow)
        call get_homo_banned_spin(tmp_det2,banned_a2,spin_hp,sign_hp,idx_hp,1,$N_int)
        if (banned_a2) cycle
        call i_h_j_double_alpha_beta_homo(tmp_det,tmp_det2,$N_int,hij_hp,spin_hp,sign_hp,idx_hp)
        v_t(k_a) = v_t(k_a) + hij_hp * u_t(l_a)
      enddo
    enddo
  enddo
  !$OMP END DO 

  !$OMP DO SCHEDULE(dynamic,64)
  do k_a=istart+ishift,iend,istep


    ! Single and double alpha excitations
    ! ===================================
    
    
    ! Initial determinant is at k_a in alpha-major representation
    ! -----------------------------------------------------------------------
    
    krow = psi_bilinear_matrix_rows(k_a)
    ASSERT (krow <= N_det_alpha_unique)

    kcol = psi_bilinear_matrix_columns(k_a)
    ASSERT (kcol <= N_det_beta_unique)
    
    tmp_det(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, krow)
    tmp_det(1:$N_int,2) = psi_det_beta_unique (1:$N_int, kcol)
    call get_homo_banned_ab(tmp_det,banned_ab1,spin_hp,sign_hp,idx_hp,$N_int)
    if (banned_ab1) cycle

    
    ! Initial determinant is at k_b in beta-major representation
    ! ----------------------------------------------------------------------
    
    k_b = psi_bilinear_matrix_order_transp_reverse(k_a)
    ASSERT (k_b <= N_det)

    spindet(1:$N_int) = tmp_det(1:$N_int,1)
    
    ! Loop inside the beta column to gather all the connected alphas
    lcol = psi_bilinear_matrix_columns(k_a)
    l_a = psi_bilinear_matrix_columns_loc(lcol)
    do i=1,N_det_alpha_unique
      if (l_a > N_det) exit
      lcol = psi_bilinear_matrix_columns(l_a)
      if (lcol /= kcol) exit
      lrow = psi_bilinear_matrix_rows(l_a)
      ASSERT (lrow <= N_det_alpha_unique)

      buffer(1:$N_int,i) = psi_det_alpha_unique(1:$N_int, lrow)
      idx(i) = l_a
      l_a = l_a+1
    enddo
    i = i-1
    
    !call get_all_spin_singles_and_doubles_$N_int(                    &
    !    buffer, idx, spindet, i,                                     &
    !    singles_a, doubles, n_singles_a, n_doubles )
    call get_all_spin_singles_and_doubles(                           &
        buffer, idx, spindet, $N_int, i,                             &
        singles_a, doubles, n_singles_a, n_doubles )

    ! Compute Hij for all alpha singles
    ! ----------------------------------

    tmp_det2(1:$N_int,2) = psi_det_beta_unique (1:$N_int, kcol)
    do i=1,n_singles_a
      l_a = singles_a(i)
      ASSERT (l_a <= N_det)

      lrow = psi_bilinear_matrix_rows(l_a)
      ASSERT (lrow <= N_det_alpha_unique)

      tmp_det2(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, lrow)
      call get_homo_banned_spin(tmp_det2,banned_a2,spin_hp,sign_hp,idx_hp,1,$N_int)
      if (banned_a2) cycle
      call i_h_j_mono_spin_homo(tmp_det,tmp_det2,$N_int,1, hij_hp,spin_hp,sign_hp,idx_hp)

      v_t(k_a) = v_t(k_a) + hij_hp * u_t(l_a)
      ! single => sij = 0 
    enddo
    
    ! Compute Hij for all alpha doubles
    ! ----------------------------------
    
    do i=1,n_doubles
      l_a = doubles(i)
      ASSERT (l_a <= N_det)

      lrow = psi_bilinear_matrix_rows(l_a)
      ASSERT (lrow <= N_det_alpha_unique)

      call get_homo_banned_single_spin(psi_det_alpha_unique(1,lrow),banned_a2,spin_hp,sign_hp,idx_hp,1,$N_int)
      if (banned_a2) cycle
      call i_h_j_double_spin_homo(tmp_det(1,1),psi_det_alpha_unique(1,lrow),$N_int,1,hij_hp,spin_hp,sign_hp,idx_hp)
      v_t(k_a) = v_t(k_a) + hij_hp * u_t(l_a)
      ! same spin => sij = 0
    enddo
    

    ! Single and double beta excitations
    ! ==================================

    
    ! Initial determinant is at k_a in alpha-major representation
    ! -----------------------------------------------------------------------
    
    krow = psi_bilinear_matrix_rows(k_a)
    kcol = psi_bilinear_matrix_columns(k_a)
    
    tmp_det(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, krow)
    tmp_det(1:$N_int,2) = psi_det_beta_unique (1:$N_int, kcol)
    
    !! should already be done from top of loop?
    call get_homo_banned_ab(tmp_det,banned_ab1,spin_hp,sign_hp,idx_hp,$N_int)
    if (banned_ab1) then
      print*,'WARNING: problem with excitation checking in',irp_here
      cycle
    endif
    
    spindet(1:$N_int) = tmp_det(1:$N_int,2)
    
    ! Initial determinant is at k_b in beta-major representation
    ! -----------------------------------------------------------------------

    k_b = psi_bilinear_matrix_order_transp_reverse(k_a) 
    ASSERT (k_b <= N_det)
    
    ! Loop inside the alpha row to gather all the connected betas
    lrow = psi_bilinear_matrix_transp_rows(k_b)
    l_b = psi_bilinear_matrix_transp_rows_loc(lrow)
    do i=1,N_det_beta_unique
      if (l_b > N_det) exit
      lrow = psi_bilinear_matrix_transp_rows(l_b)
      if (lrow /= krow) exit
      lcol = psi_bilinear_matrix_transp_columns(l_b)
      ASSERT (lcol <= N_det_beta_unique)

      buffer(1:$N_int,i) = psi_det_beta_unique(1:$N_int, lcol)
      idx(i) = l_b
      l_b = l_b+1
    enddo
    i = i-1
  
    !call get_all_spin_singles_and_doubles_$N_int(                    &
    !    buffer, idx, spindet, i,                                     &
    !    singles_b, doubles, n_singles_b, n_doubles )
    call get_all_spin_singles_and_doubles(                           &
        buffer, idx, spindet, $N_int, i,                             &
        singles_b, doubles, n_singles_b, n_doubles )
    
    ! Compute Hij for all beta singles
    ! ----------------------------------
    
    tmp_det2(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, krow)
    do i=1,n_singles_b
      l_b = singles_b(i)
      ASSERT (l_b <= N_det)

      lcol = psi_bilinear_matrix_transp_columns(l_b)
      ASSERT (lcol <= N_det_beta_unique)

      tmp_det2(1:$N_int,2) = psi_det_beta_unique (1:$N_int, lcol)
      call get_homo_banned_spin(tmp_det2,banned_b2,spin_hp,sign_hp,idx_hp,2,$N_int)
      if (banned_b2) cycle
      call i_h_j_mono_spin_homo(tmp_det,tmp_det2,$N_int,2, hij_hp,spin_hp,sign_hp,idx_hp)
      l_a = psi_bilinear_matrix_transp_order(l_b)
      ASSERT (l_a <= N_det)
      v_t(k_a) = v_t(k_a) + hij_hp * u_t(l_a)
      ! single => sij = 0 
    enddo
    
    ! Compute Hij for all beta doubles
    ! ----------------------------------
    
    do i=1,n_doubles
      l_b = doubles(i)
      ASSERT (l_b <= N_det)

      lcol = psi_bilinear_matrix_transp_columns(l_b)
      ASSERT (lcol <= N_det_beta_unique)

      call get_homo_banned_single_spin(psi_det_beta_unique(1,lcol),banned_b2,spin_hp,sign_hp,idx_hp,2,$N_int)
      if (banned_b2) cycle
      call i_h_j_double_spin_homo(tmp_det(1,2),psi_det_beta_unique(1,lcol),$N_int,2,hij_hp,spin_hp,sign_hp,idx_hp)
      l_a = psi_bilinear_matrix_transp_order(l_b)
      ASSERT (l_a <= N_det)

      v_t(k_a) = v_t(k_a) + hij_hp * u_t(l_a)
      ! same spin => sij = 0 
    enddo


    ! Diagonal contribution
    ! =====================

    
    ! Initial determinant is at k_a in alpha-major representation
    ! -----------------------------------------------------------------------
    
    krow = psi_bilinear_matrix_rows(k_a)
    ASSERT (krow <= N_det_alpha_unique)

    kcol = psi_bilinear_matrix_columns(k_a)
    ASSERT (kcol <= N_det_beta_unique)
    
    tmp_det(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, krow)
    tmp_det(1:$N_int,2) = psi_det_beta_unique (1:$N_int, kcol)

    call get_homo_banned_ab(tmp_det,banned_ab1,spin_hp,sign_hp,idx_hp,$N_int)
    if (banned_ab1) then
      print*,'WARNING: problem with excitation checking in:',irp_here
      cycle
    endif
    
    double precision, external :: diag_H_mat_elem, diag_S_mat_elem
    hii = diag_h_mat_elem(tmp_det,$N_int)

!    do ii=1,N_hp
!      if(exc_is_banned_ab1(ii)) then
!        hii_hp(ii)=0.d0
!      else
        tmp_det2=tmp_det
        na=elec_num_tab(spin_hp)
        nb=elec_num_tab(iand(spin_hp,1)+1)
        hii_hp=hii
        if (sign_hp>0) then
          call ac_operator(idx_hp,spin_hp,tmp_det2,hii_hp,$N_int,na,nb)
        else
          call a_operator(idx_hp,spin_hp,tmp_det2,hii_hp,$N_int,na,nb)
        endif
 !     endif
      v_t(k_a) = v_t(k_a) + hii_hp * u_t(k_a)
!    enddo
    

  end do
  !$OMP END DO 
  deallocate(buffer, singles_a, singles_b, doubles, idx)
  !$OMP END PARALLEL
  deallocate(idx0)
end

SUBST [ N_int ]

1;;
2;;
3;;
4;;
N_int;;

END_TEMPLATE



subroutine i_h_j_double_spin_homo(key_i,key_j,Nint,ispin,hij_hp,spin_hp,sign_hp,idx_hp)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! todo: maybe make new get_double_excitation_spin?
  !       the 4 index ordering is already done in there, so we could avoid duplicating that work
  ! Returns <i|H|j> where i and j are determinants differing by a same-spin double excitation
  END_DOC
  integer, intent(in)            :: Nint,ispin
  integer(bit_kind), intent(in)  :: key_i(Nint), key_j(Nint)
  complex*16, intent(out)  :: hij_hp
  integer, intent(in) :: spin_hp, idx_hp
  double precision, intent(in) :: sign_hp
  complex*16 :: hij0 
  double precision :: phase_hp
  integer                        :: exc(0:2,2)
  double precision               :: phase
  complex*16, external     :: get_mo_bielec_integral
  integer :: i1,i2,i3,i4,j2,j3,ii

  PROVIDE big_array_exchange_integrals mo_bielec_integrals_in_map

  call get_double_excitation_spin(key_i,key_j,exc,phase,Nint)
  hij0 = phase*(get_mo_bielec_integral(                               &
      exc(1,1),                                                      &
      exc(2,1),                                                      &
      exc(1,2),                                                      &
      exc(2,2), mo_integrals_map) -                                  &
      get_mo_bielec_integral(                                        &
      exc(1,1),                                                      &
      exc(2,1),                                                      &
      exc(2,2),                                                      &
      exc(1,2), mo_integrals_map) )

  ASSERT (exc(1,1) < exc(2,1))
  ASSERT (exc(1,2) < exc(2,2))
  i1=min(exc(1,1),exc(1,2))
  j2=max(exc(1,1),exc(1,2))
  j3=min(exc(2,1),exc(2,2))
  i4=max(exc(2,1),exc(2,2))
  i2=min(j2,j3)
  i3=max(j2,j3)

      if (ispin.eq.spin_hp) then
        if ((idx_hp.lt.i1).or.(idx_hp.gt.i4)) then
          phase_hp=1.d0
        else if ((idx_hp.lt.i2).or.(idx_hp.gt.i3)) then
          phase_hp=-1.d0
        else
          phase_hp=1.d0
        endif
      else
        phase_hp=1.d0
      endif
    hij_hp = hij0 * phase_hp
end

subroutine i_h_j_mono_spin_homo(key_i,key_j,Nint,spin,hij_hp,spin_hp,sign_hp,idx_hp)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! todo: change this to use normal version of get_mono_excitation_from_fock
  !       all info needed is in phase and hij, h/p part can happen after getting hij the normal way
  ! Returns <i|H|j> where i and j are determinants differing by a single excitation
  END_DOC
  integer, intent(in)            :: Nint, spin
  integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
  complex*16, intent(out)  :: hij_hp
  integer, intent(in) :: spin_hp, idx_hp
  double precision, intent(in) :: sign_hp
  !double precision :: phase_hp(N_hp)
  complex*16 :: hij0
  
  integer                        :: exc(0:2,2)
  double precision               :: phase

  PROVIDE big_array_exchange_integrals mo_bielec_integrals_in_map
 
  call get_mono_excitation_spin(key_i(1,spin),key_j(1,spin),exc,phase,Nint)

  call get_mono_excitation_from_fock_homo(key_i,key_j,exc(1,1),exc(1,2),spin,phase,hij_hp,spin_hp,sign_hp,idx_hp)
end

subroutine get_mono_excitation_from_fock_homo(det_1,det_2,h,p,spin,phase,hij_hp,spin_hp,sign_hp,idx_hp)
  use bitmasks
  implicit none
  integer,intent(in) :: h,p,spin
  double precision, intent(in)  :: phase
  integer(bit_kind), intent(in) :: det_1(N_int,2), det_2(N_int,2)
  complex*16, intent(out) :: hij_hp
  integer, intent(in) :: spin_hp, idx_hp
  double precision, intent(in) :: sign_hp
  double precision :: phase_hp
  complex*16 :: hij0
  integer :: low,high

  integer(bit_kind) :: differences(N_int,2)
  integer(bit_kind) :: hole(N_int,2)
  integer(bit_kind) :: partcl(N_int,2)
  integer :: occ_hole(N_int*bit_kind_size,2)
  integer :: occ_partcl(N_int*bit_kind_size,2)
  integer :: n_occ_ab_hole(2),n_occ_ab_partcl(2)
  integer :: i0,i,ii
  do i = 1, N_int
    differences(i,1) = xor(det_1(i,1),ref_closed_shell_bitmask(i,1))
    differences(i,2) = xor(det_1(i,2),ref_closed_shell_bitmask(i,2))
    hole(i,1) = iand(differences(i,1),ref_closed_shell_bitmask(i,1))
    hole(i,2) = iand(differences(i,2),ref_closed_shell_bitmask(i,2))
    partcl(i,1) = iand(differences(i,1),det_1(i,1))
    partcl(i,2) = iand(differences(i,2),det_1(i,2))
  enddo
  call bitstring_to_list_ab(hole, occ_hole, n_occ_ab_hole, N_int)
  call bitstring_to_list_ab(partcl, occ_partcl, n_occ_ab_partcl, N_int)
  hij0 = fock_operator_closed_shell_ref_bitmask(h,p)
  ! holes :: direct terms
  do i0 = 1, n_occ_ab_hole(1)
    i = occ_hole(i0,1)
    hij0 -= big_array_coulomb_integrals(i,h,p) ! get_mo_bielec_integral_schwartz(h,i,p,i,mo_integrals_map)
  enddo
  do i0 = 1, n_occ_ab_hole(2)
    i = occ_hole(i0,2)
    hij0 -= big_array_coulomb_integrals(i,h,p) !get_mo_bielec_integral_schwartz(h,i,p,i,mo_integrals_map)
  enddo
 
  ! holes :: exchange terms
  do i0 = 1, n_occ_ab_hole(spin)
    i = occ_hole(i0,spin)
    hij0 += big_array_exchange_integrals(i,h,p) ! get_mo_bielec_integral_schwartz(h,i,i,p,mo_integrals_map)
  enddo
 
  ! particles :: direct terms
  do i0 = 1, n_occ_ab_partcl(1)
    i = occ_partcl(i0,1)
    hij0 += big_array_coulomb_integrals(i,h,p)!get_mo_bielec_integral_schwartz(h,i,p,i,mo_integrals_map)
  enddo
  do i0 = 1, n_occ_ab_partcl(2)
    i = occ_partcl(i0,2)
    hij0 += big_array_coulomb_integrals(i,h,p) !get_mo_bielec_integral_schwartz(h,i,p,i,mo_integrals_map)
  enddo
 
  ! particles :: exchange terms
  do i0 = 1, n_occ_ab_partcl(spin)
    i = occ_partcl(i0,spin)
    hij0 -= big_array_exchange_integrals(i,h,p)!get_mo_bielec_integral_schwartz(h,i,i,p,mo_integrals_map)
  enddo

  low=min(h,p)
  high=max(h,p)

!!  do ii=1,N_hp
!!    if (.not.allowed_hp(ii)) then
!!      phase_hp(ii) = 0.d0
!!      cycle
!!    else if (spin_hp(ii).ne.spin) then
!!      phase_hp(ii) = 1.d0
!!    else
!!      if ((low.lt.idx_hp(ii)).and.(high.gt.idx_hp(ii))) then
!!        phase_hp(ii) = -1.d0
!!      else
!!        phase_hp(ii) = 1.d0
!!      endif
!!    endif
!!  enddo
!!
!!  do ii=1,N_hp
!!    if (allowed_hp(ii)) then
!!      hij_hp(ii) = hij + sign_hp(ii) * big_array_coulomb_integrals(idx_hp(ii),h,p)
!!      if (spin.eq.spin_hp(ii)) then
!!        hij_hp(ii) = hij_hp(ii) - sign_hp(ii) * big_array_exchange_integrals(idx_hp(ii),h,p)
!!      endif
!!    else
!!      hij_hp(ii) = 0.d0
!!    endif
!!  enddo
!!
!!  do ii=1,N_hp
!!    hij_hp(ii) = hij_hp(ii) * phase_hp(ii) * phase
!!  enddo

!  do ii=1,N_hp
!    if (.not.allowed_hp(ii)) then
!      phase_hp(ii) = 0.d0
!      hij_hp(ii) = 0.d0
!      cycle
    if (spin.eq.spin_hp) then
      hij_hp = hij0 + sign_hp *(big_array_coulomb_integrals(idx_hp,h,p) - big_array_exchange_integrals(idx_hp,h,p))
      if ((low.lt.idx_hp).and.(high.gt.idx_hp)) then
        phase_hp = -1.d0
      else
        phase_hp = 1.d0
      endif
    else
      phase_hp = 1.d0
      hij_hp = hij0 + sign_hp * big_array_coulomb_integrals(idx_hp,h,p)
    endif
    hij_hp = hij_hp * phase * phase_hp
!  enddo

end


subroutine i_h_j_double_alpha_beta_homo(key_i,key_j,Nint,hij_hp,spin_hp,sign_hp,idx_hp)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns <i|H|j> where i and j are determinants differing by an opposite-spin double excitation
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
  complex*16, intent(out)  :: hij_hp
  complex*16 :: hij0
  integer, intent(in) :: spin_hp, idx_hp
  double precision, intent(in) :: sign_hp
  double precision :: phase_hp
  integer :: i

  integer :: lowhigh(2,2)
  integer                        :: exc(0:2,2,2)
  double precision               :: phase, phase2
  complex*16, external     :: get_mo_bielec_integral

  PROVIDE big_array_exchange_integrals mo_bielec_integrals_in_map

  call get_mono_excitation_spin(key_i(1,1),key_j(1,1),exc(0,1,1),phase,Nint)
  call get_mono_excitation_spin(key_i(1,2),key_j(1,2),exc(0,1,2),phase2,Nint)
  phase = phase*phase2

  if (exc(1,1,1) == exc(1,2,2)) then
    hij0 =  big_array_exchange_integrals(exc(1,1,1),exc(1,1,2),exc(1,2,1))
  else if (exc(1,2,1) == exc(1,1,2)) then
    hij0 =  big_array_exchange_integrals(exc(1,2,1),exc(1,1,1),exc(1,2,2))
  else
    hij0 = get_mo_bielec_integral(                              &
        exc(1,1,1),                                                  &
        exc(1,1,2),                                                  &
        exc(1,2,1),                                                  &
        exc(1,2,2) ,mo_integrals_map)
  endif
  
  !todo: clean this up
  ! if new particle/hole is between p/h of single exc of same spin, then parity changes, otherwise stays the same
  ! value of Hij for double excitation is unchanged (new p/h is not one of the indices involved in the excitation)
  
  lowhigh(1,1)=min(exc(1,1,1),exc(1,2,1))
  lowhigh(2,1)=max(exc(1,1,1),exc(1,2,1))
  lowhigh(1,2)=min(exc(1,1,2),exc(1,2,2))
  lowhigh(2,2)=max(exc(1,1,2),exc(1,2,2))
!  do i=1,N_hp
!    if (.not.allowed_hp(i)) then
!      phase_hp(i)=0.d0
    if ((idx_hp.gt.lowhigh(1,spin_hp)).and.(idx_hp.lt.lowhigh(2,spin_hp))) then
      phase_hp=-1.d0
    else
      phase_hp=1.d0
    endif
    hij_hp=hij0*phase*phase_hp
!  enddo
end
