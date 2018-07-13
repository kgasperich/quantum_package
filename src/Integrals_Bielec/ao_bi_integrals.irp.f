double precision function ao_bielec_integral(i,j,k,l)
  implicit none
  BEGIN_DOC
  !  integral of the AO basis <ik|jl> or (ij|kl)
  !     i(r1) j(r1) 1/r12 k(r2) l(r2)
  END_DOC

  integer,intent(in)             :: i,j,k,l
  integer                        :: p,q,r,s
  double precision               :: I_center(3),J_center(3),K_center(3),L_center(3)
  integer                        :: num_i,num_j,num_k,num_l,dim1,I_power(3),J_power(3),K_power(3),L_power(3)
  double precision               :: integral
  include 'Utils/constants.include.F'
  double precision               :: P_new(0:max_dim,3),P_center(3),fact_p,pp
  double precision               :: Q_new(0:max_dim,3),Q_center(3),fact_q,qq
  integer                        :: iorder_p(3), iorder_q(3)
  double precision               :: ao_bielec_integral_schwartz_accel
  
   if (ao_prim_num(i) * ao_prim_num(j) * ao_prim_num(k) * ao_prim_num(l) > 1024 ) then
     ao_bielec_integral = ao_bielec_integral_schwartz_accel(i,j,k,l)
     return
   endif

  dim1 = n_pt_max_integrals
  
  num_i = ao_nucl(i)
  num_j = ao_nucl(j)
  num_k = ao_nucl(k)
  num_l = ao_nucl(l)
  ao_bielec_integral = 0.d0
  
  if (num_i /= num_j .or. num_k /= num_l .or. num_j /= num_k)then
    do p = 1, 3
      I_power(p) = ao_power(i,p)
      J_power(p) = ao_power(j,p)
      K_power(p) = ao_power(k,p)
      L_power(p) = ao_power(l,p)
      I_center(p) = nucl_coord(num_i,p)
      J_center(p) = nucl_coord(num_j,p)
      K_center(p) = nucl_coord(num_k,p)
      L_center(p) = nucl_coord(num_l,p)
    enddo
    
    double precision               :: coef1, coef2, coef3, coef4
    double precision               :: p_inv,q_inv
    double precision               :: general_primitive_integral

    do p = 1, ao_prim_num(i)
      coef1 = ao_coef_normalized_ordered_transp(p,i)
      do q = 1, ao_prim_num(j)
        coef2 = coef1*ao_coef_normalized_ordered_transp(q,j)
        call give_explicit_poly_and_gaussian(P_new,P_center,pp,fact_p,iorder_p,&
            ao_expo_ordered_transp(p,i),ao_expo_ordered_transp(q,j),                 &
            I_power,J_power,I_center,J_center,dim1)
        p_inv = 1.d0/pp
        do r = 1, ao_prim_num(k)
          coef3 = coef2*ao_coef_normalized_ordered_transp(r,k)
          do s = 1, ao_prim_num(l)
            coef4 = coef3*ao_coef_normalized_ordered_transp(s,l)
            call give_explicit_poly_and_gaussian(Q_new,Q_center,qq,fact_q,iorder_q,&
                ao_expo_ordered_transp(r,k),ao_expo_ordered_transp(s,l),             &
                K_power,L_power,K_center,L_center,dim1)
            q_inv = 1.d0/qq
            integral = general_primitive_integral(dim1,              &
                P_new,P_center,fact_p,pp,p_inv,iorder_p,             &
                Q_new,Q_center,fact_q,qq,q_inv,iorder_q)
            ao_bielec_integral = ao_bielec_integral +  coef4 * integral
          enddo ! s
        enddo  ! r
      enddo   ! q
    enddo    ! p
    
  else
    
    do p = 1, 3
      I_power(p) = ao_power(i,p)
      J_power(p) = ao_power(j,p)
      K_power(p) = ao_power(k,p)
      L_power(p) = ao_power(l,p)
    enddo
    double  precision              :: ERI

    do p = 1, ao_prim_num(i)
      coef1 = ao_coef_normalized_ordered_transp(p,i)
      do q = 1, ao_prim_num(j)
        coef2 = coef1*ao_coef_normalized_ordered_transp(q,j)
        do r = 1, ao_prim_num(k)
          coef3 = coef2*ao_coef_normalized_ordered_transp(r,k)
          do s = 1, ao_prim_num(l)
            coef4 = coef3*ao_coef_normalized_ordered_transp(s,l)
            integral = ERI(                                          &
                ao_expo_ordered_transp(p,i),ao_expo_ordered_transp(q,j),ao_expo_ordered_transp(r,k),ao_expo_ordered_transp(s,l),&
                I_power(1),J_power(1),K_power(1),L_power(1),         &
                I_power(2),J_power(2),K_power(2),L_power(2),         &
                I_power(3),J_power(3),K_power(3),L_power(3))
            ao_bielec_integral = ao_bielec_integral + coef4 * integral
          enddo ! s
        enddo  ! r
      enddo   ! q
    enddo    ! p
    
  endif
  
end

double precision function ao_bielec_integral_schwartz_accel(i,j,k,l)
  implicit none
  BEGIN_DOC
  !  integral of the AO basis <ik|jl> or (ij|kl)
  !     i(r1) j(r1) 1/r12 k(r2) l(r2)
  END_DOC
  integer,intent(in)             :: i,j,k,l
  integer                        :: p,q,r,s
  double precision               :: I_center(3),J_center(3),K_center(3),L_center(3)
  integer                        :: num_i,num_j,num_k,num_l,dim1,I_power(3),J_power(3),K_power(3),L_power(3)
  double precision               :: integral
  include 'Utils/constants.include.F'
  double precision               :: P_new(0:max_dim,3),P_center(3),fact_p,pp
  double precision               :: Q_new(0:max_dim,3),Q_center(3),fact_q,qq
  integer                        :: iorder_p(3), iorder_q(3)
  double precision, allocatable  :: schwartz_kl(:,:)
  double precision               :: schwartz_ij
  
  dim1 = n_pt_max_integrals
  
  num_i = ao_nucl(i)
  num_j = ao_nucl(j)
  num_k = ao_nucl(k)
  num_l = ao_nucl(l)
  ao_bielec_integral_schwartz_accel = 0.d0
  double precision               :: thr
  thr = ao_integrals_threshold*ao_integrals_threshold
  
  allocate(schwartz_kl(0:ao_prim_num(l),0:ao_prim_num(k)))


  if (num_i /= num_j .or. num_k /= num_l .or. num_j /= num_k)then
    do p = 1, 3
      I_power(p) = ao_power(i,p)
      J_power(p) = ao_power(j,p)
      K_power(p) = ao_power(k,p)
      L_power(p) = ao_power(l,p)
      I_center(p) = nucl_coord(num_i,p)
      J_center(p) = nucl_coord(num_j,p)
      K_center(p) = nucl_coord(num_k,p)
      L_center(p) = nucl_coord(num_l,p)
    enddo
    
    schwartz_kl(0,0) = 0.d0
    do r = 1, ao_prim_num(k)
      coef1 = ao_coef_normalized_ordered_transp(r,k)*ao_coef_normalized_ordered_transp(r,k)
      schwartz_kl(0,r) = 0.d0
      do s = 1, ao_prim_num(l)
        coef2 = coef1 * ao_coef_normalized_ordered_transp(s,l) * ao_coef_normalized_ordered_transp(s,l)
        call give_explicit_poly_and_gaussian(Q_new,Q_center,qq,fact_q,iorder_q,&
            ao_expo_ordered_transp(r,k),ao_expo_ordered_transp(s,l),                 &
            K_power,L_power,K_center,L_center,dim1)
        q_inv = 1.d0/qq
        schwartz_kl(s,r) = general_primitive_integral(dim1,          &
            Q_new,Q_center,fact_q,qq,q_inv,iorder_q,                 &
            Q_new,Q_center,fact_q,qq,q_inv,iorder_q)                 &
            * coef2 
        schwartz_kl(0,r) = max(schwartz_kl(0,r),schwartz_kl(s,r))
      enddo
      schwartz_kl(0,0) = max(schwartz_kl(0,r),schwartz_kl(0,0))
    enddo

    do p = 1, ao_prim_num(i)
      double precision               :: coef1
      coef1 = ao_coef_normalized_ordered_transp(p,i)
      do q = 1, ao_prim_num(j)
        double precision               :: coef2
        coef2 = coef1*ao_coef_normalized_ordered_transp(q,j)
        double precision               :: p_inv,q_inv
        call give_explicit_poly_and_gaussian(P_new,P_center,pp,fact_p,iorder_p,&
            ao_expo_ordered_transp(p,i),ao_expo_ordered_transp(q,j),                 &
            I_power,J_power,I_center,J_center,dim1)
        p_inv = 1.d0/pp
        schwartz_ij = general_primitive_integral(dim1,               &
            P_new,P_center,fact_p,pp,p_inv,iorder_p,                 &
            P_new,P_center,fact_p,pp,p_inv,iorder_p) *               &
            coef2*coef2
        if (schwartz_kl(0,0)*schwartz_ij < thr) then
           cycle
        endif
        do r = 1, ao_prim_num(k)
          if (schwartz_kl(0,r)*schwartz_ij < thr) then
             cycle
          endif
          double precision               :: coef3
          coef3 = coef2*ao_coef_normalized_ordered_transp(r,k)
          do s = 1, ao_prim_num(l)
            double precision               :: coef4
            if (schwartz_kl(s,r)*schwartz_ij < thr) then
               cycle
            endif
            coef4 = coef3*ao_coef_normalized_ordered_transp(s,l)
            double precision               :: general_primitive_integral
            call give_explicit_poly_and_gaussian(Q_new,Q_center,qq,fact_q,iorder_q,&
                ao_expo_ordered_transp(r,k),ao_expo_ordered_transp(s,l),             &
                K_power,L_power,K_center,L_center,dim1)
            q_inv = 1.d0/qq
            integral = general_primitive_integral(dim1,              &
                P_new,P_center,fact_p,pp,p_inv,iorder_p,             &
                Q_new,Q_center,fact_q,qq,q_inv,iorder_q)
            ao_bielec_integral_schwartz_accel = ao_bielec_integral_schwartz_accel + coef4 * integral
          enddo ! s
        enddo  ! r
      enddo   ! q
    enddo    ! p
    
  else
    
    do p = 1, 3
      I_power(p) = ao_power(i,p)
      J_power(p) = ao_power(j,p)
      K_power(p) = ao_power(k,p)
      L_power(p) = ao_power(l,p)
    enddo
    double  precision              :: ERI

    schwartz_kl(0,0) = 0.d0
    do r = 1, ao_prim_num(k)
      coef1 = ao_coef_normalized_ordered_transp(r,k)*ao_coef_normalized_ordered_transp(r,k)
      schwartz_kl(0,r) = 0.d0
      do s = 1, ao_prim_num(l)
        coef2 = coef1*ao_coef_normalized_ordered_transp(s,l)*ao_coef_normalized_ordered_transp(s,l)
        schwartz_kl(s,r) = ERI(                                      &
            ao_expo_ordered_transp(r,k),ao_expo_ordered_transp(s,l),ao_expo_ordered_transp(r,k),ao_expo_ordered_transp(s,l),&
            K_power(1),L_power(1),K_power(1),L_power(1),             &
            K_power(2),L_power(2),K_power(2),L_power(2),             &
            K_power(3),L_power(3),K_power(3),L_power(3)) * &
            coef2
        schwartz_kl(0,r) = max(schwartz_kl(0,r),schwartz_kl(s,r))
      enddo
      schwartz_kl(0,0) = max(schwartz_kl(0,r),schwartz_kl(0,0))
    enddo

    do p = 1, ao_prim_num(i)
      coef1 = ao_coef_normalized_ordered_transp(p,i)
      do q = 1, ao_prim_num(j)
        coef2 = coef1*ao_coef_normalized_ordered_transp(q,j)
        schwartz_ij = ERI(                                          &
                ao_expo_ordered_transp(p,i),ao_expo_ordered_transp(q,j),ao_expo_ordered_transp(p,i),ao_expo_ordered_transp(q,j),&
                I_power(1),J_power(1),I_power(1),J_power(1),         &
                I_power(2),J_power(2),I_power(2),J_power(2),         &
                I_power(3),J_power(3),I_power(3),J_power(3))*coef2*coef2
        if (schwartz_kl(0,0)*schwartz_ij < thr) then
           cycle
        endif
        do r = 1, ao_prim_num(k)
          if (schwartz_kl(0,r)*schwartz_ij < thr) then
             cycle
          endif
          coef3 = coef2*ao_coef_normalized_ordered_transp(r,k)
          do s = 1, ao_prim_num(l)
            if (schwartz_kl(s,r)*schwartz_ij < thr) then
               cycle
            endif
            coef4 = coef3*ao_coef_normalized_ordered_transp(s,l)
            integral = ERI(                                          &
                ao_expo_ordered_transp(p,i),ao_expo_ordered_transp(q,j),ao_expo_ordered_transp(r,k),ao_expo_ordered_transp(s,l),&
                I_power(1),J_power(1),K_power(1),L_power(1),         &
                I_power(2),J_power(2),K_power(2),L_power(2),         &
                I_power(3),J_power(3),K_power(3),L_power(3))
            ao_bielec_integral_schwartz_accel = ao_bielec_integral_schwartz_accel +  coef4 * integral
          enddo ! s
        enddo  ! r
      enddo   ! q
    enddo    ! p
    
  endif
  deallocate (schwartz_kl)
  
end


integer function ao_l4(i,j,k,l)
  implicit none
  BEGIN_DOC
! Computes the product of l values of i,j,k,and l
  END_DOC
  integer, intent(in) :: i,j,k,l
  ao_l4 = ao_l(i)*ao_l(j)*ao_l(k)*ao_l(l)
end



subroutine compute_ao_bielec_integrals(j,k,l,sze,buffer_value)
  implicit none
  use map_module
  
  BEGIN_DOC
  ! Compute AO 1/r12 integrals for all i and fixed j,k,l
  END_DOC
  
  include 'Utils/constants.include.F'
  integer, intent(in)            :: j,k,l,sze
  real(integral_kind), intent(out) :: buffer_value(sze)
  double precision               :: ao_bielec_integral
  
  integer                        :: i
  
  if (ao_overlap_abs(j,l) < thresh) then
    buffer_value = 0._integral_kind
    return
  endif
  if (ao_bielec_integral_schwartz(j,l) < thresh ) then
    buffer_value = 0._integral_kind
    return
  endif
  
  do i = 1, ao_num
    if (ao_overlap_abs(i,k)*ao_overlap_abs(j,l) < thresh) then
      buffer_value(i) = 0._integral_kind
      cycle
    endif
    if (ao_bielec_integral_schwartz(i,k)*ao_bielec_integral_schwartz(j,l) < thresh ) then
      buffer_value(i) = 0._integral_kind
      cycle
    endif
    !DIR$ FORCEINLINE
    buffer_value(i) = ao_bielec_integral(i,k,j,l)
  enddo
  
end

BEGIN_PROVIDER [ logical, ao_bielec_integrals_in_map ]
  implicit none
  use f77_zmq
  use map_module
  BEGIN_DOC
  !  Map of Atomic integrals
  !     i(r1) j(r2) 1/r12 k(r1) l(r2)
  END_DOC
  
  PROVIDE read_ao_integrals disk_access_ao_integrals
  if (read_ao_integrals) then
    print*,'Reading the AO integrals'
      call map_load_from_disk(trim(ezfio_filename)//'/work/ao_ints',ao_integrals_map)
      print*, 'AO integrals provided'
      ao_bielec_integrals_in_map = .True.
      return
  endif
  
END_PROVIDER
 
BEGIN_PROVIDER [ double precision, ao_bielec_integral_schwartz,(ao_num,ao_num)  ]
  implicit none
  BEGIN_DOC
  !  Needed to compute Schwartz inequalities
  END_DOC
  
  integer                        :: i,k
  double precision               :: ao_bielec_integral,cpu_1,cpu_2, wall_1, wall_2
  
  ao_bielec_integral_schwartz(1,1) = ao_bielec_integral(1,1,1,1)
  !$OMP PARALLEL DO PRIVATE(i,k)                                     &
      !$OMP DEFAULT(NONE)                                            &
      !$OMP SHARED (ao_num,ao_bielec_integral_schwartz)              &
      !$OMP SCHEDULE(dynamic)
  do i=1,ao_num
    do k=1,i
      ao_bielec_integral_schwartz(i,k) = dsqrt(ao_bielec_integral(i,k,i,k))
      ao_bielec_integral_schwartz(k,i) = ao_bielec_integral_schwartz(i,k)
    enddo
  enddo
  !$OMP END PARALLEL DO
  
END_PROVIDER



subroutine compute_ao_integrals_jl(j,l,n_integrals,buffer_i,buffer_value)
  implicit none
  use map_module
  BEGIN_DOC
  !  Parallel client for AO integrals
  END_DOC
  
  integer, intent(in)            :: j,l
  integer,intent(out)            :: n_integrals
  integer(key_kind),intent(out)  :: buffer_i(ao_num*ao_num)
  real(integral_kind),intent(out) :: buffer_value(ao_num*ao_num)

  integer                        :: i,k
  double precision               :: ao_bielec_integral,cpu_1,cpu_2, wall_1, wall_2
  double precision               :: integral, wall_0
  double precision               :: thr
  integer                        :: kk, m, j1, i1

  thr = ao_integrals_threshold
  
  n_integrals = 0
  
  j1 = j+ishft(l*l-l,-1)
  do k = 1, ao_num           ! r1
    i1 = ishft(k*k-k,-1)
    if (i1 > j1) then
      exit
    endif
    do i = 1, k
      i1 += 1
      if (i1 > j1) then
        exit
      endif
      if (ao_overlap_abs(i,k)*ao_overlap_abs(j,l) < thr) then
        cycle
      endif
      if (ao_bielec_integral_schwartz(i,k)*ao_bielec_integral_schwartz(j,l) < thr ) then
        cycle
      endif
      !DIR$ FORCEINLINE
      integral = ao_bielec_integral(i,k,j,l)  ! i,k : r1    j,l : r2
      if (abs(integral) < thr) then
        cycle
      endif
      n_integrals += 1
      !DIR$ FORCEINLINE
      call bielec_integrals_index(i,j,k,l,buffer_i(n_integrals))
      buffer_value(n_integrals) = integral
    enddo
  enddo
    
end
