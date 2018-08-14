 BEGIN_PROVIDER [ complex*16, Fock_matrix_mo, (mo_tot_num,mo_tot_num) ]
&BEGIN_PROVIDER [ double precision, Fock_matrix_diag_mo, (mo_tot_num)]
   implicit none
   BEGIN_DOC
   ! Fock matrix on the MO basis.
   ! For open shells, the ROHF Fock Matrix is
   !
   !  |   F-K    |  F + K/2  |    F     |
   !  |---------------------------------|
   !  | F + K/2  |     F     |  F - K/2 |
   !  |---------------------------------|
   !  |    F     |  F - K/2  |  F + K   |
   !
   ! F = 1/2 (Fa + Fb)
   !
   ! K = Fb - Fa
   !
   END_DOC
   integer                        :: i,j,n
   if (elec_alpha_num == elec_beta_num) then
     Fock_matrix_mo = Fock_matrix_mo_alpha
   else
     
     do j=1,elec_beta_num
       ! F-K
       do i=1,elec_beta_num
         Fock_matrix_mo(i,j) = 0.5d0*(Fock_matrix_mo_alpha(i,j)+Fock_matrix_mo_beta(i,j))&
             - (Fock_matrix_mo_beta(i,j) - Fock_matrix_mo_alpha(i,j))
       enddo
       ! F+K/2
       do i=elec_beta_num+1,elec_alpha_num
         Fock_matrix_mo(i,j) = 0.5d0*(Fock_matrix_mo_alpha(i,j)+Fock_matrix_mo_beta(i,j))&
             + 0.5d0*(Fock_matrix_mo_beta(i,j) - Fock_matrix_mo_alpha(i,j))
       enddo
       ! F
       do i=elec_alpha_num+1, mo_tot_num
         Fock_matrix_mo(i,j) = 0.5d0*(Fock_matrix_mo_alpha(i,j)+Fock_matrix_mo_beta(i,j))
       enddo
     enddo

     do j=elec_beta_num+1,elec_alpha_num
       ! F+K/2
       do i=1,elec_beta_num
         Fock_matrix_mo(i,j) = 0.5d0*(Fock_matrix_mo_alpha(i,j)+Fock_matrix_mo_beta(i,j))&
             + 0.5d0*(Fock_matrix_mo_beta(i,j) - Fock_matrix_mo_alpha(i,j))
       enddo
       ! F
       do i=elec_beta_num+1,elec_alpha_num
         Fock_matrix_mo(i,j) = 0.5d0*(Fock_matrix_mo_alpha(i,j)+Fock_matrix_mo_beta(i,j))
       enddo
       ! F-K/2
       do i=elec_alpha_num+1, mo_tot_num
         Fock_matrix_mo(i,j) = 0.5d0*(Fock_matrix_mo_alpha(i,j)+Fock_matrix_mo_beta(i,j))&
             - 0.5d0*(Fock_matrix_mo_beta(i,j) - Fock_matrix_mo_alpha(i,j))
       enddo
     enddo

     do j=elec_alpha_num+1, mo_tot_num
       ! F
       do i=1,elec_beta_num
         Fock_matrix_mo(i,j) = 0.5d0*(Fock_matrix_mo_alpha(i,j)+Fock_matrix_mo_beta(i,j))
       enddo
       ! F-K/2
       do i=elec_beta_num+1,elec_alpha_num
         Fock_matrix_mo(i,j) = 0.5d0*(Fock_matrix_mo_alpha(i,j)+Fock_matrix_mo_beta(i,j))&
             - 0.5d0*(Fock_matrix_mo_beta(i,j) - Fock_matrix_mo_alpha(i,j))
       enddo
       ! F+K
       do i=elec_alpha_num+1,mo_tot_num
         Fock_matrix_mo(i,j) = 0.5d0*(Fock_matrix_mo_alpha(i,j)+Fock_matrix_mo_beta(i,j)) &
             + (Fock_matrix_mo_beta(i,j) - Fock_matrix_mo_alpha(i,j))
       enddo
     enddo
     
   endif

   do i = 1, mo_tot_num
     Fock_matrix_diag_mo(i) = real(Fock_matrix_mo(i,i))
     if (dabs(imag(Fock_matrix_mo(i,i))) .gt. 1.0d-12) then
       stop 'diagonal elements of Fock matrix should be real'
     endif
   enddo
END_PROVIDER
 
 
 
 BEGIN_PROVIDER [ complex*16, Fock_matrix_ao_alpha, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ complex*16, Fock_matrix_ao_beta,  (ao_num, ao_num) ]
 implicit none
 BEGIN_DOC
 ! Alpha Fock matrix in AO basis set
 END_DOC
 
 integer                        :: i,j
 do j=1,ao_num
   do i=1,ao_num
     Fock_matrix_ao_alpha(i,j) = ao_mono_elec_integral(i,j) + ao_bi_elec_integral_alpha(i,j)
     Fock_matrix_ao_beta (i,j) = ao_mono_elec_integral(i,j) + ao_bi_elec_integral_beta (i,j)
   enddo
 enddo

END_PROVIDER


 BEGIN_PROVIDER [ complex*16, ao_bi_elec_integral_alpha, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ complex*16, ao_bi_elec_integral_beta ,  (ao_num, ao_num) ]
 use map_module
 implicit none
 BEGIN_DOC
 ! Alpha Fock matrix in AO basis set
 END_DOC
 
 integer                        :: i,j,k,l,k1,r,s
 integer                        :: i0,j0,k0,l0
 integer                        :: ii(4), jj(4), kk(4), ll(4), k2
 integer*8                      :: p,q
 complex*16                     :: integral, c0, c1, c2 
 double precision               :: local_threshold
 complex*16                     :: ao_bielec_integral
 complex*16, allocatable        :: ao_bi_elec_integral_alpha_tmp(:,:)
 complex*16, allocatable        :: ao_bi_elec_integral_beta_tmp(:,:)

 ao_bi_elec_integral_alpha = (0.d0,0.d0)
 ao_bi_elec_integral_beta  = (0.d0,0.d0)

   PROVIDE ao_bielec_integrals_in_map 
           
   integer(omp_lock_kind) :: lck(ao_num)
   integer*8                      :: i8
   integer(cache_map_size_kind)   :: n_elements_max, n_elements
   integer(key_kind), allocatable :: keys(:)
   double precision, allocatable  :: values(:)
   complex*16, parameter    :: i_sign(4) = (/(0.d0,-1.d0),(0.d0,-1.d0),(0.d0,1.d0),(0.d0,1.d0)/)
   logical :: imag_part
   integer(key_kind) :: pair_key

   !$OMP PARALLEL DEFAULT(NONE)                                      &
       !$OMP PRIVATE(i,j,l,k1,k,integral,ii,jj,kk,ll,i8,keys,values,n_elements_max, &
       !$OMP  n_elements,ao_bi_elec_integral_alpha_tmp,ao_bi_elec_integral_beta_tmp,pair_key, &
       !$OMP  c0,c1,c2)&
       !$OMP SHARED(ao_num,HF_density_matrix_ao_alpha,HF_density_matrix_ao_beta,i_sign, &
       !$OMP  ao_integrals_map, ao_bi_elec_integral_alpha, ao_bi_elec_integral_beta) 

   call get_cache_map_n_elements_max(ao_integrals_map,n_elements_max)
   allocate(keys(n_elements_max), values(n_elements_max))
   allocate(ao_bi_elec_integral_alpha_tmp(ao_num,ao_num), &
            ao_bi_elec_integral_beta_tmp(ao_num,ao_num))
   ao_bi_elec_integral_alpha_tmp = (0.d0,0.d0)
   ao_bi_elec_integral_beta_tmp  = (0.d0,0.d0)

   !$OMP DO SCHEDULE(dynamic,64)
   !DIR$ NOVECTOR
    do i8=0_8,ao_integrals_map%map_size
      n_elements = n_elements_max
      call get_cache_map(ao_integrals_map,i8,keys,values,n_elements)
      do k1=1,n_elements
 !       call bielec_integrals_index_reverse(kk,ii,ll,jj,keys(k1))
        call bielec_integrals_index_reverse_4fold(ii,jj,kk,ll,keys(k1),pair_key)
        if (keys(k1) <= pair_key) then
          !values(k1) is real part
          do k2=1,4
            if (ii(k2)==0) then
              cycle
            endif
            i = ii(k2)
            j = jj(k2)
            k = kk(k2)
            l = ll(k2)
            integral = values(k1)
            
            c0 = (HF_density_matrix_ao_alpha(l,j)+HF_density_matrix_ao_beta(l,j)) * integral
            c1 = HF_density_matrix_ao_alpha(k,j) * integral
            c2 = HF_density_matrix_ao_beta(k,j) * integral


            ao_bi_elec_integral_alpha_tmp(i,k) += c0
            ao_bi_elec_integral_beta_tmp (i,k) += c0
            ao_bi_elec_integral_alpha_tmp(i,l) -= c1
            ao_bi_elec_integral_beta_tmp (i,l) -= c2
          enddo
        else
          !values(k1) is imaginary part of <kl|ij>,<lk|ji>
          !       and -1*imaginary part of <ij|kl>,<ji|lk>
          do k2=1,4
            if (ii(k2)==0) then
              cycle
            endif
            i = ii(k2)
            j = jj(k2)
            k = kk(k2)
            l = ll(k2)
            integral = i_sign(k2)*values(k1)
            
            c0 = (HF_density_matrix_ao_alpha(l,j)+HF_density_matrix_ao_beta(l,j)) * integral
            c1 = HF_density_matrix_ao_alpha(k,j) * integral
            c2 = HF_density_matrix_ao_beta(k,j) * integral


            ao_bi_elec_integral_alpha_tmp(i,k) += c0
            ao_bi_elec_integral_beta_tmp (i,k) += c0
            ao_bi_elec_integral_alpha_tmp(i,l) -= c1
            ao_bi_elec_integral_beta_tmp (i,l) -= c2
          enddo
        endif
      enddo
    enddo
   !$OMP END DO NOWAIT
   !$OMP CRITICAL
   ao_bi_elec_integral_alpha += ao_bi_elec_integral_alpha_tmp
   !$OMP END CRITICAL
   !$OMP CRITICAL
   ao_bi_elec_integral_beta  += ao_bi_elec_integral_beta_tmp
   !$OMP END CRITICAL
   deallocate(keys,values,ao_bi_elec_integral_alpha_tmp,ao_bi_elec_integral_beta_tmp)
   !$OMP END PARALLEL

END_PROVIDER

BEGIN_PROVIDER [ complex*16, Fock_matrix_mo_alpha, (mo_tot_num,mo_tot_num) ]
   implicit none
   BEGIN_DOC
   ! Fock matrix on the MO basis
   END_DOC
   call complex_ao_to_mo(Fock_matrix_ao_alpha,size(Fock_matrix_ao_alpha,1), &
                 Fock_matrix_mo_alpha,size(Fock_matrix_mo_alpha,1))
END_PROVIDER
 
 
BEGIN_PROVIDER [ complex*16, Fock_matrix_mo_beta, (mo_tot_num,mo_tot_num) ]
   implicit none
   BEGIN_DOC
   ! Fock matrix on the MO basis
   END_DOC
   call complex_ao_to_mo(Fock_matrix_ao_beta,size(Fock_matrix_ao_beta,1), &
                 Fock_matrix_mo_beta,size(Fock_matrix_mo_beta,1))
END_PROVIDER
 
BEGIN_PROVIDER [ double precision, HF_energy ]
 implicit none
 BEGIN_DOC
 ! Hartree-Fock energy
 END_DOC
 complex*16 :: HF_energy_tmp
 HF_energy_tmp = nuclear_repulsion
 
 integer                        :: i,j
 do j=1,ao_num
   do i=1,ao_num
     HF_energy_tmp += 0.5d0 * (                                          &
         (ao_mono_elec_integral(i,j) + Fock_matrix_ao_alpha(i,j) ) *  HF_density_matrix_ao_alpha(i,j) +&
         (ao_mono_elec_integral(i,j) + Fock_matrix_ao_beta (i,j) ) *  HF_density_matrix_ao_beta (i,j) )
   enddo
 enddo

 if (abs(imag(HF_energy_tmp)) .gt. 1.0d-12) then
   stop 'HF energy should be real'
 endif
 HF_energy = real(HF_energy_tmp)
  
END_PROVIDER
 
 BEGIN_PROVIDER [ double precision, one_elec_energy ]
&BEGIN_PROVIDER [ double precision, two_elec_energy ]
 implicit none
 BEGIN_DOC
 ! Hartree-Fock energy
 END_DOC
 complex*16 :: e1_tmp, e2_tmp
 e1_tmp = (0.d0,0.d0) 
 e2_tmp = (0.d0,0.d0) 
 integer                        :: i,j
 do j=1,ao_num
   do i=1,ao_num
     e1_tmp +=  (                                          &
         (ao_mono_elec_integral(i,j) ) *  HF_density_matrix_ao_alpha(i,j) +&
         (ao_mono_elec_integral(i,j) ) *  HF_density_matrix_ao_beta (i,j) )
     e2_tmp += 0.5 * (                                          &
         (ao_bi_elec_integral_alpha(i,j) ) *  HF_density_matrix_ao_alpha(i,j) +&
         (ao_bi_elec_integral_beta(i,j) ) *  HF_density_matrix_ao_beta (i,j) )
   enddo
 enddo

 if (abs(imag(e1_tmp)) .gt. 1.0d-12) then
   stop '1-elec energy should be real'
 endif
 if (abs(imag(e2_tmp)) .gt. 1.0d-12) then
   stop '2-elec energy should be real'
 endif
 one_elec_energy = real(e1_tmp)
 two_elec_energy = real(e2_tmp)
  
END_PROVIDER


BEGIN_PROVIDER [ complex*16, Fock_matrix_ao, (ao_num, ao_num) ]
 implicit none
 BEGIN_DOC
 ! Fock matrix in AO basis set
 END_DOC
 
 if ( (elec_alpha_num == elec_beta_num).and. &
      (level_shift == 0.) ) &
 then
   integer                        :: i,j
   do j=1,ao_num
     do i=1,ao_num
       Fock_matrix_ao(i,j) = Fock_matrix_ao_alpha(i,j)
     enddo
   enddo
!   call zlacpy('X',ao_num,ao_num, &
!         Fock_matrix_ao_alpha, size(Fock_matrix_ao_alpha,1), &
!         Fock_matrix_ao, size(Fock_matrix_ao,1))
 else
   call mo_to_ao(Fock_matrix_mo,size(Fock_matrix_mo,1), &
      Fock_matrix_ao,size(Fock_matrix_ao,1)) 
 endif
END_PROVIDER


