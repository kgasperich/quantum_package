

 BEGIN_PROVIDER [ integer, homo_idx ]
&BEGIN_PROVIDER [ integer, homo_idx_int ]
&BEGIN_PROVIDER [ integer, homo_idx_bit ]
&BEGIN_PROVIDER [ integer, homo_spin ]
&BEGIN_PROVIDER [ double precision, homo_sign ]
  implicit none
  BEGIN_DOC
  ! description of particles/holes to be used in spectral density calculation
  ! green_idx: orbital index of particle/hole
  ! green_idx_{int,bit}: location of idx within determinant bitstring
  ! green_spin: 1(alpha) or 2(beta)
  ! green_sign: 1(particle) or -1(hole)
  END_DOC
  integer :: s1,s2,i1,i2
  integer :: i

  integer :: idx_homo_lumo(2), spin_homo_lumo(2)
  ! needs psi_det, mo_tot_num, N_int, mo_bielec_integral_jj, mo_mono_elec_integral_diag
  call get_homo_lumo(psi_det(1:N_int,1:2,1),N_int,mo_tot_num,idx_homo_lumo,spin_homo_lumo)

  homo_idx=idx_homo_lumo(1)
  homo_spin=spin_homo_lumo(1)
  homo_sign=-1.d0
  call get_orb_int_bit(homo_idx,homo_idx_int,homo_idx_bit)
  print*,homo_idx,homo_idx_int,homo_idx_bit

END_PROVIDER


 BEGIN_PROVIDER [ double precision, homo_phase, (N_det) ]
&BEGIN_PROVIDER [ logical, homo_allowed, (N_det) ]
  implicit none
  BEGIN_DOC
  ! for each det in psi, compute phase for each particle/hole excitation
  ! each element should be +/-1 or 0
  END_DOC
  integer :: i
  double precision :: phase_tmp
  logical :: allowed_tmp
  PROVIDE psi_det homo_idx
  
  do i=1,N_det
    call get_phase_homo(homo_idx_int,homo_idx_bit,homo_spin,homo_sign,psi_det(1,1,i),phase_tmp,N_int,allowed_tmp)
    homo_phase(i) = phase_tmp
    homo_allowed(i) = allowed_tmp
  enddo
END_PROVIDER

subroutine get_phase_homo(h_idx_int,h_idx_bit,h_spin,h_sign,det_in,h_det_phase,nint,is_allowed)
  use bitmasks
  implicit none
  integer, intent(in) :: nint
  integer, intent(in) :: h_idx_int, h_idx_bit,h_spin
  double precision, intent(in) :: h_sign
  integer(bit_kind), intent(in) :: det_in(nint,2)
  double precision, intent(out) :: h_det_phase
  logical, intent(out) :: is_allowed

  integer(bit_kind) :: tmp_spindet(nint), pm(nint,2)
  double precision, parameter :: phase_dble(0:1) = (/1.d0,-1.d0/)

  integer :: i
  logical :: is_filled

  tmp_spindet(1:nint) = det_in(1:nint,h_spin)
  call spinorb_is_filled_int_bit(tmp_spindet,h_idx_int,h_idx_bit,nint,is_filled)

  is_allowed = ((h_sign<0).eqv.(is_filled))
  if (is_allowed) then
    call get_mask_phase_new(det_in,pm,nint)
    h_det_phase = phase_dble(popcnt(iand(ibset(0_bit_kind,h_idx_bit),pm(h_idx_int,h_spin))))
  else
    h_det_phase=0.d0
  endif
end



BEGIN_PROVIDER [ complex*16, u1_lanczos_homo, (N_det) ]
  implicit none
  BEGIN_DOC
  ! initial lanczos vectors
  ! must be normalized
  END_DOC

  integer :: i
  
    do i=1,N_det
      u1_lanczos_homo(i)=homo_phase(i)*psi_coef(i,1)
    enddo
    call normalize_complex(u1_lanczos_homo,N_det)

END_PROVIDER

 BEGIN_PROVIDER [ double precision, alpha_lanczos_homo, (n_lanczos_iter) ]
&BEGIN_PROVIDER [ double precision, beta_lanczos_homo, (n_lanczos_iter) ]
&BEGIN_PROVIDER [ complex*16, un_lanczos_homo, (N_det) ]
&BEGIN_PROVIDER [ complex*16, vn_lanczos_homo, (N_det) ]
&BEGIN_PROVIDER [ double precision, lanczos_eigvals_homo, (n_lanczos_iter) ]
  implicit none
  BEGIN_DOC
  ! for each particle/hole:
  ! provide alpha and beta for tridiagonal form of H
  ! un, vn lanczos vectors from latest iteration
  ! lanczos_eigvals: eigenvalues of tridiagonal form of H
  END_DOC
!  complex*16, allocatable :: work(:,:)
!  double precision, allocatable :: alpha_tmp(:),beta_tmp(:)
!  double precision, allocatable :: alpha_tmp_vec(:,:), beta_tmp_vec(:,:)
  complex*16, allocatable :: work(:)
  double precision :: alpha_tmp,beta_tmp
  double precision, allocatable :: alpha_tmp_vec(:), beta_tmp_vec(:)
  integer :: i,j
  integer :: n_lanc_new_tmp, n_lanc_old_tmp
  call ezfio_get_green_n_lanczos_iter(n_lanc_new_tmp)
  call ezfio_get_green_n_lanczos_complete(n_lanc_old_tmp)
 
  if ((n_lanczos_complete).gt.0) then
!    allocate(alpha_tmp_vec(n_lanczos_complete,n_green_vec),beta_tmp_vec(n_lanczos_complete,n_green_vec))
    allocate(alpha_tmp_vec(n_lanczos_complete),beta_tmp_vec(n_lanczos_complete))
    logical :: has_un_lanczos_homo, has_vn_lanczos_homo
    call ezfio_has_green_un_lanczos_homo(has_un_lanczos_homo)
    call ezfio_has_green_vn_lanczos_homo(has_vn_lanczos_homo)
    if (has_un_lanczos_homo.and.has_vn_lanczos_homo) then
      call ezfio_get_green_un_lanczos_homo(un_lanczos_homo)
      call ezfio_get_green_vn_lanczos_homo(vn_lanczos_homo)
!      if (lanczos_debug_print) then
!        print*,'uu,vv read from disk'
!        do i=1,n_lanczos_debug
!          write(6,'(4(E25.15))')un_lanczos_homo(i),vn_lanczos_homo(i)
!        enddo
!      endif
    else
      print*,'problem reading lanczos vectors for restart'
      stop
    endif
    logical :: has_alpha_lanczos, has_beta_lanczos
    call ezfio_has_green_alpha_lanczos_homo(has_alpha_lanczos)
    call ezfio_has_green_beta_lanczos_homo(has_beta_lanczos)
    if (has_alpha_lanczos.and.has_beta_lanczos) then
      call ezfio_set_green_n_lanczos_iter(n_lanc_old_tmp)
      call ezfio_get_green_alpha_lanczos_homo(alpha_tmp_vec)
      call ezfio_get_green_beta_lanczos_homo(beta_tmp_vec)
      call ezfio_set_green_n_lanczos_iter(n_lanc_new_tmp)
        do i=1,n_lanczos_complete
          alpha_lanczos_homo(i)=alpha_tmp_vec(i)
          beta_lanczos_homo(i)=beta_tmp_vec(i)
        enddo
    else
      print*,'problem reading lanczos alpha, beta for restart'
      stop
    endif
    deallocate(alpha_tmp_vec,beta_tmp_vec)
  else
    call write_time(6)
    print*,'no saved lanczos vectors. starting lanczos'
    PROVIDE u1_lanczos_homo
    un_lanczos_homo=u1_lanczos_homo
    !allocate(work(N_det,n_green_vec),alpha_tmp(n_green_vec),beta_tmp(n_green_vec))
    allocate(work(N_det))
    call lanczos_h_init_homo(un_lanczos_homo,vn_lanczos_homo,work,N_det,alpha_tmp,beta_tmp,&
                           homo_spin,homo_sign,homo_idx)
    alpha_lanczos_homo(1)=alpha_tmp
    beta_lanczos_homo(1)=beta_tmp
    n_lanczos_complete=1
    !deallocate(work,alpha_tmp,beta_tmp)
    deallocate(work)
  endif

  !allocate(work(N_det,n_green_vec),alpha_tmp(n_green_vec),beta_tmp(n_green_vec))
  allocate(work(N_det))
  do i=n_lanczos_complete+1,n_lanczos_iter
    call write_time(6)
    print*,'starting lanczos iteration',i
    call lanczos_h_step_homo(un_lanczos_homo,vn_lanczos_homo,work,N_det,alpha_tmp,beta_tmp,&
                           homo_spin,homo_sign,homo_idx)
    alpha_lanczos_homo(i)=alpha_tmp
    beta_lanczos_homo(i)=beta_tmp
    n_lanczos_complete=n_lanczos_complete+1
  enddo
  !deallocate(work,alpha_tmp,beta_tmp)
  deallocate(work)

  call ezfio_set_green_alpha_lanczos_homo(alpha_lanczos_homo)
  call ezfio_set_green_beta_lanczos_homo(beta_lanczos_homo)
  call ezfio_set_green_un_lanczos_homo(un_lanczos_homo)
  call ezfio_set_green_vn_lanczos_homo(vn_lanczos_homo)
  call ezfio_set_green_n_lanczos_complete(n_lanczos_complete)

  call diag_lanczos_vals_homo(alpha_lanczos_homo, beta_lanczos_homo, n_lanczos_complete, lanczos_eigvals_homo,&
                            n_lanczos_iter)
  call ezfio_set_green_lanczos_eigvals_homo(lanczos_eigvals_homo)

END_PROVIDER

!BEGIN_PROVIDER [ double precision, delta_omega ]
!  implicit none
!  BEGIN_DOC
!  ! step size between frequency points for spectral density calculation
!  ! calculated from min, max, and number of steps
!  END_DOC
!  delta_omega=(omega_max-omega_min)/n_omega
!END_PROVIDER
!
!BEGIN_PROVIDER [ double precision, omega_list, (n_omega) ]
!  implicit none
!  BEGIN_DOC
!  ! list of frequencies at which to compute spectral density
!  END_DOC
!
!  integer :: i
!  double precision :: omega_i
!  PROVIDE delta_omega
!  do i=1,n_omega
!    omega_list(i) = omega_min + (i-1)*delta_omega
!  enddo
!
!END_PROVIDER


BEGIN_PROVIDER [ double precision, spectral_lanczos_homo, (n_omega) ]
  implicit none
  BEGIN_DOC
  ! spectral density A(omega) calculated from lanczos alpha/beta
  ! calculated for n_omega points between omega_min and omega_max
  END_DOC

  integer :: i,j
  double precision :: omega_i
  complex*16 :: z_i
  double precision :: spec_lanc
  PROVIDE delta_omega alpha_lanczos_homo beta_lanczos_homo omega_list
  do i=1,n_omega
    omega_i = omega_list(i)
    z_i = dcmplx(omega_i,gf_epsilon)
    spectral_lanczos_homo(i) = spec_lanc(n_lanczos_iter,alpha_lanczos_homo,beta_lanczos_homo,z_i)
  enddo

END_PROVIDER

!double precision function spec_lanc(n_lanc_iter,alpha,beta,z)
!  include 'constants.include.F'
!  implicit none
!  BEGIN_DOC
!  ! input:
!  !   alpha, beta: from tridiagonal form of H (obtain via lanczos)
!  !                beta and alpha same size (beta(1) is not used)
!  !   n_lanc_iter: size of alpha, beta
!  !             z: omega + i*epsilon
!  !                omega is frequency for which spectral density is to be computed
!  !                epsilon is magnitude of infinitesimal imaginary term
!  ! output:
!  !     spec_lanc: spectral density A(omega)
!  !
!  ! uses inv_pi=(1.d0/pi) from constants 
!  END_DOC
!  integer, intent(in) :: n_lanc_iter
!  double precision, intent(in) :: alpha(n_lanc_iter), beta(n_lanc_iter)
!  complex*16, intent(in) :: z
!
!  complex*16 bigAj2,bigAj1,bigAj0
!  complex*16 bigBj2,bigBj1,bigBj0
!  integer :: j
!  ! init for j=1
!  ! bigAj2 is A(j-2)
!  ! bigAj1 is A(j-1)
!  ! etc.
!
!  bigAj2=1.d0        ! A(-1)
!  bigAj1=0.d0        ! A(0)
!  bigAj0=1.d0        ! A(1)
!
!  bigBj2=0.d0        ! B(-1)
!  bigBj1=1.d0        ! B(0)
!  bigBj0=z-alpha(1)  ! B(1)
!
!  do j=2,n_lanc_iter
!    bigAj2=bigAj1
!    bigAj1=bigAj0
!    bigAj0=(z-alpha(j))*bigAj1 - beta(j)**2*bigAj2
!    
!    bigBj2=bigBj1
!    bigBj1=bigBj0
!    bigBj0=(z-alpha(j))*bigBj1 - beta(j)**2*bigBj2
!  enddo
!  spec_lanc=-imag(bigAj0/bigBj0)*inv_pi
!end


subroutine lanczos_h_init_homo(uu,vv,work,sze,alpha_i,beta_i,spin_hp,sign_hp,idx_hp)
  implicit none
  integer, intent(in) :: sze
  complex*16, intent(in)    :: uu(sze)
  complex*16, intent(out)   :: vv(sze)
  complex*16 :: work(sze)
  double precision, intent(out) :: alpha_i, beta_i
  integer, intent(in) :: spin_hp, idx_hp
  double precision, intent(in) :: sign_hp

  double precision, external :: dznrm2
  complex*16, external :: u_dot_v_complex
  integer :: i,j

  BEGIN_DOC
  ! initial step for lanczos tridiagonalization of H for multiple holes/particles
  ! uu is array of initial vectors u1 (creation/annihilation operator applied to psi)
  ! output vv is array of lanczos v1 (one for each hole/particle)
  END_DOC

  print *,'starting lanczos'
  print *,'sze = ',sze

  ! |uu> is |u(1)>

  ! |w(1)> = H|u(1)>
  ! |work> is now |w(1)>
  call compute_hu_homo(uu,work,sze,spin_hp,sign_hp,idx_hp)

  ! alpha(n+1) = <u(n+1)|w(n+1)>
!  do i=1,ng
    alpha_i=real(u_dot_v_complex(uu,work,sze))
!  enddo

!  do j=1,ng
    do i=1,sze
      vv(i)=work(i)-alpha_i*uu(i)
      write(6,'(7(E25.15))')uu(i),vv(i),work(i),alpha_i
    enddo
!  enddo
  
  beta_i=0.d0
  ! |vv> is |v(1)>
  ! |uu> is |u(1)>
end

subroutine lanczos_h_step_homo(uu,vv,work,sze,alpha_i,beta_i,spin_hp,sign_hp,idx_hp)
  implicit none
  integer, intent(in) :: sze
  complex*16, intent(inout) :: uu(sze),vv(sze)
  complex*16, intent(out) :: work(sze)
  double precision, intent(out) :: alpha_i, beta_i
  integer, intent(in) :: spin_hp, idx_hp
  double precision, intent(in) :: sign_hp

  double precision, external :: dznrm2
  complex*16, external :: u_dot_v_complex
  integer :: i,j
  complex*16 :: tmp_c16
  BEGIN_DOC
  ! lanczos tridiagonalization of H
  ! n_lanc_iter is number of lanczos iterations
  ! u1 is initial lanczos vector
  ! u1 should be normalized
  END_DOC

  ! |vv> is |v(n)>
  ! |uu> is |u(n)>

  ! compute beta(n+1)
 ! do j=1,ng
    beta_i=dznrm2(sze,vv,1)
  ! |vv> is now |u(n+1)>
    call zdscal(sze,(1.d0/beta_i),vv,1)
 ! enddo

  ! |w(n+1)> = H|u(n+1)>
  ! |work> is now |w(n+1)>
  call compute_hu_homo(vv,work,sze,spin_hp,sign_hp,idx_hp)

  ! alpha(n+1) = <u(n+1)|w(n+1)>
!  do i=1,ng
    alpha_i=real(u_dot_v_complex(vv,work,sze))
!  enddo

!  do j=1,ng
    do i=1,sze
      tmp_c16=work(i)-alpha_i*vv(i)-beta_i*uu(i)
      uu(i)=vv(i)
      vv(i)=tmp_c16
    enddo
!  enddo
  ! |vv> is |v(n+1)>
  ! |uu> is |u(n+1)>
end




subroutine compute_hu_homo(vec1,vec2,h_size,spin_hp,sign_hp,idx_hp)
  implicit none
  integer, intent(in)     :: h_size
  complex*16, intent(in)  :: vec1(h_size)
  complex*16, intent(out) :: vec2(h_size)
  integer, intent(in) :: spin_hp, idx_hp
  double precision, intent(in) :: sign_hp
  complex*16 :: vec1_tmp(h_size)
  integer :: i,j
  BEGIN_DOC
  ! |vec2> = H|vec1>
  !
  ! TODO: implement
  ! maybe reuse parts of H_S2_u_0_nstates_{openmp,zmq}?
  END_DOC

  vec1_tmp(1:h_size) = vec1(1:h_size)
  call h_u_0_homo_openmp(vec2,vec1_tmp,h_size,spin_hp,sign_hp,idx_hp)

!  do j=1,n_hp
    do i=1,h_size
      if (cdabs(vec1_tmp(i) - vec1(i)).gt.1.d-6) then
        print*,'ERROR: vec1 was changed by h_u_0_openmp'
      endif
    enddo
!  enddo
end

subroutine diag_lanczos_vals_homo(alpha, beta, nlanc, vals, sze)
  implicit none
  BEGIN_DOC
  ! diagonalization of tridiagonal form of H
  ! this returns eigenvalues in vals
  END_DOC
  integer, intent(in) :: nlanc,sze
  !double precision, intent(in) :: alpha(ng,sze), beta(sze)
  double precision, intent(in) :: alpha(sze), beta(sze)
  double precision, intent(out) :: vals(sze)
  double precision :: work(1), beta_tmp(nlanc-1), vecs(1)
  integer :: i,info,ig
  
!  do ig=1,ng
    vals(1)=alpha(1)
    do i=2,nlanc
      vals(i)=alpha(i)
      beta_tmp(i-1)=beta(i)
    enddo
  
    call dstev('N', nlanc, vals, beta_tmp, vecs, 1, work, info)
    if (info.gt.0) then
      print *,'WARNING: diagonalization of tridiagonal form of H did not converge'
    else if (info.lt.0) then
      print *,'WARNING: argument to dstev had illegal value'
    endif
 ! enddo
end
