BEGIN_PROVIDER [ complex*16, u1_lanczos, (N_det) ]
  implicit none
  BEGIN_DOC
  ! initial lanczos vector
  ! must be normalized
  END_DOC

  integer :: i
  
  do i=1,N_det
    u1_lanczos(i)=1.d0/(dble(i))**0.1d0
  enddo
  call normalize_complex(u1_lanczos,N_det)

END_PROVIDER

subroutine init_u1_lanczos(u_in,sze)
  implicit none
  integer, intent(in) :: sze
  complex*16, intent(inout) :: u_in(sze)
  integer :: i
  do i=1,sze
    u_in(i)=1.d0/(dble(i))**0.1d0
!    u_in(i)=1.d0
  enddo
  call normalize_complex(u_in,sze)
end

 BEGIN_PROVIDER [ double precision, alpha_lanczos, (n_lanczos_iter) ]
&BEGIN_PROVIDER [ double precision, beta_lanczos, (n_lanczos_iter) ]
&BEGIN_PROVIDER [ complex*16, un_lanczos, (N_det) ]
&BEGIN_PROVIDER [ complex*16, vn_lanczos, (N_det) ]
&BEGIN_PROVIDER [ double precision, lanczos_eigvals, (n_lanczos_iter) ]
  implicit none
  BEGIN_DOC
  ! provide alpha and beta for tridiagonal form of H
  END_DOC
  PROVIDE lanczos_debug_print n_lanczos_debug
  complex*16, allocatable :: work(:)
  double precision :: alpha_tmp,beta_tmp
  double precision, allocatable :: alpha_tmp_vec(:), beta_tmp_vec(:)
  integer :: i
  integer :: n_lanc_new_tmp, n_lanc_old_tmp
  call ezfio_get_green_n_lanczos_iter(n_lanc_new_tmp)
  call ezfio_get_green_n_lanczos_complete(n_lanc_old_tmp)
 
  if ((n_lanczos_complete).gt.0) then
    allocate(alpha_tmp_vec(n_lanczos_complete),beta_tmp_vec(n_lanczos_complete))
    logical :: has_un_lanczos, has_vn_lanczos
    call ezfio_has_green_un_lanczos(has_un_lanczos)
    call ezfio_has_green_vn_lanczos(has_vn_lanczos)
    if (has_un_lanczos.and.has_vn_lanczos) then
      call ezfio_get_green_un_lanczos(un_lanczos)
      call ezfio_get_green_vn_lanczos(vn_lanczos)
      if (lanczos_debug_print) then
        print*,'uu,vv read from disk'
        do i=1,n_lanczos_debug
          write(6,'(4(E25.15))')un_lanczos(i),vn_lanczos(i)
        enddo
      endif
    else
      print*,'problem reading lanczos vectors for restart'
      stop
    endif
    logical :: has_alpha_lanczos, has_beta_lanczos
    call ezfio_has_green_alpha_lanczos(has_alpha_lanczos)
    call ezfio_has_green_beta_lanczos(has_beta_lanczos)
    if (has_alpha_lanczos.and.has_beta_lanczos) then
      call ezfio_set_green_n_lanczos_iter(n_lanc_old_tmp)
      call ezfio_get_green_alpha_lanczos(alpha_tmp_vec)
      call ezfio_get_green_beta_lanczos(beta_tmp_vec)
      call ezfio_set_green_n_lanczos_iter(n_lanc_new_tmp)
      do i=1,n_lanczos_complete
        alpha_lanczos(i)=alpha_tmp_vec(i)
        beta_lanczos(i)=beta_tmp_vec(i)
      enddo
    else
      print*,'problem reading lanczos alpha, beta for restart'
      stop
    endif
    deallocate(alpha_tmp_vec,beta_tmp_vec)
  else
    call write_time(6)
    print*,'no saved lanczos vectors. starting lanczos'
    call init_u1_lanczos(un_lanczos,N_det)
    allocate(work(N_det))
    call lanczos_h_init(un_lanczos,vn_lanczos,work,N_det,alpha_tmp,beta_tmp)
    alpha_lanczos(1)=alpha_tmp
    beta_lanczos(1)=beta_tmp
    n_lanczos_complete=1
    deallocate(work)
  endif

  allocate(work(N_det))
  do i=n_lanczos_complete+1,n_lanczos_iter
    call write_time(6)
    print*,'starting lanczos iteration',i
    call lanczos_h_step(un_lanczos,vn_lanczos,work,N_det,alpha_tmp,beta_tmp)
    alpha_lanczos(i)=alpha_tmp
    beta_lanczos(i)=beta_tmp
    n_lanczos_complete=n_lanczos_complete+1
  enddo
  deallocate(work)

  call ezfio_set_green_alpha_lanczos(alpha_lanczos)
  call ezfio_set_green_beta_lanczos(beta_lanczos)
  call ezfio_set_green_un_lanczos(un_lanczos)
  call ezfio_set_green_vn_lanczos(vn_lanczos)
  call ezfio_set_green_n_lanczos_complete(n_lanczos_complete)

  call diag_lanczos_vals(alpha_lanczos, beta_lanczos, n_lanczos_complete, lanczos_eigvals, n_lanczos_iter)
  call ezfio_set_green_lanczos_eigvals(lanczos_eigvals)

END_PROVIDER

BEGIN_PROVIDER [ double precision, delta_omega ]
  implicit none
  BEGIN_DOC
  ! step size between frequency points for spectral density calculation
  ! calculated from min, max, and number of steps
  END_DOC
  delta_omega=(omega_max-omega_min)/n_omega
END_PROVIDER

BEGIN_PROVIDER [ double precision, omega_list, (n_omega) ]
  implicit none
  BEGIN_DOC
  ! test
  ! 
  END_DOC

  integer :: i
  double precision :: omega_i
  PROVIDE delta_omega
  do i=1,n_omega
    omega_list(i) = omega_min + (i-1)*delta_omega
  enddo

END_PROVIDER


BEGIN_PROVIDER [ double precision, spectral_lanczos, (n_omega) ]
  implicit none
  BEGIN_DOC
  ! spectral density A(omega) calculated from lanczos alpha/beta
  ! calculated for n_omega points between omega_min and omega_max
  END_DOC

  integer :: i
  double precision :: omega_i
  complex*16 :: z_i
  double precision :: spec_lanc
  PROVIDE delta_omega alpha_lanczos beta_lanczos omega_list
  do i=1,n_omega
    omega_i = omega_list(i)
    z_i = dcmplx(omega_i,gf_epsilon)
    spectral_lanczos(i) = spec_lanc(n_lanczos_iter,alpha_lanczos,beta_lanczos,z_i)
  enddo

END_PROVIDER

double precision function spec_lanc(n_lanc_iter,alpha,beta,z)
  include 'constants.include.F'
  implicit none
  BEGIN_DOC
  ! input:
  !   alpha, beta: from tridiagonal form of H (obtain via lanczos)
  !                beta and alpha same size (beta(1) is not used)
  !   n_lanc_iter: size of alpha, beta
  !             z: omega + i*epsilon
  !                omega is frequency for which spectral density is to be computed
  !                epsilon is magnitude of infinitesimal imaginary term
  ! output:
  !     spec_lanc: spectral density A(omega)
  !
  ! uses inv_pi=(1.d0/pi) from constants 
  END_DOC
  integer, intent(in) :: n_lanc_iter
  double precision, intent(in) :: alpha(n_lanc_iter), beta(n_lanc_iter)
  complex*16, intent(in) :: z

  complex*16 bigAj2,bigAj1,bigAj0
  complex*16 bigBj2,bigBj1,bigBj0
  integer :: j
  ! init for j=1
  ! bigAj2 is A(j-2)
  ! bigAj1 is A(j-1)
  ! etc.

  bigAj2=1.d0        ! A(-1)
  bigAj1=0.d0        ! A(0)
  bigAj0=1.d0        ! A(1)

  bigBj2=0.d0        ! B(-1)
  bigBj1=1.d0        ! B(0)
  bigBj0=z-alpha(1)  ! B(1)

  do j=2,n_lanc_iter
    bigAj2=bigAj1
    bigAj1=bigAj0
    bigAj0=(z-alpha(j))*bigAj1 - beta(j)**2*bigAj2
    
    bigBj2=bigBj1
    bigBj1=bigBj0
    bigBj0=(z-alpha(j))*bigBj1 - beta(j)**2*bigBj2
  enddo
  spec_lanc=-imag(bigAj0/bigBj0)*inv_pi
end


subroutine lanczos_h_init(uu,vv,work,sze,alpha_i,beta_i)
  implicit none
  integer, intent(in) :: sze
  complex*16, intent(inout) :: uu(sze)
  complex*16, intent(out)   :: vv(sze)
  complex*16 :: work(sze)
  double precision, intent(out) :: alpha_i, beta_i

  double precision, external :: dznrm2
  complex*16, external :: u_dot_v_complex
  integer :: i

  BEGIN_DOC
  ! lanczos tridiagonalization of H
  ! n_lanc_iter is number of lanczos iterations
  ! u1 is initial lanczos vector
  ! u1 should be normalized
  END_DOC

  print *,'starting lanczos'
  print *,'sze = ',sze
  ! exit if u1 is not normalized
!  beta_norm = dznrm2(h_size,u1,1)
!  if (dabs(beta_norm-1.d0) .gt. 1.d-6) then
!    print *, 'Error: initial Lanczos vector is not normalized'
!    stop -1
!  endif

  ! |uu> is |u(1)>

  ! |w(1)> = H|u(1)>
  ! |work> is now |w(1)>
  call compute_hu(uu,work,sze)

  ! alpha(n+1) = <u(n+1)|w(n+1)>
  alpha_i=real(u_dot_v_complex(uu,work,sze))

  do i=1,sze
    vv(i)=work(i)-alpha_i*uu(i)
  enddo
  beta_i=0.d0
  if (lanczos_debug_print) then
    print*,'init uu,vv,work'
    do i=1,n_lanczos_debug
      write(6,'(6(E25.15))')uu(i),vv(i),work(i)
    enddo
  endif
  ! |vv> is |v(1)>
  ! |uu> is |u(1)>
end

subroutine lanczos_h_step(uu,vv,work,sze,alpha_i,beta_i)
  implicit none
  integer, intent(in) :: sze
  complex*16, intent(inout) :: uu(sze),vv(sze)
  complex*16, intent(out) :: work(sze)
  double precision, intent(out) :: alpha_i, beta_i

  double precision, external :: dznrm2
  complex*16, external :: u_dot_v_complex
  integer :: i
  complex*16 :: tmp_c16
  BEGIN_DOC
  ! lanczos tridiagonalization of H
  ! n_lanc_iter is number of lanczos iterations
  ! u1 is initial lanczos vector
  ! u1 should be normalized
  END_DOC

  ! exit if u1 is not normalized
!  beta_norm = dznrm2(h_size,u1,1)
!  if (dabs(beta_norm-1.d0) .gt. 1.d-6) then
!    print *, 'Error: initial Lanczos vector is not normalized'
!    stop -1
!  endif

  ! |vv> is |v(n)>
  ! |uu> is |u(n)>

  ! compute beta(n+1)
  beta_i=dznrm2(sze,vv,1)
  if (lanczos_debug_print) then
    print*,'uu,vv in'
    do i=1,n_lanczos_debug
      write(6,'(4(E25.15))')uu(i),vv(i)
    enddo
  endif
  ! |vv> is now |u(n+1)>
  call zdscal(sze,(1.d0/beta_i),vv,1)

  ! |w(n+1)> = H|u(n+1)>
  ! |work> is now |w(n+1)>
  call compute_hu(vv,work,sze)

  if (lanczos_debug_print) then
    print*,'vv,work'
    do i=1,n_lanczos_debug
      write(6,'(4(E25.15))')vv(i),work(i)
    enddo
  endif

  ! alpha(n+1) = <u(n+1)|w(n+1)>
  alpha_i=real(u_dot_v_complex(vv,work,sze))

  do i=1,sze
    tmp_c16=work(i)-alpha_i*vv(i)-beta_i*uu(i)
    uu(i)=vv(i)
    vv(i)=tmp_c16
  enddo
  ! |vv> is |v(n+1)>
  ! |uu> is |u(n+1)>
end



subroutine lanczos_h(n_lanc_iter,alpha,beta,u1)
  implicit none
  integer, intent(in) :: n_lanc_iter
  double precision, intent(out) :: alpha(n_lanc_iter), beta(n_lanc_iter)
  complex*16, intent(in) :: u1(N_det)
  integer :: h_size
  double precision :: beta_norm, beta_norm_inv
  complex*16, allocatable :: vec1(:), vec2(:), vec3(:)
  complex*16 :: vec_tmp
  double precision, external :: dznrm2
  complex*16, external :: u_dot_v_complex

  integer :: i,j,l
  h_size=N_det
  BEGIN_DOC
  ! lanczos tridiagonalization of H
  ! n_lanc_iter is number of lanczos iterations
  ! u1 is initial lanczos vector
  ! u1 should be normalized
  END_DOC

  print *,'starting lanczos'
  print *,'h_size = ',h_size
!  print *,'initial vector:'
!  do i=1,h_size
!    print *,u1(i)
!  enddo
  ! exit if u1 is not normalized
  beta_norm = dznrm2(h_size,u1,1)
  if (dabs(beta_norm-1.d0) .gt. 1.d-6) then
    print *, 'Error: initial Lanczos vector is not normalized'
    stop -1
  endif

  allocate(vec1(h_size),  &
           vec2(h_size),  &
           vec3(h_size))  
 
  do i=1,h_size
    vec1(i)=u1(i)
  enddo

  ! |w1> = H|u1>
  ! |vec2> = H|vec1>
  call compute_hu(vec1,vec2,h_size)!! TODO: not implemented

  ! alpha(1) = <u1|H|u1> = <u1|w1>
  !          = <vec1|vec2>
  alpha(1)=real(u_dot_v_complex(vec1,vec2,h_size)) 

  ! |v1> = |w1> - alpha(1)*|u1>
  ! |vec3> = |vec2> - alpha(1)*|vec1>
  do i=1,h_size
    vec3(i)=vec2(i)-alpha(1)*vec1(i)
  enddo
  do j=2,n_lanc_iter
    call write_time(6)
    print *,'starting lanczos iteration:',j
    !! vec1 is |u(j-1)>
    !! vec3 is |v(j-1)>

    ! beta(j) = sqrt(<v(j-1)|v(j-1)>)
    beta_norm=dznrm2(h_size,vec3,1)

    ! TODO: check for beta=0?
    beta_norm_inv=1.d0/beta_norm

    ! normalize |v(j-1)> to form |u(j)>
    call zdscal(h_size,beta_norm_inv,vec3,1)
    !! vec3 is |u(j)>

    ! |w(j)> = H|u(j)>
    call compute_hu(vec3,vec2,h_size)!! TODO: not implemented
    !! vec2 is |w(j)>

    alpha(j)=real(u_dot_v_complex(vec2,vec3,h_size)) 
    beta(j)=beta_norm

    ! |v(j)> = |w(j)> - alpha(j)*|u(j)> - beta(j)*|u(j-1)>
    do l=1,h_size
      vec_tmp=vec2(l)-alpha(j)*vec3(l)-beta(j)*vec1(l)
      vec1(l)=vec3(l)
      vec3(l)=vec_tmp
    enddo
    !! vec1 is |u(j)>
    !! vec3 is |v(j)>
  enddo

end


subroutine compute_hu(vec1,vec2,h_size)
  implicit none
  integer, intent(in)     :: h_size
  complex*16, intent(in)  :: vec1(h_size)
  complex*16, intent(out) :: vec2(h_size)
  complex*16 :: vec1_tmp(h_size)
  integer :: i
  BEGIN_DOC
  ! |vec2> = H|vec1>
  !
  ! TODO: implement
  ! maybe reuse parts of H_S2_u_0_nstates_{openmp,zmq}?
  END_DOC

  vec1_tmp(1:h_size) = vec1(1:h_size)
  call h_u_0_openmp(vec2,vec1_tmp,h_size)

  do i=1,h_size
    if (cdabs(vec1_tmp(i) - vec1(i)).gt.1.d-6) then
      print*,'ERROR: vec1 was changed by h_u_0_openmp'
    endif
  enddo
end

subroutine compute_hu2(vec1,vec2,h_size)
  implicit none
  integer, intent(in)     :: h_size
  complex*16, intent(in)  :: vec1(h_size)
  complex*16, intent(out) :: vec2(h_size)
  complex*16, allocatable :: u_tmp(:,:), s_tmp(:,:),v_tmp(:,:)
  integer :: i
  BEGIN_DOC
  ! |vec2> = H|vec1>
  !
  ! TODO: implement
  ! maybe reuse parts of H_S2_u_0_nstates_{openmp,zmq}?
  END_DOC

  allocate(u_tmp(1,h_size),s_tmp(1,h_size),v_tmp(1,h_size))

  u_tmp(1,1:h_size) = vec1(1:h_size)
  call h_s2_u_0_nstates_openmp(v_tmp,s_tmp,u_tmp,1,h_size)

  do i=1,h_size
    if (cdabs(u_tmp(1,i) - vec1(i)).gt.1.d-6) then
      print*,'ERROR: vec1 was changed by h_u_0_openmp'
    endif
  enddo
  vec2(1:h_size)=v_tmp(1,1:h_size)
  deallocate(u_tmp,v_tmp,s_tmp)
end



subroutine diag_lanczos_vals_vecs(alpha, beta, nlanc, vals, vecs, sze)
  implicit none
  BEGIN_DOC
  ! diagonalization of tridiagonal form of H
  ! this returns eigenvalues and eigenvectors in vals,vecs
  END_DOC
  integer, intent(in) :: nlanc,sze
  double precision, intent(in) :: alpha(sze), beta(sze)
  double precision, intent(out) :: vals(sze), vecs(sze,sze)
  double precision :: work(2*nlanc-2), beta_tmp(nlanc-1)
  integer :: i,info
  
  vals(1)=alpha(1)
  do i=2,nlanc
    vals(i)=alpha(i)
    beta_tmp(i-1)=beta(i)
  enddo

  call dstev('V', nlanc, vals, beta_tmp, vecs, sze, work, info)
  if (info.gt.0) then
    print *,'WARNING: diagonalization of tridiagonal form of H did not converge'
  else if (info.lt.0) then
    print *,'WARNING: argument to dstev had illegal value'
  endif
end

subroutine diag_lanczos_vals(alpha, beta, nlanc, vals, sze)
  implicit none
  BEGIN_DOC
  ! diagonalization of tridiagonal form of H
  ! this returns eigenvalues in vals
  END_DOC
  integer, intent(in) :: nlanc,sze
  double precision, intent(in) :: alpha(sze), beta(sze)
  double precision, intent(out) :: vals(sze)
  double precision :: work(1), beta_tmp(nlanc-1), vecs(1)
  integer :: i,info
  
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
end
