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

 BEGIN_PROVIDER [ double precision, alpha_lanczos, (n_lanczos_iter) ]
&BEGIN_PROVIDER [ double precision, beta_lanczos, (n_lanczos_iter) ]
  implicit none
  BEGIN_DOC
  ! provide alpha and beta for tridiagonal form of H
  END_DOC

  PROVIDE u1_lanczos
  call lanczos_h(n_lanczos_iter, alpha_lanczos, beta_lanczos, u1_lanczos) 

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

