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
  END_DOC
  delta_omega=(omega_max-omega_min)/n_omega
END_PROVIDER

BEGIN_PROVIDER [ double precision, spectral_lanczos(n_omega) ]
  implicit none
  BEGIN_DOC
  ! spectral density A(omega) calculated from lanczos
  END_DOC

  integer :: i
  double precision :: omega_i
  complex*16 :: z_i
  
  do i=1,n_omega
    omega_i = omega_min + (i-1)*delta_omega
    z_i = dcmplx(omega_i,gf_epsilon)
    spectral_lanczos(i) = spec_lanc(n_lanczos_iter,alpha_lanczos,beta_lanczos,z_i)
  enddo

END_PROVIDER

double precision function spec_lanc(n_lanc_iter,alpha,beta,z)
  include 'constants.include.F'
  implicit none
  integer, intent(in) :: n_lanc_iter
  double precision, intent(in) :: alpha(n_lanc_iter), beta(n_lanc_iter)
  complex*16, intent(in) :: z

  complex*16 bigAj2,bigAj1,bigAj0
  complex*16 bigBj2,bigBj1,bigBj0
  
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
  double precision :: vec1norm2, beta_norm, beta_norm_inv
  complex*16, allocatable :: vec1(:), vec2(:), vec3(:)
  complex*16 :: vec_tmp
  double precision, external :: dznrm2

  integer :: i,j,l

  BEGIN_DOC
  ! lanczos
  ! n_lanc_iter is number of lanczos iterations
  !
  END_DOC
  h_size = N_det

  allocate(vec1(h_size),  &
           vec2(h_size),  &
           vec3(h_size))  
  
  do i=1,h_size
    vec1(i)=u1(i)
  enddo
!!  changed so that this is done externally
!!  ! construct initial lanczos vector
!!  ! TODO: why is it chosen this way?
!!  ! TODO: for more flexibility, maybe make u1 intent(in) and define this elsewhere
!!  vec1norm2=0.d0
!!  do i=1,h_size
!!    vec1(i)=1.d0/(dble(i))**0.1d0
!!  !  vec1norm2=vec1norm2+vec1(i)**2
!!    vec1norm2=vec1norm2+cdabs(vec1(i))**2
!!  enddo
!!
!!  ! normalize |vec1>, copy to |u1> and keep copy |vec1>
!!  do i=1,h_size
!!    vec1(i)=vec1(i)/dsqrt(vec1norm2)
!!    u1(i)=vec1(i)
!!  enddo

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
      vec3(l)=vec_temp
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
  BEGIN_DOC
  ! TODO: implement
  ! maybe reuse parts of H_S2_u_0_nstates_{openmp,zmq}?
  END_DOC

end



subroutine h_s2_u_0(e_0,u_0,n,keys_tmp,Nint,N_st,sze)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes e_0 = <u_0|H|u_0>/<u_0|u_0>
  !
  ! n : number of determinants
  !
  END_DOC
  integer, intent(in)            :: n,Nint, N_st, sze
  double precision, intent(out)  :: e_0(N_st)
  complex*16, intent(inout) :: u_0(sze,N_st)
  integer(bit_kind),intent(in)   :: keys_tmp(Nint,2,n)
  
  complex*16, allocatable  :: v_0(:,:), s_0(:,:), u_1(:,:)
  double precision               :: u_dot_u_complex,diag_H_mat_elem
  complex*16               :: u_dot_v_complex
  integer                        :: i,j

  if ((sze > 100000).and.distributed_davidson) then
    allocate (v_0(sze,N_states_diag),s_0(sze,N_states_diag), u_1(sze,N_states_diag))
    u_1(1:sze,1:N_states) = u_0(1:sze,1:N_states) 
    u_1(1:sze,N_states+1:N_states_diag) = 0.d0
    call H_S2_u_0_nstates_zmq(v_0,s_0,u_1,N_states_diag,sze)
    deallocate(u_1)
  else
    allocate (v_0(sze,N_st),s_0(sze,N_st))
    call H_S2_u_0_nstates_openmp(v_0,s_0,u_0,N_st,sze)
  endif
  double precision :: norm
  do i=1,N_st
    norm = u_dot_u_complex(u_0(1,i),n)
    if (norm /= 0.d0) then
      e_0(i) = real(u_dot_v_complex(v_0(1,i),u_0(1,i),n))
    else
      e_0(i) = 0.d0
    endif
  enddo
  deallocate (s_0, v_0)
end

subroutine u_0_H_u_0(e_0,u_0,n,keys_tmp,Nint,N_st,sze)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes e_0 = <u_0|H|u_0>/<u_0|u_0>
  !
  ! n : number of determinants
  !
  END_DOC
  integer, intent(in)            :: n,Nint, N_st, sze
  double precision, intent(out)  :: e_0(N_st)
  complex*16, intent(inout) :: u_0(sze,N_st)
  integer(bit_kind),intent(in)   :: keys_tmp(Nint,2,n)
  
  complex*16, allocatable  :: v_0(:,:), s_0(:,:), u_1(:,:)
  double precision               :: u_dot_u_complex,diag_H_mat_elem
  complex*16               :: u_dot_v_complex
  integer                        :: i,j

  if ((sze > 100000).and.distributed_davidson) then
    allocate (v_0(sze,N_states_diag),s_0(sze,N_states_diag), u_1(sze,N_states_diag))
    u_1(1:sze,1:N_states) = u_0(1:sze,1:N_states) 
    u_1(1:sze,N_states+1:N_states_diag) = 0.d0
    call H_S2_u_0_nstates_zmq(v_0,s_0,u_1,N_states_diag,sze)
    deallocate(u_1)
  else
    allocate (v_0(sze,N_st),s_0(sze,N_st))
    call H_S2_u_0_nstates_openmp(v_0,s_0,u_0,N_st,sze)
  endif
  double precision :: norm
  do i=1,N_st
    norm = u_dot_u_complex(u_0(1,i),n)
    if (norm /= 0.d0) then
      e_0(i) = real(u_dot_v_complex(v_0(1,i),u_0(1,i),n))
    else
      e_0(i) = 0.d0
    endif
  enddo
  deallocate (s_0, v_0)
end

