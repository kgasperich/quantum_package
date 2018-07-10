program rotate_psi
  implicit none
  integer :: i,j
  double precision :: theta1,costheta,sintheta
  double precision :: x1,y1,x2,y2
  print *,  'theta?'
  read(*,*) theta1
  costheta = dcos(theta1)
  sintheta = dsin(theta1)
  do i=1,psi_det_size
    do j=1,N_states
      x1 = real(psi_coef(i,j))
      y1 = imag(psi_coef(i,j))
      x2 = costheta*x1 - sintheta*y1
      y2 = sintheta*x1 + costheta*y1
      psi_coef(i,j) = dcmplx(x2,y2)
    enddo
  enddo
  call save_wavefunction

end
