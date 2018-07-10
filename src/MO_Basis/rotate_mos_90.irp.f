program rotate_mos_90
  implicit none
  integer :: i,j
  double precision :: theta1,costheta,sintheta
  double precision :: x1,y1,x2,y2
  print*,'rotating MOs by 90 degrees'
  theta1 = dacos(-1.d0)/2.d0
  costheta = dcos(theta1)
  sintheta = dsin(theta1)
  do i=1,ao_num
    do j=1,mo_num
      x1 = real(mo_coef(i,j))
      y1 = imag(mo_coef(i,j))
      x2 = costheta*x1 - sintheta*y1
      y2 = sintheta*x1 + costheta*y1
      mo_coef(i,j) = dcmplx(x2,y2)
    enddo
  enddo
  call save_mos

end
