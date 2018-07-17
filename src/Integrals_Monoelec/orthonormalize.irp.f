subroutine orthonormalize_mos
  implicit none
  integer :: m,p,s
  m = size(mo_coef,1)
  p = size(mo_overlap,1)
  call ortho_lowdin_complex(mo_overlap,p,mo_tot_num,mo_coef,m,ao_num)
  mo_label = 'Orthonormalized'
  SOFT_TOUCH mo_coef mo_label 
end


