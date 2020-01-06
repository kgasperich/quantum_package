program intdump
  implicit none
  BEGIN_DOC
! Performs timing benchmarks on integral access
  END_DOC
  integer                        :: i,j,k,l,m
  integer*8                      :: ii,jj
  double precision               :: r, cpu
  integer*8                      :: cpu0, cpu1, count_rate, count_max
  complex*16                     :: c,get_mo_bielec_integral
  
  PROVIDE mo_bielec_integrals_in_map
  print *,  mo_tot_num, 'MOs'
  print *,  'num to print?'
  read *,m

  do l=1,m
    do k=1,m
      do j=1,m
        do i=1,m
          c = get_mo_bielec_integral(i,j,k,l,mo_integrals_map)
          write(6,'(4(I5),2(E25.15))')i,j,k,l,c
        enddo
      enddo
    enddo
  enddo

end
