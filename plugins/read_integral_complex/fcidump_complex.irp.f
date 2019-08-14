program fcidump_complex
  implicit none

  integer :: iunit
  integer :: getunitandopen
  integer :: i,j,k,l
!  integer :: ii(8), jj(8), kk(8),ll(8)
!  integer*8 :: m
  character*(2), allocatable :: A(:)

  iunit = getunitandopen('fcidump.txt','w')
! print system info (mos, symm, spin, Nelec)
  write(iunit,*) '&FCI NORB=', mo_tot_num, ', NELEC=', elec_num, &
   ', MS2=', (elec_alpha_num-elec_beta_num), ','
  allocate (A(mo_tot_num))
  A = '1,'
  write(iunit,*) 'ORBSYM=', (A(i), i=1,mo_tot_num) 
  write(iunit,*) 'ISYM=0,'
  write(iunit,*) '/'
  deallocate(A)

  complex*16 :: integral,integral1,integral2

  PROVIDE mo_bielec_integrals_in_map
  integer :: ki,kj,kk,kl, ii,ij,ik,il, kjkl2,kikk2,jl2,ik2
!  integer :: ik,jl
  complex*16 :: get_mo_bielec_integral

  do kl=1, num_kpts
    do kj=1, kl
      call idx2_tri_int(kj,kl,kjkl2)
      do kk=1,kl
        ki=kconserv(kl,kk,kj)
        if ((kl == kj) .and. (ki > kk)) cycle
        call idx2_tri_int(ki,kk,kikk2)
        if (kikk2 > kjkl2) cycle
        do il=1,mo_tot_num_per_kpt
          l=il+(kl-1)*mo_tot_num_per_kpt
          do ij=1,mo_tot_num_per_kpt
            j=ij+(kj-1)*mo_tot_num_per_kpt
            if (j>l) exit
            call idx2_tri_int(j,l,jl2)
            do ik=1,mo_tot_num_per_kpt
              k=ik+(kk-1)*mo_tot_num_per_kpt
              if (k>l) exit
              do ii=1,mo_tot_num_per_kpt
                i=ii+(ki-1)*mo_tot_num_per_kpt
                if ((j==l) .and. (i>k)) exit
                call idx2_tri_int(i,k,ik2)
                if (ik2 > jl2) exit
                integral = get_mo_bielec_integral(i,j,k,l,mo_integrals_map)
                if (cdabs(integral) > mo_integrals_threshold) then
!                  write (iunit,'(4(I6,X),2(E25.15,X))') i,j,k,l, real(integral),  imag(integral)
                  write (iunit,'(2(E25.15,X),4(I6,X))') real(integral), imag(integral),i,k,j,l
                endif
              enddo !ii
            enddo !ik
          enddo !ij
        enddo !il
      enddo !kk
    enddo !kj
  enddo !kl

!! from real version of fcidump (for reference)
!
!  do l=1,mo_tot_num
!   do k=1,mo_tot_num
!    do j=l,mo_tot_num
!     do i=k,mo_tot_num
!       if (i>=j) then
!          integral = get_mo_bielec_integral(i,j,k,l,mo_integrals_map)
!          if (dabs(integral) > mo_integrals_threshold) then 
!            print *, integral, i,k,j,l
!          endif
!       end if
!     enddo
!    enddo
!   enddo
!  enddo

  do j=1,mo_tot_num
    do i=j,mo_tot_num
      integral1 = mo_kinetic_integral(i,j)
      integral2 = mo_nucl_elec_integral(i,j)
      integral = integral1 + integral2
      if (cdabs(integral) > mo_integrals_threshold) then
        write(iunit,*) real(integral), imag(integral), i, j, 0, 0
      endif
    enddo
  enddo

  write(iunit,*) 0.d0, 0.d0, 0, 0, 0, 0
  
  close(iunit)

end
