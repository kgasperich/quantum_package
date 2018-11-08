program print_integrals
  implicit none
  
  integer :: iunit
  integer :: getunitandopen

  integer ::i,j,k,l
  !double precision :: integral
  complex*16 :: integral

  iunit = getunitandopen('kinetic_mo','w')
  do i=1,mo_tot_num
    do j=1,i
      integral = mo_kinetic_integral(i,j)
      if (cdabs(integral) > mo_integrals_threshold) then
        write(iunit,*) i,j, real(integral), imag(integral)
      endif
    enddo
  enddo
  close(iunit)

  iunit = getunitandopen('overlap_mo','w')
  do i=1,mo_tot_num
    do j=1,i
      integral = mo_overlap(i,j)
      if (cdabs(integral) > mo_integrals_threshold) then
        write(iunit,*) i,j, real(integral), imag(integral)
      endif
    enddo
  enddo
  close(iunit)

  iunit = getunitandopen('nuclear_mo','w')
  do i=1,mo_tot_num
    do j=1,i
      integral = mo_nucl_elec_integral(i,j)
      if (cdabs(integral) > mo_integrals_threshold) then
        write(iunit,*) i,j, real(integral), imag(integral)
      endif
    enddo
  enddo
  close(iunit)

  if (use_df_mo) then
  iunit = getunitandopen('df_mo','w')
  do i=1,mo_num_per_kpt
    do j=1,mo_num_per_kpt
      do k=1,df_num
        do l=1,num_kpt_pairs
          integral = df_mo_integral_array(i,j,k,l)
          if (cdabs(integral) > mo_integrals_threshold) then
            write(iunit,*) i,j,k,l, real(integral), imag(integral)
          endif
        enddo
      enddo
    enddo
  enddo
  close(iunit)
  endif


  PROVIDE mo_bielec_integrals_in_map
  iunit = getunitandopen('bielec_mo_complex','w')
  integer :: ki,kj,kk,kl, ii,ij,ik,il, kjkl2,kikk2,jl2,ik2
  integer*8 :: intcount
!  integer :: ik,jl
  complex*16 :: integral12
  complex*16 :: get_mo_bielec_integral

  intcount=0
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
                intcount += 1
                integral12 = get_mo_bielec_integral(i,j,k,l,mo_integrals_map)
                if (cdabs(integral12) > mo_integrals_threshold) then
                  write (iunit,'(4(I6,X),2(E25.15,X))') i,j,k,l, real(integral12),  imag(integral12)
                endif
              enddo !ii
            enddo !ik
          enddo !ij
        enddo !il
      enddo !kk
    enddo !kj
  enddo !kl

  close(iunit)
  print*,'number of nonzero ints by kpt symmetry: ',intcount
end
