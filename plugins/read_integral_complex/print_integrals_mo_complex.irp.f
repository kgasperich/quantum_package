program print_integrals
  implicit none
  
  integer :: iunit
  integer :: getunitandopen

  integer ::i,j,k,l
  double precision :: integral

  iunit = getunitandopen('kinetic_mo','w')
  do i=1,mo_tot_num
    do j=1,mo_tot_num
      write(iunit,*) i,j, mo_kinetic_integral(i,j)
    enddo
  enddo
  close(iunit)

  iunit = getunitandopen('overlap_mo','w')
  do i=1,mo_tot_num
    do j=1,mo_tot_num
      write(iunit,*) i,j, mo_overlap(i,j)
    enddo
  enddo
  close(iunit)

  iunit = getunitandopen('nuclear_mo','w')
  do i=1,mo_tot_num
    do j=1,mo_tot_num
      write(iunit,*) i,j, mo_nucl_elec_integral(i,j)
    enddo
  enddo
  close(iunit)


  PROVIDE mo_bielec_integrals_in_map
  iunit = getunitandopen('bielec_mo_complex','w')
  integer :: ik,jl
  complex*16 :: integral12
  complex*16 :: get_mo_bielec_integral
  do l=1,mo_tot_num
    do j=1,l
      call idx2_tri_int(j,l,jl)
      do k=1,l
        do i=1,l
          if (j==l .and. i>k) then
            exit
          endif
          call idx2_tri_int(i,k,ik)
          if (ik > jl) then
            exit
          endif
          integral12 = get_mo_bielec_integral(i,j,k,l,mo_integrals_map)
          if (cdabs(integral12) > mo_integrals_threshold) then
            write (iunit,'(4(I6,X),2(E25.15,X))') i,j,k,l, real(integral12),  imag(integral12)
          endif
        enddo
      enddo
    enddo
  enddo

  close(iunit)
end
