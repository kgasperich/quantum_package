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
  integer(key_kind) :: i0,i1a,i1b,i2a,i2b
  !double precision :: get_mo_bielec_integral
  complex*16 :: integral12
  complex*16 :: get_mo_bielec_integral
  do l=1,mo_tot_num
    print*,'l=',l
    do k=1,mo_tot_num
      do j=l,mo_tot_num
        do i=max(j,k),mo_tot_num
          if (i.ne.j .or. l.le.k) then
            !ijkl is unique
    !        call bielec_integrals_index(i,j,k,l,i0)
            call bielec_integrals_index_2fold(i,j,k,l,i1a)
            call bielec_integrals_index_2fold(k,l,i,j,i1b)
            call bielec_integrals_index_2fold(i,l,k,j,i2a)
            call bielec_integrals_index_2fold(k,j,i,l,i2b)
            ! (i1a.eq.i1b) and (i1a.lt.i1b) will be separate when reading 
            ! and in full complex code
            print*,'getting integral from map'
            integral12 = get_mo_bielec_integral(i,j,k,l,mo_integrals_map)
            print*,'done getting integral from map'

            if (cdabs(integral12) > mo_integrals_threshold) then
              if (i1a.le.i1b) then
                write (iunit,'(4(I6,X),2(E25.15,X))') i,j,k,l, real(integral12),  imag(integral12)
                write (iunit,'(4(I6,X),2(E25.15,X))') k,l,i,j, real(integral12), -imag(integral12)
              else
                write (iunit,'(4(I6,X),2(E25.15,X))') k,l,i,j, real(integral12), -imag(integral12)
                write (iunit,'(4(I6,X),2(E25.15,X))') i,j,k,l, real(integral12),  imag(integral12)
              endif
              if (i2a.ne.i1a .and. i2a.ne.i1b) then
                if (i2a.le.i2b) then
                  write (iunit,'(4(I6,X),2(E25.15,X))') i,l,k,j, real(integral12), imag(integral12)
                else
                  write (iunit,'(4(I6,X),2(E25.15,X))') k,j,i,l, real(integral12), -imag(integral12)
                endif
              endif
            endif
          endif
        enddo
      enddo
    enddo
  enddo
  close(iunit)
end
