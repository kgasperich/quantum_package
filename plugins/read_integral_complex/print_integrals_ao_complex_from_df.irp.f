program print_integrals

  PROVIDE ezfio_filename
  call ezfio_set_integrals_bielec_disk_access_ao_integrals('Write')
  call run
end

subroutine run
  implicit none
  
  integer :: iunit
  integer :: getunitandopen

  integer ::i,j,k,l
  complex*16 :: integral
!  call ao_map_fill_from_df

  PROVIDE ao_bielec_integrals_in_map
  iunit = getunitandopen('bielec_ao','w')

!  integer(key_kind) :: i0,i1a,i1b,i2a,i2b

  integer :: ik,jl
  complex*16 :: integral12
  complex*16 :: get_ao_bielec_integral
  do l=1,ao_num
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
          integral12 = get_ao_bielec_integral(i,j,k,l,ao_integrals_map)
          if (cdabs(integral12) > ao_integrals_threshold) then
            write (iunit,'(4(I6,X),2(E25.15,X))') i,j,k,l, real(integral12),  imag(integral12)
          endif
        enddo
      enddo
    enddo
  enddo

  close(iunit)
end
