program print_integrals

  PROVIDE ezfio_filename
  call ezfio_set_integrals_monoelec_disk_access_ao_one_integrals('Read')
  call ezfio_set_utils_disk_access_ao_overlap_integrals('Read')
!  call ezfio_set_integrals_bielec_disk_access_ao_integrals('Read')
  call run
end

subroutine run
  implicit none
  
  integer :: iunit
  integer :: getunitandopen

  integer ::i,j,k,l
  complex*16 :: integral

  iunit = getunitandopen('kinetic_ao','w')
  do i=1,ao_num
    do j=1,i
      integral = ao_kinetic_integral(i,j)
      if (cdabs(integral) > ao_integrals_threshold) then
        write(iunit,*) i,j, real(integral), imag(integral)
      endif
    enddo
  enddo
  close(iunit)
  
  iunit = getunitandopen('overlap_ao','w')
  do i=1,ao_num
    do j=1,i
      integral = ao_overlap(i,j)
      if (cdabs(integral) > ao_integrals_threshold) then
        write(iunit,*) i,j, real(integral), imag(integral)
      endif
    enddo
  enddo
  close(iunit)
  
  iunit = getunitandopen('nuclear_ao','w')
  do i=1,ao_num
    do j=1,i
      integral = ao_nucl_elec_integral(i,j)
      if (cdabs(integral) > ao_integrals_threshold) then
        write(iunit,*) i,j, real(integral), imag(integral)
      endif
    enddo
  enddo
  close(iunit)

  if (use_df_ao) then
  iunit = getunitandopen('df_ao','w')
  do i=1,ao_num_per_kpt
    do j=1,ao_num_per_kpt
      do k=1,df_num
        do l=1,num_kpt_pairs
          integral = df_ao_integral_array(i,j,k,l)
          if (cdabs(integral) > ao_integrals_threshold) then
            write(iunit,*) i,j,k,l, real(integral), imag(integral)
          endif
        enddo
      enddo
    enddo
  enddo
  close(iunit)
  endif

  PROVIDE ao_bielec_integrals_in_map
  iunit = getunitandopen('bielec_ao','w')

!  integer(key_kind) :: i0,i1a,i1b,i2a,i2b

  integer :: ki,kj,kk,kl, ii,ij,ik,il, kjkl2,kikk2,jl2,ik2
  integer*8 :: intcount
  complex*16 :: integral12
  complex*16 :: get_ao_bielec_integral

  intcount=0
  do kl=1, num_kpts
    do kj=1, kl
      call idx2_tri_int(kj,kl,kjkl2)
      do kk=1,kl
        ki=kconserv(kl,kk,kj)
        if ((kl == kj) .and. (ki > kk)) cycle
        call idx2_tri_int(ki,kk,kikk2)
        if (kikk2 > kjkl2) cycle
        do il=1,ao_num_per_kpt
          l=il+(kl-1)*ao_num_per_kpt
          do ij=1,ao_num_per_kpt
            j=ij+(kj-1)*ao_num_per_kpt
            if (j>l) exit
            call idx2_tri_int(j,l,jl2)
            do ik=1,ao_num_per_kpt
              k=ik+(kk-1)*ao_num_per_kpt
              if (k>l) exit
              do ii=1,ao_num_per_kpt
                i=ii+(ki-1)*ao_num_per_kpt
                if ((j==l) .and. (i>k)) exit
                call idx2_tri_int(i,k,ik2)
                if (ik2 > jl2) exit
                intcount += 1
                integral12 = get_ao_bielec_integral(i,j,k,l,ao_integrals_map)
                if (cdabs(integral12) > ao_integrals_threshold) then
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
