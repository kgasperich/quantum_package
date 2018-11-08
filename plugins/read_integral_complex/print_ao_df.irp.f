program print_integrals

  PROVIDE ezfio_filename
  call ezfio_set_integrals_bielec_disk_access_df_ao_integral_array('Read')
  call run
end

subroutine run
  implicit none
  
  integer :: iunit
  integer :: getunitandopen

  integer :: kl,kj,mu,kjkl,l,j
  complex*16 :: integral


!  integer(key_kind) :: i0,i1a,i1b,i2a,i2b


  iunit = getunitandopen('df_ao_array','w')
  
  do kl=1,num_kpts
    do kj=1,kl
      call idx2_tri_int(kj,kl,kjkl)
      do mu=1,df_num
        do l=1,ao_num_per_kpt
!          write(iunit,'(I3,X,I3,X,I3,X,I3,X,A3,100(F16.7))'),kl,kj,kjkl,mu,' | ',df_ao_integral_array(l,:,mu,kjkl)
          do j=1,ao_num_per_kpt
            write(iunit,'(I3,X,I3,X,I3,X,I3,X,A3,100(F16.7))'),l,j,mu,kjkl,'  ',df_ao_integral_array(l,j,mu,kjkl)
          enddo 
        enddo
      enddo
    enddo
  enddo

  close(iunit)
end
