complex*16 function ao_bielec_integral(i,j,k,l)
  use map_module
  implicit none
  BEGIN_DOC
  !  integral of the AO basis <ik|jl> or (ij|kl)
  !     i(r1) j(r1) 1/r12 k(r2) l(r2)
  END_DOC

  integer,intent(in)             :: i,j,k,l
  complex*16                     :: get_ao_bielec_integral

  PROVIDE ao_bielec_integrals_in_map

  ao_bielec_integral = get_ao_bielec_integral(i,k,j,l,ao_integrals_map)
  return
end


BEGIN_PROVIDER [ logical, ao_bielec_integrals_in_map ]
  implicit none
  use f77_zmq
  use map_module
  BEGIN_DOC
  !  Map of Atomic integrals
  !     i(r1) j(r2) 1/r12 k(r1) l(r2)
  END_DOC
  
  PROVIDE read_ao_integrals disk_access_ao_integrals
  if (read_ao_integrals) then
    print*,'Reading the AO integrals'
    call map_load_from_disk(trim(ezfio_filename)//'/work/ao_ints',ao_integrals_map)
    print*, 'AO integrals provided'
    ao_bielec_integrals_in_map = .True.
    return
  else if (read_df_integral_array) then
    call ao_map_fill_from_df
  else
    print*,'ERROR: complex AO integrals must be provided'
    stop
  endif
  
END_PROVIDER
 
subroutine ao_map_fill_from_df

  n_integrals = 0
  do l = 1, ao_tot_num
    do j = 1, l
      call idx2_tri_int(j,l,jl) 
      do k = 1, l
        do i = 1, l
          if ((j.eq.l).and.(i.gt.k)) exit
          call idx2_tri_int(i,k,ik)
          if (ik.gt.jl) exit
          integral = (0.d0,0.d0)
          do mu = 1, df_tot_num
            integral += df_integral_array(i,k,mu) * df_integral_array(j,l,mu)
          enddo
          if (cdabs(integral) < ao_integrals_threshold) then
            cycle
          endif
          n_integrals += 1
          tmp_re = real(integral)
          tmp_im = imag(integral)
          call mo_bielec_integrals_index(i,j,k,l,tmp_idx1)
          call mo_bielec_integrals_index(k,l,i,j,tmp_idx2)
          if (tmp_idx1.eq.tmp_idx2) then
            ! there are mo_num^2 of these:
            ! is it worth accumulating the imaginary parts somewhere 
            ! in order to verify that they are actually zero?
            buffer_i1(n_integrals) = tmp_idx1
            buffer_i2(n_integrals) = tmp_idx1
            buffer_value1(n_integrals) = tmp_re
            buffer_value2(n_integrals) = 0.d0
          else if (tmp_idx1 .lt. tmp_idx2) then
            buffer_i1(n_integrals) = tmp_idx1
            buffer_i2(n_integrals) = tmp_idx2
            buffer_value1(n_integrals) = tmp_re
            buffer_value2(n_integrals) = tmp_im
          else
            buffer_i1(n_integrals) = tmp_idx2
            buffer_i2(n_integrals) = tmp_idx1
            buffer_value1(n_integrals) = tmp_re
            buffer_value2(n_integrals) = -tmp_im
          endif


          if (n_integrals == size_buffer) then
            call insert_into_ao_integrals_map(n_integrals,buffer_i1,buffer_value1,&
                real(ao_integrals_threshold,integral_kind))
            call insert_into_ao_integrals_map(n_integrals,buffer_i2,buffer_value2,&
                real(ao_integrals_threshold,integral_kind))
            n_integrals = 0
          endif
        enddo
      enddo
    enddo
  enddo
  call insert_into_ao_integrals_map(n_integrals,buffer_i1,buffer_value1,&
      real(ao_integrals_threshold,integral_kind))
  call insert_into_ao_integrals_map(n_integrals,buffer_i2,buffer_value2,&
      real(ao_integrals_threshold,integral_kind))

end

BEGIN_PROVIDER [ double precision, ao_bielec_integral_schwartz,(ao_num,ao_num)  ]
  implicit none
  BEGIN_DOC
  !  Needed to compute Schwartz inequalities
  END_DOC
  
  integer                        :: i,k
  double precision               :: cpu_1,cpu_2, wall_1, wall_2
  complex*16                     :: ao_bielec_integral
  
!  ao_bielec_integral_schwartz(1,1) = real(ao_bielec_integral(1,1,1,1))
  !$OMP PARALLEL DO PRIVATE(i,k)                                     &
      !$OMP DEFAULT(NONE)                                            &
      !$OMP SHARED (ao_num,ao_bielec_integral_schwartz)              &
      !$OMP SCHEDULE(dynamic)
  do i=1,ao_num
    do k=1,i
      ao_bielec_integral_schwartz(i,k) = dsqrt(real(ao_bielec_integral(i,k,k,i)))
      ao_bielec_integral_schwartz(k,i) = ao_bielec_integral_schwartz(i,k)
    enddo
  enddo
  !$OMP END PARALLEL DO

END_PROVIDER



