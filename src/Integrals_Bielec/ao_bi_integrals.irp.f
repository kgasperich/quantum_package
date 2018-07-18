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
  endif
  
END_PROVIDER
 
BEGIN_PROVIDER [ double precision, ao_bielec_integral_schwartz,(ao_num,ao_num)  ]
  implicit none
  BEGIN_DOC
  !  Needed to compute Schwartz inequalities
  END_DOC
  
  integer                        :: i,k
  double precision               :: cpu_1,cpu_2, wall_1, wall_2
  complex*16                     :: ao_bielec_integral
  
  ao_bielec_integral_schwartz(1,1) = real(ao_bielec_integral(1,1,1,1))
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



