BEGIN_PROVIDER [ complex*16, ao_nucl_elec_integral, (ao_num,ao_num)]
   BEGIN_DOC
   ! interaction nuclear electron
   END_DOC
   implicit none
   double precision               :: alpha, beta, gama, delta
   integer                        :: num_A,num_B
   double precision               :: A_center(3),B_center(3),C_center(3)
   integer                        :: power_A(3),power_B(3)
   integer                        :: i,j,k,l,n_pt_in,m
   double precision               :: overlap_x,overlap_y,overlap_z,overlap,dx,NAI_pol_mult
   
   if (read_ao_one_integrals) then
    call read_one_e_integrals_complex('ao_ne_integral', ao_nucl_elec_integral,      &
            size(ao_nucl_elec_integral,1), size(ao_nucl_elec_integral,2))
     print *,  'AO N-e integrals read from disk'
   else
     print *, 'complex AO N-e integrals must be provided'
     stop
   endif
   
END_PROVIDER

 BEGIN_PROVIDER [ complex*16, ao_nucl_elec_integral_per_atom, (ao_num,ao_num,nucl_num)]
 BEGIN_DOC
! ao_nucl_elec_integral_per_atom(i,j,k) = -<AO(i)|1/|r-Rk|AO(j)> 
! where Rk is the geometry of the kth atom
 END_DOC
  print *, 'ao_nucl_elec_integral_per_atom not implemented for k-points'
  stop
END_PROVIDER



end


