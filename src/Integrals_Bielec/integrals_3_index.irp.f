 BEGIN_PROVIDER [complex*16, big_array_coulomb_integrals, (mo_tot_num,mo_tot_num, mo_tot_num)]
&BEGIN_PROVIDER [complex*16, big_array_exchange_integrals,(mo_tot_num,mo_tot_num, mo_tot_num)]
 implicit none
 integer :: i,j,k,l,kk,dk,di,k0
 complex*16 :: get_mo_bielec_integral
 complex*16 :: integral
 

  big_array_coulomb_integrals = (0.d0,0.d0)
  big_array_exchange_integrals = (0.d0,0.d0)
  do kk = 0, num_kpts-1
    k0 = kk*mo_tot_num_per_kpt
    do dk = 1, mo_tot_num_per_kpt
      k = dk + k0
      do di = dk, mo_tot_num_per_kpt
        i = di + k0
        do j = 1, mo_tot_num
          l = j
          integral = get_mo_bielec_integral(i,j,k,l,mo_integrals_map)
          big_array_coulomb_integrals(j,i,k) = integral
          big_array_coulomb_integrals(j,k,i) = conjg(integral)
          integral = get_mo_bielec_integral(i,j,l,k,mo_integrals_map)
          big_array_exchange_integrals(j,i,k) = integral
          big_array_exchange_integrals(j,k,i) = conjg(integral)
        enddo
      enddo
    enddo
  enddo

END_PROVIDER 
