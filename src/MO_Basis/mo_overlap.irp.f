
BEGIN_PROVIDER [ complex*16, mo_overlap,(mo_tot_num,mo_tot_num)]
  implicit none
  call complex_ao_to_mo(ao_overlap, size(ao_overlap,1), &
                     mo_overlap, size(mo_overlap,1))
END_PROVIDER

