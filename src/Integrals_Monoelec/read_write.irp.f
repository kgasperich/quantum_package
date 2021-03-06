 BEGIN_PROVIDER [ logical, read_ao_one_integrals ]
&BEGIN_PROVIDER [ logical, write_ao_one_integrals ]
&BEGIN_PROVIDER [ logical, read_mo_one_integrals ]
&BEGIN_PROVIDER [ logical, write_mo_one_integrals ]
   
   BEGIN_DOC
   ! One level of abstraction for disk_access_ao_one_integrals and disk_access_mo_one_integrals
   END_DOC
   implicit none
   
   if (disk_access_ao_one_integrals.EQ.'Read') then
     read_ao_one_integrals =  .True.
     write_ao_one_integrals = .False.
     
   else if  (disk_access_ao_one_integrals.EQ.'Write') then
     read_ao_one_integrals = .False.
     write_ao_one_integrals =  .True.
     
   else if (disk_access_ao_one_integrals.EQ.'None') then
     read_ao_one_integrals = .False.
     write_ao_one_integrals = .False.
     
   else
     print *, 'integrals_monoelec/disk_access_ao_one_integrals has a wrong type'
     stop 1
     
   endif
   
   if (disk_access_mo_one_integrals.EQ.'Read') then
     read_mo_one_integrals =  .True.
     write_mo_one_integrals = .False.
     
   else if  (disk_access_mo_one_integrals.EQ.'Write') then
     read_mo_one_integrals = .False.
     write_mo_one_integrals =  .True.
     
   else if (disk_access_mo_one_integrals.EQ.'None') then
     read_mo_one_integrals = .False.
     write_mo_one_integrals = .False.
     
   else
     print *, 'integrals_monoelec/disk_access_mo_one_integrals has a wrong type'
     stop 1
     
   endif
   
END_PROVIDER
