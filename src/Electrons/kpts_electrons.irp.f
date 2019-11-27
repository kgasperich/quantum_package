 BEGIN_PROVIDER [ integer, elec_alpha_num_per_kpt ]
&BEGIN_PROVIDER [ integer, elec_beta_num_per_kpt ]
&BEGIN_PROVIDER [ integer, elec_alpha_num_extra ]
&BEGIN_PROVIDER [ integer, elec_beta_num_extra ]

 implicit none
 BEGIN_DOC
 !  Numbers of alpha ("up") , beta ("down") and total electrons per k-point
 END_DOC
 PROVIDE ezfio_filename

 elec_alpha_num_per_kpt = elec_alpha_num / num_kpts
 elec_beta_num_per_kpt = elec_beta_num / num_kpts

 elec_alpha_num_extra = mod(elec_alpha_num,num_kpts)
 elec_beta_num_extra = mod(elec_beta_num,num_kpts)

END_PROVIDER

 BEGIN_PROVIDER [ integer, elec_num_per_kpt ]
&BEGIN_PROVIDER [ integer, elec_num_per_kpt_tab, (2)]
&BEGIN_PROVIDER [ integer, elec_num_extra ]
&BEGIN_PROVIDER [ integer, elec_num_extra_tab, (2)]

 implicit none
 BEGIN_DOC
 !  Numbers of alpha ("up") , beta ("down") and total electrons
 END_DOC
 PROVIDE ezfio_filename

 elec_num_per_kpt_tab(1) = elec_alpha_num_per_kpt
 elec_num_per_kpt_tab(2) = elec_beta_num_per_kpt
 elec_num_per_kpt = elec_alpha_num_per_kpt+elec_beta_num_per_kpt

 elec_num_extra_tab(1) = elec_alpha_num_extra
 elec_num_extra_tab(2) = elec_beta_num_extra
 elec_num_extra = elec_alpha_num_extra + elec_beta_num_extra


END_PROVIDER

