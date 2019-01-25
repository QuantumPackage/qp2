 BEGIN_PROVIDER [ integer, elec_num ]
&BEGIN_PROVIDER [ integer, elec_num_tab, (2)]

 implicit none
 BEGIN_DOC
 !  Numbers of alpha ("up") , beta ("down") and total electrons
 END_DOC
 PROVIDE ezfio_filename

 elec_num_tab(1) = elec_alpha_num
 elec_num_tab(2) = elec_beta_num
 elec_num = elec_alpha_num+elec_beta_num

END_PROVIDER

