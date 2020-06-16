BEGIN_PROVIDER [ integer, ao_num_per_kpt ]
 implicit none
 BEGIN_DOC
 ! number of aos per kpt.
 END_DOC
 ao_num_per_kpt = ao_num/kpt_num
END_PROVIDER
