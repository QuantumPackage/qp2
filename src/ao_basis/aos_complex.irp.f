BEGIN_PROVIDER [ integer, ao_kpt_num ]
 implicit none
 BEGIN_DOC
 ! number of aos per kpt.
 END_DOC
 ao_kpt_num = ao_num/kpt_num
END_PROVIDER
