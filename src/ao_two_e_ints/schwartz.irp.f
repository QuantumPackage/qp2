BEGIN_TEMPLATE 

BEGIN_PROVIDER [ double precision , ao_two_e_integral$_erf_schwartz, (ao_num, ao_num) ]
 implicit none
 integer :: i,j
 double precision ::  get_ao_two_e_integral$_erf
 double precision ::  get_ao$_erf_integ_chol
 if(do_ao_cholesky)then
  do i = 1, ao_num
   do j = 1, ao_num
    ao_two_e_integral$_erf_schwartz(j,i) = get_ao$_erf_integ_chol(i,j,i,i)
   enddo
  enddo
 else
  PROVIDE ao_two_e_integrals$_erf_in_map
  do i = 1, ao_num
   do j = 1, ao_num
    ao_two_e_integral$_erf_schwartz(j,i) = get_ao_two_e_integral$_erf(i, i, j, j) 
   enddo
  enddo
 endif

END_PROVIDER 

SUBST [ _erf ]
  
;;
_erf;;
_cgtos;;

END_TEMPLATE 
