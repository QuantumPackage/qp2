BEGIN_PROVIDER [ double precision, ao_cart_integrals_threshold]
 implicit none
 ao_cart_integrals_threshold = ao_integrals_threshold

END_PROVIDER 

BEGIN_PROVIDER [ double precision, ao_cart_cholesky_threshold]
 implicit none
 ao_cart_cholesky_threshold = ao_cholesky_threshold

END_PROVIDER 

BEGIN_PROVIDER [ logical, do_ao_cart_cholesky]
 implicit none
 do_ao_cart_cholesky = do_ao_cholesky

END_PROVIDER 
