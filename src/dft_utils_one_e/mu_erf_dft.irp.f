BEGIN_PROVIDER [double precision, mu_erf_dft]
 implicit none
 BEGIN_DOC
! range separation parameter used in RS-DFT. It is set to mu_erf in order to be consistent with the two electrons integrals erf
 END_DOC
 mu_erf_dft = mu_erf

END_PROVIDER
