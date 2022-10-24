use bitmasks
use f77_zmq


! ---

BEGIN_PROVIDER [ integer, nthreads_davidson ]
 implicit none
 BEGIN_DOC
 ! Number of threads for Davidson
 END_DOC
 nthreads_davidson = nproc
 character*(32) :: env
 call getenv('QP_NTHREADS_DAVIDSON',env)
 if (trim(env) /= '') then
   read(env,*) nthreads_davidson
   call write_int(6,nthreads_davidson,'Target number of threads for <Psi|H|Psi>')
 endif
END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, threshold_davidson_pt2 ]
 implicit none
 BEGIN_DOC
 ! Threshold of Davidson's algorithm, using PT2 as a guide
 END_DOC
 threshold_davidson_pt2 = threshold_davidson

END_PROVIDER

! ---

