BEGIN_PROVIDER [ integer, nthreads_pt2 ]
 implicit none
 BEGIN_DOC
 ! Number of threads for Davidson
 END_DOC
 nthreads_pt2 = nproc
 character*(32) :: env
 call getenv('QP_NTHREADS_PT2',env)
 if (trim(env) /= '') then
   read(env,*) nthreads_pt2
   call write_int(6,nthreads_pt2,'Target number of threads for PT2')
 endif
END_PROVIDER

