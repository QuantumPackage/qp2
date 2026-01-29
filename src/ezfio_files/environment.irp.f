BEGIN_PROVIDER [ integer, nthreads_pt2 ]
 implicit none
 BEGIN_DOC
 ! Number of threads for Davidson
 END_DOC
 nthreads_pt2 = nproc
 character*(32) :: env
 call getenv('QP_NTHREADS_PT2',env)
 if (trim(env) /= '') then
   call lock_io()
   read(env,*) nthreads_pt2
   call unlock_io()
   call write_int(6,nthreads_pt2,'Target number of threads for PT2')
 endif
END_PROVIDER

BEGIN_PROVIDER [ integer, nthreads_davidson ]
 implicit none
 BEGIN_DOC
 ! Number of threads for Davidson
 END_DOC
 nthreads_davidson = nproc
 character*(32) :: env
 call getenv('QP_NTHREADS_DAVIDSON',env)
 if (trim(env) /= '') then
   call lock_io
   read(env,*) nthreads_davidson
   call unlock_io
   call write_int(6,nthreads_davidson,'Target number of threads for <Psi|H|Psi>')
 endif
END_PROVIDER

BEGIN_PROVIDER [ integer, nproc ]
  use omp_lib
  implicit none
  BEGIN_DOC
  ! Number of default OpenMP threads
  END_DOC

  nproc = 1
  !$OMP PARALLEL
  !$OMP MASTER
  !$ nproc = omp_get_num_threads()
  !$OMP END MASTER
  !$OMP END PARALLEL
END_PROVIDER


BEGIN_PROVIDER [ integer, nproc_max ]
  implicit none
  BEGIN_DOC
! Max number of threads
  END_DOC
  nproc_max = nproc
  nproc_max = max(nproc_max, nthreads_pt2)
  nproc_max = max(nproc_max, nthreads_davidson)
  nproc_max += 1 ! For collector
END_PROVIDER

