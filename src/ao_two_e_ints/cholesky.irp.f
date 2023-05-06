BEGIN_PROVIDER [ integer, cholesky_ao_num_guess ]
 implicit none
 BEGIN_DOC
 ! Number of Cholesky vectors in AO basis
 END_DOC

 integer :: i,j,k,l
 double precision :: xnorm0, x, integral
 double precision, external :: ao_two_e_integral

 cholesky_ao_num_guess = 0
 xnorm0 = 0.d0
 x = 0.d0
 do j=1,ao_num
   do i=1,ao_num
     integral = ao_two_e_integral(i,i,j,j)
     if (integral > ao_integrals_threshold) then
       cholesky_ao_num_guess += 1
     else
       x += integral
     endif
   enddo
 enddo
 print *, 'Cholesky decomposition of AO integrals'
 print *, '--------------------------------------'
 print *, ''
 print *, 'Estimated Error: ', x
 print *, 'Guess size: ', cholesky_ao_num_guess, '(', 100.d0*dble(cholesky_ao_num_guess)/dble(ao_num*ao_num), ' %)'

END_PROVIDER

 BEGIN_PROVIDER [ integer, cholesky_ao_num ]
&BEGIN_PROVIDER [ double precision, cholesky_ao, (ao_num, ao_num, cholesky_ao_num_guess) ]
 use mmap_module
 implicit none
 BEGIN_DOC
 ! Cholesky vectors in AO basis: (ik|a):
 ! <ij|kl> = (ik|jl) = sum_a (ik|a).(a|jl)
 END_DOC

 type(c_ptr) :: ptr
 integer :: fd, i,j,k,l, rank
 double precision, pointer :: ao_integrals(:,:,:,:)
 double precision, external :: ao_two_e_integral

 ! Store AO integrals in a memory mapped file
 call mmap(trim(ezfio_work_dir)//'ao_integrals', &
   (/ int(ao_num,8), int(ao_num,8), int(ao_num,8), int(ao_num,8) /), &
   8, fd, .False., ptr)
 call c_f_pointer(ptr, ao_integrals, (/ao_num, ao_num, ao_num, ao_num/))

 double precision :: integral
 logical, external :: ao_two_e_integral_zero
 !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,l, integral) SCHEDULE(dynamic)
 do l=1,ao_num
  do j=1,l
   do k=1,ao_num
    do i=1,k
     if (ao_two_e_integral_zero(i,j,k,l)) cycle
     integral = ao_two_e_integral(i,k,j,l)
     ao_integrals(i,k,j,l) = integral
     ao_integrals(k,i,j,l) = integral
     ao_integrals(i,k,l,j) = integral
     ao_integrals(k,i,l,j) = integral
    enddo
   enddo
  enddo
 enddo
 !$OMP END PARALLEL DO

 ! Call Lapack
 cholesky_ao_num = cholesky_ao_num_guess
 call pivoted_cholesky(ao_integrals, cholesky_ao_num, ao_integrals_threshold, ao_num*ao_num, cholesky_ao)
 print *, 'Rank: ', cholesky_ao_num, '(', 100.d0*dble(cholesky_ao_num)/dble(ao_num*ao_num), ' %)'

 ! Remove mmap
 double precision, external :: getUnitAndOpen
 call munmap( &
   (/ int(ao_num,8), int(ao_num,8), int(ao_num,8), int(ao_num,8) /), &
   8, fd, ptr)
 open(unit=99,file=trim(ezfio_work_dir)//'ao_integrals')
 close(99, status='delete')

END_PROVIDER

BEGIN_PROVIDER [ double precision, cholesky_ao_transp, (cholesky_ao_num, ao_num, ao_num) ]
 implicit none
 BEGIN_DOC
! Transposed of the Cholesky vectors in AO basis set
 END_DOC
 integer :: i,j,k
 do j=1,ao_num
  do i=1,ao_num
   do k=1,ao_num
    cholesky_ao_transp(k,i,j) = cholesky_ao(i,j,k)
   enddo
  enddo
 enddo
END_PROVIDER

