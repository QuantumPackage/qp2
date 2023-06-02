BEGIN_PROVIDER [ integer, cholesky_ao_num_guess ]
 implicit none
 BEGIN_DOC
 ! Number of Cholesky vectors in AO basis
 END_DOC

 cholesky_ao_num_guess = ao_num*ao_num / 2
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
 integer :: fd, i,j,k,l,m,rank
 double precision, pointer :: ao_integrals(:,:,:,:)
 double precision, external :: ao_two_e_integral

 ! Store AO integrals in a memory mapped file
 call mmap(trim(ezfio_work_dir)//'ao_integrals', &
   (/ int(ao_num,8), int(ao_num,8), int(ao_num,8), int(ao_num,8) /), &
   8, fd, .False., ptr)
 call c_f_pointer(ptr, ao_integrals, (/ao_num, ao_num, ao_num, ao_num/))

 print*, 'Providing the AO integrals (Cholesky)'
 call wall_time(wall_1)
 call cpu_time(cpu_1)

 ao_integrals = 0.d0

 double precision :: integral, cpu_1, cpu_2, wall_1, wall_2
 logical, external :: ao_two_e_integral_zero
  double precision, external :: get_ao_two_e_integral

 if (read_ao_two_e_integrals) then
   PROVIDE ao_two_e_integrals_in_map

   !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,l, integral, wall_2)
   do m=0,9
     do l=1+m,ao_num,10
       !$OMP DO SCHEDULE(dynamic)
       do j=1,l
         do k=1,ao_num
           do i=1,min(k,j)
             if (ao_two_e_integral_zero(i,j,k,l)) cycle
             integral = get_ao_two_e_integral(i,j,k,l, ao_integrals_map)
             ao_integrals(i,k,j,l) = integral
             ao_integrals(k,i,j,l) = integral
             ao_integrals(i,k,l,j) = integral
             ao_integrals(k,i,l,j) = integral
             ao_integrals(j,l,i,k) = integral
             ao_integrals(j,l,k,i) = integral
             ao_integrals(l,j,i,k) = integral
             ao_integrals(l,j,k,i) = integral
           enddo
         enddo
       enddo
       !$OMP END DO NOWAIT
     enddo
     !$OMP MASTER
     call wall_time(wall_2)
     print '(I10,'' %  in'', 4X, F10.2, '' s.'')', (m+1) * 10, wall_2-wall_1
     !$OMP END MASTER
   enddo
   !$OMP END PARALLEL

 else

   !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,l, integral, wall_2)
   do m=0,9
     do l=1+m,ao_num,10
       !$OMP DO SCHEDULE(dynamic)
       do j=1,l
         do k=1,ao_num
           do i=1,min(k,j)
             if (ao_two_e_integral_zero(i,j,k,l)) cycle
             integral = ao_two_e_integral(i,k,j,l)
             ao_integrals(i,k,j,l) = integral
             ao_integrals(k,i,j,l) = integral
             ao_integrals(i,k,l,j) = integral
             ao_integrals(k,i,l,j) = integral
             ao_integrals(j,l,i,k) = integral
             ao_integrals(j,l,k,i) = integral
             ao_integrals(l,j,i,k) = integral
             ao_integrals(l,j,k,i) = integral
           enddo
         enddo
       enddo
       !$OMP END DO NOWAIT
     enddo
     !$OMP MASTER
     call wall_time(wall_2)
     print '(I10,'' %  in'', 4X, F10.2, '' s.'')', (m+1) * 10, wall_2-wall_1
     !$OMP END MASTER
   enddo
   !$OMP END PARALLEL

   call wall_time(wall_2)
   call cpu_time(cpu_2)
   print*, 'AO integrals provided:'
   print*, ' cpu  time :',cpu_2 - cpu_1, 's'
   print*, ' wall time :',wall_2 - wall_1, 's  ( x ', (cpu_2-cpu_1)/(wall_2-wall_1+tiny(1.d0)), ' )'

 endif

 ! Call Lapack
 cholesky_ao_num = cholesky_ao_num_guess
 call pivoted_cholesky(ao_integrals, cholesky_ao_num, ao_cholesky_threshold, ao_num*ao_num, cholesky_ao)
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

