 use gpu
 BEGIN_PROVIDER [double precision, big_array_coulomb_integrals, (mo_num,mo_num, mo_num)]
&BEGIN_PROVIDER [double precision, big_array_exchange_integrals,(mo_num,mo_num, mo_num)]
 implicit none
 BEGIN_DOC
 ! big_array_coulomb_integrals(j,i,k)  = <ij|kj> = (ik|jj)
 !
 ! big_array_exchange_integrals(j,i,k) = <ij|jk> = (ij|kj)
 END_DOC
 integer :: i,j,k,l,a
 double precision :: get_two_e_integral
 double precision :: integral

 if (do_mo_cholesky) then

    type(gpu_double2) :: buffer_jj
    type(gpu_double3) :: buffer

    call gpu_allocate(buffer_jj,cholesky_mo_num,mo_num)
    call gpu_allocate(buffer,mo_num,mo_num,mo_num)
    do j=1,mo_num
      buffer_jj%f(:,j) = cholesky_mo_transp_d(0)%f(:,j,j)
    enddo

    call gpu_dgemm(blas_handle, 'T','N', mo_num*mo_num,mo_num,cholesky_mo_num, 1.d0, &
        cholesky_mo_transp_d(0)%f(1,1,1), cholesky_mo_num, &
        buffer_jj%f(1,1), cholesky_mo_num, 0.d0, &
        buffer%f(1,1,1), mo_num*mo_num)

    call gpu_synchronize()
    do j = 1, mo_num
      do k = 1, mo_num
        do i = 1, mo_num
          big_array_coulomb_integrals(j,i,k) = buffer%f(i,k,j)
        enddo
      enddo
    enddo
    call gpu_deallocate(buffer_jj)

    do j = 1, mo_num
      call gpu_dgemm(blas_handle,'T','N',mo_num,mo_num,cholesky_mo_num, 1.d0, &
        cholesky_mo_transp_d(0)%f(1,1,j), cholesky_mo_num, &
        cholesky_mo_transp_d(0)%f(1,1,j), cholesky_mo_num, 0.d0, &
        buffer%f(1,1,j), mo_num)
    enddo

    call gpu_synchronize()
    do j = 1, mo_num
      do k=1,mo_num
        do i=1,mo_num
          big_array_exchange_integrals(j,i,k) = buffer%f(i,k,j)
       enddo
     enddo
    enddo

    call gpu_deallocate(buffer)

 else

   do k = 1, mo_num
     do i = 1, mo_num
       do j = 1, mo_num
         l = j
         integral = get_two_e_integral(i,j,k,l,mo_integrals_map)
         big_array_coulomb_integrals(j,i,k) = integral
         l = j
         integral = get_two_e_integral(i,j,l,k,mo_integrals_map)
         big_array_exchange_integrals(j,i,k) = integral
       enddo
     enddo
   enddo

 endif

END_PROVIDER

