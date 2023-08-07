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

 if (do_ao_cholesky) then

    double precision, allocatable :: buffer_jj(:,:), buffer(:,:,:)
    allocate(buffer_jj(cholesky_mo_num,mo_num), buffer(mo_num,mo_num,mo_num))
    do j=1,mo_num
      buffer_jj(:,j) = cholesky_mo_transp(:,j,j)
    enddo

    call dgemm('T','N', mo_num*mo_num,mo_num,cholesky_mo_num, 1.d0, &
        cholesky_mo_transp, cholesky_mo_num, &
        buffer_jj, cholesky_mo_num, 0.d0, &
        buffer, mo_num*mo_num)

    do k = 1, mo_num
      do i = 1, mo_num
        do j = 1, mo_num
          big_array_coulomb_integrals(j,i,k) = buffer(i,k,j)
        enddo
      enddo
    enddo
    deallocate(buffer_jj)

    allocate(buffer_jj(mo_num,mo_num))

    do j = 1, mo_num

      call dgemm('T','N',mo_num,mo_num,cholesky_mo_num, 1.d0, &
        cholesky_mo_transp(1,1,j), cholesky_mo_num, &
        cholesky_mo_transp(1,1,j), cholesky_mo_num, 0.d0, &
        buffer_jj, mo_num)

      do k=1,mo_num
        do i=1,mo_num
          big_array_exchange_integrals(j,i,k) = buffer_jj(i,k)
       enddo
     enddo
    enddo

    deallocate(buffer_jj)

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

