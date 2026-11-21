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
 PROVIDE mo_integrals_map

 if (do_mo_cholesky) then


    double precision, allocatable :: buffer(:,:,:)
    allocate(buffer(mo_num,mo_num,mo_num))

    if (gpu_num > 0) then
      type(gpu_double2) :: buffer_jj_d
      type(gpu_double3) :: buffer_d

      call gpu_allocate(buffer_jj_d,cholesky_mo_num,mo_num)
      call gpu_allocate(buffer_d,mo_num,mo_num,mo_num)
      !$OMP PARALLEL DO PRIVATE(j)
      do j=1,mo_num
        call gpu_copy(cholesky_mo_transp_d(0)%f(1,j,j), buffer_jj_d%f(1,j), cholesky_mo_num)
      enddo
      !$OMP END PARALLEL DO

      call gpu_dgemm(blas_handle, 'T','N', mo_num*mo_num,mo_num,cholesky_mo_num, 1.d0, &
          cholesky_mo_transp_d(0)%f(1,1,1), cholesky_mo_num, &
          buffer_jj_d%f(1,1), cholesky_mo_num, 0.d0, &
          buffer_d%f(1,1,1), mo_num*mo_num)

      call gpu_download(buffer_d,buffer)
      call gpu_deallocate(buffer_jj_d)

    else
      double precision, allocatable :: buffer_jj(:,:)
      allocate(buffer_jj(cholesky_mo_num,mo_num))
      !$OMP PARALLEL DO PRIVATE(j)
      do j=1,mo_num
        buffer_jj(1:cholesky_mo_num,j) = cholesky_mo_transp(1:cholesky_mo_num,j,j)
      enddo
      !$OMP END PARALLEL DO

      call dgemm('T','N', mo_num*mo_num,mo_num,cholesky_mo_num, 1.d0, &
          cholesky_mo_transp(1,1,1), cholesky_mo_num, &
          buffer_jj(1,1), cholesky_mo_num, 0.d0, &
          buffer(1,1,1), mo_num*mo_num)

      deallocate(buffer_jj)

    endif

    do j = 1, mo_num
      do k = 1, mo_num
        do i = 1, mo_num
          big_array_coulomb_integrals(j,i,k) = buffer(i,k,j)
        enddo
      enddo
    enddo

    if (gpu_num > 0) then
      do j = 1, mo_num
        call gpu_dgemm(blas_handle,'T','N',mo_num,mo_num,cholesky_mo_num, 1.d0, &
          cholesky_mo_transp_d(0)%f(1,1,j), cholesky_mo_num, &
          cholesky_mo_transp_d(0)%f(1,1,j), cholesky_mo_num, 0.d0, &
          buffer_d%f(1,1,j), mo_num)
      enddo

      call gpu_download(buffer_d,buffer)
      call gpu_deallocate(buffer_d)

    else

      do j = 1, mo_num
        call dgemm('T','N',mo_num,mo_num,cholesky_mo_num, 1.d0, &
          cholesky_mo_transp(1,1,j), cholesky_mo_num, &
          cholesky_mo_transp(1,1,j), cholesky_mo_num, 0.d0, &
          buffer(1,1,j), mo_num)
      enddo

    endif

    do j = 1, mo_num
      do k=1,mo_num
        do i=1,mo_num
          big_array_exchange_integrals(j,i,k) = buffer(i,k,j)
       enddo
     enddo
    enddo

 else

   real(integral_kind)            :: tmp
   integer(key_kind)              :: idx
   do k = 1, mo_num
     do i = 1, mo_num
       do j = 1, mo_num
         l = j
         call two_e_integrals_index(i,j,k,l,idx)
         call map_get(mo_integrals_map,idx,tmp)
         big_array_coulomb_integrals(j,i,k) = dble(tmp)
         l = j
         call two_e_integrals_index(i,j,l,k,idx)
         call map_get(mo_integrals_map,idx,tmp)
         big_array_exchange_integrals(j,i,k) = dble(tmp)
       enddo
     enddo
   enddo

 endif

END_PROVIDER

