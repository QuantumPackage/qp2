! F

subroutine gen_f_space(det,n1,n2,list1,list2,f)

  implicit none

  integer, intent(in)           :: n1,n2
  integer, intent(in)           :: list1(n1),list2(n2)
  integer(bit_kind), intent(in) :: det(N_int,2)
  double precision, intent(out) :: f(n1,n2)

  double precision, allocatable :: tmp_F(:,:)
  integer                       :: i1,i2,idx1,idx2

  allocate(tmp_F(mo_num,mo_num))

  call get_fock_matrix_spin(det,1,tmp_F)

  !$OMP PARALLEL &
  !$OMP SHARED(tmp_F,f,n1,n2,list1,list2) &
  !$OMP PRIVATE(idx1,idx2,i1,i2)&
  !$OMP DEFAULT(NONE)
  !$OMP DO collapse(1)
  do i2 = 1, n2
    do i1 = 1, n1
      idx2 = list2(i2)
      idx1 = list1(i1)
      f(i1,i2) = tmp_F(idx1,idx2)
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  deallocate(tmp_F)

end

! V

subroutine gen_v_space(n1,n2,n3,n4,list1,list2,list3,list4,v)

  implicit none

  integer, intent(in)           :: n1,n2,n3,n4
  integer, intent(in)           :: list1(n1),list2(n2),list3(n3),list4(n4)
  double precision, intent(out) :: v(n1,n2,n3,n4)

  integer                       :: i1,i2,i3,i4,idx1,idx2,idx3,idx4,k

  if (do_ao_cholesky) then
    double precision, allocatable :: buffer(:,:,:,:)
    double precision, allocatable :: v1(:,:,:), v2(:,:,:)
    allocate(v1(cholesky_mo_num,n1,n3), v2(cholesky_mo_num,n2,n4))
    allocate(buffer(n1,n3,n2,n4))

    call gen_v_space_chol(n1,n3,list1,list3,v1,cholesky_mo_num)
    call gen_v_space_chol(n2,n4,list2,list4,v2,cholesky_mo_num)

    call dgemm('T','N', n1*n3, n2*n4, cholesky_mo_num, 1.d0, &
         v1, cholesky_mo_num, &
         v2, cholesky_mo_num, 0.d0, buffer, n1*n3)

    deallocate(v1,v2)

    !$OMP PARALLEL DO PRIVATE(i1,i2,i3,i4)
    do i4 = 1, n4
      do i3 = 1, n3
        do i2 = 1, n2
          do i1 = 1, n1
            v(i1,i2,i3,i4) = buffer(i1,i3,i2,i4)
          enddo
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

  else
    double precision              :: get_two_e_integral

    PROVIDE mo_two_e_integrals_in_map

    !$OMP PARALLEL &
    !$OMP SHARED(n1,n2,n3,n4,list1,list2,list3,list4,v,mo_integrals_map) &
    !$OMP PRIVATE(i1,i2,i3,i4,idx1,idx2,idx3,idx4)&
    !$OMP DEFAULT(NONE)
    !$OMP DO collapse(3)
    do i4 = 1, n4
      do i3 = 1, n3
        do i2 = 1, n2
          do i1 = 1, n1
            idx4 = list4(i4)
            idx3 = list3(i3)
            idx2 = list2(i2)
            idx1 = list1(i1)
            v(i1,i2,i3,i4) = get_two_e_integral(idx1,idx2,idx3,idx4,mo_integrals_map)
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

  endif

end

subroutine gen_v_space_chol(n1,n3,list1,list3,v,ldv)

  implicit none

  integer, intent(in)           :: n1,n3,ldv
  integer, intent(in)           :: list1(n1),list3(n3)
  double precision, intent(out) :: v(ldv,n1,n3)

  integer                       :: i1,i3,idx1,idx3,k

  !$OMP PARALLEL DO PRIVATE(i1,i3,idx1,idx3,k)
  do i3=1,n3
    idx3 = list3(i3)
    do i1=1,n1
      idx1 = list1(i1)
      do k=1,cholesky_mo_num
        v(k,i1,i3) = cholesky_mo_transp(k,idx1,idx3)
      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO

end

! full

BEGIN_PROVIDER [double precision, cc_space_v, (mo_num,mo_num,mo_num,mo_num)]
  implicit none
  if (do_ao_cholesky) then
    integer                       :: i1,i2,i3,i4
    double precision, allocatable :: buffer(:,:,:)
    call set_multiple_levels_omp(.False.)
    !$OMP PARALLEL &
    !$OMP SHARED(cc_space_v,mo_num,cholesky_mo_transp,cholesky_mo_num) &
    !$OMP PRIVATE(i1,i2,i3,i4,k,buffer)&
    !$OMP DEFAULT(NONE)
    allocate(buffer(mo_num,mo_num,mo_num))
    !$OMP DO
    do i4 = 1, mo_num
      call dgemm('T','N', mo_num*mo_num, mo_num, cholesky_mo_num, 1.d0, &
           cholesky_mo_transp, cholesky_mo_num, &
           cholesky_mo_transp(1,1,i4), cholesky_mo_num, 0.d0, buffer, mo_num*mo_num)
      do i2 = 1, mo_num
        do i3 = 1, mo_num
          do i1 = 1, mo_num
            cc_space_v(i1,i2,i3,i4) = buffer(i1,i3,i2)
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO
    deallocate(buffer)
    !$OMP END PARALLEL
  else
    integer          :: i,j,k,l
    double precision :: get_two_e_integral

    PROVIDE mo_two_e_integrals_in_map

    !$OMP PARALLEL &
    !$OMP SHARED(cc_space_v,mo_num,mo_integrals_map) &
    !$OMP PRIVATE(i,j,k,l) &
    !$OMP DEFAULT(NONE)

    !$OMP DO collapse(3)
    do l = 1, mo_num
      do k = 1, mo_num
        do j = 1, mo_num
          do i = 1, mo_num
            cc_space_v(i,j,k,l) = get_two_e_integral(i,j,k,l,mo_integrals_map)
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
  endif

END_PROVIDER

! oooo

BEGIN_PROVIDER [double precision, cc_space_v_oooo, (cc_nOa, cc_nOa, cc_nOa, cc_nOa)]

  implicit none

  if (do_ao_cholesky) then

    integer :: i1, i2, i3, i4
    integer :: n1, n2, n3, n4
    
    n1 = size(cc_space_v_oooo,1)
    n2 = size(cc_space_v_oooo,2)
    n3 = size(cc_space_v_oooo,3)
    n4 = size(cc_space_v_oooo,4)

    double precision, allocatable :: buffer(:,:,:,:)
    allocate(buffer(n1,n3,n2,n4))

    call dgemm('T','N', n1*n3, n2*n4, cholesky_mo_num, 1.d0, &
         cc_space_v_oo_chol, cholesky_mo_num, &
         cc_space_v_oo_chol, cholesky_mo_num, 0.d0, buffer, n1*n3)

    !$OMP PARALLEL DO PRIVATE(i1,i2,i3,i4) COLLAPSE(2)
    do i4 = 1, n4
      do i3 = 1, n3
        do i2 = 1, n2
          do i1 = 1, n1
            cc_space_v_oooo(i1,i2,i3,i4) = buffer(i1,i3,i2,i4)
          enddo
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    deallocate(buffer)

  else
    call gen_v_space(cc_nOa,cc_nOa,cc_nOa,cc_nOa, cc_list_occ,cc_list_occ,cc_list_occ,cc_list_occ, cc_space_v_oooo)
  endif

END_PROVIDER

! vooo

BEGIN_PROVIDER [double precision, cc_space_v_vooo, (cc_nVa, cc_nOa, cc_nOa, cc_nOa)]

  implicit none

  if (do_ao_cholesky) then

    integer :: i1, i2, i3, i4
    integer :: n1, n2, n3, n4
    
    n1 = size(cc_space_v_vooo,1)
    n2 = size(cc_space_v_vooo,2)
    n3 = size(cc_space_v_vooo,3)
    n4 = size(cc_space_v_vooo,4)

    double precision, allocatable :: buffer(:,:,:,:)
    allocate(buffer(n1,n3,n2,n4))

    call dgemm('T','N', n1*n3, n2*n4, cholesky_mo_num, 1.d0, &
         cc_space_v_vo_chol, cholesky_mo_num, &
         cc_space_v_oo_chol, cholesky_mo_num, 0.d0, buffer, n1*n3)

    !$OMP PARALLEL DO PRIVATE(i1,i2,i3,i4) COLLAPSE(2)
    do i4 = 1, n4
      do i3 = 1, n3
        do i2 = 1, n2
          do i1 = 1, n1
            cc_space_v_vooo(i1,i2,i3,i4) = buffer(i1,i3,i2,i4)
          enddo
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    deallocate(buffer)

  else
    call gen_v_space(cc_nVa,cc_nOa,cc_nOa,cc_nOa, cc_list_vir,cc_list_occ,cc_list_occ,cc_list_occ, cc_space_v_vooo)
  endif

END_PROVIDER

! ovoo

BEGIN_PROVIDER [double precision, cc_space_v_ovoo, (cc_nOa, cc_nVa, cc_nOa, cc_nOa)]

  implicit none


  if (do_ao_cholesky) then

    integer :: i1, i2, i3, i4
    integer :: n1, n2, n3, n4
    
    n1 = size(cc_space_v_ovoo,1)
    n2 = size(cc_space_v_ovoo,2)
    n3 = size(cc_space_v_ovoo,3)
    n4 = size(cc_space_v_ovoo,4)

    !$OMP PARALLEL DO PRIVATE(i1,i2,i3,i4) COLLAPSE(2)
    do i4 = 1, n4
      do i3 = 1, n3
        do i2 = 1, n2
          do i1 = 1, n1
            cc_space_v_ovoo(i1,i2,i3,i4) = cc_space_v_vooo(i2,i1,i4,i3)
          enddo
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

  else
    call gen_v_space(cc_nOa,cc_nVa,cc_nOa,cc_nOa, cc_list_occ,cc_list_vir,cc_list_occ,cc_list_occ, cc_space_v_ovoo)
  endif

END_PROVIDER

! oovo

BEGIN_PROVIDER [double precision, cc_space_v_oovo, (cc_nOa, cc_nOa, cc_nVa, cc_nOa)]

  implicit none

  if (do_ao_cholesky) then

    integer :: i1, i2, i3, i4
    integer :: n1, n2, n3, n4
    
    n1 = size(cc_space_v_oovo,1)
    n2 = size(cc_space_v_oovo,2)
    n3 = size(cc_space_v_oovo,3)
    n4 = size(cc_space_v_oovo,4)

    !$OMP PARALLEL DO PRIVATE(i1,i2,i3,i4) COLLAPSE(2)
    do i4 = 1, n4
      do i3 = 1, n3
        do i2 = 1, n2
          do i1 = 1, n1
            cc_space_v_oovo(i1,i2,i3,i4) = cc_space_v_vooo(i3,i2,i1,i4)
          enddo
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

  else
    call gen_v_space(cc_nOa,cc_nOa,cc_nVa,cc_nOa, cc_list_occ,cc_list_occ,cc_list_vir,cc_list_occ, cc_space_v_oovo)
  endif

END_PROVIDER

! ooov

BEGIN_PROVIDER [double precision, cc_space_v_ooov, (cc_nOa, cc_nOa, cc_nOa, cc_nVa)]

  implicit none

  if (do_ao_cholesky) then

    integer :: i1, i2, i3, i4
    integer :: n1, n2, n3, n4
    
    n1 = size(cc_space_v_oovo,1)
    n2 = size(cc_space_v_oovo,2)
    n3 = size(cc_space_v_oovo,3)
    n4 = size(cc_space_v_oovo,4)

    !$OMP PARALLEL DO PRIVATE(i1,i2,i3,i4) COLLAPSE(2)
    do i4 = 1, n4
      do i3 = 1, n3
        do i2 = 1, n2
          do i1 = 1, n1
            cc_space_v_ooov(i1,i2,i3,i4) = cc_space_v_ovoo(i1,i4,i3,i2)
          enddo
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

  else
    call gen_v_space(cc_nOa,cc_nOa,cc_nOa,cc_nVa, cc_list_occ,cc_list_occ,cc_list_occ,cc_list_vir, cc_space_v_ooov)
  endif

END_PROVIDER

! vvoo

BEGIN_PROVIDER [double precision, cc_space_v_vvoo, (cc_nVa, cc_nVa, cc_nOa, cc_nOa)]

  implicit none

  if (do_ao_cholesky) then

    integer :: i1, i2, i3, i4
    integer :: n1, n2, n3, n4
    
    n1 = size(cc_space_v_vvoo,1)
    n2 = size(cc_space_v_vvoo,2)
    n3 = size(cc_space_v_vvoo,3)
    n4 = size(cc_space_v_vvoo,4)

    double precision, allocatable :: buffer(:,:,:,:)
    allocate(buffer(n1,n3,n2,n4))

    call dgemm('T','N', n1*n3, n2*n4, cholesky_mo_num, 1.d0, &
         cc_space_v_vo_chol, cholesky_mo_num, &
         cc_space_v_vo_chol, cholesky_mo_num, 0.d0, buffer, n1*n3)

    !$OMP PARALLEL DO PRIVATE(i1,i2,i3,i4) COLLAPSE(2)
    do i4 = 1, n4
      do i3 = 1, n3
        do i2 = 1, n2
          do i1 = 1, n1
            cc_space_v_vvoo(i1,i2,i3,i4) = buffer(i1,i3,i2,i4)
          enddo
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    deallocate(buffer)

  else
    call gen_v_space(cc_nVa,cc_nVa,cc_nOa,cc_nOa, cc_list_vir,cc_list_vir,cc_list_occ,cc_list_occ, cc_space_v_vvoo)
  endif

END_PROVIDER

! vovo

BEGIN_PROVIDER [double precision, cc_space_v_vovo, (cc_nVa, cc_nOa, cc_nVa, cc_nOa)]

  implicit none

  if (do_ao_cholesky) then

    integer :: i1, i2, i3, i4
    integer :: n1, n2, n3, n4
    
    n1 = size(cc_space_v_vovo,1)
    n2 = size(cc_space_v_vovo,2)
    n3 = size(cc_space_v_vovo,3)
    n4 = size(cc_space_v_vovo,4)

    double precision, allocatable :: buffer(:,:,:,:)
    allocate(buffer(n1,n3,n2,n4))

    call dgemm('T','N', n1*n3, n2*n4, cholesky_mo_num, 1.d0, &
         cc_space_v_vv_chol, cholesky_mo_num, &
         cc_space_v_oo_chol, cholesky_mo_num, 0.d0, buffer, n1*n3)

    !$OMP PARALLEL DO PRIVATE(i1,i2,i3,i4) COLLAPSE(2)
    do i4 = 1, n4
      do i3 = 1, n3
        do i2 = 1, n2
          do i1 = 1, n1
            cc_space_v_vovo(i1,i2,i3,i4) = buffer(i1,i3,i2,i4)
          enddo
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    deallocate(buffer)

  else
    call gen_v_space(cc_nVa,cc_nOa,cc_nVa,cc_nOa, cc_list_vir,cc_list_occ,cc_list_vir,cc_list_occ, cc_space_v_vovo)
  endif

END_PROVIDER

! voov

BEGIN_PROVIDER [double precision, cc_space_v_voov, (cc_nVa, cc_nOa, cc_nOa, cc_nVa)]

  implicit none

  if (do_ao_cholesky) then

    integer :: i1, i2, i3, i4
    integer :: n1, n2, n3, n4
    
    n1 = size(cc_space_v_voov,1)
    n2 = size(cc_space_v_voov,2)
    n3 = size(cc_space_v_voov,3)
    n4 = size(cc_space_v_voov,4)

    !$OMP PARALLEL DO PRIVATE(i1,i2,i3,i4) COLLAPSE(2)
    do i4 = 1, n4
      do i3 = 1, n3
        do i2 = 1, n2
          do i1 = 1, n1
            cc_space_v_voov(i1,i2,i3,i4) = cc_space_v_vvoo(i1,i4,i3,i2)
          enddo
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

  else
    call gen_v_space(cc_nVa,cc_nOa,cc_nOa,cc_nVa, cc_list_vir,cc_list_occ,cc_list_occ,cc_list_vir, cc_space_v_voov)
  endif

END_PROVIDER

! ovvo

BEGIN_PROVIDER [double precision, cc_space_v_ovvo, (cc_nOa, cc_nVa, cc_nVa, cc_nOa)]

  implicit none

  if (do_ao_cholesky) then

    integer :: i1, i2, i3, i4
    integer :: n1, n2, n3, n4
    
    n1 = size(cc_space_v_ovvo,1)
    n2 = size(cc_space_v_ovvo,2)
    n3 = size(cc_space_v_ovvo,3)
    n4 = size(cc_space_v_ovvo,4)

    !$OMP PARALLEL DO PRIVATE(i1,i2,i3,i4) COLLAPSE(2)
    do i4 = 1, n4
      do i3 = 1, n3
        do i2 = 1, n2
          do i1 = 1, n1
            cc_space_v_ovvo(i1,i2,i3,i4) = cc_space_v_vvoo(i3,i2,i1,i4)
          enddo
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

  else
    call gen_v_space(cc_nOa,cc_nVa,cc_nVa,cc_nOa, cc_list_occ,cc_list_vir,cc_list_vir,cc_list_occ, cc_space_v_ovvo)
  endif

END_PROVIDER

! ovov

BEGIN_PROVIDER [double precision, cc_space_v_ovov, (cc_nOa, cc_nVa, cc_nOa, cc_nVa)]

  implicit none

  if (do_ao_cholesky) then

    integer :: i1, i2, i3, i4
    integer :: n1, n2, n3, n4
    
    n1 = size(cc_space_v_ovov,1)
    n2 = size(cc_space_v_ovov,2)
    n3 = size(cc_space_v_ovov,3)
    n4 = size(cc_space_v_ovov,4)

    !$OMP PARALLEL DO PRIVATE(i1,i2,i3,i4) COLLAPSE(2)
    do i4 = 1, n4
      do i3 = 1, n3
        do i2 = 1, n2
          do i1 = 1, n1
            cc_space_v_ovov(i1,i2,i3,i4) = cc_space_v_vovo(i2,i1,i4,i3)
          enddo
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

  else
    call gen_v_space(cc_nOa,cc_nVa,cc_nOa,cc_nVa, cc_list_occ,cc_list_vir,cc_list_occ,cc_list_vir, cc_space_v_ovov)
  endif

END_PROVIDER

! oovv

BEGIN_PROVIDER [double precision, cc_space_v_oovv, (cc_nOa, cc_nOa, cc_nVa, cc_nVa)]

  implicit none

  if (do_ao_cholesky) then

    integer :: i1, i2, i3, i4
    integer :: n1, n2, n3, n4
    
    n1 = size(cc_space_v_oovv,1)
    n2 = size(cc_space_v_oovv,2)
    n3 = size(cc_space_v_oovv,3)
    n4 = size(cc_space_v_oovv,4)

    !$OMP PARALLEL DO PRIVATE(i1,i2,i3,i4) COLLAPSE(2)
    do i4 = 1, n4
      do i3 = 1, n3
        do i2 = 1, n2
          do i1 = 1, n1
            cc_space_v_oovv(i1,i2,i3,i4) = cc_space_v_vvoo(i3,i4,i1,i2)
          enddo
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

  else
    call gen_v_space(cc_nOa,cc_nOa,cc_nVa,cc_nVa, cc_list_occ,cc_list_occ,cc_list_vir,cc_list_vir, cc_space_v_oovv)
  endif

END_PROVIDER

! vvvo

BEGIN_PROVIDER [double precision, cc_space_v_vvvo, (cc_nVa, cc_nVa, cc_nVa, cc_nOa)]

  implicit none

  call gen_v_space(cc_nVa,cc_nVa,cc_nVa,cc_nOa, cc_list_vir,cc_list_vir,cc_list_vir,cc_list_occ, cc_space_v_vvvo)

END_PROVIDER

! vvov

BEGIN_PROVIDER [double precision, cc_space_v_vvov, (cc_nVa, cc_nVa, cc_nOa, cc_nVa)]

  implicit none

  call gen_v_space(cc_nVa,cc_nVa,cc_nOa,cc_nVa, cc_list_vir,cc_list_vir,cc_list_occ,cc_list_vir, cc_space_v_vvov)

END_PROVIDER

! vovv

BEGIN_PROVIDER [double precision, cc_space_v_vovv, (cc_nVa, cc_nOa, cc_nVa, cc_nVa)]

  implicit none

  call gen_v_space(cc_nVa,cc_nOa,cc_nVa,cc_nVa, cc_list_vir,cc_list_occ,cc_list_vir,cc_list_vir, cc_space_v_vovv)

END_PROVIDER

! ovvv

BEGIN_PROVIDER [double precision, cc_space_v_ovvv, (cc_nOa, cc_nVa, cc_nVa, cc_nVa)]

  implicit none

  call gen_v_space(cc_nOa,cc_nVa,cc_nVa,cc_nVa, cc_list_occ,cc_list_vir,cc_list_vir,cc_list_vir, cc_space_v_ovvv)

END_PROVIDER

! vvvv

BEGIN_PROVIDER [double precision, cc_space_v_vvvv, (cc_nVa, cc_nVa, cc_nVa, cc_nVa)]

  implicit none

  call gen_v_space(cc_nVa,cc_nVa,cc_nVa,cc_nVa, cc_list_vir,cc_list_vir,cc_list_vir,cc_list_vir, cc_space_v_vvvv)

END_PROVIDER

BEGIN_PROVIDER [double precision, cc_space_v_vv_chol, (cholesky_mo_num, cc_nVa, cc_nVa)]

  implicit none

  call gen_v_space_chol(cc_nVa, cc_nVa, cc_list_vir, cc_list_vir, cc_space_v_vv_chol, cholesky_mo_num)

END_PROVIDER

BEGIN_PROVIDER [double precision, cc_space_v_vo_chol, (cholesky_mo_num, cc_nVa, cc_nOa)]

  implicit none

  call gen_v_space_chol(cc_nVa, cc_nOa, cc_list_vir, cc_list_occ, cc_space_v_vo_chol, cholesky_mo_num)

END_PROVIDER

BEGIN_PROVIDER [double precision, cc_space_v_ov_chol, (cholesky_mo_num, cc_nOa, cc_nVa)]

  implicit none

  call gen_v_space_chol(cc_nOa, cc_nVa, cc_list_occ, cc_list_vir, cc_space_v_ov_chol, cholesky_mo_num)

END_PROVIDER

BEGIN_PROVIDER [double precision, cc_space_v_oo_chol, (cholesky_mo_num, cc_nOa, cc_nOa)]

  implicit none

  call gen_v_space_chol(cc_nOa, cc_nOa, cc_list_occ, cc_list_occ, cc_space_v_oo_chol, cholesky_mo_num)

END_PROVIDER

! ppqq

BEGIN_PROVIDER [double precision, cc_space_v_ppqq, (cc_n_mo, cc_n_mo)]

  implicit none

  BEGIN_DOC
  ! <pp|qq> integrals for general MOs (excepted core and deleted ones)
  END_DOC

  integer                       :: p,q
  double precision, allocatable :: tmp_v(:,:,:,:)

  allocate(tmp_v(cc_n_mo,cc_n_mo,cc_n_mo,cc_n_mo))

  call gen_v_space(cc_n_mo,cc_n_mo,cc_n_mo,cc_n_mo, cc_list_gen,cc_list_gen,cc_list_gen,cc_list_gen, tmp_v)

  do q = 1, cc_n_mo
    do p = 1, cc_n_mo
      cc_space_v_ppqq(p,q) = tmp_v(p,p,q,q)
    enddo
  enddo

  deallocate(tmp_v)

END_PROVIDER

! aaii

BEGIN_PROVIDER [double precision, cc_space_v_aaii, (cc_nVa,cc_nOa)]

  implicit none

  BEGIN_DOC
  ! <aa|ii> integrals
  ! a: virtual MO
  ! i: occupied MO
  END_DOC

  integer :: a,i

  do i = 1, cc_nOa
    do a = 1, cc_nVa
      cc_space_v_aaii(a,i) = cc_space_v_vvoo(a,a,i,i)
    enddo
  enddo

  FREE cc_space_v_vvoo

END_PROVIDER

! iiaa

BEGIN_PROVIDER [double precision, cc_space_v_iiaa, (cc_nOa,cc_nVa)]

  implicit none

  BEGIN_DOC
  ! <ii|aa> integrals
  ! a: virtual MO
  ! i: occupied MO
  END_DOC

  integer :: a,i

  do a = 1, cc_nVa
    do i = 1, cc_nOa
      cc_space_v_iiaa(i,a) = cc_space_v_oovv(i,i,a,a)
    enddo
  enddo

  FREE cc_space_v_oovv

END_PROVIDER

! iijj

BEGIN_PROVIDER [double precision, cc_space_v_iijj, (cc_nOa,cc_nOa)]

  implicit none

  BEGIN_DOC
  ! <ii|jj> integrals
  ! i,j: occupied MO
  END_DOC

  integer :: i,j

  do j = 1, cc_nOa
    do i = 1, cc_nOa
      cc_space_v_iijj(i,j) = cc_space_v_oooo(i,i,j,j)
    enddo
  enddo

  FREE cc_space_v_oooo

END_PROVIDER

! aabb

BEGIN_PROVIDER [double precision, cc_space_v_aabb, (cc_nVa,cc_nVa)]

  implicit none

  BEGIN_DOC
  ! <aa|bb> integrals
  ! a,b: virtual MO
  END_DOC

  integer :: a,b

  do b = 1, cc_nVa
    do a = 1, cc_nVa
      cc_space_v_aabb(a,b) = cc_space_v_vvvv(a,a,b,b)
    enddo
  enddo

  FREE cc_space_v_vvvv

END_PROVIDER

! iaia

BEGIN_PROVIDER [double precision, cc_space_v_iaia, (cc_nOa,cc_nVa)]

  implicit none

  BEGIN_DOC
  ! <ia|ia> integrals
  ! a: virtual MO
  ! i: occupied MO
  END_DOC

  integer :: a,i

  do a = 1, cc_nVa
    do i = 1, cc_nOa
      cc_space_v_iaia(i,a) = cc_space_v_ovov(i,a,i,a)
    enddo
  enddo

  FREE cc_space_v_ovov

END_PROVIDER

! iaai

BEGIN_PROVIDER [double precision, cc_space_v_iaai, (cc_nOa,cc_nVa)]

  implicit none

  BEGIN_DOC
  ! <ia|ai> integrals
  ! a: virtual MO
  ! i: inactive MO
  END_DOC

  integer :: a,i

  do a = 1, cc_nVa
    do i = 1, cc_nOa
      cc_space_v_iaai(i,a) = cc_space_v_ovvo(i,a,a,i)
    enddo
  enddo

  FREE cc_space_v_ovvo

END_PROVIDER

! aiia

BEGIN_PROVIDER [double precision, cc_space_v_aiia, (cc_nVa,cc_nOa)]

  implicit none

  BEGIN_DOC
  ! <ai|ia> integrals
  ! a: virtual MO
  ! i: inactive MO
  END_DOC

  integer :: a,i

  do i = 1, cc_nOa
    do a = 1, cc_nVa
      cc_space_v_aiia(a,i) = cc_space_v_voov(a,i,i,a)
    enddo
  enddo

  FREE cc_space_v_voov

END_PROVIDER

! oovv

BEGIN_PROVIDER [double precision, cc_space_w_oovv, (cc_nOa, cc_nOa, cc_nVa, cc_nVa)]

  implicit none

  double precision, allocatable :: tmp_v(:,:,:,:)
  integer :: i,j,a,b

  allocate(tmp_v(cc_nOa,cc_nOa,cc_nVa,cc_nVa))

  call gen_v_space(cc_nOa,cc_nOa,cc_nVa,cc_nVa, cc_list_occ,cc_list_occ,cc_list_vir,cc_list_vir, tmp_v)

  !$OMP PARALLEL &
  !$OMP SHARED(cc_nVa,cc_nOa,tmp_v,cc_space_w_oovv) &
  !$OMP PRIVATE(i,j,a,b)&
  !$OMP DEFAULT(NONE)
  !$OMP DO
  do b = 1, cc_nVa
    do a = 1, cc_nVa
      do j = 1, cc_nOa
        do i = 1, cc_nOa
          cc_space_w_oovv(i,j,a,b) = 2d0 * tmp_v(i,j,a,b) - tmp_v(j,i,a,b)
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  deallocate(tmp_v)

END_PROVIDER

! vvoo

BEGIN_PROVIDER [double precision, cc_space_w_vvoo, (cc_nVa, cc_nVa, cc_nOa, cc_nOa)]

  implicit none

  double precision, allocatable :: tmp_v(:,:,:,:)
  integer :: i,j,a,b

  allocate(tmp_v(cc_nVa,cc_nVa,cc_nOa,cc_nOa))

  call gen_v_space(cc_nVa,cc_nVa,cc_nOa,cc_nOa, cc_list_vir,cc_list_vir,cc_list_occ,cc_list_occ, tmp_v)

  !$OMP PARALLEL &
  !$OMP SHARED(cc_nVa,cc_nOa,tmp_v,cc_space_w_vvoo) &
  !$OMP PRIVATE(i,j,a,b)&
  !$OMP DEFAULT(NONE)
  !$OMP DO
  do j = 1, cc_nOa
    do i = 1, cc_nOa
      do b = 1, cc_nVa
        do a = 1, cc_nVa
          cc_space_w_vvoo(a,b,i,j) = 2d0 * tmp_v(a,b,i,j) - tmp_v(b,a,i,j)
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  deallocate(tmp_v)

END_PROVIDER

! F_oo

BEGIN_PROVIDER [double precision, cc_space_f_oo, (cc_nOa, cc_nOa)]

  implicit none

  call gen_f_space(psi_det(1,1,cc_ref), cc_nOa,cc_nOa, cc_list_occ,cc_list_occ, cc_space_f_oo)

END_PROVIDER

! F_ov

BEGIN_PROVIDER [double precision, cc_space_f_ov, (cc_nOa, cc_nVa)]

  implicit none

  call gen_f_space(psi_det(1,1,cc_ref), cc_nOa,cc_nVa, cc_list_occ,cc_list_vir, cc_space_f_ov)

END_PROVIDER

! F_vo

BEGIN_PROVIDER [double precision, cc_space_f_vo, (cc_nVa, cc_nOa)]

  implicit none

  call gen_f_space(psi_det(1,1,cc_ref), cc_nVa,cc_nOa, cc_list_vir,cc_list_occ, cc_space_f_vo)

END_PROVIDER

! F_vv

BEGIN_PROVIDER [double precision, cc_space_f_vv, (cc_nVa, cc_nVa)]

  implicit none

  call gen_f_space(psi_det(1,1,cc_ref), cc_nVa,cc_nVa, cc_list_vir,cc_list_vir, cc_space_f_vv)

END_PROVIDER

! F_o

BEGIN_PROVIDER [double precision, cc_space_f_o, (cc_nOa)]

  implicit none

  integer :: i

  do i = 1, cc_nOa
    cc_space_f_o(i) = cc_space_f_oo(i,i)
  enddo

END_PROVIDER

! F_v

BEGIN_PROVIDER [double precision, cc_space_f_v, (cc_nVa)]

  implicit none

  integer :: i

  do i = 1, cc_nVa
    cc_space_f_v(i) = cc_space_f_vv(i,i)
  enddo

END_PROVIDER

! Shift

subroutine shift_idx_spin(s,n_S,shift)

  implicit none

  BEGIN_DOC
  ! Shift for the partitionning alpha/beta of the spin orbitals
  ! n_S(1): number of spin alpha in the correspondong list
  ! n_S(2): number of spin beta in the correspondong list
  END_DOC

  integer, intent(in)  :: s, n_S(2)
  integer, intent(out) :: shift

  if (s == 1) then
    shift = 0
  else
    shift = n_S(1)
  endif

end

! F

subroutine gen_f_spin(det, n1,n2, n1_S,n2_S, list1,list2, dim1,dim2, f)

  implicit none

  BEGIN_DOC
  ! Compute the Fock matrix corresponding to two lists of spin orbitals.
  ! Ex: occ/occ, occ/vir,...
  END_DOC

  integer(bit_kind), intent(in) :: det(N_int,2)
  integer, intent(in)           :: n1,n2, n1_S(2), n2_S(2)
  integer, intent(in)           :: list1(n1,2), list2(n2,2)
  integer, intent(in)           :: dim1, dim2

  double precision, intent(out) :: f(dim1, dim2)

  double precision, allocatable :: tmp_F(:,:)
  integer                       :: i,j, idx_i,idx_j,i_shift,j_shift
  integer                       :: tmp_i,tmp_j
  integer                       :: si,sj,s
  PROVIDE big_array_exchange_integrals big_array_coulomb_integrals

  allocate(tmp_F(mo_num,mo_num))

  do sj = 1, 2
    call shift_idx_spin(sj,n2_S,j_shift)
    do si = 1, 2
      call shift_idx_spin(si,n1_S,i_shift)
      s = si + sj

      if (s == 2 .or. s == 4) then
        call get_fock_matrix_spin(det,sj,tmp_F)
      else
        do j = 1, mo_num
          do i = 1, mo_num
            tmp_F(i,j) = 0d0
          enddo
        enddo
      endif

      do tmp_j = 1, n2_S(sj)
        j = list2(tmp_j,sj)
        idx_j = tmp_j + j_shift
        do tmp_i = 1, n1_S(si)
          i = list1(tmp_i,si)
          idx_i = tmp_i + i_shift
          f(idx_i,idx_j) = tmp_F(i,j)
        enddo
      enddo

    enddo
  enddo

  deallocate(tmp_F)

end

! Get F

subroutine get_fock_matrix_spin(det,s,f)

  implicit none

  BEGIN_DOC
  ! Fock matrix alpha or beta of an arbitrary det
  END_DOC

  integer(bit_kind), intent(in) :: det(N_int,2)
  integer, intent(in)           :: s

  double precision, intent(out) :: f(mo_num,mo_num)

  integer                       :: p,q,i,s1,s2
  integer(bit_kind)             :: res(N_int,2)
  logical                       :: ok
  double precision              :: mo_two_e_integral

  if (s == 1) then
    s1 = 1
    s2 = 2
  else
    s1 = 2
    s2 = 1
  endif

  PROVIDE big_array_coulomb_integrals big_array_exchange_integrals

  !$OMP PARALLEL &
  !$OMP SHARED(f,mo_num,s1,s2,N_int,det,mo_one_e_integrals,big_array_coulomb_integrals,big_array_exchange_integrals) &
  !$OMP PRIVATE(p,q,ok,i,res)&
  !$OMP DEFAULT(NONE)
  !$OMP DO collapse(1)
  do q = 1, mo_num
    do p = 1, mo_num
      f(p,q) = mo_one_e_integrals(p,q)
      do i = 1, mo_num
        call apply_hole(det, s1, i, res, ok, N_int)
        if (ok) then
!          f(p,q) = f(p,q) + mo_two_e_integral(p,i,q,i) - mo_two_e_integral(p,i,i,q)
          f(p,q) = f(p,q) + big_array_coulomb_integrals(i,p,q) - big_array_exchange_integrals(i,p,q)
        endif
      enddo
      do i = 1, mo_num
        call apply_hole(det, s2, i, res, ok, N_int)
        if (ok) then
          f(p,q) = f(p,q) + big_array_coulomb_integrals(i,p,q)
        endif
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

end

! V

subroutine gen_v_spin(n1,n2,n3,n4, n1_S,n2_S,n3_S,n4_S, list1,list2,list3,list4, dim1,dim2,dim3,dim4, v)

  implicit none

   BEGIN_DOC
  ! Compute the bi electronic integrals corresponding to four lists of spin orbitals.
  ! Ex: occ/occ/occ/occ, occ/vir/occ/vir, ...
  END_DOC

  integer, intent(in)           :: n1,n2,n3,n4,n1_S(2),n2_S(2),n3_S(2),n4_S(2)
  integer, intent(in)           :: list1(n1,2), list2(n2,2), list3(n3,2), list4(n4,2)
  integer, intent(in)           :: dim1, dim2, dim3, dim4
  double precision, intent(out) :: v(dim1,dim2,dim3,dim4)

  double precision              :: mo_two_e_integral
  integer                       :: i,j,k,l,idx_i,idx_j,idx_k,idx_l
  integer                       :: i_shift,j_shift,k_shift,l_shift
  integer                       :: tmp_i,tmp_j,tmp_k,tmp_l
  integer                       :: si,sj,sk,sl,s

  PROVIDE cc_space_v

  !$OMP PARALLEL &
  !$OMP SHARED(cc_space_v,n1_S,n2_S,n3_S,n4_S,list1,list2,list3,list4,v) &
  !$OMP PRIVATE(s,si,sj,sk,sl,i_shift,j_shift,k_shift,l_shift, &
  !$OMP i,j,k,l,idx_i,idx_j,idx_k,idx_l,&
  !$OMP tmp_i,tmp_j,tmp_k,tmp_l)&
  !$OMP DEFAULT(NONE)

  do sl = 1, 2
    call shift_idx_spin(sl,n4_S,l_shift)
    do sk = 1, 2
      call shift_idx_spin(sk,n3_S,k_shift)
      do sj = 1, 2
        call shift_idx_spin(sj,n2_S,j_shift)
        do si = 1, 2
          call shift_idx_spin(si,n1_S,i_shift)

          s = si+sj+sk+sl
          ! <aa||aa> or <bb||bb>
          if (s == 4 .or. s == 8) then
            !$OMP DO collapse(3)
            do tmp_l = 1, n4_S(sl)
              do tmp_k = 1, n3_S(sk)
                do tmp_j = 1, n2_S(sj)
                  do tmp_i = 1, n1_S(si)
                    l = list4(tmp_l,sl)
                    idx_l = tmp_l + l_shift
                    k = list3(tmp_k,sk)
                    idx_k = tmp_k + k_shift
                    j = list2(tmp_j,sj)
                    idx_j = tmp_j + j_shift
                    i = list1(tmp_i,si)
                    idx_i = tmp_i + i_shift
                       !v(idx_i,idx_j,idx_k,idx_l) = mo_two_e_integral(i,j,k,l) - mo_two_e_integral(j,i,k,l)
                       v(idx_i,idx_j,idx_k,idx_l) = cc_space_v(i,j,k,l) - cc_space_v(j,i,k,l)
                  enddo
                enddo
              enddo
            enddo
            !$OMP END DO

          ! <ab||ab> or <ba||ba>
          elseif (si == sk .and. sj == sl) then
            !$OMP DO collapse(3)
            do tmp_l = 1, n4_S(sl)
              do tmp_k = 1, n3_S(sk)
                do tmp_j = 1, n2_S(sj)
                  do tmp_i = 1, n1_S(si)
                    l = list4(tmp_l,sl)
                    idx_l = tmp_l + l_shift
                    k = list3(tmp_k,sk)
                    idx_k = tmp_k + k_shift
                    j = list2(tmp_j,sj)
                    idx_j = tmp_j + j_shift
                    i = list1(tmp_i,si)
                    idx_i = tmp_i + i_shift
                       !v(idx_i,idx_j,idx_k,idx_l) = mo_two_e_integral(i,j,k,l)
                       v(idx_i,idx_j,idx_k,idx_l) = cc_space_v(i,j,k,l)
                  enddo
                enddo
              enddo
            enddo
            !$OMP END DO

          ! <ab||ba> or <ba||ab>
          elseif (si == sl .and. sj == sk) then
            !$OMP DO collapse(3)
            do tmp_l = 1, n4_S(sl)
              do tmp_k = 1, n3_S(sk)
                do tmp_j = 1, n2_S(sj)
                  do tmp_i = 1, n1_S(si)
                    l = list4(tmp_l,sl)
                    idx_l = tmp_l + l_shift
                    k = list3(tmp_k,sk)
                    idx_k = tmp_k + k_shift
                    j = list2(tmp_j,sj)
                    idx_j = tmp_j + j_shift
                    i = list1(tmp_i,si)
                    idx_i = tmp_i + i_shift
                       !v(idx_i,idx_j,idx_k,idx_l) = - mo_two_e_integral(j,i,k,l)
                       v(idx_i,idx_j,idx_k,idx_l) = - cc_space_v(j,i,k,l)
                  enddo
                enddo
              enddo
            enddo
            !$OMP END DO
          else
             !$OMP DO collapse(3)
            do tmp_l = 1, n4_S(sl)
              do tmp_k = 1, n3_S(sk)
                do tmp_j = 1, n2_S(sj)
                  do tmp_i = 1, n1_S(si)
                    l = list4(tmp_l,sl)
                    idx_l = tmp_l + l_shift
                    k = list3(tmp_k,sk)
                    idx_k = tmp_k + k_shift
                    j = list2(tmp_j,sj)
                    idx_j = tmp_j + j_shift
                    i = list1(tmp_i,si)
                    idx_i = tmp_i + i_shift
                       v(idx_i,idx_j,idx_k,idx_l) = 0d0
                  enddo
                enddo
              enddo
            enddo
            !$OMP END DO
          endif

        enddo
      enddo
    enddo
  enddo
  !$OMP END PARALLEL

end

! V_3idx

subroutine gen_v_spin_3idx(n1,n2,n3,n4, idx_l, n1_S,n2_S,n3_S,n4_S, list1,list2,list3,list4, dim1,dim2,dim3, v_l)

  implicit none

   BEGIN_DOC
  ! Compute the bi electronic integrals corresponding to four lists of spin orbitals.
  ! Ex: occ/occ/occ/occ, occ/vir/occ/vir, ...
  END_DOC

  integer, intent(in)           :: n1,n2,n3,n4,idx_l,n1_S(2),n2_S(2),n3_S(2),n4_S(2)
  integer, intent(in)           :: list1(n1,2), list2(n2,2), list3(n3,2), list4(n4,2)
  integer, intent(in)           :: dim1, dim2, dim3
  double precision, intent(out) :: v_l(dim1,dim2,dim3)

  double precision              :: mo_two_e_integral
  integer                       :: i,j,k,l,idx_i,idx_j,idx_k
  integer                       :: i_shift,j_shift,k_shift,l_shift
  integer                       :: tmp_i,tmp_j,tmp_k,tmp_l
  integer                       :: si,sj,sk,sl,s

  PROVIDE cc_space_v

  if (idx_l <= n4_S(1)) then
    sl = 1
  else
    sl = 2
  endif
  call shift_idx_spin(sl,n4_S,l_shift)
  tmp_l = idx_l - l_shift
  l = list4(tmp_l,sl)

  !$OMP PARALLEL &
  !$OMP SHARED(l,sl,idx_l,cc_space_v,n1_S,n2_S,n3_S,n4_S,list1,list2,list3,list4,v_l) &
  !$OMP PRIVATE(s,si,sj,sk,i_shift,j_shift,k_shift, &
  !$OMP i,j,k,idx_i,idx_j,idx_k,&
  !$OMP tmp_i,tmp_j,tmp_k)&
  !$OMP DEFAULT(NONE)

  do sk = 1, 2
    call shift_idx_spin(sk,n3_S,k_shift)
    do sj = 1, 2
      call shift_idx_spin(sj,n2_S,j_shift)
      do si = 1, 2
        call shift_idx_spin(si,n1_S,i_shift)

        s = si+sj+sk+sl
        ! <aa||aa> or <bb||bb>
        if (s == 4 .or. s == 8) then
          !$OMP DO collapse(2)
          do tmp_k = 1, n3_S(sk)
            do tmp_j = 1, n2_S(sj)
              do tmp_i = 1, n1_S(si)
                k = list3(tmp_k,sk)
                idx_k = tmp_k + k_shift
                j = list2(tmp_j,sj)
                idx_j = tmp_j + j_shift
                i = list1(tmp_i,si)
                idx_i = tmp_i + i_shift
                   !v(idx_i,idx_j,idx_k,idx_l) = mo_two_e_integral(i,j,k,l) - mo_two_e_integral(j,i,k,l)
                   v_l(idx_i,idx_j,idx_k) = cc_space_v(i,j,k,l) - cc_space_v(j,i,k,l)
              enddo
            enddo
          enddo
          !$OMP END DO

        ! <ab||ab> or <ba||ba>
        elseif (si == sk .and. sj == sl) then
          !$OMP DO collapse(2)
          do tmp_k = 1, n3_S(sk)
            do tmp_j = 1, n2_S(sj)
              do tmp_i = 1, n1_S(si)
                k = list3(tmp_k,sk)
                idx_k = tmp_k + k_shift
                j = list2(tmp_j,sj)
                idx_j = tmp_j + j_shift
                i = list1(tmp_i,si)
                idx_i = tmp_i + i_shift
                   !v(idx_i,idx_j,idx_k,idx_l) = mo_two_e_integral(i,j,k,l)
                   v_l(idx_i,idx_j,idx_k) = cc_space_v(i,j,k,l)
              enddo
            enddo
          enddo
          !$OMP END DO

        ! <ab||ba> or <ba||ab>
        elseif (si == sl .and. sj == sk) then
          !$OMP DO collapse(2)
          do tmp_k = 1, n3_S(sk)
            do tmp_j = 1, n2_S(sj)
              do tmp_i = 1, n1_S(si)
                k = list3(tmp_k,sk)
                idx_k = tmp_k + k_shift
                j = list2(tmp_j,sj)
                idx_j = tmp_j + j_shift
                i = list1(tmp_i,si)
                idx_i = tmp_i + i_shift
                   !v(idx_i,idx_j,idx_k,idx_l) = - mo_two_e_integral(j,i,k,l)
                   v_l(idx_i,idx_j,idx_k) = - cc_space_v(j,i,k,l)
              enddo
            enddo
          enddo
          !$OMP END DO
        else
          !$OMP DO collapse(2)
          do tmp_k = 1, n3_S(sk)
            do tmp_j = 1, n2_S(sj)
              do tmp_i = 1, n1_S(si)
                k = list3(tmp_k,sk)
                idx_k = tmp_k + k_shift
                j = list2(tmp_j,sj)
                idx_j = tmp_j + j_shift
                i = list1(tmp_i,si)
                idx_i = tmp_i + i_shift
                   v_l(idx_i,idx_j,idx_k) = 0d0
              enddo
            enddo
          enddo
          !$OMP END DO
        endif

      enddo
    enddo
  enddo
  !$OMP END PARALLEL

end

! V_3idx_ij_l

subroutine gen_v_spin_3idx_ij_l(n1,n2,n3,n4, idx_k, n1_S,n2_S,n3_S,n4_S, list1,list2,list3,list4, dim1,dim2,dim3, v_k)

  implicit none

   BEGIN_DOC
  ! Compute the bi electronic integrals corresponding to four lists of spin orbitals.
  ! Ex: occ/occ/occ/occ, occ/vir/occ/vir, ...
  END_DOC

  integer, intent(in)           :: n1,n2,n3,n4,idx_k,n1_S(2),n2_S(2),n3_S(2),n4_S(2)
  integer, intent(in)           :: list1(n1,2), list2(n2,2), list3(n3,2), list4(n4,2)
  integer, intent(in)           :: dim1, dim2, dim3
  double precision, intent(out) :: v_k(dim1,dim2,dim3)

  double precision              :: mo_two_e_integral
  integer                       :: i,j,k,l,idx_i,idx_j,idx_l
  integer                       :: i_shift,j_shift,k_shift,l_shift
  integer                       :: tmp_i,tmp_j,tmp_k,tmp_l
  integer                       :: si,sj,sk,sl,s

  PROVIDE cc_space_v

  if (idx_k <= n3_S(1)) then
    sk = 1
  else
    sk = 2
  endif
  call shift_idx_spin(sk,n3_S,k_shift)
  tmp_k = idx_k - k_shift
  k = list3(tmp_k,sk)

  !$OMP PARALLEL &
  !$OMP SHARED(k,sk,idx_k,cc_space_v,n1_S,n2_S,n3_S,n4_S,list1,list2,list3,list4,v_k) &
  !$OMP PRIVATE(s,si,sj,sl,i_shift,j_shift,l_shift, &
  !$OMP i,j,l,idx_i,idx_j,idx_l,&
  !$OMP tmp_i,tmp_j,tmp_l)&
  !$OMP DEFAULT(NONE)

  do sl = 1, 2
    call shift_idx_spin(sl,n4_S,l_shift)
    do sj = 1, 2
      call shift_idx_spin(sj,n2_S,j_shift)
      do si = 1, 2
        call shift_idx_spin(si,n1_S,i_shift)

        s = si+sj+sk+sl
        ! <aa||aa> or <bb||bb>
        if (s == 4 .or. s == 8) then
          !$OMP DO collapse(2)
          do tmp_l = 1, n4_S(sl)
            do tmp_j = 1, n2_S(sj)
              do tmp_i = 1, n1_S(si)
                l = list4(tmp_l,sl)
                idx_l = tmp_l + l_shift
                j = list2(tmp_j,sj)
                idx_j = tmp_j + j_shift
                i = list1(tmp_i,si)
                idx_i = tmp_i + i_shift
                   !v(idx_i,idx_j,idx_k,idx_l) = mo_two_e_integral(i,j,k,l) - mo_two_e_integral(j,i,k,l)
                   v_k(idx_i,idx_j,idx_l) = cc_space_v(i,j,k,l) - cc_space_v(j,i,k,l)
              enddo
            enddo
          enddo
          !$OMP END DO

        ! <ab||ab> or <ba||ba>
        elseif (si == sk .and. sj == sl) then
          !$OMP DO collapse(2)
          do tmp_l = 1, n4_S(sl)
            do tmp_j = 1, n2_S(sj)
              do tmp_i = 1, n1_S(si)
                l = list4(tmp_l,sl)
                idx_l = tmp_l + l_shift
                j = list2(tmp_j,sj)
                idx_j = tmp_j + j_shift
                i = list1(tmp_i,si)
                idx_i = tmp_i + i_shift
                   !v(idx_i,idx_j,idx_k,idx_l) = mo_two_e_integral(i,j,k,l)
                   v_k(idx_i,idx_j,idx_l) = cc_space_v(i,j,k,l)
              enddo
            enddo
          enddo
          !$OMP END DO

        ! <ab||ba> or <ba||ab>
        elseif (si == sl .and. sj == sk) then
          !$OMP DO collapse(2)
          do tmp_l = 1, n4_S(sl)
            do tmp_j = 1, n2_S(sj)
              do tmp_i = 1, n1_S(si)
                l = list4(tmp_l,sl)
                idx_l = tmp_l + l_shift
                j = list2(tmp_j,sj)
                idx_j = tmp_j + j_shift
                i = list1(tmp_i,si)
                idx_i = tmp_i + i_shift
                   !v(idx_i,idx_j,idx_k,idx_l) = - mo_two_e_integral(j,i,k,l)
                   v_k(idx_i,idx_j,idx_l) = - cc_space_v(j,i,k,l)
              enddo
            enddo
          enddo
          !$OMP END DO
        else
          !$OMP DO collapse(2)
          do tmp_l = 1, n4_S(sl)
            do tmp_j = 1, n2_S(sj)
              do tmp_i = 1, n1_S(si)
                l = list4(tmp_l,sl)
                idx_l = tmp_l + l_shift
                j = list2(tmp_j,sj)
                idx_j = tmp_j + j_shift
                i = list1(tmp_i,si)
                idx_i = tmp_i + i_shift
                   v_k(idx_i,idx_j,idx_l) = 0d0
              enddo
            enddo
          enddo
          !$OMP END DO
        endif

      enddo
    enddo
  enddo
  !$OMP END PARALLEL

end

! V_3idx_i_kl

subroutine gen_v_spin_3idx_i_kl(n1,n2,n3,n4, idx_j, n1_S,n2_S,n3_S,n4_S, list1,list2,list3,list4, dim1,dim2,dim3, v_j)

  implicit none

   BEGIN_DOC
  ! Compute the bi electronic integrals corresponding to four lists of spin orbitals.
  ! Ex: occ/occ/occ/occ, occ/vir/occ/vir, ...
  END_DOC

  integer, intent(in)           :: n1,n2,n3,n4,idx_j,n1_S(2),n2_S(2),n3_S(2),n4_S(2)
  integer, intent(in)           :: list1(n1,2), list2(n2,2), list3(n3,2), list4(n4,2)
  integer, intent(in)           :: dim1, dim2, dim3
  double precision, intent(out) :: v_j(dim1,dim2,dim3)

  double precision              :: mo_two_e_integral
  integer                       :: i,j,k,l,idx_i,idx_k,idx_l
  integer                       :: i_shift,j_shift,k_shift,l_shift
  integer                       :: tmp_i,tmp_j,tmp_k,tmp_l
  integer                       :: si,sj,sk,sl,s

  PROVIDE cc_space_v

  if (idx_j <= n2_S(1)) then
    sj = 1
  else
    sj = 2
  endif
  call shift_idx_spin(sj,n2_S,j_shift)
  tmp_j = idx_j - j_shift
  j = list2(tmp_j,sj)

  !$OMP PARALLEL &
  !$OMP SHARED(j,sj,idx_j,cc_space_v,n1_S,n2_S,n3_S,n4_S,list1,list2,list3,list4,v_j) &
  !$OMP PRIVATE(s,si,sk,sl,i_shift,l_shift,k_shift, &
  !$OMP i,k,l,idx_i,idx_k,idx_l,&
  !$OMP tmp_i,tmp_k,tmp_l)&
  !$OMP DEFAULT(NONE)

  do sl = 1, 2
    call shift_idx_spin(sl,n4_S,l_shift)
    do sk = 1, 2
      call shift_idx_spin(sk,n3_S,k_shift)
      do si = 1, 2
        call shift_idx_spin(si,n1_S,i_shift)

        s = si+sj+sk+sl
        ! <aa||aa> or <bb||bb>
        if (s == 4 .or. s == 8) then
          !$OMP DO collapse(2)
          do tmp_l = 1, n4_S(sl)
            do tmp_k = 1, n3_S(sk)
              do tmp_i = 1, n1_S(si)
                l = list4(tmp_l,sl)
                idx_l = tmp_l + l_shift
                k = list3(tmp_k,sk)
                idx_k = tmp_k + k_shift
                i = list1(tmp_i,si)
                idx_i = tmp_i + i_shift
                   !v(idx_i,idx_j,idx_k,idx_l) = mo_two_e_integral(i,j,k,l) - mo_two_e_integral(j,i,k,l)
                   v_j(idx_i,idx_k,idx_l) = cc_space_v(i,j,k,l) - cc_space_v(j,i,k,l)
              enddo
            enddo
          enddo
          !$OMP END DO

        ! <ab||ab> or <ba||ba>
        elseif (si == sk .and. sj == sl) then
          !$OMP DO collapse(2)
          do tmp_l = 1, n4_S(sl)
            do tmp_k = 1, n3_S(sk)
              do tmp_i = 1, n1_S(si)
                l = list4(tmp_l,sl)
                idx_l = tmp_l + l_shift
                k = list3(tmp_k,sk)
                idx_k = tmp_k + k_shift
                i = list1(tmp_i,si)
                idx_i = tmp_i + i_shift
                   !v(idx_i,idx_j,idx_k,idx_l) = mo_two_e_integral(i,j,k,l)
                   v_j(idx_i,idx_k,idx_l) = cc_space_v(i,j,k,l)
              enddo
            enddo
          enddo
          !$OMP END DO

        ! <ab||ba> or <ba||ab>
        elseif (si == sl .and. sj == sk) then
          !$OMP DO collapse(2)
          do tmp_l = 1, n4_S(sl)
            do tmp_k = 1, n3_S(sk)
              do tmp_i = 1, n1_S(si)
                l = list4(tmp_l,sl)
                idx_l = tmp_l + l_shift
                k = list3(tmp_k,sk)
                idx_k = tmp_k + k_shift
                i = list1(tmp_i,si)
                idx_i = tmp_i + i_shift
                   !v(idx_i,idx_j,idx_k,idx_l) = - mo_two_e_integral(j,i,k,l)
                   v_j(idx_i,idx_k,idx_l) = - cc_space_v(j,i,k,l)
              enddo
            enddo
          enddo
          !$OMP END DO
        else
          !$OMP DO collapse(2)
          do tmp_l = 1, n4_S(sl)
            do tmp_k = 1, n3_S(sk)
              do tmp_i = 1, n1_S(si)
                l = list4(tmp_l,sl)
                idx_l = tmp_l + l_shift
                k = list3(tmp_k,sk)
                idx_k = tmp_k + k_shift
                i = list1(tmp_i,si)
                idx_i = tmp_i + i_shift
                   v_j(idx_i,idx_k,idx_l) = 0d0
              enddo
            enddo
          enddo
          !$OMP END DO
        endif

      enddo
    enddo
  enddo
  !$OMP END PARALLEL

end
