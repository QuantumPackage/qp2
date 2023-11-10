program tc_bi_ortho

  BEGIN_DOC
  ! TODO
  END_DOC
  implicit none

  my_grid_becke = .True.
  PROVIDE tc_grid1_a tc_grid1_r
  my_n_pt_r_grid = tc_grid1_r
  my_n_pt_a_grid = tc_grid1_a
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  call ERI_dump()
  call KMat_tilde_dump()
  call LMat_tilde_dump()

end

! ---

subroutine KMat_tilde_dump()

  implicit none
  integer              :: i, j, k, l
  integer              :: isym, ms2, st, iii
  character(16)        :: corb
  double precision     :: t1, t2
  integer, allocatable :: orbsym(:)

  print *, ' generating FCIDUMP'
  call wall_time(t1)

  PROVIDE mo_bi_ortho_tc_two_e_chemist
  PROVIDE mo_bi_ortho_tc_one_e

  isym = 1
  ms2  = elec_alpha_num - elec_beta_num
  st   = 0
  iii  = 0

  allocate(orbsym(mo_num))
  orbsym(1:mo_num) = 1

  open(33, file='FCIDUMP', action='write')  

    write(33,'("&",a)') 'FCI'
    write(33,'(1x,a,"=",i0,",")') 'NORB', mo_num
    write(33,'(1x,a,"=",i0,",")') 'NELEC', elec_num
    write(33,'(1x,a,"=",i0,",")') 'MS2', ms2
    write(33,'(1x,a,"=",i0,",")') 'ISYM', isym
    write(corb,'(i0)') mo_num
    write(33,'(1x,a,"=",'//corb//'(i0,","))') 'ORBSYM', orbsym
    write(33,'(1x,a,"=",i0,",")') 'ST', st
    write(33,'(1x,a,"=",i0,",")') 'III', iii
    write(33,'(1x,a,"=",i0,",")') 'OCC', (elec_num-ms2)/2+ms2
    write(33,'(1x,a,"=",i0,",")') 'CLOSED', 2*elec_alpha_num
    write(33,'(1x,"/")')

    do l = 1, mo_num
      do k = 1, mo_num
        do j = 1, mo_num
          do i = 1, mo_num
            ! TCHint convention
            write(33, '(ES15.7, 4X, 4(I4, 2X))') mo_bi_ortho_tc_two_e_chemist(j,i,l,k), i, j, k, l
          enddo
        enddo
      enddo
    enddo

    do j = 1, mo_num
      do i = 1, mo_num
        ! TCHint convention
        write(33, '(ES15.7, 4X, 4(I4, 2X))') mo_bi_ortho_tc_one_e(i,j), i, j, 0, 0
      enddo
    enddo

  close(33)

  deallocate(orbsym)

  call wall_time(t2)
  print *, ' end after (min)', (t2-t1)/60.d0

  return
end subroutine KMat_tilde_dump

! ---

subroutine ERI_dump()

  implicit none
  integer                       :: i, j, k, l
  double precision, allocatable :: a1(:,:,:,:), a2(:,:,:,:)

  PROVIDE mo_r_coef mo_l_coef

  allocate(a2(ao_num,ao_num,ao_num,mo_num))

  call dgemm( 'T', 'N', ao_num*ao_num*ao_num, mo_num, ao_num, 1.d0   &
            , ao_two_e_coul(1,1,1,1), ao_num, mo_l_coef(1,1), ao_num &
            , 0.d0 , a2(1,1,1,1), ao_num*ao_num*ao_num)

  allocate(a1(ao_num,ao_num,mo_num,mo_num))

  call dgemm( 'T', 'N', ao_num*ao_num*mo_num, mo_num, ao_num, 1.d0 &
            , a2(1,1,1,1), ao_num, mo_r_coef(1,1), ao_num          &
            , 0.d0, a1(1,1,1,1), ao_num*ao_num*mo_num)

  deallocate(a2)
  allocate(a2(ao_num,mo_num,mo_num,mo_num))

  call dgemm( 'T', 'N', ao_num*mo_num*mo_num, mo_num, ao_num, 1.d0 &
            , a1(1,1,1,1), ao_num, mo_l_coef(1,1), ao_num          &
            , 0.d0, a2(1,1,1,1), ao_num*mo_num*mo_num)

  deallocate(a1)
  allocate(a1(mo_num,mo_num,mo_num,mo_num))

  call dgemm( 'T', 'N', mo_num*mo_num*mo_num, mo_num, ao_num, 1.d0 &
            , a2(1,1,1,1), ao_num, mo_r_coef(1,1), ao_num          &
            , 0.d0, a1(1,1,1,1), mo_num*mo_num*mo_num)

  deallocate(a2)

  open(33, file='ERI.dat', action='write')  
    do l = 1, mo_num
      do k = 1, mo_num
        do j = 1, mo_num
          do i = 1, mo_num
            write(33, '(4(I4, 2X), 4X, ES15.7)') i, j, k, l, a1(i,j,k,l)
          enddo
        enddo
      enddo
    enddo
  close(33)

  deallocate(a1)

  return
end subroutine ERI_dump

! ---

subroutine LMat_tilde_dump()

  implicit none
  integer          :: i, j, k, l, m, n
  double precision :: integral
  double precision :: t1, t2

  print *, ' generating TCDUMP'
  call wall_time(t1)

  PROVIDE mo_l_coef mo_r_coef

  open(33, file='TCDUMP', action='write')  
    write(33, '(4X, I4)') mo_num
    do n = 1, mo_num
      do m = 1, mo_num
        do l = 1, mo_num
          do k = 1, mo_num
            do j = 1, mo_num
              do i = 1, mo_num
                ! < i j k | -L | l m n > with a BI-ORTHONORMAL MOLECULAR ORBITALS
                call give_integrals_3_body_bi_ort(i, j, k, l, m, n, integral)
                !write(33, '(6(I4, 2X), 4X, E15.7)') i, j, k, l, m, n, integral
                ! TCHint convention
                if(dabs(integral).gt.1d-10) then
                  write(33, '(ES15.7, 4X, 6(I4, 2X))') -integral/3.d0, i, j, k, l, m, n
                  !write(33, '(ES15.7, 4X, 6(I4, 2X))') -integral/3.d0, l, m, n, i, j, k
                endif
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  close(33)

  call wall_time(t2)
  print *, ' end after (min)', (t2-t1)/60.d0

  return
end subroutine LMat_tilde_dump

! ---
