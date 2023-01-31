
! ---

BEGIN_PROVIDER [double precision, ao_two_e_tc_tot, (ao_num, ao_num, ao_num, ao_num) ]

  BEGIN_DOC
  !
  ! ao_two_e_tc_tot(k,i,l,j) = (ki|V^TC(r_12)|lj) = <lk| V^TC(r_12) |ji> where V^TC(r_12) is the total TC operator 
  !
  ! including both hermitian and non hermitian parts. THIS IS IN CHEMIST NOTATION. 
  !
  ! WARNING :: non hermitian ! acts on "the right functions" (i,j)
  !
  END_DOC

  integer                    :: i, j, k, l
  double precision           :: integral_sym, integral_nsym
  double precision, external :: get_ao_tc_sym_two_e_pot

  provide j1b_type

  if(j1b_type .eq. 3) then

    do j = 1, ao_num
      do l = 1, ao_num
        do i = 1, ao_num
          do k = 1, ao_num
            ao_two_e_tc_tot(k,i,l,j) = ao_tc_int_chemist(k,i,l,j)
            !write(222,*) ao_two_e_tc_tot(k,i,l,j) 
          enddo
        enddo
      enddo
    enddo

  else

    PROVIDE ao_tc_sym_two_e_pot_in_map

    do j = 1, ao_num
      do l = 1, ao_num
        do i = 1, ao_num
          do k = 1, ao_num

            integral_sym  = get_ao_tc_sym_two_e_pot(i, j, k, l, ao_tc_sym_two_e_pot_map)
            ! ao_non_hermit_term_chemist(k,i,l,j) = < k l | [erf( mu r12) - 1] d/d_r12 | i j > on the AO basis
            integral_nsym = ao_non_hermit_term_chemist(k,i,l,j)

            !print *, ' sym     integ = ', integral_sym
            !print *, ' non-sym integ = ', integral_nsym

            ao_two_e_tc_tot(k,i,l,j) = integral_sym + integral_nsym 
            !write(111,*) ao_two_e_tc_tot(k,i,l,j) 
          enddo
        enddo
      enddo
    enddo

  endif

END_PROVIDER 

! ---

double precision function bi_ortho_mo_ints(l, k, j, i)

  BEGIN_DOC
  !
  ! <mo^L_k mo^L_l | V^TC(r_12) | mo^R_i mo^R_j>
  !
  ! WARNING :: very naive, super slow, only used to DEBUG.
  !
  END_DOC

  implicit none
  integer, intent(in) :: i, j, k, l
  integer             :: m, n, p, q

  bi_ortho_mo_ints = 0.d0
  do m = 1, ao_num
    do p = 1, ao_num
      do n = 1, ao_num
        do q = 1, ao_num
          !                                   p1h1p2h2   l1                  l2              r1               r2
          bi_ortho_mo_ints += ao_two_e_tc_tot(n,q,m,p) * mo_l_coef(m,l) * mo_l_coef(n,k) * mo_r_coef(p,j) * mo_r_coef(q,i)
        enddo
      enddo
    enddo
  enddo

end function bi_ortho_mo_ints

! ---

! TODO :: transform into DEGEMM

BEGIN_PROVIDER [double precision, mo_bi_ortho_tc_two_e_chemist, (mo_num, mo_num, mo_num, mo_num)]

  BEGIN_DOC
  !
  ! mo_bi_ortho_tc_two_e_chemist(k,i,l,j) = <k l|V(r_12)|i j> where i,j are right MOs and k,l are left MOs
  !
  END_DOC

  implicit none
  integer                       :: i, j, k, l, m, n, p, q
  double precision, allocatable :: mo_tmp_1(:,:,:,:), mo_tmp_2(:,:,:,:)

  allocate(mo_tmp_1(mo_num,ao_num,ao_num,ao_num))
  mo_tmp_1 = 0.d0

  do m = 1, ao_num
    do p = 1, ao_num
      do n = 1, ao_num
        do q = 1, ao_num
          do k = 1, mo_num
            !       (k n|p m)    = sum_q c_qk * (q n|p m)
            mo_tmp_1(k,n,p,m) += mo_l_coef_transp(k,q) * ao_two_e_tc_tot(q,n,p,m)
          enddo
        enddo
      enddo
    enddo
  enddo

  allocate(mo_tmp_2(mo_num,mo_num,ao_num,ao_num))
  mo_tmp_2 = 0.d0

  do m = 1, ao_num
    do p = 1, ao_num
      do n = 1, ao_num
        do i = 1, mo_num
          do k = 1, mo_num
            !       (k i|p m) = sum_n c_ni * (k n|p m)
            mo_tmp_2(k,i,p,m) += mo_r_coef_transp(i,n) * mo_tmp_1(k,n,p,m)
          enddo
        enddo
      enddo
    enddo
  enddo
  deallocate(mo_tmp_1)

  allocate(mo_tmp_1(mo_num,mo_num,mo_num,ao_num))
  mo_tmp_1 = 0.d0
  do m = 1, ao_num
    do p = 1, ao_num
      do l = 1, mo_num
        do i = 1, mo_num
          do k = 1, mo_num
            mo_tmp_1(k,i,l,m) += mo_l_coef_transp(l,p) * mo_tmp_2(k,i,p,m)
          enddo
        enddo
      enddo
    enddo
  enddo
  deallocate(mo_tmp_2)

  mo_bi_ortho_tc_two_e_chemist = 0.d0 
  do m = 1, ao_num
    do j = 1, mo_num
      do l = 1, mo_num
        do i = 1, mo_num
          do k = 1, mo_num
            mo_bi_ortho_tc_two_e_chemist(k,i,l,j) += mo_r_coef_transp(j,m) * mo_tmp_1(k,i,l,m)
          enddo
        enddo
      enddo
    enddo
  enddo
  deallocate(mo_tmp_1)

END_PROVIDER 

! ---

BEGIN_PROVIDER [double precision, mo_bi_ortho_tc_two_e, (mo_num, mo_num, mo_num, mo_num)]

  BEGIN_DOC
  !
  ! mo_bi_ortho_tc_two_e(k,l,i,j) = <k l| V(r_12) |i j> where i,j are right MOs and k,l are left MOs
  !
  ! the potential V(r_12) contains ALL TWO-E CONTRIBUTION OF THE TC-HAMILTONIAN
  !
  END_DOC

  implicit none
  integer :: i, j, k, l

  do j = 1, mo_num
    do i = 1, mo_num
      do l = 1, mo_num
        do k = 1, mo_num
           !              < k l | V12 | i j >                          (k i|l j)
           mo_bi_ortho_tc_two_e(k,l,i,j) = mo_bi_ortho_tc_two_e_chemist(k,i,l,j)
        enddo
      enddo
    enddo
  enddo

END_PROVIDER 

! ---


 BEGIN_PROVIDER [ double precision, mo_bi_ortho_tc_two_e_jj, (mo_num,mo_num) ]
&BEGIN_PROVIDER [ double precision, mo_bi_ortho_tc_two_e_jj_exchange, (mo_num,mo_num) ]
&BEGIN_PROVIDER [ double precision, mo_bi_ortho_tc_two_e_jj_anti, (mo_num,mo_num) ]
  implicit none
  BEGIN_DOC
  ! mo_bi_ortho_tc_two_e_jj(i,j) = J_ij = <ji|W-K|ji>
  ! mo_bi_ortho_tc_two_e_jj_exchange(i,j) = K_ij = <ij|W-K|ji>
  ! mo_bi_ortho_tc_two_e_jj_anti(i,j) = J_ij - K_ij
  END_DOC

  integer                        :: i,j
  double precision               :: get_two_e_integral

  mo_bi_ortho_tc_two_e_jj = 0.d0
  mo_bi_ortho_tc_two_e_jj_exchange = 0.d0

  do i=1,mo_num
    do j=1,mo_num
      mo_bi_ortho_tc_two_e_jj(i,j) = mo_bi_ortho_tc_two_e(j,i,j,i)
      mo_bi_ortho_tc_two_e_jj_exchange(i,j) = mo_bi_ortho_tc_two_e(i,j,j,i)
      mo_bi_ortho_tc_two_e_jj_anti(i,j) = mo_bi_ortho_tc_two_e_jj(i,j) - mo_bi_ortho_tc_two_e_jj_exchange(i,j)
    enddo
  enddo

END_PROVIDER

 BEGIN_PROVIDER [double precision, tc_2e_3idx_coulomb_integrals, (mo_num,mo_num, mo_num)]
&BEGIN_PROVIDER [double precision, tc_2e_3idx_exchange_integrals,(mo_num,mo_num, mo_num)]
 implicit none
 BEGIN_DOC
 ! tc_2e_3idx_coulomb_integrals(j,k,i)  = <jk|ji> 
 !
 ! tc_2e_3idx_exchange_integrals(j,k,i) = <kj|ji> 
 END_DOC
 integer :: i,j,k,l
 double precision :: get_two_e_integral
 double precision :: integral

 do i = 1, mo_num
  do k = 1, mo_num
   do j = 1, mo_num
     tc_2e_3idx_coulomb_integrals(j, k,i) = mo_bi_ortho_tc_two_e(j ,k ,j ,i ) 
     tc_2e_3idx_exchange_integrals(j,k,i) = mo_bi_ortho_tc_two_e(k ,j ,j ,i ) 
   enddo
  enddo
 enddo

END_PROVIDER
