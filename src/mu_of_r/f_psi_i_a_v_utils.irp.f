subroutine give_f_ii_val_ab(r1,r2,f_ii_val_ab,two_bod_dens)
 implicit none
 BEGIN_DOC
! contribution from purely inactive orbitals to f_{\Psi^B}(r_1,r_2) for a CAS wave function 
 END_DOC
 double precision, intent(in) :: r1(3),r2(3)
 double precision, intent(out):: f_ii_val_ab,two_bod_dens
 integer :: i,j,m,n,i_m,i_n
 integer :: i_i,i_j
 double precision, allocatable :: mos_array_inact_r1(:),mos_array_inact_r2(:)
 double precision, allocatable :: mos_array_basis_r1(:),mos_array_basis_r2(:)
 double precision, allocatable :: mos_array_r1(:) , mos_array_r2(:) 
 double precision :: get_two_e_integral
 ! You get all orbitals in r_1 and r_2
 allocate(mos_array_r1(mo_num) , mos_array_r2(mo_num) )
 call give_all_mos_at_r(r1,mos_array_r1) 
 call give_all_mos_at_r(r2,mos_array_r2) 
 ! You extract the inactive orbitals 
 allocate(mos_array_inact_r1(n_inact_orb) , mos_array_inact_r2(n_inact_orb) )
 do i_m = 1, n_inact_orb
  mos_array_inact_r1(i_m) = mos_array_r1(list_inact(i_m))
 enddo
 do i_m = 1, n_inact_orb
  mos_array_inact_r2(i_m) = mos_array_r2(list_inact(i_m))
 enddo

 ! You extract the orbitals belonging to the space \mathcal{B}
 allocate(mos_array_basis_r1(n_basis_orb) , mos_array_basis_r2(n_basis_orb) )
 do i_m = 1, n_basis_orb 
  mos_array_basis_r1(i_m) = mos_array_r1(list_basis(i_m))
  mos_array_basis_r2(i_m) = mos_array_r2(list_basis(i_m))
 enddo

 f_ii_val_ab = 0.d0
 two_bod_dens = 0.d0
 ! You browse all OCCUPIED ALPHA electrons in the \mathcal{A} space
 do m = 1, n_inact_orb ! electron 1 
 ! You browse all OCCUPIED BETA  electrons in the \mathcal{A} space
  do n = 1, n_inact_orb ! electron 2 
   ! two_bod_dens(r_1,r_2) = n_alpha(r_1) * n_beta(r_2)
   two_bod_dens += mos_array_inact_r1(m) * mos_array_inact_r1(m) * mos_array_inact_r2(n) * mos_array_inact_r2(n) 
   ! You browse all COUPLE OF ORBITALS in the \mathacal{B} space 
   do i = 1, n_basis_orb
    do j = 1, n_basis_orb
     !                              2 1 2 1
     f_ii_val_ab +=  two_e_int_ii_f(j,i,n,m)   * mos_array_inact_r1(m) * mos_array_basis_r1(i)  & 
                                               * mos_array_inact_r2(n) * mos_array_basis_r2(j)    
    enddo
   enddo
  enddo
 enddo
end


subroutine give_f_ia_val_ab(r1,r2,f_ia_val_ab,two_bod_dens,istate)
 BEGIN_DOC
! contribution from inactive and active orbitals to f_{\Psi^B}(r_1,r_2) for the "istate" state of a CAS wave function 
 END_DOC
 implicit none
 integer, intent(in) :: istate
 double precision, intent(in) :: r1(3),r2(3)
 double precision, intent(out):: f_ia_val_ab,two_bod_dens
 integer :: i,orb_i,a,orb_a,n,m,b
 double precision :: rho
 double precision, allocatable :: mos_array_r1(:) , mos_array_r2(:) 
 double precision, allocatable :: mos_array_inact_r1(:),mos_array_inact_r2(:)
 double precision, allocatable :: mos_array_basis_r1(:),mos_array_basis_r2(:)
 double precision, allocatable :: mos_array_act_r1(:),mos_array_act_r2(:)
 double precision, allocatable :: integrals_array(:,:),rho_tilde(:,:),v_tilde(:,:)

 f_ia_val_ab = 0.d0
 two_bod_dens = 0.d0
 ! You get all orbitals in r_1 and r_2
 allocate(mos_array_r1(mo_num) , mos_array_r2(mo_num) )
 call give_all_mos_at_r(r1,mos_array_r1) 
 call give_all_mos_at_r(r2,mos_array_r2) 

 ! You extract the inactive orbitals 
 allocate( mos_array_inact_r1(n_inact_orb) , mos_array_inact_r2(n_inact_orb) )
 do i = 1, n_inact_orb
  mos_array_inact_r1(i) = mos_array_r1(list_inact(i))
 enddo
 do i= 1, n_inact_orb
  mos_array_inact_r2(i) = mos_array_r2(list_inact(i))
 enddo

 ! You extract the active orbitals 
 allocate( mos_array_act_r1(n_basis_orb) , mos_array_act_r2(n_basis_orb) )
 do i= 1, n_act_orb
  mos_array_act_r1(i) = mos_array_r1(list_act(i))
 enddo
 do i= 1, n_act_orb
  mos_array_act_r2(i) = mos_array_r2(list_act(i))
 enddo

 ! You extract the orbitals belonging to the space \mathcal{B}
 allocate( mos_array_basis_r1(n_basis_orb) , mos_array_basis_r2(n_basis_orb) )
 do i= 1, n_basis_orb
  mos_array_basis_r1(i) = mos_array_r1(list_basis(i))
 enddo
 do i= 1, n_basis_orb
  mos_array_basis_r2(i) = mos_array_r2(list_basis(i))
 enddo

 ! Contracted density : intermediate quantity
 ! rho_tilde(i,a) = \sum_b rho(b,a) * phi_i(1) * phi_j(2)
 allocate(rho_tilde(n_inact_orb,n_act_orb))
 two_bod_dens = 0.d0
 do a = 1, n_act_orb 
  do i = 1, n_inact_orb
   rho_tilde(i,a) = 0.d0
   do b = 1, n_act_orb 
    rho = one_e_act_dm_beta_mo_for_dft(b,a,istate) + one_e_act_dm_alpha_mo_for_dft(b,a,istate) 
    two_bod_dens += mos_array_inact_r1(i) * mos_array_inact_r1(i) * mos_array_act_r2(a) * mos_array_act_r2(b) * rho
    rho_tilde(i,a) += rho * mos_array_inact_r1(i) * mos_array_act_r2(b)
   enddo
  enddo
 enddo

 ! Contracted two-e integrals : intermediate quantity
 ! v_tilde(i,a) = \sum_{m,n} phi_m(1) * phi_n(2) < i a | m n >
 allocate( v_tilde(n_act_orb,n_act_orb)   )
 allocate( integrals_array(mo_num,mo_num) )
 v_tilde = 0.d0
 do a = 1, n_act_orb 
  orb_a = list_act(a)
  do i = 1, n_inact_orb
   v_tilde(i,a) = 0.d0
   orb_i = list_inact(i)
!   call get_mo_two_e_integrals_ij(orb_i,orb_a,mo_num,integrals_array,mo_integrals_map) 
   do m = 1, n_basis_orb
    do n = 1, n_basis_orb
!     v_tilde(i,a) += integrals_array(n,m) * mos_array_basis_r2(n) * mos_array_basis_r1(m)
     v_tilde(i,a) += two_e_int_ia_f(n,m,i,a) * mos_array_basis_r2(n) * mos_array_basis_r1(m)
    enddo
   enddo
  enddo
 enddo

 do a = 1, n_act_orb
  do i = 1, n_inact_orb
   f_ia_val_ab += v_tilde(i,a) * rho_tilde(i,a)
  enddo
 enddo
end


subroutine give_f_aa_val_ab(r1,r2,f_aa_val_ab,two_bod_dens,istate)
 BEGIN_DOC
! contribution from purely active orbitals to f_{\Psi^B}(r_1,r_2) for the "istate" state of a CAS wave function 
 END_DOC
 implicit none
 integer, intent(in) :: istate
 double precision, intent(in) :: r1(3),r2(3)
 double precision, intent(out):: f_aa_val_ab,two_bod_dens
 integer :: i,orb_i,a,orb_a,n,m,b,c,d
 double precision :: rho
 double precision, allocatable :: mos_array_r1(:) , mos_array_r2(:) 
 double precision, allocatable :: mos_array_basis_r1(:),mos_array_basis_r2(:)
 double precision, allocatable :: mos_array_act_r1(:),mos_array_act_r2(:)
 double precision, allocatable :: integrals_array(:,:),rho_tilde(:,:),v_tilde(:,:)

 f_aa_val_ab = 0.d0
 two_bod_dens = 0.d0
 ! You get all orbitals in r_1 and r_2
 allocate(mos_array_r1(mo_num) , mos_array_r2(mo_num) )
 call give_all_mos_at_r(r1,mos_array_r1) 
 call give_all_mos_at_r(r2,mos_array_r2) 

 ! You extract the active orbitals 
 allocate( mos_array_act_r1(n_basis_orb) , mos_array_act_r2(n_basis_orb) )
 do i= 1, n_act_orb
  mos_array_act_r1(i) = mos_array_r1(list_act(i))
 enddo
 do i= 1, n_act_orb
  mos_array_act_r2(i) = mos_array_r2(list_act(i))
 enddo

 ! You extract the orbitals belonging to the space \mathcal{B}
 allocate( mos_array_basis_r1(n_basis_orb) , mos_array_basis_r2(n_basis_orb) )
 do i= 1, n_basis_orb
  mos_array_basis_r1(i) = mos_array_r1(list_basis(i))
 enddo
 do i= 1, n_basis_orb
  mos_array_basis_r2(i) = mos_array_r2(list_basis(i))
 enddo

 ! Contracted density : intermediate quantity
 ! rho_tilde(i,a) = \sum_b rho(b,a) * phi_i(1) * phi_j(2)
 allocate(rho_tilde(n_act_orb,n_act_orb))
 two_bod_dens = 0.d0
 rho_tilde = 0.d0
 do a = 1, n_act_orb  ! 1
  do b = 1, n_act_orb  ! 2
   do c = 1, n_act_orb  ! 1
    do d = 1, n_act_orb  ! 2
     rho = mos_array_act_r1(c) * mos_array_act_r2(d) * act_2_rdm_ab_mo(d,c,b,a,istate)
     rho_tilde(b,a) += rho 
     two_bod_dens += rho * mos_array_act_r1(a) * mos_array_act_r2(b) 
    enddo
   enddo
  enddo
 enddo

 ! Contracted two-e integrals : intermediate quantity
 ! v_tilde(i,a) = \sum_{m,n} phi_m(1) * phi_n(2) < i a | m n >
 allocate( v_tilde(n_act_orb,n_act_orb)   )
 v_tilde = 0.d0
 do a = 1, n_act_orb 
  do b = 1, n_act_orb
   v_tilde(b,a) = 0.d0
   do m = 1, n_basis_orb
    do n = 1, n_basis_orb
     v_tilde(b,a) += two_e_int_aa_f(n,m,b,a) * mos_array_basis_r2(n) * mos_array_basis_r1(m)
    enddo
   enddo
  enddo
 enddo

 do a = 1, n_act_orb
  do b = 1, n_act_orb
   f_aa_val_ab += v_tilde(b,a) * rho_tilde(b,a)
  enddo
 enddo
end

BEGIN_PROVIDER [double precision, two_e_int_aa_f, (n_basis_orb,n_basis_orb,n_act_orb,n_act_orb)]
 implicit none
 BEGIN_DOC
! list of two-electron integrals (built with the MOs belonging to the \mathcal{B} space) 
! 
! needed to compute the function f_{ii}(r_1,r_2) 
!
! two_e_int_aa_f(j,i,n,m) = < j i | n m > where all orbitals belong to "list_basis"
 END_DOC
 integer :: orb_i,orb_j,i,j,orb_m,orb_n,m,n
 double precision :: integrals_array(mo_num,mo_num),get_two_e_integral
 do orb_m = 1, n_act_orb ! electron 1 
  m = list_act(orb_m)
  do orb_n = 1, n_act_orb ! electron 2 
   n = list_act(orb_n)
   call get_mo_two_e_integrals_ij(m,n,mo_num,integrals_array,mo_integrals_map) 
   do orb_i = 1, n_basis_orb ! electron 1 
    i = list_basis(orb_i)
    do orb_j = 1, n_basis_orb ! electron 2 
     j = list_basis(orb_j)
     !              2       1     2     1
     two_e_int_aa_f(orb_j,orb_i,orb_n,orb_m) = get_two_e_integral(m,n,i,j,mo_integrals_map) 
!     two_e_int_aa_f(orb_j,orb_i,orb_n,orb_m) = integrals_array(j,i) 
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [double precision, two_e_int_ia_f, (n_basis_orb,n_basis_orb,n_inact_orb,n_act_orb)]
 implicit none
 BEGIN_DOC
! list of two-electron integrals (built with the MOs belonging to the \mathcal{B} space) 
! 
! needed to compute the function f_{ia}(r_1,r_2) 
!
! two_e_int_aa_f(j,i,n,m) = < j i | n m > where all orbitals belong to "list_basis"
 END_DOC
 integer :: orb_i,orb_j,i,j,orb_m,orb_n,m,n
 double precision :: integrals_array(mo_num,mo_num),get_two_e_integral
 do orb_m = 1, n_act_orb ! electron 1 
  m = list_act(orb_m)
  do orb_n = 1, n_inact_orb ! electron 2 
   n = list_inact(orb_n)
   call get_mo_two_e_integrals_ij(m,n,mo_num,integrals_array,mo_integrals_map) 
   do orb_i = 1, n_basis_orb ! electron 1 
    i = list_basis(orb_i)
    do orb_j = 1, n_basis_orb ! electron 2 
     j = list_basis(orb_j)
     !              2       1     2     1
!     two_e_int_ia_f(orb_j,orb_i,orb_n,orb_m) = get_two_e_integral(m,n,i,j,mo_integrals_map) 
     two_e_int_ia_f(orb_j,orb_i,orb_n,orb_m) = integrals_array(j,i) 
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [double precision, two_e_int_ii_f, (n_basis_orb,n_basis_orb,n_inact_orb,n_inact_orb)]
 implicit none
 BEGIN_DOC
! list of two-electron integrals (built with the MOs belonging to the \mathcal{B} space) 
! 
! needed to compute the function f_{ii}(r_1,r_2) 
!
! two_e_int_ii_f(j,i,n,m) = < j i | n m > where all orbitals belong to "list_basis"
 END_DOC
 integer :: orb_i,orb_j,i,j,orb_m,orb_n,m,n
 double precision :: get_two_e_integral,integrals_array(mo_num,mo_num)
 do orb_m = 1, n_inact_orb ! electron 1 
  m = list_inact(orb_m)
  do orb_n = 1, n_inact_orb ! electron 2 
   n = list_inact(orb_n)
   call get_mo_two_e_integrals_ij(m,n,mo_num,integrals_array,mo_integrals_map) 
   do orb_i = 1, n_basis_orb ! electron 1 
    i = list_basis(orb_i)
    do orb_j = 1, n_basis_orb ! electron 2 
     j = list_basis(orb_j)
     !              2       1     2     1
!     two_e_int_ii_f(orb_j,orb_i,orb_n,orb_m) = get_two_e_integral(m,n,i,j,mo_integrals_map) 
     two_e_int_ii_f(orb_j,orb_i,orb_n,orb_m) = integrals_array(j,i) 
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 

