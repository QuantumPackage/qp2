BEGIN_PROVIDER [double precision, two_e_int_hf_f, (n_basis_orb,n_basis_orb,n_max_occ_val_orb_for_hf,n_max_occ_val_orb_for_hf)]
 implicit none
 BEGIN_DOC
! list of two-electron integrals (built with the MOs belonging to the \mathcal{B} space) 
! 
! needed to compute the function f_{HF}(r_1,r_2) 
!
! two_e_int_hf_f(j,i,n,m) = < j i | n m > where all orbitals belong to "list_basis"
 END_DOC
 integer :: orb_i,orb_j,i,j,orb_m,orb_n,m,n
 double precision :: get_two_e_integral
 do orb_m = 1, n_max_occ_val_orb_for_hf! electron 1 
  m = list_valence_orb_for_hf(orb_m,1)
  do orb_n = 1, n_max_occ_val_orb_for_hf! electron 2 
   n = list_valence_orb_for_hf(orb_n,1)
   do orb_i = 1, n_basis_orb ! electron 1 
    i = list_basis(orb_i)
    do orb_j = 1, n_basis_orb ! electron 2 
     j = list_basis(orb_j)
     !              2       1     2     1
     two_e_int_hf_f(orb_j,orb_i,orb_n,orb_m) = get_two_e_integral(m,n,i,j,mo_integrals_map) 
    enddo
   enddo
  enddo
 enddo
END_PROVIDER 

subroutine f_HF_valence_ab(r1,r2,f_HF_val_ab,two_bod_dens)
 implicit none
 BEGIN_DOC
! f_HF_val_ab(r1,r2) = function f_{\Psi^B}(X_1,X_2) of Eq. (22) of J. Chem. Phys. 149, 194301 (2018) 
! 
! for alpha beta spins and an HF wave function and excluding the "core" orbitals (see Eq. 16a of Phys.Chem.Lett.2019, 10, 2931   2937)
! 
! two_bod_dens = on-top pair density of the HF wave function 
!
! < HF | wee_{\alpha\beta} | HF > =  \int (r1,r2) f_HF_ab(r1,r2) excluding all contributions from "core" "electrons" 
 END_DOC
 double precision, intent(in) :: r1(3), r2(3)
 double precision, intent(out):: f_HF_val_ab,two_bod_dens
 integer :: i,j,m,n,i_m,i_n
 integer :: i_i,i_j
 double precision :: mo_two_e_integral
 double precision, allocatable :: mos_array_r1(:)
 double precision, allocatable :: mos_array_r2(:)
 double precision, allocatable :: mos_array_valence_r1(:),mos_array_valence_r2(:)
 double precision, allocatable :: mos_array_valence_hf_r1(:),mos_array_valence_hf_r2(:)
 double precision :: get_two_e_integral
 allocate(mos_array_valence_r1(n_basis_orb) , mos_array_valence_r2(n_basis_orb), mos_array_r1(mo_num), mos_array_r2(mo_num))
 allocate(mos_array_valence_hf_r1(n_occ_val_orb_for_hf(1)) , mos_array_valence_hf_r2(n_occ_val_orb_for_hf(2)) )
 ! You get all orbitals in r_1 and r_2
 call give_all_mos_at_r(r1,mos_array_r1) 
 call give_all_mos_at_r(r2,mos_array_r2) 
 ! You extract the occupied ALPHA/BETA orbitals belonging to the space \mathcal{A}
 do i_m = 1, n_occ_val_orb_for_hf(1)
  mos_array_valence_hf_r1(i_m) = mos_array_r1(list_valence_orb_for_hf(i_m,1))
 enddo
 do i_m = 1, n_occ_val_orb_for_hf(2)
  mos_array_valence_hf_r2(i_m) = mos_array_r2(list_valence_orb_for_hf(i_m,2))
 enddo

 ! You extract the orbitals belonging to the space \mathcal{B}
 do i_m = 1, n_basis_orb 
  mos_array_valence_r1(i_m) = mos_array_r1(list_basis(i_m))
  mos_array_valence_r2(i_m) = mos_array_r2(list_basis(i_m))
 enddo


 f_HF_val_ab = 0.d0
 two_bod_dens = 0.d0
 ! You browse all OCCUPIED ALPHA electrons in the \mathcal{A} space
 do m = 1, n_occ_val_orb_for_hf(1)! electron 1 
 ! You browse all OCCUPIED BETA  electrons in the \mathcal{A} space
  do n = 1, n_occ_val_orb_for_hf(2)! electron 2 
   ! two_bod_dens(r_1,r_2) = n_alpha(r_1) * n_beta(r_2)
   two_bod_dens += mos_array_valence_hf_r1(m) * mos_array_valence_hf_r1(m) * mos_array_valence_hf_r2(n) * mos_array_valence_hf_r2(n) 
   ! You browse all COUPLE OF ORBITALS in the \mathacal{B} space 
   do i = 1, n_basis_orb
    do j = 1, n_basis_orb
     !                                             2 1 2 1
     f_HF_val_ab +=  two_e_int_hf_f(j,i,n,m)  & 
     * mos_array_valence_r1(i) * mos_array_valence_hf_r1(m)  & 
     * mos_array_valence_r2(j) * mos_array_valence_hf_r2(n)    
    enddo
   enddo
  enddo
 enddo
end


subroutine integral_f_HF_valence_ab(r1,int_f_HF_val_ab)
 implicit none
 BEGIN_DOC
! in_f_HF_val_ab(r_1) = \int dr_2 f_{\Psi^B}(r_1,r_2) 
! 
! where f_{\Psi^B}(r_1,r_2) is defined by Eq. (22) of J. Chem. Phys. 149, 194301 (2018) 
!
! for alpha beta spins and an HF wave function and excluding the "core" orbitals (see Eq. 16a of Phys.Chem.Lett.2019, 10, 2931   2937)
!
! Such function can be used to test if the f_HF_val_ab(r_1,r_2) is correctly built. 
! 
! < HF | wee_{\alpha\beta} | HF > =  \int (r1) int_f_HF_val_ab(r_1)
 END_DOC
 double precision, intent(in) :: r1(3)
 double precision, intent(out):: int_f_HF_val_ab
 integer :: i,j,m,n,i_m,i_n
 integer :: i_i,i_j
 double precision :: mo_two_e_integral
 double precision :: mos_array_r1(mo_num)
 double precision, allocatable  :: mos_array_valence_r1(:)
 double precision, allocatable  :: mos_array_valence_hf_r1(:)
 double precision :: get_two_e_integral
 call give_all_mos_at_r(r1,mos_array_r1) 
 allocate(mos_array_valence_r1( n_basis_orb ))
 allocate(mos_array_valence_hf_r1( n_occ_val_orb_for_hf(1) ) )
 do i_m = 1, n_occ_val_orb_for_hf(1)
  mos_array_valence_hf_r1(i_m) = mos_array_r1(list_valence_orb_for_hf(i_m,1))
 enddo
 do i_m = 1, n_basis_orb 
  mos_array_valence_r1(i_m) = mos_array_r1(list_basis(i_m))
 enddo

 int_f_HF_val_ab = 0.d0
 ! You browse all OCCUPIED ALPHA electrons in the \mathcal{A} space
 do m = 1, n_occ_val_orb_for_hf(1)! electron 1 
 ! You browse all OCCUPIED BETA  electrons in the \mathcal{A} space
  do n = 1, n_occ_val_orb_for_hf(2)! electron 2 
   ! You browse all ORBITALS in the \mathacal{B} space 
   do i = 1, n_basis_orb
   ! due to integration in real-space and the use of orthonormal MOs, a Kronecker delta_jn shoes up 
    j = n
     !                                             2 1 2 1
     int_f_HF_val_ab +=  two_e_int_hf_f(j,i,n,m)  & 
     * mos_array_valence_r1(i) * mos_array_valence_hf_r1(m) 
   enddo
  enddo
 enddo
end
