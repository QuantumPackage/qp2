subroutine give_all_perm_for_three_e(n,l,k,m,j,i,idx_list,phase)
 implicit none
 BEGIN_DOC
 ! returns all the list of permutting indices for the antimmetrization of 
 !
 ! (k^dagger l^dagger n^dagger m j i)  <nlk|L|mji>   when all indices have the same spins
 !
 ! idx_list(:,i) == list of the 6 indices corresponding the permutation "i"
 !
 ! phase(i) == phase of the permutation "i" 
 !
 ! there are in total 6 permutations with different indices
 END_DOC
 integer, intent(in)  :: n,l,k,m,j,i
 integer, intent(out) :: idx_list(6,6) 
 double precision :: phase(6)
 integer :: list(6)
 !!! CYCLIC PERMUTATIONS 
 phase(1:3) = 1.d0
 !!! IDENTITY PERMUTATION 
 list = (/n,l,k,m,j,i/)
 idx_list(:,1) = list(:)
 !!! FIRST CYCLIC PERMUTATION 
 list = (/n,l,k,j,i,m/)
 idx_list(:,2) = list(:)
 !!! FIRST CYCLIC PERMUTATION 
 list = (/n,l,k,i,m,j/)
 idx_list(:,3) = list(:)

 !!! NON CYCLIC PERMUTATIONS 
 phase(1:3) = -1.d0
 !!! PARTICLE 1 is FIXED
 list = (/n,l,k,j,m,i/)
 idx_list(:,4) = list(:)
 !!! PARTICLE 2 is FIXED
 list = (/n,l,k,i,j,m/)
 idx_list(:,5) = list(:)
 !!! PARTICLE 3 is FIXED
 list = (/n,l,k,m,i,j/)
 idx_list(:,6) = list(:)

end

double precision function sym_3_e_int_from_6_idx_tensor(n,l,k,m,j,i)
 implicit none
 BEGIN_DOC
 ! returns all good combinations of permutations of integrals with the good signs 
 !
 ! for a given (k^dagger l^dagger n^dagger m j i)  <nlk|L|mji> when all indices have the same spins
 END_DOC
 integer, intent(in)  :: n,l,k,m,j,i
 sym_3_e_int_from_6_idx_tensor = three_body_ints_bi_ort(n,l,k,m,j,i) & ! direct 
                               + three_body_ints_bi_ort(n,l,k,j,i,m) & ! 1st cyclic permutation  
                               + three_body_ints_bi_ort(n,l,k,i,m,j) & ! 2nd cyclic permutation  
                               - three_body_ints_bi_ort(n,l,k,j,m,i) & ! elec 1 is kept fixed 
                               - three_body_ints_bi_ort(n,l,k,i,j,m) & ! elec 2 is kept fixed
                               - three_body_ints_bi_ort(n,l,k,m,i,j)   ! elec 3 is kept fixed

end

double precision function direct_sym_3_e_int(n,l,k,m,j,i)
 implicit none
 BEGIN_DOC
 ! returns all good combinations of permutations of integrals with the good signs 
 !
 ! for a given (k^dagger l^dagger n^dagger m j i)  <nlk|L|mji> when all indices have the same spins
 END_DOC
 integer, intent(in)  :: n,l,k,m,j,i
 double precision :: integral
 direct_sym_3_e_int = 0.d0
 call give_integrals_3_body_bi_ort(n,l,k,m,j,i,integral)   ! direct 
 direct_sym_3_e_int += integral
 call give_integrals_3_body_bi_ort(n,l,k,j,i,m,integral)   ! 1st cyclic permutation  
 direct_sym_3_e_int += integral
 call give_integrals_3_body_bi_ort(n,l,k,i,m,j,integral)   ! 2nd cyclic permutation  
 direct_sym_3_e_int += integral
 call give_integrals_3_body_bi_ort(n,l,k,j,m,i,integral)   ! elec 1 is kept fixed 
 direct_sym_3_e_int += -integral
 call give_integrals_3_body_bi_ort(n,l,k,i,j,m,integral)   ! elec 2 is kept fixed
 direct_sym_3_e_int += -integral
 call give_integrals_3_body_bi_ort(n,l,k,m,i,j,integral)   ! elec 3 is kept fixed
 direct_sym_3_e_int += -integral

end

double precision function three_e_diag_parrallel_spin(m,j,i)
 implicit none
 integer, intent(in) :: i,j,m
  three_e_diag_parrallel_spin = three_e_3_idx_direct_bi_ort(m,j,i)  ! direct
  three_e_diag_parrallel_spin += three_e_3_idx_cycle_1_bi_ort(m,j,i) + three_e_3_idx_cycle_2_bi_ort(m,j,i) & ! two cyclic permutations 
  - three_e_3_idx_exch23_bi_ort(m,j,i) - three_e_3_idx_exch13_bi_ort(m,j,i)  & ! two first exchange 
  - three_e_3_idx_exch12_bi_ort(m,j,i) ! last exchange 
end

double precision function three_e_single_parrallel_spin(m,j,k,i)
 implicit none
 integer, intent(in) :: i,k,j,m
  three_e_single_parrallel_spin = three_e_4_idx_direct_bi_ort(m,j,k,i)  ! direct
  three_e_single_parrallel_spin += three_e_4_idx_cycle_1_bi_ort(m,j,k,i) + three_e_4_idx_cycle_2_bi_ort(m,j,k,i) & ! two cyclic permutations 
  - three_e_4_idx_exch23_bi_ort(m,j,k,i) - three_e_4_idx_exch13_bi_ort(m,j,k,i)  & ! two first exchange 
  - three_e_4_idx_exch12_bi_ort(m,j,k,i) ! last exchange 
end

double precision function three_e_double_parrallel_spin(m,l,j,k,i)
 implicit none
 integer, intent(in) :: i,k,j,m,l
  three_e_double_parrallel_spin = three_e_5_idx_direct_bi_ort(m,l,j,k,i)  ! direct
  three_e_double_parrallel_spin += three_e_5_idx_cycle_1_bi_ort(m,l,j,k,i) + three_e_5_idx_cycle_2_bi_ort(m,l,j,k,i) & ! two cyclic permutations 
  - three_e_5_idx_exch23_bi_ort(m,l,j,k,i) - three_e_5_idx_exch13_bi_ort(m,l,j,k,i)  & ! two first exchange 
  - three_e_5_idx_exch12_bi_ort(m,l,j,k,i) ! last exchange 
end
