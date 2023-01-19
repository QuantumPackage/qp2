  use bitmasks
BEGIN_PROVIDER [integer(bit_kind), ref_closed_shell_bitmask, (N_int,2)]
 implicit none
 integer :: i,i0
 integer :: n_occ_ab(2)
 integer :: occ(N_int*bit_kind_size,2)
 call bitstring_to_list_ab(ref_bitmask, occ, n_occ_ab, N_int)
 ! do the closed shell determinant
 do i = 1, N_int
  ref_closed_shell_bitmask(i,1) = ref_bitmask(i,1)
  ref_closed_shell_bitmask(i,2) = ref_bitmask(i,2)
 enddo
 do i0 = elec_beta_num+1, elec_alpha_num
   i=occ(i0,1)
   call clear_bit_to_integer(i,ref_closed_shell_bitmask(1,1),N_int)
 enddo


END_PROVIDER


BEGIN_PROVIDER [double precision, fock_operator_closed_shell_ref_bitmask, (mo_num, mo_num) ]
 implicit none
 integer :: i0,j0,i,j,k0,k
 integer :: n_occ_ab(2)
 integer :: occ(N_int*bit_kind_size,2)
 integer :: n_occ_ab_virt(2)
 integer :: occ_virt(N_int*bit_kind_size,2)
 integer(bit_kind) :: key_test(N_int)
 integer(bit_kind) :: key_virt(N_int,2)
 fock_operator_closed_shell_ref_bitmask = 0.d0
 call bitstring_to_list_ab(ref_closed_shell_bitmask, occ, n_occ_ab, N_int)
 do i = 1, N_int
  key_virt(i,1) = full_ijkl_bitmask(i)
  key_virt(i,2) = full_ijkl_bitmask(i)
  key_virt(i,1) = xor(key_virt(i,1),ref_closed_shell_bitmask(i,1))
  key_virt(i,2) = xor(key_virt(i,2),ref_closed_shell_bitmask(i,2))
 enddo
 double precision, allocatable :: array_coulomb(:),array_exchange(:)
 allocate (array_coulomb(mo_num),array_exchange(mo_num))
 call bitstring_to_list_ab(key_virt, occ_virt, n_occ_ab_virt, N_int)
 ! docc ---> virt single excitations
 do i0 = 1,  n_occ_ab(1)
  i=occ(i0,1)
  do j0 = 1, n_occ_ab_virt(1)
   j = occ_virt(j0,1)
   call get_mo_two_e_integrals_coulomb_ii(i,j,mo_num,array_coulomb,mo_integrals_map)
   call get_mo_two_e_integrals_exch_ii(i,j,mo_num,array_exchange,mo_integrals_map)
   double precision :: accu
   accu = 0.d0
   do k0 = 1, n_occ_ab(1)
    k = occ(k0,1)
    accu += 2.d0 * array_coulomb(k) - array_exchange(k)
   enddo
   fock_operator_closed_shell_ref_bitmask(i,j) = accu + mo_one_e_integrals(i,j)
   fock_operator_closed_shell_ref_bitmask(j,i) = accu + mo_one_e_integrals(i,j)
  enddo
 enddo

 ! virt ---> virt single excitations
 do i0 = 1,  n_occ_ab_virt(1)
  i=occ_virt(i0,1)
  do j0 = 1, n_occ_ab_virt(1)
   j = occ_virt(j0,1)
   call get_mo_two_e_integrals_coulomb_ii(i,j,mo_num,array_coulomb,mo_integrals_map)
   call get_mo_two_e_integrals_exch_ii(i,j,mo_num,array_exchange,mo_integrals_map)
   accu = 0.d0
   do k0 = 1, n_occ_ab(1)
    k = occ(k0,1)
    accu += 2.d0 * array_coulomb(k) - array_exchange(k)
   enddo
   fock_operator_closed_shell_ref_bitmask(i,j) = accu+ mo_one_e_integrals(i,j)
   fock_operator_closed_shell_ref_bitmask(j,i) = accu+ mo_one_e_integrals(i,j)
  enddo
 enddo

 ! docc ---> docc single excitations
 do i0 = 1,  n_occ_ab(1)
  i=occ(i0,1)
  do j0 = 1, n_occ_ab(1)
   j = occ(j0,1)
   call get_mo_two_e_integrals_coulomb_ii(i,j,mo_num,array_coulomb,mo_integrals_map)
   call get_mo_two_e_integrals_exch_ii(i,j,mo_num,array_exchange,mo_integrals_map)
   accu = 0.d0
   do k0 = 1, n_occ_ab(1)
    k = occ(k0,1)
    accu += 2.d0 * array_coulomb(k) - array_exchange(k)
   enddo
   fock_operator_closed_shell_ref_bitmask(i,j) = accu+ mo_one_e_integrals(i,j)
   fock_operator_closed_shell_ref_bitmask(j,i) = accu+ mo_one_e_integrals(i,j)
  enddo
 enddo
 deallocate(array_coulomb,array_exchange)

END_PROVIDER

subroutine get_single_excitation_from_fock(det_1,det_2,h,p,spin,phase,hij)
 use bitmasks
 implicit none
 integer,intent(in) :: h,p,spin
 double precision, intent(in)  :: phase
 integer(bit_kind), intent(in) :: det_1(N_int,2), det_2(N_int,2)
 double precision, intent(out) :: hij
 integer(bit_kind) :: differences(N_int,2)
 integer(bit_kind) :: hole(N_int,2)
 integer(bit_kind) :: partcl(N_int,2)
 integer :: occ_hole(N_int*bit_kind_size,2)
 integer :: occ_partcl(N_int*bit_kind_size,2)
 integer :: n_occ_ab_hole(2),n_occ_ab_partcl(2)
 integer :: i0,i
 double precision :: buffer_c(mo_num),buffer_x(mo_num)
 do i=1, mo_num
   buffer_c(i) = big_array_coulomb_integrals(i,h,p)
   buffer_x(i) = big_array_exchange_integrals(i,h,p)
 enddo
 do i = 1, N_int
  differences(i,1) = xor(det_1(i,1),ref_closed_shell_bitmask(i,1))
  differences(i,2) = xor(det_1(i,2),ref_closed_shell_bitmask(i,2))
  hole(i,1) = iand(differences(i,1),ref_closed_shell_bitmask(i,1))
  hole(i,2) = iand(differences(i,2),ref_closed_shell_bitmask(i,2))
  partcl(i,1) = iand(differences(i,1),det_1(i,1))
  partcl(i,2) = iand(differences(i,2),det_1(i,2))
 enddo
 call bitstring_to_list_ab(hole, occ_hole, n_occ_ab_hole, N_int)
 call bitstring_to_list_ab(partcl, occ_partcl, n_occ_ab_partcl, N_int)
 hij = fock_operator_closed_shell_ref_bitmask(h,p)
 ! holes :: direct terms
 do i0 = 1, n_occ_ab_hole(1)
  i = occ_hole(i0,1)
  hij -= buffer_c(i)
 enddo
 do i0 = 1, n_occ_ab_hole(2)
  i = occ_hole(i0,2)
  hij -= buffer_c(i)
 enddo

 ! holes :: exchange terms
 do i0 = 1, n_occ_ab_hole(spin)
  i = occ_hole(i0,spin)
  hij += buffer_x(i)
 enddo

 ! particles :: direct terms
 do i0 = 1, n_occ_ab_partcl(1)
  i = occ_partcl(i0,1)
  hij += buffer_c(i)
 enddo
 do i0 = 1, n_occ_ab_partcl(2)
  i = occ_partcl(i0,2)
  hij += buffer_c(i)
 enddo

 ! particles :: exchange terms
 do i0 = 1, n_occ_ab_partcl(spin)
  i = occ_partcl(i0,spin)
  hij -= buffer_x(i)
 enddo
 hij = hij * phase

end

