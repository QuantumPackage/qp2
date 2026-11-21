  use bitmasks
subroutine single_excitation_wee(det_1,det_2,h,p,spin,phase,hij)
 use bitmasks
 implicit none
  BEGIN_DOC
  ! For a determinant det_1 that differs from the closed-shell
  ! reference (ref_closed_shell_bitmask), computes the
  ! two-electron part of the Hamiltonian matrix element
  ! <det_1|H|det_2>, where det_2 is the single excitation (h, p)
  ! of spin `spin` from det_1:
  ! hij = fock_wee_closed_shell(h, p)
  !       - Coulomb integrals of the holes
  !       + exchange integrals of the same-spin holes
  !       + Coulomb integrals of the particles
  !       - exchange integrals of the same-spin particles,
  ! times phase.
  ! Holes are the spin orbitals occupied in the reference but
  ! not in det_1, particles the opposite. det_2 is a single
  ! excitation from det_1 (cf. its use in i_H_j_s2); it is not
  ! used because the hole and particle are already given as
  ! parameters (author-corrected 2026-10-31).
  END_DOC
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
 hij = fock_wee_closed_shell(h,p)
 ! holes :: direct terms
 do i0 = 1, n_occ_ab_hole(1)
  i = occ_hole(i0,1)
  hij -= big_array_coulomb_integrals(i,h,p) ! get_mo_two_e_integral_schwartz(h,i,p,i,mo_integrals_map)
 enddo
 do i0 = 1, n_occ_ab_hole(2)
  i = occ_hole(i0,2)
  hij -= big_array_coulomb_integrals(i,h,p) !get_mo_two_e_integral_schwartz(h,i,p,i,mo_integrals_map)
 enddo

 ! holes :: exchange terms
 do i0 = 1, n_occ_ab_hole(spin)
  i = occ_hole(i0,spin)
  hij += big_array_exchange_integrals(i,h,p) ! get_mo_two_e_integral_schwartz(h,i,i,p,mo_integrals_map)
 enddo

 ! particles :: direct terms
 do i0 = 1, n_occ_ab_partcl(1)
  i = occ_partcl(i0,1)
  hij += big_array_coulomb_integrals(i,h,p)!get_mo_two_e_integral_schwartz(h,i,p,i,mo_integrals_map)
 enddo
 do i0 = 1, n_occ_ab_partcl(2)
  i = occ_partcl(i0,2)
  hij += big_array_coulomb_integrals(i,h,p) !get_mo_two_e_integral_schwartz(h,i,p,i,mo_integrals_map)
 enddo

 ! particles :: exchange terms
 do i0 = 1, n_occ_ab_partcl(spin)
  i = occ_partcl(i0,spin)
  hij -= big_array_exchange_integrals(i,h,p)!get_mo_two_e_integral_schwartz(h,i,i,p,mo_integrals_map)
 enddo
 hij = hij * phase

end


BEGIN_PROVIDER [double precision, fock_wee_closed_shell, (mo_num, mo_num) ]
 implicit none
  BEGIN_DOC
  ! Two-electron part of the Fock matrix of the closed-shell
  ! reference, for orbitals i and j in the union of the
  ! occupied and virtual orbitals:
  ! fock_wee_closed_shell(i,j) = sum_k [ 2 (i j | k k) - (i k | j k) ]
  ! with k running over the occupied orbitals.
  ! Only the alpha orbitals are used to build the occupied
  ! and virtual lists, as the reference is closed-shell.
  ! Entries for other pairs of orbitals are not set.
  END_DOC
 integer :: i0,j0,i,j,k0,k
 integer :: n_occ_ab(2)
 integer :: occ(N_int*bit_kind_size,2)
 integer :: n_occ_ab_virt(2)
 integer :: occ_virt(N_int*bit_kind_size,2)
 integer(bit_kind) :: key_test(N_int)
 integer(bit_kind) :: key_virt(N_int,2)

 call bitstring_to_list_ab(ref_closed_shell_bitmask, occ, n_occ_ab, N_int)
 do i = 1, N_int
  key_virt(i,1) = full_ijkl_bitmask(i)
  key_virt(i,2) = full_ijkl_bitmask(i)
  key_virt(i,1) = xor(key_virt(i,1),ref_closed_shell_bitmask(i,1))
  key_virt(i,2) = xor(key_virt(i,2),ref_closed_shell_bitmask(i,2))
 enddo
 call bitstring_to_list_ab(key_virt, occ_virt, n_occ_ab_virt, N_int)
 ! docc ---> virt single excitations
 do i0 = 1,  n_occ_ab(1)
  i=occ(i0,1)
  do j0 = 1, n_occ_ab_virt(1)
   j = occ_virt(j0,1)
   double precision :: accu
   accu = 0.d0
   do k0 = 1, n_occ_ab(1)
    k = occ(k0,1)
    accu += 2.d0 * big_array_coulomb_integrals(k,i,j) - big_array_exchange_integrals(k,i,j)
   enddo
   fock_wee_closed_shell(i,j) = accu
   fock_wee_closed_shell(j,i) = accu
  enddo
 enddo

 ! virt ---> virt single excitations
 do i0 = 1,  n_occ_ab_virt(1)
  i=occ_virt(i0,1)
  do j0 = 1, n_occ_ab_virt(1)
   j = occ_virt(j0,1)
   accu = 0.d0
   do k0 = 1, n_occ_ab(1)
    k = occ(k0,1)
    accu += 2.d0 * big_array_coulomb_integrals(k,i,j) - big_array_exchange_integrals(k,i,j)
   enddo
   fock_wee_closed_shell(i,j) = accu
   fock_wee_closed_shell(j,i) = accu
  enddo
 enddo

 ! docc ---> docc single excitations
 do i0 = 1,  n_occ_ab(1)
  i=occ(i0,1)
  do j0 = 1, n_occ_ab(1)
   j = occ(j0,1)
   accu = 0.d0
   do k0 = 1, n_occ_ab(1)
    k = occ(k0,1)
    accu += 2.d0 * big_array_coulomb_integrals(k,i,j) - big_array_exchange_integrals(k,i,j)
   enddo
   fock_wee_closed_shell(i,j) = accu
   fock_wee_closed_shell(j,i) = accu
  enddo
 enddo

END_PROVIDER


