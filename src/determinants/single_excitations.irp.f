  use bitmasks
BEGIN_PROVIDER [integer(bit_kind), ref_closed_shell_bitmask, (N_int,2)]
 implicit none
 integer :: i,i0,k
 integer :: n_occ_ab(2)
 integer :: occ(N_int*bit_kind_size,2)
 call bitstring_to_list_ab(ref_bitmask, occ, n_occ_ab, N_int)
 ! do the closed shell determinant
 do i = 1, N_int
  ref_closed_shell_bitmask(i,1) = ref_bitmask(i,1)
  ref_closed_shell_bitmask(i,2) = ref_bitmask(i,2)
 enddo
 if (is_complex) then
   !todo: check this
   do k=1,kpt_num
     call bitstring_to_list_ab(ref_bitmask_kpts(1,1,k),occ,n_occ_ab,N_int)
     do i0=elec_beta_num_kpts(k)+1,elec_alpha_num_kpts(k)
       i=occ(i0,1)
       call clear_bit_to_integer(i,ref_closed_shell_bitmask(1,1),N_int)
     enddo
   enddo
 else
   do i0 = elec_beta_num+1, elec_alpha_num
     i=occ(i0,1)
     call clear_bit_to_integer(i,ref_closed_shell_bitmask(1,1),N_int)
   enddo
 endif
END_PROVIDER

BEGIN_PROVIDER [double precision, fock_op_cshell_ref_bitmask, (mo_num, mo_num) ]
 implicit none
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
   fock_op_cshell_ref_bitmask(i,j) = accu + mo_one_e_integrals(i,j)
   fock_op_cshell_ref_bitmask(j,i) = accu + mo_one_e_integrals(i,j)
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
   fock_op_cshell_ref_bitmask(i,j) = accu+ mo_one_e_integrals(i,j)
   fock_op_cshell_ref_bitmask(j,i) = accu+ mo_one_e_integrals(i,j)
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
   fock_op_cshell_ref_bitmask(i,j) = accu+ mo_one_e_integrals(i,j)
   fock_op_cshell_ref_bitmask(j,i) = accu+ mo_one_e_integrals(i,j)
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
 hij = fock_op_cshell_ref_bitmask(h,p)
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

!============================================!
!                                            !
!                 complex                    !
!                                            !
!============================================!

BEGIN_PROVIDER [complex*16, fock_op_cshell_ref_bitmask_cplx, (mo_num, mo_num) ]
 implicit none
 integer :: i0,j0,i,j,k0,k
 integer :: n_occ_ab(2)
 integer :: occ(N_int*bit_kind_size,2)
 integer :: n_occ_ab_virt(2)
 integer :: occ_virt(N_int*bit_kind_size,2)
 integer(bit_kind) :: key_test(N_int)
 integer(bit_kind) :: key_virt(N_int,2)
 complex*16 :: accu

 call bitstring_to_list_ab(ref_closed_shell_bitmask, occ, n_occ_ab, N_int)
 do i = 1, N_int
  key_virt(i,1) = full_ijkl_bitmask(i)
  key_virt(i,2) = full_ijkl_bitmask(i)
  key_virt(i,1) = xor(key_virt(i,1),ref_closed_shell_bitmask(i,1))
  key_virt(i,2) = xor(key_virt(i,2),ref_closed_shell_bitmask(i,2))
 enddo
 complex*16, allocatable :: array_coulomb(:),array_exchange(:)
 allocate (array_coulomb(mo_num),array_exchange(mo_num))
 call bitstring_to_list_ab(key_virt, occ_virt, n_occ_ab_virt, N_int)
 ! docc ---> virt single excitations
 do i0 = 1,  n_occ_ab(1)
  i=occ(i0,1)
  do j0 = 1, n_occ_ab_virt(1)
   j = occ_virt(j0,1)
   ! <ia|ja>
   call get_mo_two_e_integrals_coulomb_ii_complex(i,j,mo_num,array_coulomb,mo_integrals_map,mo_integrals_map_2)
   ! <ia|aj>
   call get_mo_two_e_integrals_exch_ii_complex(i,j,mo_num,array_exchange,mo_integrals_map,mo_integrals_map_2)
   accu = (0.d0,0.d0)
   do k0 = 1, n_occ_ab(1)
    k = occ(k0,1)
    accu += 2.d0 * array_coulomb(k) - array_exchange(k)
   enddo
   fock_op_cshell_ref_bitmask_cplx(i,j) = accu + mo_one_e_integrals_complex(i,j)
   !fock_op_cshell_ref_bitmask_cplx(j,i) = dconjg(accu) + mo_one_e_integrals_complex(j,i)
   fock_op_cshell_ref_bitmask_cplx(j,i) = dconjg(fock_op_cshell_ref_bitmask_cplx(i,j))
  enddo
 enddo

 ! virt ---> virt single excitations
 do i0 = 1,  n_occ_ab_virt(1)
  i=occ_virt(i0,1)
  do j0 = 1, n_occ_ab_virt(1)
   j = occ_virt(j0,1)
   call get_mo_two_e_integrals_coulomb_ii_complex(i,j,mo_num,array_coulomb,mo_integrals_map,mo_integrals_map_2)
   call get_mo_two_e_integrals_exch_ii_complex(i,j,mo_num,array_exchange,mo_integrals_map,mo_integrals_map_2)
   accu = (0.d0,0.d0)
   do k0 = 1, n_occ_ab(1)
    k = occ(k0,1)
    accu += 2.d0 * array_coulomb(k) - array_exchange(k)
   enddo
   fock_op_cshell_ref_bitmask_cplx(i,j) = accu+ mo_one_e_integrals_complex(i,j)
   fock_op_cshell_ref_bitmask_cplx(j,i) = dconjg(accu)+ mo_one_e_integrals_complex(j,i)
  enddo
 enddo

 ! docc ---> docc single excitations
 do i0 = 1,  n_occ_ab(1)
  i=occ(i0,1)
  do j0 = 1, n_occ_ab(1)
   j = occ(j0,1)
   call get_mo_two_e_integrals_coulomb_ii_complex(i,j,mo_num,array_coulomb,mo_integrals_map,mo_integrals_map_2)
   call get_mo_two_e_integrals_exch_ii_complex(i,j,mo_num,array_exchange,mo_integrals_map,mo_integrals_map_2)
   accu = (0.d0,0.d0)
   do k0 = 1, n_occ_ab(1)
    k = occ(k0,1)
    accu += 2.d0 * array_coulomb(k) - array_exchange(k)
   enddo
   fock_op_cshell_ref_bitmask_cplx(i,j) = accu+ mo_one_e_integrals_complex(i,j)
   fock_op_cshell_ref_bitmask_cplx(j,i) = dconjg(accu)+ mo_one_e_integrals_complex(j,i)
  enddo
 enddo
 deallocate(array_coulomb,array_exchange)

END_PROVIDER

subroutine get_single_excitation_from_fock_complex(det_1,det_2,h,p,spin,phase,hij)
 use bitmasks
 implicit none
 integer,intent(in) :: h,p,spin
 double precision, intent(in)  :: phase
 integer(bit_kind), intent(in) :: det_1(N_int,2), det_2(N_int,2)
 complex*16, intent(out) :: hij
 integer(bit_kind) :: differences(N_int,2)
 integer(bit_kind) :: hole(N_int,2)
 integer(bit_kind) :: partcl(N_int,2)
 integer :: occ_hole(N_int*bit_kind_size,2)
 integer :: occ_partcl(N_int*bit_kind_size,2)
 integer :: n_occ_ab_hole(2),n_occ_ab_partcl(2)
 integer :: i0,i
 complex*16 :: buffer_c(mo_num),buffer_x(mo_num)
 do i=1, mo_num
   buffer_c(i) = big_array_coulomb_integrals_complex(i,h,p)
   buffer_x(i) = big_array_exchange_integrals_complex(i,h,p)
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
 hij = fock_op_cshell_ref_bitmask_cplx(h,p)
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

!============================================!
!                                            !
!                    kpts                    !
!                                            !
!============================================!

BEGIN_PROVIDER [integer(bit_kind), ref_closed_shell_bitmask_kpts, (N_int,2,kpt_num)]
 implicit none
 integer :: i,k
 do k = 1, kpt_num
   do i = 1, N_int
     ref_closed_shell_bitmask_kpts(i,1,k) = iand(ref_closed_shell_bitmask(i,1),kpts_bitmask(i,k))
     ref_closed_shell_bitmask_kpts(i,2,k) = iand(ref_closed_shell_bitmask(i,2),kpts_bitmask(i,k))
   enddo
 enddo
END_PROVIDER

BEGIN_PROVIDER [complex*16, fock_op_cshell_ref_bitmask_kpts, (mo_num_per_kpt, mo_num_per_kpt,kpt_num) ]
 implicit none
 integer :: i0,j0,i,j,k0,k,kblock,kvirt
 integer :: i_i, i_j, i_k, kocc
 integer :: n_occ_ab(2,kpt_num)
 integer :: occ(N_int*bit_kind_size,2,kpt_num)
 integer :: n_occ_ab_virt(2)
 integer :: occ_virt(N_int*bit_kind_size,2)
 integer(bit_kind) :: key_test(N_int)
 integer(bit_kind) :: key_virt(N_int,2)
 complex*16 :: accu
 complex*16, allocatable :: array_coulomb(:),array_exchange(:)

 do kblock = 1,kpt_num
   call bitstring_to_list_ab(ref_closed_shell_bitmask_kpts(1,1,kblock), &
                            occ(1,1,kblock), n_occ_ab(1,kblock), N_int)
 enddo
 allocate (array_coulomb(mo_num_per_kpt),array_exchange(mo_num_per_kpt))
 do kblock = 1,kpt_num
   ! get virt orbs for this kpt
   do i = 1, N_int
     key_virt(i,1) = iand(full_ijkl_bitmask(i),kpts_bitmask(i,kblock))
     key_virt(i,2) = iand(full_ijkl_bitmask(i),kpts_bitmask(i,kblock))
     key_virt(i,1) = xor(key_virt(i,1),ref_closed_shell_bitmask_kpts(i,1,kblock))
     key_virt(i,2) = xor(key_virt(i,2),ref_closed_shell_bitmask_kpts(i,2,kblock))
   enddo
   call bitstring_to_list_ab(key_virt, occ_virt, n_occ_ab_virt, N_int)
   ! docc ---> virt single excitations
   do i0 = 1,  n_occ_ab(1,kblock)
    i=occ(i0,1,kblock)
    i_i = mod(i-1,mo_num_per_kpt)+1
    do j0 = 1, n_occ_ab_virt(1)
     j = occ_virt(j0,1)
     i_j = mod(j-1,mo_num_per_kpt)+1
     accu = (0.d0,0.d0)
     do kocc = 1,kpt_num
      ! <ia|ja>
      array_coulomb(1:mo_num_per_kpt) = big_array_coulomb_integrals_kpts(1:mo_num_per_kpt,kocc,i_i,i_j,kblock)
      ! <ia|aj>
      array_exchange(1:mo_num_per_kpt) = big_array_exchange_integrals_kpts(1:mo_num_per_kpt,kocc,i_i,i_j,kblock)
      do k0 = 1, n_occ_ab(1,kocc)
       k = occ(k0,1,kocc)
       i_k = mod(k-1,mo_num_per_kpt)+1
       accu += 2.d0 * array_coulomb(i_k) - array_exchange(i_k)
      enddo
     enddo
     fock_op_cshell_ref_bitmask_kpts(i_i,i_j,kblock) = accu + mo_one_e_integrals_kpts(i_i,i_j,kblock)
     !fock_op_cshell_ref_bitmask_cplx(j,i) = dconjg(accu) + mo_one_e_integrals_complex(j,i)
     fock_op_cshell_ref_bitmask_kpts(i_j,i_i,kblock) = dconjg(fock_op_cshell_ref_bitmask_kpts(i_i,i_j,kblock))
    enddo
   enddo
  
   ! virt ---> virt single excitations
   do i0 = 1,  n_occ_ab_virt(1)
    i=occ_virt(i0,1)
    i_i = mod(i-1,mo_num_per_kpt)+1
    do j0 = 1, n_occ_ab_virt(1)
     j = occ_virt(j0,1)
     i_j = mod(j-1,mo_num_per_kpt)+1
     accu = (0.d0,0.d0)
     do kocc = 1,kpt_num
      array_coulomb(1:mo_num_per_kpt) = big_array_coulomb_integrals_kpts(1:mo_num_per_kpt,kocc,i_i,i_j,kblock)
      array_exchange(1:mo_num_per_kpt) = big_array_exchange_integrals_kpts(1:mo_num_per_kpt,kocc,i_i,i_j,kblock)
      do k0 = 1, n_occ_ab(1,kocc)
       k = occ(k0,1,kocc)
       i_k = mod(k-1,mo_num_per_kpt)+1
       accu += 2.d0 * array_coulomb(i_k) - array_exchange(i_k)
      enddo
     enddo
     fock_op_cshell_ref_bitmask_kpts(i_i,i_j,kblock) = accu + mo_one_e_integrals_kpts(i_i,i_j,kblock)
     fock_op_cshell_ref_bitmask_kpts(i_j,i_i,kblock) = dconjg(fock_op_cshell_ref_bitmask_kpts(i_i,i_j,kblock))
    enddo
   enddo
  
   ! docc ---> docc single excitations
   do i0 = 1,  n_occ_ab(1,kblock)
    i=occ(i0,1,kblock)
    i_i = mod(i-1,mo_num_per_kpt)+1
    do j0 = 1, n_occ_ab(1,kblock)
     j = occ(j0,1,kblock)
     i_j = mod(j-1,mo_num_per_kpt)+1
     accu = (0.d0,0.d0)
     do kocc = 1,kpt_num
      array_coulomb(1:mo_num_per_kpt) = big_array_coulomb_integrals_kpts(1:mo_num_per_kpt,kocc,i_i,i_j,kblock)
      array_exchange(1:mo_num_per_kpt) = big_array_exchange_integrals_kpts(1:mo_num_per_kpt,kocc,i_i,i_j,kblock)
      do k0 = 1, n_occ_ab(1,kocc)
       k = occ(k0,1,kocc)
       i_k = mod(k-1,mo_num_per_kpt)+1
       accu += 2.d0 * array_coulomb(i_k) - array_exchange(i_k)
      enddo
     enddo
     fock_op_cshell_ref_bitmask_kpts(i_i,i_j,kblock) = accu + mo_one_e_integrals_kpts(i_i,i_j,kblock)
     fock_op_cshell_ref_bitmask_kpts(i_j,i_i,kblock) = dconjg(fock_op_cshell_ref_bitmask_kpts(i_i,i_j,kblock))
    enddo
   enddo
  enddo
 deallocate(array_coulomb,array_exchange)

END_PROVIDER

subroutine get_single_excitation_from_fock_kpts(det_1,det_2,ih,ip,spin,phase,hij)
  use bitmasks
  !called by i_h_j{,_s2,_single_spin}_complex
  ! ih, ip are indices in total mo list (not per kpt)
  implicit none
  integer,intent(in) :: ih,ip,spin
  double precision, intent(in)  :: phase
  integer(bit_kind), intent(in) :: det_1(N_int,2), det_2(N_int,2)
  complex*16, intent(out) :: hij
  integer(bit_kind) :: differences(N_int,2)
  integer(bit_kind) :: hole(N_int,2)
  integer(bit_kind) :: partcl(N_int,2)
  integer :: occ_hole(N_int*bit_kind_size,2)
  integer :: occ_partcl(N_int*bit_kind_size,2)
  integer :: n_occ_ab_hole(2),n_occ_ab_partcl(2)
  integer :: i0,i,h,p
  integer :: ki,khp,kh
  complex*16 :: buffer_c(mo_num_per_kpt),buffer_x(mo_num_per_kpt)

  call get_kpt_idx_mo(ip,khp,p)
  call get_kpt_idx_mo(ih,kh,h)
  ASSERT (kh==khp)
  !todo: omp kpts
  hij = fock_op_cshell_ref_bitmask_kpts(h,p,khp)
  do ki=1,kpt_num
    do i=1, mo_num_per_kpt
      !<hi|pi>
      buffer_c(i) = big_array_coulomb_integrals_kpts(i,ki,h,p,khp)
      !<hi|ip>
      buffer_x(i) = big_array_exchange_integrals_kpts(i,ki,h,p,khp)
    enddo
    do i = 1, N_int
      !holes in ref, not in det1
      !part in det1, not in ref
      differences(i,1) = iand(xor(det_1(i,1),ref_closed_shell_bitmask(i,1)),kpts_bitmask(i,ki))
      differences(i,2) = iand(xor(det_1(i,2),ref_closed_shell_bitmask(i,2)),kpts_bitmask(i,ki))
      !differences(i,1) = xor(det_1(i,1),ref_closed_shell_bitmask_kpts(i,1,ki))
      !differences(i,2) = xor(det_1(i,2),ref_closed_shell_bitmask_kpts(i,2,ki))
      hole(i,1) = iand(differences(i,1),ref_closed_shell_bitmask_kpts(i,1,ki))
      hole(i,2) = iand(differences(i,2),ref_closed_shell_bitmask_kpts(i,2,ki))
      partcl(i,1) = iand(differences(i,1),det_1(i,1))
      partcl(i,2) = iand(differences(i,2),det_1(i,2))
    enddo
    call bitstring_to_list_ab(hole, occ_hole, n_occ_ab_hole, N_int)
    call bitstring_to_list_ab(partcl, occ_partcl, n_occ_ab_partcl, N_int)
    ! holes :: direct terms
    do i0 = 1, n_occ_ab_hole(1)
      i = occ_hole(i0,1) - (ki-1)*mo_num_per_kpt
      hij -= buffer_c(i)
    enddo
    do i0 = 1, n_occ_ab_hole(2)
      i = occ_hole(i0,2) - (ki-1)*mo_num_per_kpt
      hij -= buffer_c(i)
    enddo
   
    ! holes :: exchange terms
    do i0 = 1, n_occ_ab_hole(spin)
      i = occ_hole(i0,spin) - (ki-1)*mo_num_per_kpt
      hij += buffer_x(i)
    enddo
   
    ! particles :: direct terms
    do i0 = 1, n_occ_ab_partcl(1)
      i = occ_partcl(i0,1) - (ki-1)*mo_num_per_kpt
      hij += buffer_c(i)
    enddo
    do i0 = 1, n_occ_ab_partcl(2)
      i = occ_partcl(i0,2) - (ki-1)*mo_num_per_kpt
      hij += buffer_c(i)
    enddo
   
    ! particles :: exchange terms
    do i0 = 1, n_occ_ab_partcl(spin)
      i = occ_partcl(i0,spin) - (ki-1)*mo_num_per_kpt
      hij -= buffer_x(i)
    enddo
  enddo
  hij = hij * phase

end

