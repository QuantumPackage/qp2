
subroutine H_apply_cis_kpts_monoexc(key_in, hole_1,particl_1,fock_diag_tmp,i_generator,iproc_in  )
  use omp_lib
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Generate all single excitations of key_in using the bit masks of holes and
  ! particles.
  ! Assume N_int is already provided.
  END_DOC
  integer,parameter              :: size_max = 8192
  
  integer          ,intent(in)   :: i_generator
  integer(bit_kind),intent(in)   :: key_in(N_int,2)
  integer(bit_kind),intent(in)   :: hole_1(N_int,2), particl_1(N_int,2)
  integer, intent(in)            :: iproc_in
  double precision, intent(in)   :: fock_diag_tmp(2,mo_num+1)
  integer(bit_kind),allocatable  :: keys_out(:,:,:)
  integer(bit_kind),allocatable  :: hole_save(:,:)
  integer(bit_kind),allocatable  :: key(:,:),hole(:,:), particle(:,:)
  integer(bit_kind),allocatable  :: hole_tmp(:,:), particle_tmp(:,:)
  integer(bit_kind),allocatable  :: hole_2(:,:), particl_2(:,:)
  integer                        :: ii,i,jj,j,k,ispin,l
  integer,allocatable            :: occ_particle(:,:), occ_hole(:,:)
  integer,allocatable            :: occ_particle_tmp(:,:), occ_hole_tmp(:,:)
  integer,allocatable            :: ib_jb_pairs(:,:)
  integer                        :: kk,pp,other_spin,key_idx
  integer                        :: N_elec_in_key_hole_1(2),N_elec_in_key_part_1(2)
  integer                        :: N_elec_in_key_hole_2(2),N_elec_in_key_part_2(2)
  logical                        :: is_a_two_holes_two_particles
  integer(bit_kind), allocatable :: key_union_hole_part(:)

  integer, allocatable           :: ia_ja_pairs(:,:,:)
  logical, allocatable           :: array_pairs(:,:)
  double precision               :: diag_H_mat_elem
  integer                        :: iproc

  integer(bit_kind)              :: key_mask(N_int, 2)

  logical :: check_double_excitation
  logical :: is_a_2h1p
  logical :: is_a_2h
  logical :: is_a_1h1p
  logical :: is_a_1h2p
  logical :: is_a_1h
  logical :: is_a_1p
  logical :: is_a_2p
  logical :: yes_no

  do k=1,N_int
    key_mask(k,1) = 0_bit_kind
    key_mask(k,2) = 0_bit_kind
  enddo

  iproc = iproc_in

  check_double_excitation = .True.
  


  

  
!$ iproc = omp_get_thread_num()
  allocate (keys_out(N_int,2,size_max), hole_save(N_int,2),          &
      key(N_int,2),hole(N_int,2), particle(N_int,2), hole_tmp(N_int,2),&
      particle_tmp(N_int,2), occ_particle(N_int*bit_kind_size,2),    &
      occ_hole(N_int*bit_kind_size,2), occ_particle_tmp(N_int*bit_kind_size,2),&
      occ_hole_tmp(N_int*bit_kind_size,2),key_union_hole_part(N_int))
  
  !!!! First couple hole particle
  do j = 1, N_int
    hole(j,1) = iand(hole_1(j,1),key_in(j,1))
    hole(j,2) = iand(hole_1(j,2),key_in(j,2))
    particle(j,1) = iand(xor(particl_1(j,1),key_in(j,1)),particl_1(j,1))
    particle(j,2) = iand(xor(particl_1(j,2),key_in(j,2)),particl_1(j,2))
  enddo

  call bitstring_to_list_ab(particle,occ_particle,N_elec_in_key_part_1,N_int)
  call bitstring_to_list_ab(hole,occ_hole,N_elec_in_key_hole_1,N_int)
  allocate (ia_ja_pairs(2,0:(elec_alpha_num)*mo_num,2))

  do ispin=1,2
    i=0
    do ii=N_elec_in_key_hole_1(ispin),1,-1             ! hole
      i_a = occ_hole(ii,ispin)
      do jj=1,N_elec_in_key_part_1(ispin)                            !particule
        j_a = occ_particle(jj,ispin)
        i += 1
        ia_ja_pairs(1,i,ispin) = i_a
        ia_ja_pairs(2,i,ispin) = j_a
      enddo
    enddo
    ia_ja_pairs(1,0,ispin) = i
  enddo

  key_idx = 0

  integer                        :: i_a,j_a,i_b,j_b,k_a,l_a,k_b,l_b
  integer(bit_kind)              :: test(N_int,2)
  double precision               :: accu
  accu = 0.d0
  do ispin=1,2
    other_spin = iand(ispin,1)+1
    
    do ii=1,ia_ja_pairs(1,0,ispin)
      i_a = ia_ja_pairs(1,ii,ispin)
      j_a = ia_ja_pairs(2,ii,ispin)
      hole = key_in
      k = shiftr(i_a-1,bit_kind_shift)+1
      j = i_a-shiftl(k-1,bit_kind_shift)-1
  
      hole(k,ispin) = ibclr(hole(k,ispin),j)
      k_a = shiftr(j_a-1,bit_kind_shift)+1
      l_a = j_a-shiftl(k_a-1,bit_kind_shift)-1
  
      hole(k_a,ispin) = ibset(hole(k_a,ispin),l_a)
      
      
      
      
      
      
      
      
      
      
      
      
      
      key_idx += 1
      do k=1,N_int
        keys_out(k,1,key_idx) = hole(k,1)
        keys_out(k,2,key_idx) = hole(k,2)
      enddo
      if (key_idx == size_max) then
        call fill_H_apply_buffer_no_selection(key_idx,keys_out,N_int,iproc)
        key_idx = 0
      endif
    enddo  ! ii
    
  enddo   ! ispin
  call fill_H_apply_buffer_no_selection(key_idx,keys_out,N_int,iproc)
  
  deallocate (ia_ja_pairs, &
      keys_out, hole_save,          &
      key,hole, particle, hole_tmp,&
      particle_tmp, occ_particle,    &
      occ_hole, occ_particle_tmp,&
      occ_hole_tmp,key_union_hole_part)
  
  

end

subroutine H_apply_cis_kpts()
  implicit none
  use omp_lib
  use bitmasks
  BEGIN_DOC
  ! Calls H_apply on the |HF| determinant and selects all connected single and double
  ! excitations (of the same symmetry). Auto-generated by the ``generate_h_apply`` script.
  END_DOC

  

  integer                        :: i_generator
  double precision               :: wall_0, wall_1
  integer(bit_kind), allocatable :: mask(:,:,:)
  integer(bit_kind), allocatable :: mask_kpts(:,:,:,:)
  integer                        :: kk
  integer                        :: ispin, k
  integer                        :: iproc
  double precision, allocatable  :: fock_diag_tmp(:,:)

  
  if (is_complex) then
    PROVIDE H_apply_buffer_allocated mo_two_e_integrals_in_map psi_det_generators psi_coef_generators_complex
  else
  PROVIDE H_apply_buffer_allocated mo_two_e_integrals_in_map psi_det_generators psi_coef_generators
  endif

  call wall_time(wall_0)

  iproc = 0
  !allocate( mask(N_int,2,6), fock_diag_tmp(2,mo_num+1) )
  allocate( mask_kpts(N_int,2,6,kpt_num), fock_diag_tmp(2,mo_num+1) )
  do i_generator=1,N_det_generators

    ! Compute diagonal of the Fock matrix
    !call build_fock_tmp(fock_diag_tmp,psi_det_generators(1,1,i_generator),N_int)
    fock_diag_tmp=0.d0

    ! Create bit masks for holes and particles
    do kk=1,kpt_num
    do ispin=1,2
      do k=1,N_int
        mask_kpts(k,ispin,s_hole,kk) =                                      &
            iand(generators_bitmask_kpts(k,ispin,s_hole,kk),  &
            psi_det_generators(k,ispin,i_generator) )
        mask_kpts(k,ispin,s_part,kk) =                                      &
            iand(generators_bitmask_kpts(k,ispin,s_part,kk),  &
            not(psi_det_generators(k,ispin,i_generator)) )
       ! mask_kpts(k,ispin,d_hole1,kk) =                                      &
       !     iand(generators_bitmask_kpts(k,ispin,d_hole1,kk),  &
       !     psi_det_generators(k,ispin,i_generator) )
       ! mask_kpts(k,ispin,d_part1,kk) =                                      &
       !     iand(generators_bitmask_kpts(k,ispin,d_part1,kk),  &
       !     not(psi_det_generators(k,ispin,i_generator)) )
       ! mask_kpts(k,ispin,d_hole2,kk) =                                      &
       !     iand(generators_bitmask_kpts(k,ispin,d_hole2,kk),  &
       !     psi_det_generators(k,ispin,i_generator) )
       ! mask_kpts(k,ispin,d_part2,kk) =                                      &
       !     iand(generators_bitmask_kpts(k,ispin,d_part2,kk),  &
       !     not(psi_det_generators(k,ispin,i_generator)) )
      enddo
    enddo
    enddo
    !if(.False.)then
    ! call H_apply_cis_kpts_diexc(psi_det_generators(1,1,i_generator),      &
    !     psi_det_generators(1,1,1),                                   &
    !     mask(1,1,d_hole1), mask(1,1,d_part1),                        &
    !     mask(1,1,d_hole2), mask(1,1,d_part2),                        &
    !     fock_diag_tmp, i_generator, iproc )
    !endif
    if(.True.)then
      do kk=1,kpt_num
     call H_apply_cis_kpts_monoexc(psi_det_generators(1,1,i_generator),    &
         mask_kpts(1,1,s_hole,kk), mask_kpts(1,1,s_part,kk ),                        &
         fock_diag_tmp, i_generator, iproc )
      enddo
    endif
    call wall_time(wall_1)
    
    if (wall_1 - wall_0 > 2.d0) then
        write(6,*)  &
       100.*float(i_generator)/float(N_det_generators), '% in ', wall_1-wall_0, 's'
        wall_0 = wall_1
    endif
  enddo

  !deallocate( mask, fock_diag_tmp )
  deallocate( mask_kpts, fock_diag_tmp )

  call copy_H_apply_buffer_to_wf
  if (s2_eig) then
    call make_s2_eigenfunction
  endif
  if (is_complex) then
    SOFT_TOUCH psi_det psi_coef_complex N_det
  else
    SOFT_TOUCH psi_det psi_coef N_det
  endif

  
  ! Sort H_jj to find the N_states lowest states
  integer                        :: i
  integer, allocatable           :: iorder(:)
  double precision, allocatable  :: H_jj(:)
  double precision, external     :: diag_h_mat_elem
  allocate(H_jj(N_det),iorder(N_det))
  !$OMP PARALLEL DEFAULT(NONE)                                       &
      !$OMP SHARED(psi_det,N_int,H_jj,iorder,N_det)                  &
      !$OMP PRIVATE(i)
  !$OMP DO
  do i = 1, N_det
    H_jj(i) = diag_h_mat_elem(psi_det(1,1,i),N_int)
    iorder(i) = i
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  call dsort(H_jj,iorder,N_det)
  if (is_complex) then
    do k=1,N_states
      psi_coef_complex(iorder(k),k) = (1.d0,0.d0)
    enddo
  else
    do k=1,N_states
      psi_coef(iorder(k),k) = 1.d0
    enddo
  endif
  deallocate(H_jj,iorder)
    

end




subroutine H_apply_cis_sym_kpts_monoexc(key_in, hole_1,particl_1,fock_diag_tmp,i_generator,iproc_in  )
  use omp_lib
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Generate all single excitations of key_in using the bit masks of holes and
  ! particles.
  ! Assume N_int is already provided.
  END_DOC
  integer,parameter              :: size_max = 8192
  
  integer          ,intent(in)   :: i_generator
  integer(bit_kind),intent(in)   :: key_in(N_int,2)
  integer(bit_kind),intent(in)   :: hole_1(N_int,2), particl_1(N_int,2)
  integer, intent(in)            :: iproc_in
  double precision, intent(in)   :: fock_diag_tmp(2,mo_num+1)
  integer(bit_kind),allocatable  :: keys_out(:,:,:)
  integer(bit_kind),allocatable  :: hole_save(:,:)
  integer(bit_kind),allocatable  :: key(:,:),hole(:,:), particle(:,:)
  integer(bit_kind),allocatable  :: hole_tmp(:,:), particle_tmp(:,:)
  integer(bit_kind),allocatable  :: hole_2(:,:), particl_2(:,:)
  integer                        :: ii,i,jj,j,k,ispin,l
  integer,allocatable            :: occ_particle(:,:), occ_hole(:,:)
  integer,allocatable            :: occ_particle_tmp(:,:), occ_hole_tmp(:,:)
  integer,allocatable            :: ib_jb_pairs(:,:)
  integer                        :: kk,pp,other_spin,key_idx
  integer                        :: N_elec_in_key_hole_1(2),N_elec_in_key_part_1(2)
  integer                        :: N_elec_in_key_hole_2(2),N_elec_in_key_part_2(2)
  logical                        :: is_a_two_holes_two_particles
  integer(bit_kind), allocatable :: key_union_hole_part(:)

  integer, allocatable           :: ia_ja_pairs(:,:,:)
  logical, allocatable           :: array_pairs(:,:)
  double precision               :: diag_H_mat_elem
  integer                        :: iproc

  integer(bit_kind)              :: key_mask(N_int, 2)

  logical :: check_double_excitation
  logical :: is_a_2h1p
  logical :: is_a_2h
  logical :: is_a_1h1p
  logical :: is_a_1h2p
  logical :: is_a_1h
  logical :: is_a_1p
  logical :: is_a_2p
  logical :: yes_no

  do k=1,N_int
    key_mask(k,1) = 0_bit_kind
    key_mask(k,2) = 0_bit_kind
  enddo

  iproc = iproc_in

  check_double_excitation = .True.
  


  

  
!$ iproc = omp_get_thread_num()
  allocate (keys_out(N_int,2,size_max), hole_save(N_int,2),          &
      key(N_int,2),hole(N_int,2), particle(N_int,2), hole_tmp(N_int,2),&
      particle_tmp(N_int,2), occ_particle(N_int*bit_kind_size,2),    &
      occ_hole(N_int*bit_kind_size,2), occ_particle_tmp(N_int*bit_kind_size,2),&
      occ_hole_tmp(N_int*bit_kind_size,2),key_union_hole_part(N_int))
  
  !!!! First couple hole particle
  do j = 1, N_int
    hole(j,1) = iand(hole_1(j,1),key_in(j,1))
    hole(j,2) = iand(hole_1(j,2),key_in(j,2))
    particle(j,1) = iand(xor(particl_1(j,1),key_in(j,1)),particl_1(j,1))
    particle(j,2) = iand(xor(particl_1(j,2),key_in(j,2)),particl_1(j,2))
  enddo

  call bitstring_to_list_ab(particle,occ_particle,N_elec_in_key_part_1,N_int)
  call bitstring_to_list_ab(hole,occ_hole,N_elec_in_key_hole_1,N_int)
  allocate (ia_ja_pairs(2,0:(elec_alpha_num)*mo_num,2))

  do ispin=1,2
    i=0
    do ii=N_elec_in_key_hole_1(ispin),1,-1             ! hole
      i_a = occ_hole(ii,ispin)
      do jj=1,N_elec_in_key_part_1(ispin)                            !particule
        j_a = occ_particle(jj,ispin)
        i += 1
        ia_ja_pairs(1,i,ispin) = i_a
        ia_ja_pairs(2,i,ispin) = j_a
      enddo
    enddo
    ia_ja_pairs(1,0,ispin) = i
  enddo

  key_idx = 0

  integer                        :: i_a,j_a,i_b,j_b,k_a,l_a,k_b,l_b
  integer(bit_kind)              :: test(N_int,2)
  double precision               :: accu
  accu = 0.d0
  do ispin=1,2
    other_spin = iand(ispin,1)+1
    
    do ii=1,ia_ja_pairs(1,0,ispin)
      i_a = ia_ja_pairs(1,ii,ispin)
      j_a = ia_ja_pairs(2,ii,ispin)
      hole = key_in
      k = shiftr(i_a-1,bit_kind_shift)+1
      j = i_a-shiftl(k-1,bit_kind_shift)-1
  
      hole(k,ispin) = ibclr(hole(k,ispin),j)
      k_a = shiftr(j_a-1,bit_kind_shift)+1
      l_a = j_a-shiftl(k_a-1,bit_kind_shift)-1
  
      hole(k_a,ispin) = ibset(hole(k_a,ispin),l_a)
      
      
      
      
      
      
      
      
      
      
      
      
      
     call connected_to_hf(hole,yes_no)
     if (.not.yes_no) cycle
    
      key_idx += 1
      do k=1,N_int
        keys_out(k,1,key_idx) = hole(k,1)
        keys_out(k,2,key_idx) = hole(k,2)
      enddo
      if (key_idx == size_max) then
        call fill_H_apply_buffer_no_selection(key_idx,keys_out,N_int,iproc)
        key_idx = 0
      endif
    enddo  ! ii
    
  enddo   ! ispin
  call fill_H_apply_buffer_no_selection(key_idx,keys_out,N_int,iproc)
  
  deallocate (ia_ja_pairs, &
      keys_out, hole_save,          &
      key,hole, particle, hole_tmp,&
      particle_tmp, occ_particle,    &
      occ_hole, occ_particle_tmp,&
      occ_hole_tmp,key_union_hole_part)
  
  

end

subroutine H_apply_cis_sym_kpts()
  implicit none
  use omp_lib
  use bitmasks
  BEGIN_DOC
  ! Calls H_apply on the |HF| determinant and selects all connected single and double
  ! excitations (of the same symmetry). Auto-generated by the ``generate_h_apply`` script.
  END_DOC

  

  integer                        :: i_generator
  double precision               :: wall_0, wall_1
  integer(bit_kind), allocatable :: mask(:,:,:)
  integer(bit_kind), allocatable :: mask_kpts(:,:,:,:)
  integer                        :: kk
  integer                        :: ispin, k
  integer                        :: iproc
  double precision, allocatable  :: fock_diag_tmp(:,:)

  
  if (is_complex) then
    PROVIDE H_apply_buffer_allocated mo_two_e_integrals_in_map psi_det_generators psi_coef_generators_complex
  else
  PROVIDE H_apply_buffer_allocated mo_two_e_integrals_in_map psi_det_generators psi_coef_generators
  endif

  call wall_time(wall_0)

  iproc = 0
  !allocate( mask(N_int,2,6), fock_diag_tmp(2,mo_num+1) )
  allocate( mask_kpts(N_int,2,6,kpt_num), fock_diag_tmp(2,mo_num+1) )
  do i_generator=1,N_det_generators

    ! Compute diagonal of the Fock matrix
    !call build_fock_tmp(fock_diag_tmp,psi_det_generators(1,1,i_generator),N_int)
    fock_diag_tmp=0.d0

    ! Create bit masks for holes and particles
    do kk=1,kpt_num
    do ispin=1,2
      do k=1,N_int
        mask_kpts(k,ispin,d_hole2,kk) =                                      &
            iand(generators_bitmask_kpts(k,ispin,d_hole2,kk),  &
            psi_det_generators(k,ispin,i_generator) )
        mask_kpts(k,ispin,d_part2,kk) =                                      &
            iand(generators_bitmask_kpts(k,ispin,d_part2,kk),  &
            not(psi_det_generators(k,ispin,i_generator)) )
       ! mask_kpts(k,ispin,d_hole1,kk) =                                      &
       !     iand(generators_bitmask_kpts(k,ispin,d_hole1,kk),  &
       !     psi_det_generators(k,ispin,i_generator) )
       ! mask_kpts(k,ispin,d_part1,kk) =                                      &
       !     iand(generators_bitmask_kpts(k,ispin,d_part1,kk),  &
       !     not(psi_det_generators(k,ispin,i_generator)) )
       ! mask_kpts(k,ispin,d_hole2,kk) =                                      &
       !     iand(generators_bitmask_kpts(k,ispin,d_hole2,kk),  &
       !     psi_det_generators(k,ispin,i_generator) )
       ! mask_kpts(k,ispin,d_part2,kk) =                                      &
       !     iand(generators_bitmask_kpts(k,ispin,d_part2,kk),  &
       !     not(psi_det_generators(k,ispin,i_generator)) )
      enddo
    enddo
    enddo
    !if(.False.)then
    ! call H_apply_cis_sym_kpts_diexc(psi_det_generators(1,1,i_generator),      &
    !     psi_det_generators(1,1,1),                                   &
    !     mask(1,1,d_hole1), mask(1,1,d_part1),                        &
    !     mask(1,1,d_hole2), mask(1,1,d_part2),                        &
    !     fock_diag_tmp, i_generator, iproc )
    !endif
    if(.True.)then
      do kk=1,kpt_num
     call H_apply_cis_sym_kpts_monoexc(psi_det_generators(1,1,i_generator),    &
         mask_kpts(1,1,s_hole,kk), mask_kpts(1,1,s_part,kk ),                        &
         fock_diag_tmp, i_generator, iproc )
      enddo
    endif
    call wall_time(wall_1)
    
    if (wall_1 - wall_0 > 2.d0) then
        write(6,*)  &
       100.*float(i_generator)/float(N_det_generators), '% in ', wall_1-wall_0, 's'
        wall_0 = wall_1
    endif
  enddo

  !deallocate( mask, fock_diag_tmp )
  deallocate( mask_kpts, fock_diag_tmp )

  call copy_H_apply_buffer_to_wf
  if (s2_eig) then
    call make_s2_eigenfunction
  endif
  if (is_complex) then
    SOFT_TOUCH psi_det psi_coef_complex N_det
  else
    SOFT_TOUCH psi_det psi_coef N_det
  endif

  
  ! Sort H_jj to find the N_states lowest states
  integer                        :: i
  integer, allocatable           :: iorder(:)
  double precision, allocatable  :: H_jj(:)
  double precision, external     :: diag_h_mat_elem
  allocate(H_jj(N_det),iorder(N_det))
  !$OMP PARALLEL DEFAULT(NONE)                                       &
      !$OMP SHARED(psi_det,N_int,H_jj,iorder,N_det)                  &
      !$OMP PRIVATE(i)
  !$OMP DO
  do i = 1, N_det
    H_jj(i) = diag_h_mat_elem(psi_det(1,1,i),N_int)
    iorder(i) = i
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  call dsort(H_jj,iorder,N_det)
  if (is_complex) then
    do k=1,N_states
      psi_coef_complex(iorder(k),k) = (1.d0,0.d0)
    enddo
  else
    do k=1,N_states
      psi_coef(iorder(k),k) = 1.d0
    enddo
  endif
  deallocate(H_jj,iorder)
    

end


