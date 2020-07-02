
subroutine H_apply_cisd_kpts_diexc(key_in, key_prev, hole_1,particl_1, hole_2, particl_2, fock_diag_tmp, i_generator, iproc_in  )
  implicit none
  integer(bit_kind), intent(in)         :: key_in(N_int, 2), hole_1(N_int, 2), hole_2(N_int, 2)
  integer(bit_kind), intent(in)         :: particl_1(N_int, 2), particl_2(N_int, 2)
  integer(bit_kind)                     :: p1_mask(N_int, 2), p2_mask(N_int, 2), tmp
  integer,intent(in)                    :: i_generator,iproc_in
  integer                               :: status(N_int*bit_kind_size, 2)
  integer                               :: highest, p1,p2,sp,ni,i,mi,nt,ns,k
  double precision, intent(in)          :: fock_diag_tmp(2,mo_num+1)
  integer(bit_kind), intent(in)         :: key_prev(N_int, 2, *)
  PROVIDE N_int
  PROVIDE N_det

  

  highest = 0
  do k=1,N_int*bit_kind_size
    status(k,1) = 0
    status(k,2) = 0
  enddo
  do sp=1,2
    do ni=1,N_int
      do i=1,bit_kind_size
        if(iand(1_bit_kind,shiftr(key_in(ni, sp), (i-1))) == 0) then
          cycle
        end if
        mi = (ni-1)*bit_kind_size+i
        status(mi, sp) = int(iand(1_bit_kind,shiftr(hole_1(ni,sp),(i-1))),4)
        status(mi, sp) = status(mi, sp) + 2*int(iand(1_bit_kind,shiftr(hole_2(ni,sp),(i-1))),4)
        if(status(mi, sp) /= 0 .and. mi > highest) then
          highest = mi
        end if
      end do
    end do
  end do

  do sp=1,2
    do p1=1,highest
      if(status(p1, sp) == 0) then
        cycle
      end if
      do p2=1,highest
        if(status(p2, sp) == 0) then
          cycle
        end if
        if((status(p1, sp) == 1 .and. status(p2, sp) > 1) .or. &
            (status(p1, sp) == 2 .and. status(p2, sp) == 3) .or. &
            (status(p1, sp) == 3 .and. status(p2, sp) == 3 .and. p2 > p1)) then
          call H_apply_cisd_kpts_diexcP(key_in, sp, p1, particl_1, sp, p2, particl_2, fock_diag_tmp, i_generator, iproc_in  )
        end if
      end do
    end do
  end do
  do p1=1,highest
    if(status(p1, 1) == 0) then
      cycle
    end if
    do p2=1,highest
      if(status(p2, 2) == 0) then
        cycle
      end if
      if((status(p1, 1) == 3) .or. &
          (status(p1, 1) == 1 .and. status(p2, 2) >= 2) .or. &
          (status(p1, 1) == 2 .and. status(p2, 2) /= 2)) then

          call H_apply_cisd_kpts_diexcP(key_in, 1, p1, particl_1, 2, p2, particl_2, fock_diag_tmp, i_generator, iproc_in  )
      end if
    end do
  end do
end subroutine


subroutine H_apply_cisd_kpts_diexcP(key_in, fs1, fh1, particl_1, fs2, fh2, particl_2, fock_diag_tmp, i_generator, iproc_in  )
  implicit none
  integer(bit_kind), intent(in)         :: key_in(N_int, 2), particl_1(N_int, 2), particl_2(N_int, 2)
  double precision, intent(in)          :: fock_diag_tmp(2,mo_num+1)
  integer(bit_kind)                     :: p1_mask(N_int, 2), p2_mask(N_int, 2), key_mask(N_int, 2)
  integer,intent(in)                    :: fs1,fs2,i_generator,iproc_in, fh1,fh2
  integer(bit_kind)                     :: miniList(N_int, 2, N_det)
  integer                               :: n_minilist, n_alpha, n_beta, deg(2), i, ni, k
  
  integer(bit_kind), parameter :: one = 1_bit_kind

  do k=1,N_int
    p1_mask(k,1) = 0_bit_kind
    p1_mask(k,2) = 0_bit_kind
    p2_mask(k,1) = 0_bit_kind
    p2_mask(k,2) = 0_bit_kind
  enddo
  p1_mask(shiftr(fh1-1,bit_kind_shift) + 1, fs1) = shiftl(one,iand(fh1-1,bit_kind_size-1))
  p2_mask(shiftr(fh2-1,bit_kind_shift) + 1, fs2) = shiftl(one,iand(fh2-1,bit_kind_size-1))

  do k=1,N_int
    key_mask(k,1) = key_in(k,1)
    key_mask(k,2) = key_in(k,2)
  enddo

  key_mask(shiftr(fh1-1,bit_kind_shift) + 1, fs1) -= shiftl(one,iand(fh1-1,bit_kind_size-1))
  key_mask(shiftr(fh2-1,bit_kind_shift) + 1, fs2) -= shiftl(one,iand(fh2-1,bit_kind_size-1))


  call H_apply_cisd_kpts_diexcOrg(key_in, key_mask, p1_mask, particl_1, p2_mask, particl_2, fock_diag_tmp, i_generator, iproc_in  )
end subroutine


subroutine H_apply_cisd_kpts_diexcOrg(key_in,key_mask,hole_1,particl_1,hole_2, particl_2, fock_diag_tmp, i_generator, iproc_in  )
  use omp_lib
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Generate all double excitations of key_in using the bit masks of holes and
  ! particles.
  ! Assume N_int is already provided.
  END_DOC
  integer,parameter              :: size_max = 8192
  
  integer          ,intent(in)   :: i_generator
  integer(bit_kind),intent(in)   :: key_in(N_int,2), key_mask(N_int, 2)
  integer(bit_kind),allocatable  :: keys_out(:,:,:)
  integer(bit_kind), intent(in)  :: hole_1(N_int,2), particl_1(N_int,2)
  integer(bit_kind), intent(in)  :: hole_2(N_int,2), particl_2(N_int,2)
  integer, intent(in)            :: iproc_in
  double precision, intent(in)   :: fock_diag_tmp(2,mo_num+1)
  integer(bit_kind), allocatable :: hole_save(:,:)
  integer(bit_kind), allocatable :: key(:,:),hole(:,:), particle(:,:)
  integer(bit_kind), allocatable :: hole_tmp(:,:), particle_tmp(:,:)
  integer(bit_kind), allocatable :: key_union_hole_part(:)
  integer                        :: ii,i,jj,j,k,ispin,l
  integer, allocatable           :: occ_particle(:,:), occ_hole(:,:)
  integer, allocatable           :: occ_particle_tmp(:,:), occ_hole_tmp(:,:)
  integer                        :: kk,pp,other_spin,key_idx
  integer                        :: N_elec_in_key_hole_1(2),N_elec_in_key_part_1(2)
  integer                        :: N_elec_in_key_hole_2(2),N_elec_in_key_part_2(2)

  double precision               :: mo_two_e_integral
  logical                        :: is_a_two_holes_two_particles
  integer, allocatable           :: ia_ja_pairs(:,:,:)
  integer, allocatable           :: ib_jb_pairs(:,:)
  double precision               :: diag_H_mat_elem
  integer                        :: iproc
  integer                        :: jtest_vvvv

  logical :: check_double_excitation
  logical :: is_a_1h1p
  logical :: is_a_1h2p
  logical :: is_a_1h
  logical :: is_a_1p
  logical :: is_a_2p
  logical :: is_a_2h1p
  logical :: is_a_2h
  logical :: b_cycle
  logical :: yes_no
  check_double_excitation = .True.
  iproc = iproc_in


  

  
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
  allocate (ia_ja_pairs(2,0:(elec_alpha_num)*mo_num,2),          &
            ib_jb_pairs(2,0:(elec_alpha_num)*mo_num))

  do ispin=1,2
    i=0
    do ii=N_elec_in_key_hole_1(ispin),1,-1             ! hole
      i_a = occ_hole(ii,ispin)
      ASSERT (i_a > 0)
      ASSERT (i_a <= mo_num)

      do jj=1,N_elec_in_key_part_1(ispin)              !particle
        j_a = occ_particle(jj,ispin)
        ASSERT (j_a > 0)
        ASSERT (j_a <= mo_num)
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
  logical, allocatable           :: array_pairs(:,:)
  allocate(array_pairs(mo_num,mo_num))
  accu = 0.d0
  do ispin=1,2
    other_spin = iand(ispin,1)+1
    
    do ii=1,ia_ja_pairs(1,0,ispin)
      i_a = ia_ja_pairs(1,ii,ispin)
      ASSERT (i_a > 0)
      ASSERT (i_a <= mo_num)
      j_a = ia_ja_pairs(2,ii,ispin)
      ASSERT (j_a > 0)
      ASSERT (j_a <= mo_num)
      hole = key_in
      k = shiftr(i_a-1,bit_kind_shift)+1
      j = i_a-shiftl(k-1,bit_kind_shift)-1
      hole(k,ispin) = ibclr(hole(k,ispin),j)
      k_a = shiftr(j_a-1,bit_kind_shift)+1
      l_a = j_a-shiftl(k_a-1,bit_kind_shift)-1
      hole(k_a,ispin) = ibset(hole(k_a,ispin),l_a)

      !!!! Second couple hole particle
      do j = 1, N_int
        hole_tmp(j,1) = iand(hole_2(j,1),hole(j,1))
        hole_tmp(j,2) = iand(hole_2(j,2),hole(j,2))
        particle_tmp(j,1) = iand(xor(particl_2(j,1),hole(j,1)),particl_2(j,1))
        particle_tmp(j,2) = iand(xor(particl_2(j,2),hole(j,2)),particl_2(j,2))
      enddo

      call bitstring_to_list_ab(particle_tmp,occ_particle_tmp,N_elec_in_key_part_2,N_int)
      call bitstring_to_list_ab(hole_tmp,occ_hole_tmp,N_elec_in_key_hole_2,N_int)

      !   hole = a^(+)_j_a(ispin) a_i_a(ispin)|key_in> : single exc :: orb(i_a,ispin) --> orb(j_a,ispin)
      hole_save = hole

      ! Build array of the non-zero integrals of second excitation
      array_pairs = .True.

        if (ispin == 1) then
        integer                        :: jjj

        i=0
        do kk = 1,N_elec_in_key_hole_2(other_spin)
          i_b = occ_hole_tmp(kk,other_spin)
          ASSERT (i_b > 0)
          ASSERT (i_b <= mo_num)
          do jjj=1,N_elec_in_key_part_2(other_spin)     ! particle
            j_b = occ_particle_tmp(jjj,other_spin)
            ASSERT (j_b > 0)
            ASSERT (j_b <= mo_num)
            if (array_pairs(i_b,j_b)) then
              
              i+= 1
              ib_jb_pairs(1,i) = i_b
              ib_jb_pairs(2,i) = j_b
            endif
          enddo
        enddo
        ib_jb_pairs(1,0) = i

        do kk = 1,ib_jb_pairs(1,0)
          hole = hole_save
          i_b = ib_jb_pairs(1,kk)
          j_b = ib_jb_pairs(2,kk)
          k = shiftr(i_b-1,bit_kind_shift)+1
          j = i_b-shiftl(k-1,bit_kind_shift)-1
          hole(k,other_spin) = ibclr(hole(k,other_spin),j)
          key = hole
          k = shiftr(j_b-1,bit_kind_shift)+1
          l = j_b-shiftl(k-1,bit_kind_shift)-1
          key(k,other_spin) = ibset(key(k,other_spin),l)
          
          
          
          
          
          
          
          
          
          
          key_idx += 1
          do k=1,N_int
            keys_out(k,1,key_idx) = key(k,1)
            keys_out(k,2,key_idx) = key(k,2)
          enddo
          ASSERT (key_idx <= size_max)
          if (key_idx == size_max) then
            call fill_H_apply_buffer_no_selection(key_idx,keys_out,N_int,iproc)
            key_idx = 0
          endif
        enddo
      endif

      !   does all the single excitations of the same spin
      i=0
      do kk = 1,N_elec_in_key_hole_2(ispin)
        i_b = occ_hole_tmp(kk,ispin)
        if (i_b <= i_a.or.i_b == j_a) cycle
        ASSERT (i_b > 0)
        ASSERT (i_b <= mo_num)
        do jjj=1,N_elec_in_key_part_2(ispin)     ! particule
          j_b = occ_particle_tmp(jjj,ispin)
          ASSERT (j_b > 0)
          ASSERT (j_b <= mo_num)
          if (j_b <= j_a) cycle
          if (array_pairs(i_b,j_b)) then
            
            i+= 1
            ib_jb_pairs(1,i) = i_b
            ib_jb_pairs(2,i) = j_b
          endif
        enddo
      enddo
      ib_jb_pairs(1,0) = i

      do kk = 1,ib_jb_pairs(1,0)
        hole = hole_save
        i_b = ib_jb_pairs(1,kk)
        j_b = ib_jb_pairs(2,kk)
        k = shiftr(i_b-1,bit_kind_shift)+1
        j = i_b-shiftl(k-1,bit_kind_shift)-1
        hole(k,ispin) = ibclr(hole(k,ispin),j)
        key = hole
        k = shiftr(j_b-1,bit_kind_shift)+1
        l = j_b-shiftl(k-1,bit_kind_shift)-1
        key(k,ispin) = ibset(key(k,ispin),l)
        
        
        
        
        
        
        
        
        
        
        key_idx += 1
        do k=1,N_int
          keys_out(k,1,key_idx) = key(k,1)
          keys_out(k,2,key_idx) = key(k,2)
        enddo
        ASSERT (key_idx <= size_max)
        if (key_idx == size_max) then
          call fill_H_apply_buffer_no_selection(key_idx,keys_out,N_int,iproc)
          key_idx = 0
        endif
      enddo ! kk

    enddo  ! ii
    
  enddo   ! ispin
  call fill_h_apply_buffer_no_selection(key_idx,keys_out,N_int,iproc)
  
  deallocate (ia_ja_pairs, ib_jb_pairs,                              &
      keys_out, hole_save,                                           &
      key,hole, particle, hole_tmp,                                  &
      particle_tmp, occ_particle,                                    &
      occ_hole, occ_particle_tmp,                                    &
      occ_hole_tmp,array_pairs,key_union_hole_part)
  
  
end

subroutine H_apply_cisd_kpts_monoexc(key_in, hole_1,particl_1,fock_diag_tmp,i_generator,iproc_in  )
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

subroutine H_apply_cisd_kpts()
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
  integer                        :: ispin, k
  integer                        :: iproc
  double precision, allocatable  :: fock_diag_tmp(:,:)

  integer :: kk,kh1,kh2,kp1,kp2
  integer(bit_kind), allocatable :: mask_kpts(:,:,:,:)
  
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
    call build_fock_tmp(fock_diag_tmp,psi_det_generators(1,1,i_generator),N_int)

    ! Create bit masks for holes and particles
    do kk=1,kpt_num
    do ispin=1,2
      do k=1,N_int
        mask_kpts(k,ispin,s_hole,kk) =                                 &
            iand(generators_bitmask_kpts(k,ispin,s_hole,kk),  &
            psi_det_generators(k,ispin,i_generator) )
        mask_kpts(k,ispin,s_part,kk) =                                 &
            iand(generators_bitmask_kpts(k,ispin,s_part,kk),  &
            not(psi_det_generators(k,ispin,i_generator)) )
        mask_kpts(k,ispin,d_hole1,kk) =                                 &
            iand(generators_bitmask_kpts(k,ispin,d_hole1,kk),  &
            psi_det_generators(k,ispin,i_generator) )
        mask_kpts(k,ispin,d_part1,kk) =                                 &
            iand(generators_bitmask_kpts(k,ispin,d_part1,kk),  &
            not(psi_det_generators(k,ispin,i_generator)) )
        mask_kpts(k,ispin,d_hole2,kk) =                                 &
            iand(generators_bitmask_kpts(k,ispin,d_hole2,kk),  &
            psi_det_generators(k,ispin,i_generator) )
        mask_kpts(k,ispin,d_part2,kk) =                                 &
            iand(generators_bitmask_kpts(k,ispin,d_part2,kk),  &
            not(psi_det_generators(k,ispin,i_generator)) )
      enddo
    enddo
    enddo
    if(.True.)then
      do kh1=1,kpt_num
        do kh2=1,kh1
          do kp1=1,kpt_num
            kp2=kconserv(kh1,kh2,kp1)
            print*,'kh1h2p1p1',kh1,kh2,kp1,kp2
            print*,'size_before: ',h_apply_buffer(iproc)%n_det
     call H_apply_cisd_kpts_diexc(psi_det_generators(1,1,i_generator),      &
         psi_det_generators(1,1,1),                                   &
         mask_kpts(1,1,d_hole1,kh1), mask_kpts(1,1,d_part1,kp1),                        &
         mask_kpts(1,1,d_hole2,kh2), mask_kpts(1,1,d_part2,kp2),                        &
         fock_diag_tmp, i_generator, iproc )
            print*,'size_after: ',h_apply_buffer(iproc)%n_det
          enddo
        enddo
      enddo
    endif
    if(.True.)then
      do kk=1,kpt_num
     call H_apply_cisd_kpts_monoexc(psi_det_generators(1,1,i_generator),    &
         mask_kpts(1,1,s_hole,kk), mask_kpts(1,1,s_part,kk),                        &
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

