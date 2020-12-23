use bitmasks



subroutine alpha_callback(delta_ij_loc, i_generator, subset, csubset, iproc)
  use bitmasks
  implicit none
  integer, intent(in)            :: i_generator, subset, csubset
  double precision,intent(inout) :: delta_ij_loc(N_states,N_det,2)
  integer, intent(in)            :: iproc

  integer :: k,l

  integer(bit_kind)              :: hole_mask(N_int,2), particle_mask(N_int,2)


  call generate_singles_and_doubles(delta_ij_loc,i_generator,subset,csubset,iproc)
end subroutine


BEGIN_PROVIDER [ integer, psi_from_sorted_gen, (N_det) ]
  implicit none
  integer :: i,inpsisor

  psi_from_sorted_gen = 0

  do i=1,N_det
    psi_from_sorted_gen(psi_det_sorted_gen_order(i)) = i
    inpsisor = psi_det_sorted_gen_order(i)
    if(inpsisor <= 0) stop "idx_non_ref_from_sorted"
  end do
END_PROVIDER


subroutine generate_singles_and_doubles(delta_ij_loc, i_generator, subset, csubset, iproc)
  use bitmasks
  implicit none
  BEGIN_DOC
! TODO
  END_DOC

  double precision,intent(inout) :: delta_ij_loc(N_states,N_det,2)
  integer, intent(in)            :: i_generator, subset, csubset
  integer, intent(in)            :: iproc


  integer                         :: h1,h2,s1,s2,s3,i1,i2,ib,sp,k,i,j,nt,ii,n
  integer(bit_kind)               :: hole(N_int,2), particle(N_int,2), mask(N_int, 2), pmask(N_int, 2)
  integer(bit_kind)  :: mmask(N_int, 2)
  logical                         :: fullMatch, ok

  integer(bit_kind) :: mobMask(N_int, 2), negMask(N_int, 2)
  integer,allocatable               :: preinteresting(:), prefullinteresting(:), interesting(:), fullinteresting(:)
  integer(bit_kind), allocatable :: minilist(:, :, :), fullminilist(:, :, :)
  logical, allocatable           :: banned(:,:,:), bannedOrb(:,:)
  integer, allocatable           :: counted(:,:), countedOrb(:,:)
  integer ::                      countedGlob, siz, lsiz

  integer, allocatable   :: indexes_end(:,:), indexes(:,:)

  logical :: monoAdo, monoBdo
  integer :: maskInd

  integer(bit_kind), allocatable:: preinteresting_det(:,:,:)
  integer ,allocatable :: abuf(:), labuf(:)

  allocate(abuf(N_det*6), labuf(N_det))
  allocate(preinteresting_det(N_int,2,N_det))


  maskInd = -1

  monoAdo = .true.
  monoBdo = .true.


  ! Masks adapted for MRCC
  do k=1,N_int
    hole    (k,1) = iand(psi_det_generators(k,1,i_generator), ior(generators_bitmask(k,1,s_hole),generators_bitmask(k,1,s_part)  ) )
    hole    (k,2) = iand(psi_det_generators(k,2,i_generator), ior(generators_bitmask(k,2,s_hole),generators_bitmask(k,2,s_part)  ) )
    particle(k,1) = iand(not(psi_det_generators(k,1,i_generator)), ior(generators_bitmask(k,1,s_part),generators_bitmask(k,1,s_hole)) )
    particle(k,2) = iand(not(psi_det_generators(k,2,i_generator)), ior(generators_bitmask(k,2,s_part),generators_bitmask(k,2,s_hole)) )
  enddo

  integer                        :: N_holes(2), N_particles(2)
  integer                        :: hole_list(N_int*bit_kind_size,2)
  integer                        :: particle_list(N_int*bit_kind_size,2)

  call bitstring_to_list_ab(hole    , hole_list    , N_holes    , N_int)
  call bitstring_to_list_ab(particle, particle_list, N_particles, N_int)

  integer :: l_a, nmax
  integer, allocatable :: indices(:), exc_degree(:), iorder(:)
  allocate (indices(N_det),  &
            exc_degree(max(N_det_alpha_unique,N_det_beta_unique)))

  PROVIDE psi_det_sorted_gen_order
  !PROVIDE psi_bilinear_matrix_columns_loc psi_det_alpha_unique psi_det_beta_unique
  !PROVIDE psi_bilinear_matrix_rows psi_det_sorted_gen_order psi_bilinear_matrix_order
  !PROVIDE psi_bilinear_matrix_transp_rows_loc psi_bilinear_matrix_transp_columns
  !PROVIDE psi_bilinear_matrix_transp_order

  k=1
  do i=1,N_det_alpha_unique
    call get_excitation_degree_spin(psi_det_alpha_unique(1,i), &
      psi_det_generators(1,1,i_generator), exc_degree(i), N_int)
  enddo

  do j=1,N_det_beta_unique
    call get_excitation_degree_spin(psi_det_beta_unique(1,j), &
      psi_det_generators(1,2,i_generator), nt, N_int)
    if (nt > 2) cycle
    do l_a=psi_bilinear_matrix_columns_loc(j), psi_bilinear_matrix_columns_loc(j+1)-1
      i = psi_bilinear_matrix_rows(l_a)
      if (nt + exc_degree(i) <= 4) then
        indices(k) = psi_det_sorted_gen_order(psi_bilinear_matrix_order(l_a))
        k=k+1
      endif
    enddo
  enddo

  do i=1,N_det_beta_unique
    call get_excitation_degree_spin(psi_det_beta_unique(1,i), &
      psi_det_generators(1,2,i_generator), exc_degree(i), N_int)
  enddo

  do j=1,N_det_alpha_unique
    call get_excitation_degree_spin(psi_det_alpha_unique(1,j), &
      psi_det_generators(1,1,i_generator), nt, N_int)
    if (nt > 1) cycle
    do l_a=psi_bilinear_matrix_transp_rows_loc(j), psi_bilinear_matrix_transp_rows_loc(j+1)-1
      i = psi_bilinear_matrix_transp_columns(l_a)
      if (exc_degree(i) < 3) cycle
      if (nt + exc_degree(i) <= 4) then
        indices(k) = psi_det_sorted_gen_order(                   &
                        psi_bilinear_matrix_order(           &
                          psi_bilinear_matrix_transp_order(l_a)))
        k=k+1
      endif
    enddo
  enddo
  nmax=k-1

  allocate(iorder(nmax))
  do i=1,nmax
    iorder(i) = i
  enddo
  call isort(indices,iorder,nmax)

  allocate(preinteresting(0:N_det_selectors), prefullinteresting(0:N_det), &
            interesting(0:N_det_selectors), fullinteresting(0:N_det))
  preinteresting(0) = 0
  prefullinteresting(0) = 0

  do i=1,N_int
    negMask(i,1) = not(psi_det_generators(i,1,i_generator))
    negMask(i,2) = not(psi_det_generators(i,2,i_generator))
  end do
  if(psi_det_generators(1,1,i_generator) /= psi_det_sorted_gen(1,1,i_generator)) stop "gen <> sorted"
  do k=1,nmax
    i = indices(k)
    mobMask(1,1) = iand(negMask(1,1), psi_det_sorted_gen(1,1,i))
    mobMask(1,2) = iand(negMask(1,2), psi_det_sorted_gen(1,2,i))
    nt = popcnt(mobMask(1, 1)) + popcnt(mobMask(1, 2))
    do j=2,N_int
      mobMask(j,1) = iand(negMask(j,1), psi_det_sorted_gen(j,1,i))
      mobMask(j,2) = iand(negMask(j,2), psi_det_sorted_gen(j,2,i))
      nt = nt + popcnt(mobMask(j, 1)) + popcnt(mobMask(j, 2))
    end do

    if(nt <= 4) then
      if(i <= N_det_selectors) then
        preinteresting(0) += 1
        preinteresting(preinteresting(0)) = i
        do j=1,N_int
          preinteresting_det(j,1,preinteresting(0)) = psi_det_sorted_gen(j,1,i)
          preinteresting_det(j,2,preinteresting(0)) = psi_det_sorted_gen(j,2,i)
        enddo
      else if(nt <= 2) then
        prefullinteresting(0) += 1
        prefullinteresting(prefullinteresting(0)) = i
      end if
    end if
  end do


  allocate(minilist(N_int, 2, N_det_selectors), fullminilist(N_int, 2, N_det))
  allocate(banned(mo_num, mo_num,2), bannedOrb(mo_num, 2))
  allocate(counted(mo_num, mo_num), countedOrb(mo_num, 2))
  allocate (indexes(0:mo_num, 0:mo_num))
  allocate (indexes_end(0:mo_num, 0:mo_num))
  integer :: nb_count
  do s1=1,2
    do i1=N_holes(s1),1,-1   ! Generate low excitations first
      h1 = hole_list(i1,s1)
      call apply_hole(psi_det_generators(1,1,i_generator), s1,h1, pmask, ok, N_int)

      negMask = not(pmask)

      interesting(0) = 0
      fullinteresting(0) = 0

      do ii=1,preinteresting(0)
        select case (N_int)
          case (1)
            mobMask(1,1) = iand(negMask(1,1), preinteresting_det(1,1,ii))
            mobMask(1,2) = iand(negMask(1,2), preinteresting_det(1,2,ii))
            nt = popcnt(mobMask(1, 1)) + popcnt(mobMask(1, 2))
          case (2)
            mobMask(1:2,1) = iand(negMask(1:2,1), preinteresting_det(1:2,1,ii))
            mobMask(1:2,2) = iand(negMask(1:2,2), preinteresting_det(1:2,2,ii))
            nt = popcnt(mobMask(1, 1)) + popcnt(mobMask(1, 2)) + &
                 popcnt(mobMask(2, 1)) + popcnt(mobMask(2, 2))
          case (3)
            mobMask(1:3,1) = iand(negMask(1:3,1), preinteresting_det(1:3,1,ii))
            mobMask(1:3,2) = iand(negMask(1:3,2), preinteresting_det(1:3,2,ii))
            nt = 0
            do j=3,1,-1
              if (mobMask(j,1) /= 0_bit_kind) then
                nt = nt+ popcnt(mobMask(j, 1))
                if (nt > 4) exit
              endif
              if (mobMask(j,2) /= 0_bit_kind) then
                nt = nt+ popcnt(mobMask(j, 2))
                if (nt > 4) exit
              endif
            end do
          case (4)
            mobMask(1:4,1) = iand(negMask(1:4,1), preinteresting_det(1:4,1,ii))
            mobMask(1:4,2) = iand(negMask(1:4,2), preinteresting_det(1:4,2,ii))
            nt = 0
            do j=4,1,-1
              if (mobMask(j,1) /= 0_bit_kind) then
                nt = nt+ popcnt(mobMask(j, 1))
                if (nt > 4) exit
              endif
              if (mobMask(j,2) /= 0_bit_kind) then
                nt = nt+ popcnt(mobMask(j, 2))
                if (nt > 4) exit
              endif
            end do
          case default
            mobMask(1:N_int,1) = iand(negMask(1:N_int,1), preinteresting_det(1:N_int,1,ii))
            mobMask(1:N_int,2) = iand(negMask(1:N_int,2), preinteresting_det(1:N_int,2,ii))
            nt = 0
            do j=N_int,1,-1
              if (mobMask(j,1) /= 0_bit_kind) then
                nt = nt+ popcnt(mobMask(j, 1))
                if (nt > 4) exit
              endif
              if (mobMask(j,2) /= 0_bit_kind) then
                nt = nt+ popcnt(mobMask(j, 2))
                if (nt > 4) exit
              endif
            end do
        end select

        if(nt <= 4) then
          i = preinteresting(ii)
          interesting(0) += 1
          interesting(interesting(0)) = i
          minilist(1,1,interesting(0)) = preinteresting_det(1,1,ii)
          minilist(1,2,interesting(0)) = preinteresting_det(1,2,ii)
          do j=2,N_int
            minilist(j,1,interesting(0)) = preinteresting_det(j,1,ii)
            minilist(j,2,interesting(0)) = preinteresting_det(j,2,ii)
          enddo
          if(nt <= 2) then
            fullinteresting(0) += 1
            fullinteresting(fullinteresting(0)) = i
            fullminilist(1,1,fullinteresting(0)) = preinteresting_det(1,1,ii)
            fullminilist(1,2,fullinteresting(0)) = preinteresting_det(1,2,ii)
            do j=2,N_int
              fullminilist(j,1,fullinteresting(0)) = preinteresting_det(j,1,ii)
              fullminilist(j,2,fullinteresting(0)) = preinteresting_det(j,2,ii)
            enddo
          end if
        end if

      end do

      do ii=1,prefullinteresting(0)
        i = prefullinteresting(ii)
        nt = 0
        mobMask(1,1) = iand(negMask(1,1), psi_det_sorted_gen(1,1,i))
        mobMask(1,2) = iand(negMask(1,2), psi_det_sorted_gen(1,2,i))
        nt = popcnt(mobMask(1, 1)) + popcnt(mobMask(1, 2))
        if (nt > 2) cycle
        do j=N_int,2,-1
          mobMask(j,1) = iand(negMask(j,1), psi_det_sorted_gen(j,1,i))
          mobMask(j,2) = iand(negMask(j,2), psi_det_sorted_gen(j,2,i))
          nt = nt+ popcnt(mobMask(j, 1)) + popcnt(mobMask(j, 2))
          if (nt > 2) exit
        end do

        if(nt <= 2) then
          fullinteresting(0) += 1
          fullinteresting(fullinteresting(0)) = i
          fullminilist(1,1,fullinteresting(0)) = psi_det_sorted_gen(1,1,i)
          fullminilist(1,2,fullinteresting(0)) = psi_det_sorted_gen(1,2,i)
          do j=2,N_int
            fullminilist(j,1,fullinteresting(0)) = psi_det_sorted_gen(j,1,i)
            fullminilist(j,2,fullinteresting(0)) = psi_det_sorted_gen(j,2,i)
          enddo
        end if
      end do



      do s2=s1,2
        sp = s1

        if(s1 /= s2) sp = 3

        ib = 1
        if(s1 == s2) ib = i1+1
        monoAdo = .true.
        do i2=N_holes(s2),ib,-1   ! Generate low excitations first

          h2 = hole_list(i2,s2)
          call apply_hole(pmask, s2,h2, mask, ok, N_int)
          banned = .false.
          do j=1,mo_num
            bannedOrb(j, 1) = .true.
            bannedOrb(j, 2) = .true.
          enddo
          do s3=1,2
            do i=1,N_particles(s3)
              bannedOrb(particle_list(i,s3), s3) = .false.
            enddo
          enddo
          if(s1 /= s2) then
            if(monoBdo) then
              bannedOrb(h1,s1) = .false.
            end if
            if(monoAdo) then
              bannedOrb(h2,s2) = .false.
              monoAdo = .false.
            end if
          end if

          maskInd += 1
          if(mod(maskInd, csubset) == (subset-1)) then

            call spot_isinwf(mask, fullminilist, i_generator, fullinteresting(0), banned, fullMatch, fullinteresting)
            if(fullMatch) cycle

            call count_pq(mask, sp, minilist, i_generator, interesting(0), bannedOrb, banned, countedGlob, countedOrb, counted, interesting)
            call create_indexes(countedGlob, countedOrb, counted, indexes, siz)
            indexes_end = indexes


            if(siz > size(abuf)) stop "buffer too small in alpha_factory"
            call splash_pq(mask, sp, minilist, i_generator, interesting(0), bannedOrb, banned, indexes_end, abuf, interesting)
            call alpha_callback_mask(delta_ij_loc, i_generator, sp, mask, bannedOrb, banned, indexes, indexes_end, abuf, siz, iproc)

          end if
        enddo
        if(s1 /= s2) monoBdo = .false.
      enddo
    enddo
  enddo
end subroutine


subroutine alpha_callback_mask(delta_ij_loc, i_gen, sp, mask, bannedOrb, banned, indexes, indexes_end, rabuf, siz, iproc)
  use bitmasks
  implicit none

  double precision,intent(inout) :: delta_ij_loc(N_states,N_det,2)
  integer, intent(in) :: sp, indexes(0:mo_num, 0:mo_num), siz, iproc, i_gen
  integer, intent(in) :: indexes_end(0:mo_num, 0:mo_num), rabuf(*)
  logical, intent(in) :: bannedOrb(mo_num,2), banned(mo_num, mo_num)
  integer(bit_kind), intent(in) :: mask(N_int, 2)
  integer(bit_kind) :: alpha(N_int, 2)
  integer, allocatable :: labuf(:), abuf(:), iorder(:)
  logical :: ok
  integer :: i,j,k,s,st1,st2,st3,st4,t2
  integer :: lindex(mo_num,2), lindex_end(mo_num, 2)
  integer :: s1, s2, stamo
  integer(bit_kind), allocatable :: det_minilist(:,:,:)


  lindex = 0
  lindex_end = 0
  allocate(abuf(siz), labuf(N_det), iorder(siz), det_minilist(N_int, 2, N_det))

  do i=1,siz
    abuf(i) = psi_from_sorted_gen(rabuf(i))
  end do


  st1 = indexes_end(0,0)-1 !!
  if(st1 > 0) then
    labuf(:st1) = abuf(:st1)
    do i=1,st1
      det_minilist(:,:,i) = psi_det(:,:,labuf(i))
    end do
  end if
  st1 += 1

  if(sp == 3) then
    s1 = 1
    s2 = 2
    lindex(:, 1) = indexes(1:,0)
    lindex_end(:,1) = indexes_end(1:,0)-1
    lindex(:, 2) = indexes(0, 1:)
    lindex_end(:, 2) = indexes_end(0, 1:)-1
  else if(sp == 2) then
    s1 = 2
    s2 = 2
    lindex(:, 2) = indexes(0, 1:)
    lindex_end(:, 2) = indexes_end(0, 1:)-1
  else if(sp == 1) then
    s1 = 1
    s2 = 1
    lindex(:, 1) = indexes(1:, 0)
    lindex_end(:,1) = indexes_end(1:, 0)-1
  end if

  do i=1,mo_num
  do j=1,2
    if(lindex(i,j) > 0 .and. lindex_end(i,j) > lindex(i,j)) then
      call isort(abuf(lindex(i,j)), iorder, lindex_end(i,j)-lindex(i,j)+1)
    end if
  end do
  end do


  do i=1,mo_num
    if(bannedOrb(i,s1)) cycle
    if(lindex(i,s1) /= 0) then
      st2 = st1 + 1 + lindex_end(i,s1)-lindex(i,s1)
      labuf(st1:st2-1) = abuf(lindex(i,s1):lindex_end(i,s1))
      do j=st1,st2-1
        det_minilist(:,:,j) = psi_det(:,:,labuf(j))
      end do
    else
      st2 = st1
    end if

    if(sp == 3) then
      stamo = 1
    else
      stamo = i+1
    end if

    do j=stamo,mo_num
      if(bannedOrb(j,s2) .or. banned(i,j)) cycle
      if(lindex(j,s2) /= 0) then
        k = lindex(j,s2)
        st3 = st2
        t2 = st1
        do while(k <= lindex_end(j,s2))
          if(t2 >= st2) then
            labuf(st3) = abuf(k)
            det_minilist(:,:,st3) = psi_det(:,:,abuf(k))
            st3 += 1
            k += 1
          else if(abuf(k) > labuf(t2)) then
            t2 += 1
          else if(abuf(k) < labuf(t2)) then
            labuf(st3) = abuf(k)
            det_minilist(:,:,st3) = psi_det(:,:,abuf(k))
            st3 += 1
            k += 1
          else
            k += 1
            t2 += 1
          end if
        end do
      else
        st3 = st2
      end if

      if(indexes(i,j) /= 0) then
        st4 = st3 + 1 + indexes_end(i,j)-indexes(i,j) -1!!
        labuf(st3:st4-1) = abuf(indexes(i,j):indexes_end(i,j)-1) !!
        do k=st3, st4-1
          det_minilist(:,:,k) = psi_det(:,:,labuf(k))
        end do
      else
        st4 = st3
      end if
      !APPLY PART
      if(st4 > 1) then
        call apply_particles(mask, s1, i, s2, j, alpha, ok, N_int)
        call dress_with_alpha_buffer(N_states, N_det, N_int, delta_ij_loc, i_gen, labuf, det_minilist, st4-1, alpha, iproc)
      end if
    end do
  end do
end subroutine


subroutine create_indexes(countedGlob, countedOrb, counted, indexes, siz)
  use bitmasks
  implicit none

  integer, intent(in) :: countedGlob, countedOrb(mo_num,2), counted(mo_num, mo_num)
  integer, intent(out) :: indexes(0:mo_num, 0:mo_num), siz
  integer :: tmp, i, j

  indexes(0, 0) = countedGlob
  indexes(0, 1:) = countedOrb(:, 2)
  indexes(1:, 0) = countedOrb(:, 1)
  indexes(1:, 1:) = counted(:,:)

  siz = 1

  do i=0, mo_num
  do j=0, mo_num
    if(indexes(i,j) == 0) cycle
    tmp = indexes(i,j)
    indexes(i,j) = siz
    siz += tmp
  end do
  end do

  siz -= 1
end subroutine


subroutine count_pq(mask, sp, det, i_gen, N_sel, bannedOrb, banned, countedGlob, countedOrb, counted, interesting)
  use bitmasks
  implicit none

  integer, intent(in)            :: sp, i_gen, N_sel
  integer, intent(in)            :: interesting(0:N_sel)
  integer(bit_kind),intent(in)   :: mask(N_int, 2), det(N_int, 2, N_sel)
  logical, intent(inout)         :: bannedOrb(mo_num, 2), banned(mo_num, mo_num, 2)
  integer, intent(inout)         :: countedGlob, countedOrb(mo_num, 2), counted(mo_num, mo_num)


  integer                        :: i, s, ii, j, k, l, h(0:2,2), p(0:4,2), nt
  integer(bit_kind)              :: perMask(N_int, 2), mobMask(N_int, 2), negMask(N_int, 2)

  PROVIDE psi_selectors_coef_transp
  countedGlob = 0
  countedOrb = 0
  counted = 0

  do i=1,N_int
    negMask(i,1) = not(mask(i,1))
    negMask(i,2) = not(mask(i,2))
  end do

  do i=1, N_sel
    !if (interesting(i) < 0) then
    !  stop 'prefetch interesting(i)'
    !endif

    mobMask(1,1) = iand(negMask(1,1), det(1,1,i))
    mobMask(1,2) = iand(negMask(1,2), det(1,2,i))
    nt = popcnt(mobMask(1, 1)) + popcnt(mobMask(1, 2))


    if(nt > 4) cycle

    do j=2,N_int
      mobMask(j,1) = iand(negMask(j,1), det(j,1,i))
      mobMask(j,2) = iand(negMask(j,2), det(j,2,i))
      nt = nt + popcnt(mobMask(j, 1)) + popcnt(mobMask(j, 2))
    end do

    if(nt > 4) cycle

    if (interesting(i) == i_gen) then
      do s=1,2
      do j=1,mo_num
        if(bannedOrb(j, s)) then
          if(sp == 3 .and. s == 1) then
            banned(j, :, 1) = .true.
          else if(sp == 3 .and. s == 2) then
            banned(:, j, 1) = .true.
          else if(s == sp) then
            banned(j,:,1) = .true.
            banned(:,j,1) = .true.
          end if
        end if
      end do
      end do

      if(sp == 3) then
        do j=1,mo_num
          do k=1,mo_num
            banned(j,k,2) = banned(k,j,1)
          enddo
        enddo
      else
        do k=1,mo_num
        do l=k+1,mo_num
          banned(l,k,1) = banned(k,l,1)
        end do
        end do
      end if
    end if

    call bitstring_to_list_in_selection(mobMask(1,1), p(1,1), p(0,1), N_int)
    call bitstring_to_list_in_selection(mobMask(1,2), p(1,2), p(0,2), N_int)

    perMask(1,1) = iand(mask(1,1), not(det(1,1,i)))
    perMask(1,2) = iand(mask(1,2), not(det(1,2,i)))
    do j=2,N_int
      perMask(j,1) = iand(mask(j,1), not(det(j,1,i)))
      perMask(j,2) = iand(mask(j,2), not(det(j,2,i)))
    end do

    call bitstring_to_list_in_selection(perMask(1,1), h(1,1), h(0,1), N_int)
    call bitstring_to_list_in_selection(perMask(1,2), h(1,2), h(0,2), N_int)

    if (interesting(i) >= i_gen) then
        if(nt == 4) then
          call count_d2(counted, p, sp)
        else if(nt == 3) then
          call count_d1(countedOrb, p)
        else
          countedGlob += 1
        end if
    else
        if(nt == 4) call past_d2(banned, p, sp)
        if(nt == 3) call past_d1(bannedOrb, p)
        if(nt < 3) stop "past_d0 ?"
    end if
  end do

  do i=1,mo_num
    if(bannedOrb(i,1)) countedOrb(i,1) = 0
    if(bannedOrb(i,2)) countedOrb(i,2) = 0
    do j=1,mo_num
      if(banned(i,j,1)) counted(i,j) = 0
    end do
  end do

  if(sp /= 3) then
    countedOrb(:, mod(sp, 2)+1) = 0
  end if
end



subroutine splash_pq(mask, sp, det, i_gen, N_sel, bannedOrb, banned, indexes, abuf, interesting)
  use bitmasks
  implicit none

  integer, intent(in)            :: sp, i_gen, N_sel
  integer, intent(in)            :: interesting(0:N_sel)
  integer(bit_kind),intent(in)   :: mask(N_int, 2), det(N_int, 2, N_sel)
  logical, intent(inout)         :: bannedOrb(mo_num, 2), banned(mo_num, mo_num, 2)
  integer, intent(inout)         :: indexes(0:mo_num, 0:mo_num)
  integer, intent(inout)         :: abuf(*)
  integer                        :: i, ii, j, k, l, h(0:2,2), p(0:4,2), nt, s
  integer(bit_kind)              :: perMask(N_int, 2), mobMask(N_int, 2), negMask(N_int, 2)
  integer :: phasemask(2,N_int*bit_kind_size)

  PROVIDE psi_selectors_coef_transp
  do i=1,N_int
    negMask(i,1) = not(mask(i,1))
    negMask(i,2) = not(mask(i,2))
  end do

  do i=1, N_sel ! interesting(0)
    !i = interesting(ii)
    !if (interesting(i) < 0) then
    !  stop 'prefetch interesting(i)'
    !endif
    if(interesting(i) < i_gen) cycle


    mobMask(1,1) = iand(negMask(1,1), det(1,1,i))
    mobMask(1,2) = iand(negMask(1,2), det(1,2,i))
    nt = popcnt(mobMask(1, 1)) + popcnt(mobMask(1, 2))

    if(nt > 4) cycle

    do j=2,N_int
      mobMask(j,1) = iand(negMask(j,1), det(j,1,i))
      mobMask(j,2) = iand(negMask(j,2), det(j,2,i))
      nt = nt + popcnt(mobMask(j, 1)) + popcnt(mobMask(j, 2))
    end do

    if(nt > 4) cycle


    call bitstring_to_list_in_selection(mobMask(1,1), p(1,1), p(0,1), N_int)
    call bitstring_to_list_in_selection(mobMask(1,2), p(1,2), p(0,2), N_int)

    perMask(1,1) = iand(mask(1,1), not(det(1,1,i)))
    perMask(1,2) = iand(mask(1,2), not(det(1,2,i)))
    do j=2,N_int
      perMask(j,1) = iand(mask(j,1), not(det(j,1,i)))
      perMask(j,2) = iand(mask(j,2), not(det(j,2,i)))
    end do

    call bitstring_to_list_in_selection(perMask(1,1), h(1,1), h(0,1), N_int)
    call bitstring_to_list_in_selection(perMask(1,2), h(1,2), h(0,2), N_int)

    if (interesting(i) >= i_gen) then
       if(nt == 4) then
         call get_d2(interesting(i), det(1,1,i), banned, bannedOrb, indexes, abuf, mask, h, p, sp)
       else if(nt == 3) then
         call get_d1(interesting(i), det(1,1,i), banned, bannedOrb, indexes, abuf, mask, h, p, sp)
       else
         abuf(indexes(0,0)) = interesting(i)
         indexes(0,0) += 1
       end if
    end if
  end do
end subroutine


subroutine get_d2(i_gen, gen, banned, bannedOrb, indexes, abuf, mask, h, p, sp)
  use bitmasks
  implicit none

  integer(bit_kind), intent(in) :: mask(N_int, 2), gen(N_int, 2)
  integer, intent(inout) :: abuf(*)
  integer, intent(in) :: i_gen
  logical, intent(in) :: bannedOrb(mo_num, 2), banned(mo_num, mo_num,2)
  integer, intent(inout) :: indexes(0:mo_num, 0:mo_num)
  integer, intent(in) :: h(0:2,2), p(0:4,2), sp

  !double precision, external :: get_phase_bi
  double precision, external :: mo_two_e_integral

  integer :: i, j, tip, ma, mi, puti, putj
  integer :: h1, h2, p1, p2, i1, i2
  double precision :: hij, phase

  integer, parameter:: turn2d(2,3,4) = reshape((/0,0, 0,0, 0,0,  3,4, 0,0, 0,0,  2,4, 1,4, 0,0,  2,3, 1,3, 1,2 /), (/2,3,4/))
  integer, parameter :: turn2(2) = (/2, 1/)
  integer, parameter :: turn3(2,3) = reshape((/2,3,  1,3, 1,2/), (/2,3/))

  integer :: bant
  integer :: phasemask(2,N_int*bit_kind_size)
  bant = 1

  tip = p(0,1) * p(0,2)

  ma = sp
  if(p(0,1) > p(0,2)) ma = 1
  if(p(0,1) < p(0,2)) ma = 2
  mi = mod(ma, 2) + 1

  if(sp == 3) then
    if(ma == 2) bant = 2

    if(tip == 3) then
      puti = p(1, mi)
      do i = 1, 3
        putj = p(i, ma)
        if(banned(putj,puti,bant)) cycle
        !i1 = turn3(1,i)
        !i2 = turn3(2,i)
        !p1 = p(i1, ma)
        !p2 = p(i2, ma)
        !h1 = h(1, ma)
        !h2 = h(2, ma)

        !hij = (mo_two_e_integral(p1, p2, h1, h2) - mo_two_e_integral(p2,p1, h1, h2)) * get_phase_bi(phasemask, ma, ma, h1, p1, h2, p2)
        if(ma == 1) then
          abuf(indexes(putj, puti)) = i_gen
          indexes(putj, puti) += 1
        else
          abuf(indexes(puti, putj)) = i_gen
          indexes(puti, putj) += 1
        end if
      end do
    else
      !h1 = h(1,1)
      !h2 = h(1,2)
      do j = 1,2
        putj = p(j, 2)
        !p2 = p(turn2(j), 2)
        do i = 1,2
          puti = p(i, 1)

          if(banned(puti,putj,bant)) cycle
          !p1 = p(turn2(i), 1)

          !hij = mo_two_e_integral(p1, p2, h1, h2) * get_phase_bi(phasemask, 1, 2, h1, p1, h2, p2)

          abuf(indexes(puti, putj)) = i_gen
          indexes(puti, putj) += 1
        end do
      end do
    end if

  else
    if(tip == 0) then
      !h1 = h(1, ma)
      !h2 = h(2, ma)
      do i=1,3
      puti = p(i, ma)
      do j=i+1,4
        putj = p(j, ma)
        if(banned(puti,putj,1)) cycle

        !i1 = turn2d(1, i, j)
        !i2 = turn2d(2, i, j)
        !p1 = p(i1, ma)
        !p2 = p(i2, ma)
        !hij = (mo_two_e_integral(p1, p2, h1, h2) - mo_two_e_integral(p2,p1, h1, h2)) * get_phase_bi(phasemask, ma, ma, h1, p1, h2, p2)
        abuf(indexes(puti, putj)) = i_gen
        indexes(puti, putj) += 1
      end do
      end do
    else if(tip == 3) then
      !h1 = h(1, mi)
      !h2 = h(1, ma)
      !p1 = p(1, mi)
      do i=1,3
        puti = p(turn3(1,i), ma)
        putj = p(turn3(2,i), ma)
        if(banned(puti,putj,1)) cycle
        !p2 = p(i, ma)

        !hij = mo_two_e_integral(p1, p2, h1, h2) * get_phase_bi(phasemask, mi, ma, h1, p1, h2, p2)
        abuf(indexes(min(puti, putj), max(puti, putj))) = i_gen
        indexes(min(puti, putj), max(puti, putj)) += 1
      end do
    else ! tip == 4
      puti = p(1, sp)
      putj = p(2, sp)
      if(.not. banned(puti,putj,1)) then
        !p1 = p(1, mi)
        !p2 = p(2, mi)
        !h1 = h(1, mi)
        !h2 = h(2, mi)
        !hij = (mo_two_e_integral(p1, p2, h1, h2) - mo_two_e_integral(p2,p1, h1, h2)) * get_phase_bi(phasemask, mi, mi, h1, p1, h2, p2)

        abuf(indexes(puti, putj)) = i_gen
        indexes(puti, putj) += 1
      end if
    end if
  end if
end


subroutine get_d1(i_gen, gen, banned, bannedOrb, indexes, abuf, mask, h, p, sp)
  use bitmasks
  implicit none

  integer(bit_kind), intent(in)  :: mask(N_int, 2), gen(N_int, 2)
  integer, intent(inout)         :: abuf(*)
  integer,intent(in)             :: i_gen
  logical, intent(in)            :: bannedOrb(mo_num, 2), banned(mo_num, mo_num,2)
  integer(bit_kind)              :: det(N_int, 2)
  integer, intent(inout)         :: indexes(0:mo_num, 0:mo_num)
  integer, intent(in)            :: h(0:2,2), p(0:4,2), sp
  double precision               :: hij, tmp_row(N_states, mo_num), tmp_row2(N_states, mo_num)
  !double precision, external     :: get_phase_bi
  double precision, external     :: mo_two_e_integral
  logical                        :: ok

  logical, allocatable           :: lbanned(:,:)
  integer                        :: puti, putj, ma, mi, s1, s2, i, i1, i2, j
  integer                        :: hfix, pfix, h1, h2, p1, p2, ib

  integer, parameter             :: turn2(2) = (/2,1/)
  integer, parameter             :: turn3(2,3) = reshape((/2,3,  1,3, 1,2/), (/2,3/))

  integer                        :: bant
  integer :: phasemask(2,N_int*bit_kind_size)


  allocate (lbanned(mo_num, 2))
  lbanned = bannedOrb

  do i=1, p(0,1)
    lbanned(p(i,1), 1) = .true.
  end do
  do i=1, p(0,2)
    lbanned(p(i,2), 2) = .true.
  end do

  ma = 1
  if(p(0,2) >= 2) ma = 2
  mi = turn2(ma)

  bant = 1

  if(sp == 3) then
    !move MA
    if(ma == 2) bant = 2
    puti = p(1,mi)
    hfix = h(1,ma)
    p1 = p(1,ma)
    p2 = p(2,ma)
    if(.not. bannedOrb(puti, mi)) then
      !tmp_row = 0d0
      !do putj=1, hfix-1
      !  if(lbanned(putj, ma) .or. banned(putj, puti,bant)) cycle
      !  hij = (mo_two_e_integral(p1, p2, putj, hfix)-mo_two_e_integral(p2,p1,putj,hfix)) * get_phase_bi(phasemask, ma, ma, putj, p1, hfix, p2)
      !  tmp_row(1:N_states,putj) += hij * coefs(1:N_states)
      !end do
      !do putj=hfix+1, mo_num
      !  if(lbanned(putj, ma) .or. banned(putj, puti,bant)) cycle
      !  hij = (mo_two_e_integral(p1, p2, hfix, putj)-mo_two_e_integral(p2,p1,hfix,putj)) * get_phase_bi(phasemask, ma, ma, hfix, p1, putj, p2)
      !  tmp_row(1:N_states,putj) += hij * coefs(1:N_states)
      !end do

      if(ma == 1) then
        !mat(1:N_states,1:mo_num,puti) += tmp_row(1:N_states,1:mo_num)
        abuf(indexes(0, puti)) = i_gen
        indexes(0, puti) += 1
        !countedOrb(puti, 2) -= 1
      else
        !mat(1:N_states,puti,1:mo_num) += tmp_row(1:N_states,1:mo_num)
        abuf(indexes(puti, 0)) = i_gen
        indexes(puti, 0) += 1
        !countedOrb(puti, 1) -= 1
      end if
    end if

    !MOVE MI
    !pfix = p(1,mi)
    !tmp_row = 0d0
    !tmp_row2 = 0d0
    !do puti=1,mo_num
    !  if(lbanned(puti,mi)) cycle
      !p1 fixed
    !  putj = p1
      !if(.not. banned(putj,puti,bant)) then
      !  hij = mo_two_e_integral(p2,pfix,hfix,puti) * get_phase_bi(phasemask, ma, mi, hfix, p2, puti, pfix)
      !  tmp_row(:,puti) += hij * coefs(:)
      !end if

    !  putj = p2
      !if(.not. banned(putj,puti,bant)) then
      !  hij = mo_two_e_integral(p1,pfix,hfix,puti) * get_phase_bi(phasemask, ma, mi, hfix, p1, puti, pfix)
      !  tmp_row2(:,puti) += hij * coefs(:)
      !end if
    !end do

    if(mi == 1) then
      if(.not. bannedOrb(p1, 2)) then
        abuf(indexes(0,p1)) = i_gen
        indexes(0,p1) += 1
      end if
      if(.not. bannedOrb(p2, 2)) then
        abuf(indexes(0,p2)) = i_gen
        indexes(0,p2) += 1
      end if
    else
      if(.not. bannedOrb(p1, 1)) then
        abuf(indexes(p1,0)) = i_gen
        indexes(p1,0) += 1
      end if
      if(.not. bannedOrb(p2, 1)) then
        abuf(indexes(p2,0)) = i_gen
        indexes(p2,0) += 1
      end if
    end if
  else
    if(p(0,ma) == 3) then
      do i=1,3
        !hfix = h(1,ma)
        puti = p(i, ma)
        !p1 = p(turn3(1,i), ma)
        !p2 = p(turn3(2,i), ma)
        !tmp_row = 0d0
        !do putj=1,hfix-1
        !  if(lbanned(putj,ma) .or. banned(puti,putj,1)) cycle
        !  hij = (mo_two_e_integral(p1, p2, putj, hfix)-mo_two_e_integral(p2,p1,putj,hfix)) * get_phase_bi(phasemask, ma, ma, putj, p1, hfix, p2)
        !  tmp_row(:,putj) += hij * coefs(:)
        !end do
        !do putj=hfix+1,mo_num
        !  if(lbanned(putj,ma) .or. banned(puti,putj,1)) cycle
        !  hij = (mo_two_e_integral(p1, p2, hfix, putj)-mo_two_e_integral(p2,p1,hfix,putj)) * get_phase_bi(phasemask, ma, ma, hfix, p1, putj, p2)
        !  tmp_row(:,putj) += hij * coefs(:)
        !end do

        !mat(:, :puti-1, puti) += tmp_row(:,:puti-1)
        !mat(:, puti, puti:) += tmp_row(:,puti:)
        if(.not. bannedOrb(puti, sp)) then
          if(sp == 1) then
            abuf(indexes(puti, 0)) = i_gen
            indexes(puti, 0) += 1
          else
            abuf(indexes(0, puti)) = i_gen
            indexes(0, puti) += 1
          end if
        end if
      end do
    else
      !hfix = h(1,mi)
      !pfix = p(1,mi)
      p1 = p(1,ma)
      p2 = p(2,ma)
      !tmp_row = 0d0
      !tmp_row2 = 0d0
      !do puti=1,mo_num
      !  if(lbanned(puti,ma)) cycle
      !  putj = p2
        !if(.not. banned(puti,putj,1)) then
        !  hij = mo_two_e_integral(pfix, p1, hfix, puti) * get_phase_bi(phasemask, mi, ma, hfix, pfix, puti, p1)
        !  tmp_row(:,puti) += hij * coefs(:)
        !end if

      !  putj = p1
        !if(.not. banned(puti,putj,1)) then
        !  hij = mo_two_e_integral(pfix, p2, hfix, puti) * get_phase_bi(phasemask, mi, ma, hfix, pfix, puti, p2)
        !  tmp_row2(:,puti) += hij * coefs(:)
        !end if
      !end do
      if(.not. bannedOrb(p2, sp)) then
        if(sp == 1) then
          abuf(indexes(p2, 0)) = i_gen
          indexes(p2, 0) += 1
        else
          abuf(indexes(0, p2)) = i_gen
          indexes(0, p2) += 1
        end if
      end if
      if(.not. bannedOrb(p1, sp)) then
        if(sp == 1) then
          abuf(indexes(p1, 0)) = i_gen
          indexes(p1, 0) += 1
        else
          abuf(indexes(0, p1)) = i_gen
          indexes(0, p1) += 1
        end if
      end if
    end if
  end if

 !! MONO
  !  if(sp == 3) then
  !    s1 = 1
  !    s2 = 2
  !  else
  !    s1 = sp
  !    s2 = sp
  !  end if
!
!    do i1=1,p(0,s1)
!      ib = 1
!      if(s1 == s2) ib = i1+1
!      do i2=ib,p(0,s2)
!        p1 = p(i1,s1)
!        p2 = p(i2,s2)
 !       if(bannedOrb(p1, s1) .or. bannedOrb(p2, s2) .or. banned(p1, p2, 1)) cycle
 !       call apply_particles(mask, s1, p1, s2, p2, det, ok, N_int)
!        call i_h_j(gen, det, N_int, hij)
!        !mat(:, p1, p2) += coefs(:) * hij
!        !!!!!!!! DUPLICTATE counted(p1, p2) !!!!!!!!!!!!!!!!!!!!
!      end do
!    end do
end


subroutine past_d1(bannedOrb, p)
  use bitmasks
  implicit none

  logical, intent(inout) :: bannedOrb(mo_num, 2)
  integer, intent(in) :: p(0:4, 2)
  integer :: i,s

  do s = 1, 2
    do i = 1, p(0, s)
      bannedOrb(p(i, s), s) = .true.
    end do
  end do
end


subroutine past_d2(banned, p, sp)
  use bitmasks
  implicit none

  logical, intent(inout) :: banned(mo_num, mo_num)
  integer, intent(in) :: p(0:4, 2), sp
  integer :: i,j

  if(sp == 3) then
    do i=1,p(0,1)
      do j=1,p(0,2)
        banned(p(i,1), p(j,2)) = .true.
      end do
    end do
  else
    do i=1,p(0, sp)
      do j=1,i-1
        banned(p(j,sp), p(i,sp)) = .true.
        banned(p(i,sp), p(j,sp)) = .true.
      end do
    end do
  end if
end


subroutine count_d1(countedOrb, p)
  use bitmasks
  implicit none

  integer, intent(inout) :: countedOrb(mo_num, 2)
  integer, intent(in) :: p(0:4, 2)
  integer :: i,s

  do s = 1, 2
    do i = 1, p(0, s)
      countedOrb(p(i, s), s) += 1
    end do
  end do
end


subroutine count_d2(counted, p, sp)
  use bitmasks
  implicit none

  integer, intent(inout) :: counted(mo_num, mo_num)
  integer, intent(in) :: p(0:4, 2), sp
  integer :: i,j

  if(sp == 3) then
    do i=1,p(0,1)
      do j=1,p(0,2)
        counted(p(i,1), p(j,2)) += 1
      end do
    end do
  else
    do i=1,p(0, sp)
      do j=1,i-1
        counted(p(j,sp), p(i,sp)) += 1
      end do
    end do
  end if
end



subroutine spot_isinwf(mask, det, i_gen, N, banned, fullMatch, interesting)
  use bitmasks
  implicit none

  integer, intent(in) :: i_gen, N
  integer, intent(in) :: interesting(0:N)
  integer(bit_kind),intent(in) :: mask(N_int, 2), det(N_int, 2, N)
  logical, intent(inout) :: banned(mo_num, mo_num)
  logical, intent(out) :: fullMatch


  integer :: i, j, na, nb, list(3)
  integer(bit_kind) :: myMask(N_int, 2), negMask(N_int, 2)

  fullMatch = .false.

  do i=1,N_int
    negMask(i,1) = not(mask(i,1))
    negMask(i,2) = not(mask(i,2))
  end do

  genl : do i=1, N
    do j=1, N_int
      if(iand(det(j,1,i), mask(j,1)) /= mask(j, 1)) cycle genl
      if(iand(det(j,2,i), mask(j,2)) /= mask(j, 2)) cycle genl
    end do

    if(interesting(i) < i_gen) then
      fullMatch = .true.
      return
    end if

    do j=1, N_int
      myMask(j, 1) = iand(det(j, 1, i), negMask(j, 1))
      myMask(j, 2) = iand(det(j, 2, i), negMask(j, 2))
    end do

    call bitstring_to_list_in_selection(myMask(1,1), list(1), na, N_int)
    call bitstring_to_list_in_selection(myMask(1,2), list(na+1), nb, N_int)
    banned(list(1), list(2)) = .true.
  end do genl
end


subroutine bitstring_to_list_in_selection( string, list, n_elements, Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Gives the inidices(+1) of the bits set to 1 in the bit string
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: string(Nint)
  integer, intent(out)           :: list(Nint*bit_kind_size)
  integer, intent(out)           :: n_elements

  integer                        :: i, ishift
  integer(bit_kind)              :: l

  n_elements = 0
  ishift = 2
  do i=1,Nint
    l = string(i)
    do while (l /= 0_bit_kind)
      n_elements = n_elements+1
      list(n_elements) = ishift+popcnt(l-1_bit_kind) - popcnt(l)
      l = iand(l,l-1_bit_kind)
    enddo
    ishift = ishift + bit_kind_size
  enddo

end


