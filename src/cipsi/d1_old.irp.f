
subroutine get_d1_complex_old(gen, phasemask, bannedOrb, banned, mat, mask, h, p, sp, coefs)
  !todo: indices should be okay for complex?
  use bitmasks
  implicit none

  integer(bit_kind), intent(in)  :: mask(N_int, 2), gen(N_int, 2)
  integer(bit_kind), intent(in)  :: phasemask(N_int,2)
  logical, intent(in)            :: bannedOrb(mo_num, 2), banned(mo_num, mo_num,2)
  integer(bit_kind)              :: det(N_int, 2)
  complex*16, intent(in)   :: coefs(N_states)
  complex*16, intent(inout) :: mat(N_states, mo_num, mo_num)
  integer, intent(in)            :: h(0:2,2), p(0:4,2), sp
  double precision, external     :: get_phase_bi
  complex*16, external     :: mo_two_e_integral_complex
  logical                        :: ok

  logical, allocatable           :: lbanned(:,:)
  integer                        :: puti, putj, ma, mi, s1, s2, i, i1, i2, j
  integer                        :: hfix, pfix, h1, h2, p1, p2, ib, k, l

  integer, parameter             :: turn2(2) = (/2,1/)
  integer, parameter             :: turn3(2,3) = reshape((/2,3,  1,3, 1,2/), (/2,3/))

  integer                        :: bant
  complex*16, allocatable :: hij_cache(:,:)
  complex*16               :: hij, tmp_row(N_states, mo_num), tmp_row2(N_states, mo_num)
  PROVIDE mo_integrals_map N_int

  allocate (lbanned(mo_num, 2))
  allocate (hij_cache(mo_num,2))
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
      call get_mo_two_e_integrals_complex(hfix,p1,p2,mo_num,hij_cache(1,1),mo_integrals_map,mo_integrals_map_2)
      call get_mo_two_e_integrals_complex(hfix,p2,p1,mo_num,hij_cache(1,2),mo_integrals_map,mo_integrals_map_2)
      tmp_row = (0.d0,0.d0)
      do putj=1, hfix-1
        if(lbanned(putj, ma)) cycle
        if(banned(putj, puti,bant)) cycle
        hij = hij_cache(putj,1) - hij_cache(putj,2)
        if (hij /= (0.d0,0.d0)) then
          hij = hij * get_phase_bi(phasemask, ma, ma, putj, p1, hfix, p2, N_int)
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            tmp_row(k,putj) = tmp_row(k,putj) + hij * coefs(k)
          enddo
        endif
      end do
      do putj=hfix+1, mo_num
        if(lbanned(putj, ma)) cycle
        if(banned(putj, puti,bant)) cycle
        hij = hij_cache(putj,2) - hij_cache(putj,1)
        if (hij /= (0.d0,0.d0)) then
          hij = hij * get_phase_bi(phasemask, ma, ma, hfix, p1, putj, p2, N_int)
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            tmp_row(k,putj) = tmp_row(k,putj) + hij * coefs(k)
          enddo
        endif
      end do

      if(ma == 1) then
        mat(1:N_states,1:mo_num,puti) = mat(1:N_states,1:mo_num,puti) + tmp_row(1:N_states,1:mo_num)
      else
        do l=1,mo_num
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat(k,puti,l) = mat(k,puti,l) + tmp_row(k,l)
          enddo
        enddo
      end if
    end if

    !MOVE MI
    pfix = p(1,mi)
    tmp_row = (0.d0,0.d0)
    tmp_row2 = (0.d0,0.d0)
    call get_mo_two_e_integrals_complex(hfix,pfix,p1,mo_num,hij_cache(1,1),mo_integrals_map,mo_integrals_map_2)
    call get_mo_two_e_integrals_complex(hfix,pfix,p2,mo_num,hij_cache(1,2),mo_integrals_map,mo_integrals_map_2)
    putj = p1
    do puti=1,mo_num !HOT
      if(lbanned(puti,mi)) cycle
      !p1 fixed
      putj = p1
      if(.not. banned(putj,puti,bant)) then
        hij = hij_cache(puti,2)
        if (hij /= (0.d0,0.d0)) then
          hij = hij * get_phase_bi(phasemask, ma, mi, hfix, p2, puti, pfix, N_int)
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            tmp_row(k,puti) = tmp_row(k,puti) + hij * coefs(k)
          enddo
        endif
      end if
!    enddo
!      
      putj = p2
!    do puti=1,mo_num !HOT
      if(.not. banned(putj,puti,bant)) then
        hij = hij_cache(puti,1)
        if (hij /= (0.d0,0.d0)) then
          hij = hij * get_phase_bi(phasemask, ma, mi, hfix, p1, puti, pfix, N_int)
          do k=1,N_states
            tmp_row2(k,puti) = tmp_row2(k,puti) + hij * coefs(k)
          enddo
        endif
      end if
    end do

    if(mi == 1) then
      mat(:,:,p1) = mat(:,:,p1) + tmp_row(:,:)
      mat(:,:,p2) = mat(:,:,p2) + tmp_row2(:,:)
    else
      do l=1,mo_num
        !DIR$ LOOP COUNT AVG(4)
        do k=1,N_states
          mat(k,p1,l) = mat(k,p1,l) + tmp_row(k,l)
          mat(k,p2,l) = mat(k,p2,l) + tmp_row2(k,l)
        enddo
      enddo
    end if

  else  ! sp /= 3

    if(p(0,ma) == 3) then
      do i=1,3
        hfix = h(1,ma)
        puti = p(i, ma)
        p1 = p(turn3(1,i), ma)
        p2 = p(turn3(2,i), ma)
        call get_mo_two_e_integrals_complex(hfix,p1,p2,mo_num,hij_cache(1,1),mo_integrals_map,mo_integrals_map_2)
        call get_mo_two_e_integrals_complex(hfix,p2,p1,mo_num,hij_cache(1,2),mo_integrals_map,mo_integrals_map_2)
        tmp_row = (0.d0,0.d0)
        do putj=1,hfix-1
          if(banned(putj,puti,1)) cycle
          if(lbanned(putj,ma)) cycle
          hij = hij_cache(putj,1) - hij_cache(putj,2)
          if (hij /= (0.d0,0.d0)) then
            hij = hij * get_phase_bi(phasemask, ma, ma, putj, p1, hfix, p2, N_int)
            tmp_row(:,putj) = tmp_row(:,putj) + hij * coefs(:)
          endif
        end do
        do putj=hfix+1,mo_num
          if(banned(putj,puti,1)) cycle
          if(lbanned(putj,ma)) cycle
          hij = hij_cache(putj,2) - hij_cache(putj,1)
          if (hij /= (0.d0,0.d0)) then
            hij = hij * get_phase_bi(phasemask, ma, ma, hfix, p1, putj, p2, N_int)
            tmp_row(:,putj) = tmp_row(:,putj) + hij * coefs(:)
          endif
        end do

        mat(:, :puti-1, puti) = mat(:, :puti-1, puti) + tmp_row(:,:puti-1)
        do l=puti,mo_num
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat(k, puti, l) = mat(k, puti,l) + tmp_row(k,l)
          enddo
        enddo
      end do
    else
      hfix = h(1,mi)
      pfix = p(1,mi)
      p1 = p(1,ma)
      p2 = p(2,ma)
      tmp_row = (0.d0,0.d0)
      tmp_row2 = (0.d0,0.d0)
      call get_mo_two_e_integrals_complex(hfix,p1,pfix,mo_num,hij_cache(1,1),mo_integrals_map,mo_integrals_map_2)
      call get_mo_two_e_integrals_complex(hfix,p2,pfix,mo_num,hij_cache(1,2),mo_integrals_map,mo_integrals_map_2)
      putj = p2
      do puti=1,mo_num
        if(lbanned(puti,ma)) cycle
        putj = p2
        if(.not. banned(puti,putj,1)) then
          hij = hij_cache(puti,1)
          if (hij /= (0.d0,0.d0)) then
            hij = hij * get_phase_bi(phasemask, mi, ma, hfix, pfix, puti, p1, N_int)
            !DIR$ LOOP COUNT AVG(4)
            do k=1,N_states
              tmp_row(k,puti) = tmp_row(k,puti) + hij * coefs(k)
            enddo
          endif
        end if

        putj = p1
        if(.not. banned(puti,putj,1)) then
          hij = hij_cache(puti,2)
          if (hij /= (0.d0,0.d0)) then
            hij = hij * get_phase_bi(phasemask, mi, ma, hfix, pfix, puti, p2, N_int)
            do k=1,N_states
              tmp_row2(k,puti) = tmp_row2(k,puti) + hij * coefs(k)
            enddo
          endif
        end if
      end do
      mat(:,:p2-1,p2) = mat(:,:p2-1,p2) + tmp_row(:,:p2-1)
      do l=p2,mo_num
        !DIR$ LOOP COUNT AVG(4)
        do k=1,N_states
          mat(k,p2,l) = mat(k,p2,l) + tmp_row(k,l)
        enddo
      enddo
      mat(:,:p1-1,p1) = mat(:,:p1-1,p1) + tmp_row2(:,:p1-1)
      do l=p1,mo_num
        !DIR$ LOOP COUNT AVG(4)
        do k=1,N_states
          mat(k,p1,l) = mat(k,p1,l) + tmp_row2(k,l)
        enddo
      enddo
    end if
  end if
  deallocate(lbanned,hij_cache)

 !! MONO
    if(sp == 3) then
      s1 = 1
      s2 = 2
    else
      s1 = sp
      s2 = sp
    end if

    do i1=1,p(0,s1)
      ib = 1
      if(s1 == s2) ib = i1+1
      do i2=ib,p(0,s2)
        p1 = p(i1,s1)
        p2 = p(i2,s2)
        if(bannedOrb(p1, s1) .or. bannedOrb(p2, s2) .or. banned(p1, p2, 1)) cycle
        call apply_particles(mask, s1, p1, s2, p2, det, ok, N_int)
        ! gen is a selector; mask is ionized generator; det is alpha
        ! hij is contribution to <psi|H|alpha>
        call i_h_j_complex(gen, det, N_int, hij)
        !DIR$ LOOP COUNT AVG(4)
        do k=1,N_states
          ! take conjugate to get contribution to <alpha|H|psi> instead of <psi|H|alpha>
          mat(k, p1, p2) = mat(k, p1, p2) + coefs(k) * dconjg(hij)
        enddo
      end do
    end do
end

