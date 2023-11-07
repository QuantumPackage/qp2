subroutine get_d1_new(gen, phasemask, bannedOrb, banned, mat_l, mat_r, mask, h, p, sp, coefs)
  !todo: indices should be okay for complex?
  use bitmasks
  implicit none

  integer(bit_kind), intent(in)  :: mask(N_int, 2), gen(N_int, 2)
  integer(bit_kind), intent(in)  :: phasemask(N_int,2)
  logical, intent(in)            :: bannedOrb(mo_num, 2), banned(mo_num, mo_num,2)
  integer(bit_kind)              :: det(N_int, 2)
  double precision, intent(in)   :: coefs(N_states,2)
  double precision, intent(inout) :: mat_l(N_states, mo_num, mo_num)
  double precision, intent(inout) :: mat_r(N_states, mo_num, mo_num)
  integer, intent(in)            :: h(0:2,2), p(0:4,2), sp
  double precision, external     :: get_phase_bi
  double precision, external     :: mo_two_e_integral_complex
  logical                        :: ok

  logical, allocatable           :: lbanned(:,:)
  integer                        :: puti, putj, ma, mi, s1, s2, i, i1, i2, j
  integer                        :: hfix, pfix, h1, h2, p1, p2, ib, k, l, mm

  integer, parameter             :: turn2(2) = (/2,1/)
  integer, parameter             :: turn3(2,3) = reshape((/2,3,  1,3, 1,2/), (/2,3/))

  integer                        :: bant
  double precision, allocatable :: hij_cache(:,:)
  double precision               :: hij, tmp_rowij(N_states, mo_num), tmp_rowij2(N_states, mo_num)
  double precision, allocatable :: hji_cache(:,:)
  double precision               :: hji, tmp_rowji(N_states, mo_num), tmp_rowji2(N_states, mo_num)
!  PROVIDE mo_integrals_map N_int
!  print*,'in get_d1_new'
!  call debug_det(gen,N_int)
!  print*,'coefs',coefs(1,:)

  allocate (lbanned(mo_num, 2))
  allocate (hij_cache(mo_num,2))
  allocate (hji_cache(mo_num,2))
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
!      call get_mo_two_e_integrals_complex(hfix,p1,p2,mo_num,hij_cache(1,1),mo_integrals_map,mo_integrals_map_2)
!      call get_mo_two_e_integrals_complex(hfix,p2,p1,mo_num,hij_cache(1,2),mo_integrals_map,mo_integrals_map_2)
      do mm = 1, mo_num
       hij_cache(mm,1) = mo_bi_ortho_tc_two_e(mm,hfix,p1,p2)
       hij_cache(mm,2) = mo_bi_ortho_tc_two_e(mm,hfix,p2,p1)
       hji_cache(mm,1) = mo_bi_ortho_tc_two_e(p1,p2,mm,hfix)
       hji_cache(mm,2) = mo_bi_ortho_tc_two_e(p2,p1,mm,hfix)
      enddo
      !! <alpha|H|psi>
      tmp_rowij = 0.d0
      do putj=1, hfix-1
        if(lbanned(putj, ma)) cycle
        if(banned(putj, puti,bant)) cycle
        hij = hij_cache(putj,1) - hij_cache(putj,2)
        if (hij /= 0.d0) then
          hij = hij * get_phase_bi(phasemask, ma, ma, putj, p1, hfix, p2, N_int)
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            tmp_rowij(k,putj) = tmp_rowij(k,putj) + hij * coefs(k,2)
          enddo
        endif
      end do
      do putj=hfix+1, mo_num
        if(lbanned(putj, ma)) cycle
        if(banned(putj, puti,bant)) cycle
        hij = hij_cache(putj,2) - hij_cache(putj,1)
        if (hij /= 0.d0) then
          hij = hij * get_phase_bi(phasemask, ma, ma, hfix, p1, putj, p2, N_int)
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            tmp_rowij(k,putj) = tmp_rowij(k,putj) + hij * coefs(k,2)
          enddo
        endif
      end do

      if(ma == 1) then
        mat_r(1:N_states,1:mo_num,puti) = mat_r(1:N_states,1:mo_num,puti) + tmp_rowij(1:N_states,1:mo_num)
      else
        do l=1,mo_num
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat_r(k,puti,l) = mat_r(k,puti,l) + tmp_rowij(k,l)
          enddo
        enddo
      end if

      !! <phi|H|alpha>
      tmp_rowji = 0.d0
      do putj=1, hfix-1
        if(lbanned(putj, ma)) cycle
        if(banned(putj, puti,bant)) cycle
        hji = hji_cache(putj,1) - hji_cache(putj,2)
        if (hji /= 0.d0) then
          hji = hji * get_phase_bi(phasemask, ma, ma, putj, p1, hfix, p2, N_int)
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            tmp_rowji(k,putj) = tmp_rowji(k,putj) + hji * coefs(k,1)
          enddo
        endif
      end do
      do putj=hfix+1, mo_num
        if(lbanned(putj, ma)) cycle
        if(banned(putj, puti,bant)) cycle
        hji = hji_cache(putj,2) - hji_cache(putj,1)
        if (hji /= 0.d0) then
          hji = hji * get_phase_bi(phasemask, ma, ma, hfix, p1, putj, p2, N_int)
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            tmp_rowji(k,putj) = tmp_rowji(k,putj) + hji * coefs(k,1)
          enddo
        endif
      end do

      if(ma == 1) then
        mat_l(1:N_states,1:mo_num,puti) = mat_l(1:N_states,1:mo_num,puti) + tmp_rowji(1:N_states,1:mo_num)
      else
        do l=1,mo_num
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat_l(k,puti,l) = mat_l(k,puti,l) + tmp_rowji(k,l)
          enddo
        enddo
      end if
    end if

    !MOVE MI
    pfix = p(1,mi)
    tmp_rowij  = 0.d0
    tmp_rowij2 = 0.d0
    tmp_rowji  = 0.d0
    tmp_rowji2 = 0.d0
!    call get_mo_two_e_integrals_complex(hfix,pfix,p1,mo_num,hij_cache(1,1),mo_integrals_map,mo_integrals_map_2)
!    call get_mo_two_e_integrals_complex(hfix,pfix,p2,mo_num,hij_cache(1,2),mo_integrals_map,mo_integrals_map_2)
    do mm = 1, mo_num
     hij_cache(mm,1) = mo_bi_ortho_tc_two_e(mm,hfix,pfix,p1)
     hij_cache(mm,2) = mo_bi_ortho_tc_two_e(mm,hfix,pfix,p2)
     hji_cache(mm,1) = mo_bi_ortho_tc_two_e(pfix,p1,mm,hfix)
     hji_cache(mm,2) = mo_bi_ortho_tc_two_e(pfix,p2,mm,hfix)
    enddo
    putj = p1
    !! <alpha|H|psi>
    do puti=1,mo_num !HOT
      if(lbanned(puti,mi)) cycle
      !p1 fixed
      putj = p1
      if(.not. banned(putj,puti,bant)) then
        hij = hij_cache(puti,2)
        if (hij /= 0.d0) then
          hij = hij * get_phase_bi(phasemask, ma, mi, hfix, p2, puti, pfix, N_int)
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            tmp_rowij(k,puti) = tmp_rowij(k,puti) + hij * coefs(k,2)
          enddo
        endif
      end if
!      
      putj = p2
      if(.not. banned(putj,puti,bant)) then
        hij = hij_cache(puti,1)
        if (hij /= 0.d0) then
          hij = hij * get_phase_bi(phasemask, ma, mi, hfix, p1, puti, pfix, N_int)
          do k=1,N_states
            tmp_rowij2(k,puti) = tmp_rowij2(k,puti) + hij * coefs(k,2)
          enddo
        endif
      end if
    end do

    if(mi == 1) then
      mat_r(:,:,p1) = mat_r(:,:,p1) + tmp_rowij(:,:)
      mat_r(:,:,p2) = mat_r(:,:,p2) + tmp_rowij2(:,:)
    else
      do l=1,mo_num
        !DIR$ LOOP COUNT AVG(4)
        do k=1,N_states
          mat_r(k,p1,l) = mat_r(k,p1,l) + tmp_rowij(k,l)
          mat_r(k,p2,l) = mat_r(k,p2,l) + tmp_rowij2(k,l)
        enddo
      enddo
    end if

    putj = p1
    !! <phi|H|alpha>
    do puti=1,mo_num !HOT
      if(lbanned(puti,mi)) cycle
      !p1 fixed
      putj = p1
      if(.not. banned(putj,puti,bant)) then
        hji = hji_cache(puti,2)
        if (hji /= 0.d0) then
          hji = hji * get_phase_bi(phasemask, ma, mi, hfix, p2, puti, pfix, N_int)
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            tmp_rowji(k,puti) = tmp_rowji(k,puti) + hji * coefs(k,1)
          enddo
        endif
      end if
!      
      putj = p2
      if(.not. banned(putj,puti,bant)) then
        hji = hji_cache(puti,1)
        if (hji /= 0.d0) then
          hji = hji * get_phase_bi(phasemask, ma, mi, hfix, p1, puti, pfix, N_int)
          do k=1,N_states
            tmp_rowji2(k,puti) = tmp_rowji2(k,puti) + hji * coefs(k,1)
          enddo
        endif
      end if
    end do

    if(mi == 1) then
      mat_l(:,:,p1) = mat_l(:,:,p1) + tmp_rowji(:,:)
      mat_l(:,:,p2) = mat_l(:,:,p2) + tmp_rowji2(:,:)
    else
      do l=1,mo_num
        !DIR$ LOOP COUNT AVG(4)
        do k=1,N_states
          mat_l(k,p1,l) = mat_l(k,p1,l) + tmp_rowji(k,l)
          mat_l(k,p2,l) = mat_l(k,p2,l) + tmp_rowji2(k,l)
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
!        call get_mo_two_e_integrals_complex(hfix,p1,p2,mo_num,hij_cache(1,1),mo_integrals_map,mo_integrals_map_2)
!        call get_mo_two_e_integrals_complex(hfix,p2,p1,mo_num,hij_cache(1,2),mo_integrals_map,mo_integrals_map_2)
        do mm = 1, mo_num
         hij_cache(mm,1) = mo_bi_ortho_tc_two_e(mm,hfix,p1,p2)
         hij_cache(mm,2) = mo_bi_ortho_tc_two_e(mm,hfix,p2,p1)
         hji_cache(mm,1) = mo_bi_ortho_tc_two_e(p1,p2,mm,hfix)
         hji_cache(mm,2) = mo_bi_ortho_tc_two_e(p2,p1,mm,hfix)
        enddo
    !! <alpha|H|psi>
        tmp_rowij = 0.d0
        do putj=1,hfix-1
          if(banned(putj,puti,1)) cycle
          if(lbanned(putj,ma)) cycle
          hij = hij_cache(putj,1) - hij_cache(putj,2)
          if (hij /= 0.d0) then
            hij = hij * get_phase_bi(phasemask, ma, ma, putj, p1, hfix, p2, N_int)
            tmp_rowij(:,putj) = tmp_rowij(:,putj) + hij * coefs(:,2)
          endif
        end do
        do putj=hfix+1,mo_num
          if(banned(putj,puti,1)) cycle
          if(lbanned(putj,ma)) cycle
          hij = hij_cache(putj,2) - hij_cache(putj,1)
          if (hij /= 0.d0) then
            hij = hij * get_phase_bi(phasemask, ma, ma, hfix, p1, putj, p2, N_int)
            tmp_rowij(:,putj) = tmp_rowij(:,putj) + hij * coefs(:,2)
          endif
        end do

        mat_r(:, :puti-1, puti) = mat_r(:, :puti-1, puti) + tmp_rowij(:,:puti-1)
        do l=puti,mo_num
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat_r(k, puti, l) = mat_r(k, puti,l) + tmp_rowij(k,l)
          enddo
        enddo
    !! <phi|H|alpha>
        tmp_rowji = 0.d0
        do putj=1,hfix-1
          if(banned(putj,puti,1)) cycle
          if(lbanned(putj,ma)) cycle
          hji = hji_cache(putj,1) - hji_cache(putj,2)
          if (hji /= 0.d0) then
            hji = hji * get_phase_bi(phasemask, ma, ma, putj, p1, hfix, p2, N_int)
            tmp_rowji(:,putj) = tmp_rowji(:,putj) + hji * coefs(:,1)
          endif
        end do
        do putj=hfix+1,mo_num
          if(banned(putj,puti,1)) cycle
          if(lbanned(putj,ma)) cycle
          hji = hji_cache(putj,2) - hji_cache(putj,1)
          if (hji /= 0.d0) then
            hji = hji * get_phase_bi(phasemask, ma, ma, hfix, p1, putj, p2, N_int)
            tmp_rowji(:,putj) = tmp_rowji(:,putj) + hji * coefs(:,1)
          endif
        end do

        mat_l(:, :puti-1, puti) = mat_l(:, :puti-1, puti) + tmp_rowji(:,:puti-1)
        do l=puti,mo_num
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat_l(k, puti, l) = mat_l(k, puti,l) + tmp_rowji(k,l)
          enddo
        enddo
      end do
    else
      hfix = h(1,mi)
      pfix = p(1,mi)
      p1 = p(1,ma)
      p2 = p(2,ma)
      tmp_rowij =  0.d0
      tmp_rowij2 = 0.d0
      tmp_rowji =  0.d0
      tmp_rowji2 = 0.d0
!      call get_mo_two_e_integrals_complex(hfix,p1,pfix,mo_num,hij_cache(1,1),mo_integrals_map,mo_integrals_map_2)
!      call get_mo_two_e_integrals_complex(hfix,p2,pfix,mo_num,hij_cache(1,2),mo_integrals_map,mo_integrals_map_2)
      do mm = 1, mo_num
       hij_cache(mm,1) = mo_bi_ortho_tc_two_e(mm,hfix,p1,pfix)
       hij_cache(mm,2) = mo_bi_ortho_tc_two_e(mm,hfix,p2,pfix)
       hji_cache(mm,1) = mo_bi_ortho_tc_two_e(p1,pfix,mm,hfix)
       hji_cache(mm,2) = mo_bi_ortho_tc_two_e(p2,pfix,mm,hfix)
      enddo
      putj = p2
    !! <alpha|H|psi>
      do puti=1,mo_num
        if(lbanned(puti,ma)) cycle
        putj = p2
        if(.not. banned(puti,putj,1)) then
          hij = hij_cache(puti,1)
          if (hij /= 0.d0) then
            hij = hij * get_phase_bi(phasemask, mi, ma, hfix, pfix, puti, p1, N_int)
            !DIR$ LOOP COUNT AVG(4)
            do k=1,N_states
              tmp_rowij(k,puti) = tmp_rowij(k,puti) + hij * coefs(k,2)
            enddo
          endif
        end if

        putj = p1
        if(.not. banned(puti,putj,1)) then
          hij = hij_cache(puti,2)
          if (hij /= 0.d0) then
            hij = hij * get_phase_bi(phasemask, mi, ma, hfix, pfix, puti, p2, N_int)
            do k=1,N_states
              tmp_rowij2(k,puti) = tmp_rowij2(k,puti) + hij * coefs(k,2)
            enddo
          endif
        end if
      end do
      mat_r(:,:p2-1,p2) = mat_r(:,:p2-1,p2) + tmp_rowij(:,:p2-1)
      do l=p2,mo_num
        !DIR$ LOOP COUNT AVG(4)
        do k=1,N_states
          mat_r(k,p2,l) = mat_r(k,p2,l) + tmp_rowij(k,l)
        enddo
      enddo
      mat_r(:,:p1-1,p1) = mat_r(:,:p1-1,p1) + tmp_rowij2(:,:p1-1)
      do l=p1,mo_num
        !DIR$ LOOP COUNT AVG(4)
        do k=1,N_states
          mat_r(k,p1,l) = mat_r(k,p1,l) + tmp_rowij2(k,l)
        enddo
      enddo


    !! <phi|H|alpha>
      putj = p2
      do puti=1,mo_num
        if(lbanned(puti,ma)) cycle
        putj = p2
        if(.not. banned(puti,putj,1)) then
          hji = hji_cache(puti,1)
          if (hji /= 0.d0) then
            hji = hji * get_phase_bi(phasemask, mi, ma, hfix, pfix, puti, p1, N_int)
            !DIR$ LOOP COUNT AVG(4)
            do k=1,N_states
              tmp_rowji(k,puti) = tmp_rowji(k,puti) + hji * coefs(k,1)
            enddo
          endif
        end if

        putj = p1
        if(.not. banned(puti,putj,1)) then
          hji = hji_cache(puti,2)
          if (hji /= 0.d0) then
            hji = hji * get_phase_bi(phasemask, mi, ma, hfix, pfix, puti, p2, N_int)
            do k=1,N_states
              tmp_rowji2(k,puti) = tmp_rowji2(k,puti) + hji * coefs(k,1)
            enddo
          endif
        end if
      end do
      mat_l(:,:p2-1,p2) = mat_l(:,:p2-1,p2) + tmp_rowji(:,:p2-1)
      do l=p2,mo_num
        !DIR$ LOOP COUNT AVG(4)
        do k=1,N_states
          mat_l(k,p2,l) = mat_l(k,p2,l) + tmp_rowji(k,l)
        enddo
      enddo
      mat_l(:,:p1-1,p1) = mat_l(:,:p1-1,p1) + tmp_rowji2(:,:p1-1)
      do l=p1,mo_num
        !DIR$ LOOP COUNT AVG(4)
        do k=1,N_states
          mat_l(k,p1,l) = mat_l(k,p1,l) + tmp_rowji2(k,l)
        enddo
      enddo
    end if
  end if
  deallocate(lbanned,hij_cache, hji_cache)

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
!        call i_h_j_complex(gen, det, N_int, hij)
        call htilde_mu_mat_opt_bi_ortho_no_3e(det, gen, N_int, hij)
        call htilde_mu_mat_opt_bi_ortho_no_3e(gen, det, N_int, hji)
        !DIR$ LOOP COUNT AVG(4)
        do k=1,N_states
          ! take conjugate to get contribution to <alpha|H|psi> instead of <psi|H|alpha>
!          mat_r(k, p1, p2) = mat_r(k, p1, p2) + coefs(k,1) * dconjg(hij)
          mat_r(k, p1, p2) = mat_r(k, p1, p2) + coefs(k,2) * hij
          mat_l(k, p1, p2) = mat_l(k, p1, p2) + coefs(k,1) * hji
        enddo
      end do
    end do
end

