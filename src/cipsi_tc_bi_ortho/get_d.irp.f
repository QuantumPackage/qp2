
! ---

double precision function get_phase_bi(phasemask, s1, s2, h1, p1, h2, p2, Nint)

  use bitmasks
  implicit none

  integer,           intent(in) :: Nint
  integer,           intent(in) :: s1, s2, h1, h2, p1, p2
  integer(bit_kind), intent(in) :: phasemask(Nint,2)

  double precision, save        :: res(0:1) = (/1d0, -1d0/)

  integer                       :: np
  integer                       :: h1_int, h2_int
  integer                       :: p1_int, p2_int
  integer                       :: h1_bit, h2_bit
  integer                       :: p1_bit, p2_bit
  logical                       :: change

  h1_int = shiftr(h1-1,bit_kind_shift)+1
  h1_bit = h1 - shiftl(h1_int-1,bit_kind_shift)-1

  h2_int = shiftr(h2-1,bit_kind_shift)+1
  h2_bit = h2 - shiftl(h2_int-1,bit_kind_shift)-1

  p1_int = shiftr(p1-1,bit_kind_shift)+1
  p1_bit = p1 - shiftl(p1_int-1,bit_kind_shift)-1

  p2_int = shiftr(p2-1,bit_kind_shift)+1
  p2_bit = p2 - shiftl(p2_int-1,bit_kind_shift)-1

  ! Put the phasemask bits at position 0, and add them all
  h1_bit = int( shiftr( phasemask(h1_int,s1), h1_bit ) )
  p1_bit = int( shiftr( phasemask(p1_int,s1), p1_bit ) )
  h2_bit = int( shiftr( phasemask(h2_int,s2), h2_bit ) )
  p2_bit = int( shiftr( phasemask(p2_int,s2), p2_bit ) )

  np = h1_bit + p1_bit + h2_bit + p2_bit

  if(p1 < h1) np = np + 1
  if(p2 < h2) np = np + 1

  if(s1 == s2 .and. max(h1, p1) > min(h2, p2)) np = np + 1
  get_phase_bi = res(iand(np,1))

end function get_phase_bi

! ---

subroutine get_d3_htc(gen, bannedOrb, banned, mat_m, mat_p, mask, p, sp, rcoefs, lcoefs)

  use bitmasks
  implicit none

  integer(bit_kind), intent(in)   :: gen(N_int, 2), mask(N_int, 2)
  integer,           intent(in)   :: p(0:4,2), sp
  logical,           intent(in)   :: bannedOrb(mo_num, 2), banned(mo_num, mo_num,2)
  double precision,  intent(in)   :: rcoefs(N_states), lcoefs(N_states)
  double precision, intent(inout) :: mat_m(N_states, mo_num, mo_num), mat_p(N_states, mo_num, mo_num)

  integer(bit_kind)               :: det(N_int, 2)
  integer                         :: k, h1, h2, p1, p2, puti, putj
  double precision                :: i_h_alpha, alpha_h_i 
  logical                         :: ok

  if(sp == 3) then ! AB

    h1 = p(1,1)
    h2 = p(1,2)
    do p1 = 1, mo_num
      if(bannedOrb(p1, 1)) cycle
      do p2 = 1, mo_num
        if(bannedOrb(p2,2)) cycle
        if(banned(p1, p2, 1)) cycle ! rentable?

        call apply_particles(mask, 1, p1, 2, p2, det, ok, N_int)
        call htilde_mu_mat_bi_ortho_tot(gen, det, N_int, i_h_alpha)
        call htilde_mu_mat_bi_ortho_tot(det,gen, N_int, alpha_h_i)
!        call hji_hij_mu_mat_tot(gen, det, N_int,i_h_alpha , alpha_h_i)
        if( dabs(alpha_h_i) .gt. 0.d0) then
          !DIR$ LOOP COUNT AVG(4)
          do k = 1, N_states
            mat_p(k, p1, p2) = mat_p(k, p1, p2) + rcoefs(k) * alpha_h_i 
          enddo
        endif
        if( dabs(i_h_alpha) .gt. 0.d0) then
          !DIR$ LOOP COUNT AVG(4)
          do k = 1, N_states
            mat_m(k, p1, p2) = mat_m(k, p1, p2) + lcoefs(k) * i_h_alpha 
          enddo
        endif

      enddo
    enddo

  else ! AA BB

    p1 = p(1,sp)
    p2 = p(2,sp)
    do puti = 1, mo_num
      if(bannedOrb(puti, sp)) cycle
      do putj = puti+1, mo_num
        if(bannedOrb(putj, sp)) cycle
        if(banned(puti, putj, 1)) cycle ! rentable?

        call apply_particles(mask, sp, puti, sp, putj, det, ok, N_int)
!        call hji_hij_mu_mat_tot(gen, det, N_int, i_h_alpha, alpha_h_i)
        call htilde_mu_mat_bi_ortho_tot(gen, det, N_int, i_h_alpha)
        call htilde_mu_mat_bi_ortho_tot( det,gen, N_int, alpha_h_i)
        if( dabs(alpha_h_i) .gt. 0.d0) then
          !DIR$ LOOP COUNT AVG(4)
          do k = 1, N_states
            mat_p(k, puti, putj) = mat_p(k, puti, putj) + rcoefs(k) * alpha_h_i 
          enddo
        endif
        if( dabs(i_h_alpha) .gt. 0.d0) then
          !DIR$ LOOP COUNT AVG(4)
          do k = 1, N_states
            mat_m(k, puti, putj) = mat_m(k, puti, putj) + lcoefs(k) * i_h_alpha 
          enddo
        endif

      enddo
    enddo

  endif

end subroutine get_d3_htc

! ---

subroutine get_d3_h(gen, bannedOrb, banned, mat, mask, p, sp, coefs)

  use bitmasks
  implicit none

  integer(bit_kind), intent(in)   :: gen(N_int, 2), mask(N_int, 2)
  integer,           intent(in)   :: p(0:4,2), sp
  logical,           intent(in)   :: bannedOrb(mo_num, 2), banned(mo_num, mo_num,2)
  double precision,  intent(in)   :: coefs(N_states)
  double precision, intent(inout) :: mat(N_states, mo_num, mo_num)

  integer(bit_kind)               :: det(N_int, 2)
  integer                         :: k, h1, h2, p1, p2, puti, putj
  double precision                :: hij
  logical                         :: ok

  if(sp == 3) then ! AB

    h1 = p(1,1)
    h2 = p(1,2)
    do p1 = 1, mo_num
      if(bannedOrb(p1, 1)) cycle
      do p2 = 1, mo_num
        if(bannedOrb(p2,2)) cycle
        if(banned(p1, p2, 1)) cycle ! rentable?

        call apply_particles(mask, 1, p1, 2, p2, det, ok, N_int)
        call i_h_j(gen, det, N_int, hij)
        if (hij == 0.d0) cycle
        !DIR$ LOOP COUNT AVG(4)
        do k = 1, N_states
          mat(k, p1, p2) = mat(k, p1, p2) + coefs(k) * hij
        enddo

      enddo
    enddo

  else ! AA BB

    p1 = p(1,sp)
    p2 = p(2,sp)
    do puti = 1, mo_num
      if(bannedOrb(puti, sp)) cycle
      do putj = puti+1, mo_num
        if(bannedOrb(putj, sp)) cycle
        if(banned(puti, putj, 1)) cycle ! rentable?

        call apply_particles(mask, sp, puti, sp, putj, det, ok, N_int)
        call i_h_j(gen, det, N_int, hij)
        !DIR$ LOOP COUNT AVG(4)
        do k = 1, N_states
          mat(k, puti, putj) = mat(k, puti, putj) + coefs(k) * hij
        enddo

      enddo
    enddo

  endif

end subroutine get_d3_h

! ---

subroutine get_d2(gen, phasemask, bannedOrb, banned, mat_p, mat_m, mask, h, p, sp, coefs)

  use bitmasks
  implicit none

  integer(bit_kind), intent(in) :: mask(N_int, 2), gen(N_int, 2)
  integer(bit_kind), intent(in) :: phasemask(N_int,2)
  logical, intent(in) :: bannedOrb(mo_num, 2), banned(mo_num, mo_num,2)
  double precision, intent(in) :: coefs(N_states,2)
  double precision, intent(inout) :: mat_p(N_states, mo_num, mo_num), mat_m(N_states, mo_num, mo_num)
  integer, intent(in) :: h(0:2,2), p(0:4,2), sp

  double precision, external :: get_phase_bi

  integer :: i, j, k, tip, ma, mi, puti, putj
  integer :: h1, h2, p1, p2, i1, i2
  double precision :: hij, hji, phase

  integer, parameter:: turn2d(2,3,4) = reshape((/0,0, 0,0, 0,0,  3,4, 0,0, 0,0,  2,4, 1,4, 0,0,  2,3, 1,3, 1,2 /), (/2,3,4/))
  integer, parameter :: turn2(2) = (/2, 1/)
  integer, parameter :: turn3(2,3) = reshape((/2,3,  1,3, 1,2/), (/2,3/))

  integer :: bant
  bant = 1

  tip = p(0,1) * p(0,2)

  ma = sp
  print*,'in get d2'
  if(p(0,1) > p(0,2)) ma = 1
  if(p(0,1) < p(0,2)) ma = 2
  mi = mod(ma, 2) + 1

  if(sp == 3) then
    print*,'in sp == 3'
    if(ma == 2) bant = 2

    if(tip == 3) then
      puti = p(1, mi)
      if(bannedOrb(puti, mi)) return
      h1 = h(1, ma)
      h2 = h(2, ma)

      do i = 1, 3
        putj = p(i, ma)
        if(banned(putj,puti,bant)) cycle
        i1 = turn3(1,i)
        i2 = turn3(2,i)
        p1 = p(i1, ma)
        p2 = p(i2, ma)

        ! <p1 p2|1/r12|h1 h2> --> <p1 p2| w_ee^h + t^nh | h1 h2> --> < p2 p1 | H^tilde| h1 h2 >
        ! 
        !                      <p1 p2 | h1 h2>        -            <p2 p1 | h1 h2 >
        ! < p2 p1 | H^tilde^dag| h1 h2 > = < h1 h2 | w_ee^h + t^nh | p1 p2 >
        hji = mo_bi_ortho_tc_two_e(p1, p2, h1, h2) - mo_bi_ortho_tc_two_e(p2, p1, h1, h2)
        hij = mo_bi_ortho_tc_two_e(h1, h2, p1, p2) - mo_bi_ortho_tc_two_e(h2, h1, p1, p2)
        if (hij == 0.d0) cycle

        hij = hij * get_phase_bi(phasemask, ma, ma, h1, p1, h2, p2, N_int)
        hji = hji * get_phase_bi(phasemask, ma, ma, h1, p1, h2, p2, N_int)

        if(ma == 1) then
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat_p(k, putj, puti) = mat_p(k, putj, puti) + coefs(k,1) * hij
            mat_m(k, putj, puti) = mat_m(k, putj, puti) + coefs(k,2) * hji
          enddo
        else
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat_p(k, puti, putj) = mat_p(k, puti, putj) + coefs(k,1) * hij
            mat_m(k, puti, putj) = mat_m(k, puti, putj) + coefs(k,2) * hji
          enddo
        end if
      end do
    else
      h1 = h(1,1)
      h2 = h(1,2)
      do j = 1,2
        putj = p(j, 2)
        if(bannedOrb(putj, 2)) cycle
        p2 = p(turn2(j), 2)
        do i = 1,2
          puti = p(i, 1)

          if(banned(puti,putj,bant) .or. bannedOrb(puti,1)) cycle
          p1 = p(turn2(i), 1)

          hji = mo_bi_ortho_tc_two_e(p1, p2, h1, h2)
          hij = mo_bi_ortho_tc_two_e(h1, h2, p1, p2)
          if (hij /= 0.d0) then
            hij = hij * get_phase_bi(phasemask, 1, 2, h1, p1, h2, p2, N_int)
            hji = hji * get_phase_bi(phasemask, 1, 2, h1, p1, h2, p2, N_int)
            !DIR$ LOOP COUNT AVG(4)
            do k=1,N_states
              mat_p(k, puti, putj) = mat_p(k, puti, putj) + coefs(k,1) * hij
              mat_m(k, puti, putj) = mat_m(k, puti, putj) + coefs(k,2) * hji
            enddo
          endif
        end do
      end do
    end if

  else
    print*,'NOT in sp == 3'
    if(tip == 0) then
      h1 = h(1, ma)
      h2 = h(2, ma)
      do i=1,3
        puti = p(i, ma)
        if(bannedOrb(puti,ma)) cycle
        do j=i+1,4
          putj = p(j, ma)
          if(bannedOrb(putj,ma)) cycle
          if(banned(puti,putj,1)) cycle

          i1 = turn2d(1, i, j)
          i2 = turn2d(2, i, j)
          p1 = p(i1, ma)
          p2 = p(i2, ma)
          hij = mo_bi_ortho_tc_two_e(p1, p2, h1, h2) - mo_bi_ortho_tc_two_e(p2,p1, h1, h2)
          hji = mo_bi_ortho_tc_two_e(h1, h2, p1, p2) - mo_bi_ortho_tc_two_e(h2,h1, p1, p2)
          if (hij == 0.d0) cycle

          hij = hij * get_phase_bi(phasemask, ma, ma, h1, p1, h2, p2, N_int)
          hji = hji * get_phase_bi(phasemask, ma, ma, h1, p1, h2, p2, N_int)
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat_p(k, puti, putj) = mat_p(k, puti, putj) +coefs(k,1) * hij
            mat_m(k, puti, putj) = mat_m(k, puti, putj) +coefs(k,2) * hji
          enddo
        end do
      end do
    else if(tip == 3) then
      h1 = h(1, mi)
      h2 = h(1, ma)
      p1 = p(1, mi)
      do i=1,3
        puti = p(turn3(1,i), ma)
        if(bannedOrb(puti,ma)) cycle
        putj = p(turn3(2,i), ma)
        if(bannedOrb(putj,ma)) cycle
        if(banned(puti,putj,1)) cycle
        p2 = p(i, ma)

        hij = mo_bi_ortho_tc_two_e(p1, p2, h1, h2)
        hji = mo_bi_ortho_tc_two_e(h1, h2, p1, p2)
        if (hij == 0.d0) cycle

        hij = hij * get_phase_bi(phasemask, mi, ma, h1, p1, h2, p2, N_int)
        hji = hji * get_phase_bi(phasemask, mi, ma, h1, p1, h2, p2, N_int)
        if (puti < putj) then
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat_p(k, puti, putj) = mat_p(k, puti, putj) + coefs(k,1) * hij
            mat_m(k, puti, putj) = mat_m(k, puti, putj) + coefs(k,2) * hji
          enddo
        else
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat_p(k, putj, puti) = mat_p(k, putj, puti) + coefs(k,1) * hij
            mat_m(k, putj, puti) = mat_m(k, putj, puti) + coefs(k,2) * hji
          enddo
        endif
      end do
    else ! tip == 4
      puti = p(1, sp)
      putj = p(2, sp)
      if(.not. banned(puti,putj,1)) then
        p1 = p(1, mi)
        p2 = p(2, mi)
        h1 = h(1, mi)
        h2 = h(2, mi)
        hij = (mo_bi_ortho_tc_two_e(p1, p2, h1, h2) - mo_bi_ortho_tc_two_e(p2,p1, h1, h2))
        hji = (mo_bi_ortho_tc_two_e(h1, h2, p1, p2) - mo_bi_ortho_tc_two_e(h2,h1, p1, p2))
        if (hij /= 0.d0) then
          hij = hij * get_phase_bi(phasemask, mi, mi, h1, p1, h2, p2, N_int)
          hji = hji * get_phase_bi(phasemask, mi, mi, h1, p1, h2, p2, N_int)
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat_p(k, puti, putj) = mat_p(k, puti, putj) + coefs(k,1) * hij
            mat_m(k, puti, putj) = mat_m(k, puti, putj) + coefs(k,2) * hji
          enddo
        end if
      end if
    end if
  end if

end subroutine get_d2

! ---

subroutine get_d1(gen, phasemask, bannedOrb, banned, mat_p, mat_m, mask, h, p, sp, coefs)

  use bitmasks
  implicit none

  integer(bit_kind), intent(in)  :: mask(N_int, 2), gen(N_int, 2)
  integer(bit_kind), intent(in)  :: phasemask(N_int,2)
  logical, intent(in)            :: bannedOrb(mo_num, 2), banned(mo_num, mo_num,2)
  integer(bit_kind)              :: det(N_int, 2)
  double precision, intent(in)   :: coefs(N_states,2)
  double precision, intent(inout) :: mat_p(N_states, mo_num, mo_num), mat_m(N_states, mo_num, mo_num)
  integer, intent(in)            :: h(0:2,2), p(0:4,2), sp
  double precision, external     :: get_phase_bi
  logical                        :: ok

  logical, allocatable           :: lbanned(:,:)
  integer                        :: puti, putj, ma, mi, s1, s2, i, i1, i2, j
  integer                        :: hfix, pfix, h1, h2, p1, p2, ib, k, l
  integer :: mm

  integer, parameter             :: turn2(2) = (/2,1/)
  integer, parameter             :: turn3(2,3) = reshape((/2,3,  1,3, 1,2/), (/2,3/))

  integer                        :: bant
  double precision, allocatable  :: hij_cache(:,:)
  double precision               :: hij, tmp_row_ij(N_states, mo_num), tmp_row_ij2(N_states, mo_num)
  double precision, allocatable  :: hji_cache(:,:)
  double precision               :: hji, tmp_row_ji(N_states, mo_num), tmp_row_ji2(N_states, mo_num)

  PROVIDE mo_integrals_map N_int

  allocate (lbanned(mo_num, 2))
  allocate (hij_cache(mo_num,2))
  allocate (hji_cache(mo_num,2))
  lbanned = bannedOrb
  print*,'in get d1'
  call debug_det(gen, N_int)

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
    print*,'in sp == 3'
    !move MA
    if(ma == 2) bant = 2
    puti = p(1,mi)
    hfix = h(1,ma)
    p1 = p(1,ma)
    p2 = p(2,ma)
    print*,puti, hfix,p1,p2
    if(.not. bannedOrb(puti, mi)) then
!      print*,'not banned'
      do mm = 1, mo_num
       hji_cache(mm,1) = mo_bi_ortho_tc_two_e(p1,p2,mm,hfix)
       hji_cache(mm,2) = mo_bi_ortho_tc_two_e(p2,p1,mm,hfix)
       hij_cache(mm,1) = mo_bi_ortho_tc_two_e(mm,hfix,p1,p2)
       hij_cache(mm,2) = mo_bi_ortho_tc_two_e(mm,hfix,p2,p1)
      enddo
!      call get_mo_bi_ortho_tc_two_es(hfix,p1,p2,mo_num,hij_cache(1,1),mo_integrals_map)
!      call get_mo_bi_ortho_tc_two_es(hfix,p2,p1,mo_num,hij_cache(1,2),mo_integrals_map)
      tmp_row_ij = 0d0
      tmp_row_ji = 0d0
      do putj=1, hfix-1
        if(lbanned(putj, ma)) cycle
        if(banned(putj, puti,bant)) cycle
        hij = hij_cache(putj,1) - hij_cache(putj,2)
        if (hij /= 0.d0) then
          hij = hij * get_phase_bi(phasemask, ma, ma, putj, p1, hfix, p2, N_int)
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            tmp_row_ij(k,putj) = tmp_row_ij(k,putj) + hij * coefs(k,2)
          enddo
        endif
        hji = hji_cache(putj,1) - hji_cache(putj,2)
        if (hji /= 0.d0) then
          hji = hji * get_phase_bi(phasemask, ma, ma, putj, p1, hfix, p2, N_int)
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            tmp_row_ji(k,putj) = tmp_row_ji(k,putj) + hji * coefs(k,1)
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
            tmp_row_ij(k,putj) = tmp_row_ij(k,putj) + hij * coefs(k,2)
          enddo
        endif
        hji = hji_cache(putj,2) - hji_cache(putj,1)
        if (hji /= 0.d0) then
          hji = hji * get_phase_bi(phasemask, ma, ma, hfix, p1, putj, p2, N_int)
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            tmp_row_ji(k,putj) = tmp_row_ji(k,putj) + hji * coefs(k,1)
          enddo
        endif
      end do

      if(ma == 1) then
        mat_p(1:N_states,1:mo_num,puti) = mat_p(1:N_states,1:mo_num,puti) + tmp_row_ij(1:N_states,1:mo_num)
        mat_m(1:N_states,1:mo_num,puti) = mat_m(1:N_states,1:mo_num,puti) + tmp_row_ji(1:N_states,1:mo_num)
      else
        do l=1,mo_num
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat_p(k,puti,l) = mat_p(k,puti,l) + tmp_row_ij(k,l)
            mat_m(k,puti,l) = mat_m(k,puti,l) + tmp_row_ji(k,l)
          enddo
        enddo
      end if
    end if

    !MOVE MI
    pfix = p(1,mi)
    tmp_row_ij = 0d0
    tmp_row_ij2 = 0d0
    tmp_row_ji = 0d0
    tmp_row_ji2 = 0d0
!    call get_mo_bi_ortho_tc_two_es(hfix,pfix,p1,mo_num,hij_cache(1,1),mo_integrals_map)
!    call get_mo_bi_ortho_tc_two_es(hfix,pfix,p2,mo_num,hij_cache(1,2),mo_integrals_map)
    do mm = 1, mo_num
     hij_cache(mm,1) = mo_bi_ortho_tc_two_e(mm,hfix,pfix,p1)
     hij_cache(mm,2) = mo_bi_ortho_tc_two_e(mm,hfix,pfix,p2)
     hji_cache(mm,1) = mo_bi_ortho_tc_two_e(pfix,p1,mm,hfix)
     hji_cache(mm,2) = mo_bi_ortho_tc_two_e(pfix,p2,mm,hfix)
    enddo
    putj = p1
    do puti = 1, mo_num !HOT

      if(lbanned(puti,mi)) cycle
      !p1 fixed
      putj = p1
      if(.not. banned(putj,puti,bant)) then
        hij = hij_cache(puti,2)
        if (hij /= 0.d0) then
          hij = hij * get_phase_bi(phasemask, ma, mi, hfix, p2, puti, pfix, N_int)
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            tmp_row_ij(k,puti) = tmp_row_ij(k,puti) + hij * coefs(k,2)
          enddo
        endif
        hji = hji_cache(puti,2)
        if (hji /= 0.d0) then
          hji = hji * get_phase_bi(phasemask, ma, mi, hfix, p2, puti, pfix, N_int)
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            tmp_row_ji(k,puti) = tmp_row_ji(k,puti) + hji * coefs(k,1)
          enddo
        endif
      endif

      putj = p2
      if(.not. banned(putj,puti,bant)) then
        hij = hij_cache(puti,1)
        if (hij /= 0.d0) then
          hij = hij * get_phase_bi(phasemask, ma, mi, hfix, p1, puti, pfix, N_int)
          do k=1,N_states
            tmp_row_ij2(k,puti) = tmp_row_ij2(k,puti) + hij * coefs(k,2)
          enddo
        endif
        hji = hji_cache(puti,1)
        if (hji /= 0.d0) then
          hji = hji * get_phase_bi(phasemask, ma, mi, hfix, p1, puti, pfix, N_int)
          do k=1,N_states
            tmp_row_ji2(k,puti) = tmp_row_ji2(k,puti) + hji * coefs(k,1)
          enddo
        endif
      endif

    enddo

    if(mi == 1) then
      mat_p(:,:,p1) = mat_p(:,:,p1) + tmp_row_ij(:,:)
      mat_p(:,:,p2) = mat_p(:,:,p2) + tmp_row_ij2(:,:)
      mat_m(:,:,p1) = mat_m(:,:,p1) + tmp_row_ji(:,:)
      mat_m(:,:,p2) = mat_m(:,:,p2) + tmp_row_ji2(:,:)
    else
      do l=1,mo_num
        !DIR$ LOOP COUNT AVG(4)
        do k=1,N_states
          mat_p(k,p1,l) = mat_p(k,p1,l) + tmp_row_ij(k,l)
          mat_p(k,p2,l) = mat_p(k,p2,l) + tmp_row_ij2(k,l)
          mat_m(k,p1,l) = mat_m(k,p1,l) + tmp_row_ji(k,l)
          mat_m(k,p2,l) = mat_m(k,p2,l) + tmp_row_ji2(k,l)
        enddo
      enddo
    end if

  else  ! sp /= 3
    print*,'not in sp == 3'

    if(p(0,ma) == 3) then
      do i=1,3
        hfix = h(1,ma)
        puti = p(i, ma)
        p1 = p(turn3(1,i), ma)
        p2 = p(turn3(2,i), ma)
!        call get_mo_bi_ortho_tc_two_es(hfix,p1,p2,mo_num,hij_cache(1,1),mo_integrals_map)
!        call get_mo_bi_ortho_tc_two_es(hfix,p2,p1,mo_num,hij_cache(1,2),mo_integrals_map)
        do mm = 1, mo_num
         hji_cache(mm,1) = mo_bi_ortho_tc_two_e(p1,p2,mm,hfix)
         hji_cache(mm,2) = mo_bi_ortho_tc_two_e(p2,p1,mm,hfix)
         hij_cache(mm,1) = mo_bi_ortho_tc_two_e(mm,hfix,p1,p2)
         hij_cache(mm,2) = mo_bi_ortho_tc_two_e(mm,hfix,p2,p1)
        enddo
        tmp_row_ij = 0d0
        tmp_row_ji = 0d0
        do putj=1,hfix-1
          if(banned(putj,puti,1)) cycle
          if(lbanned(putj,ma)) cycle
          hij = hij_cache(putj,1) - hij_cache(putj,2)
          if (hij /= 0.d0) then
            hij = hij * get_phase_bi(phasemask, ma, ma, putj, p1, hfix, p2, N_int)
            tmp_row_ij(:,putj) = tmp_row_ij(:,putj) + hij * coefs(:,1)
          endif
          hji = hji_cache(putj,1) - hji_cache(putj,2)
          if (hji /= 0.d0) then
            hji = hji * get_phase_bi(phasemask, ma, ma, putj, p1, hfix, p2, N_int)
            tmp_row_ji(:,putj) = tmp_row_ji(:,putj) + hji * coefs(:,2)
          endif
        end do
        do putj=hfix+1,mo_num
          if(banned(putj,puti,1)) cycle
          if(lbanned(putj,ma)) cycle
          hij = hij_cache(putj,2) - hij_cache(putj,1)
          if (hij /= 0.d0) then
            hij = hij * get_phase_bi(phasemask, ma, ma, hfix, p1, putj, p2, N_int)
            tmp_row_ij(:,putj) = tmp_row_ij(:,putj) + hij * coefs(:,1)
          endif
          hji = hji_cache(putj,2) - hji_cache(putj,1)
          if (hji /= 0.d0) then
            hji = hji * get_phase_bi(phasemask, ma, ma, hfix, p1, putj, p2, N_int)
            tmp_row_ji(:,putj) = tmp_row_ji(:,putj) + hji * coefs(:,2)
          endif
        end do

        mat_p(:, :puti-1, puti) = mat_p(:, :puti-1, puti) + tmp_row_ij(:,:puti-1)
        mat_m(:, :puti-1, puti) = mat_m(:, :puti-1, puti) + tmp_row_ji(:,:puti-1)
        do l=puti,mo_num
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat_p(k, puti, l) = mat_p(k, puti,l) + tmp_row_ij(k,l)
            mat_m(k, puti, l) = mat_m(k, puti,l) + tmp_row_ji(k,l)
          enddo
        enddo
      end do
    else
      hfix = h(1,mi)
      pfix = p(1,mi)
      p1 = p(1,ma)
      p2 = p(2,ma)
      tmp_row_ij = 0d0
      tmp_row_ij2 = 0d0
      tmp_row_ji = 0d0
      tmp_row_ji2 = 0d0
!      call get_mo_bi_ortho_tc_two_es(hfix,p1,pfix,mo_num,hij_cache(1,1),mo_integrals_map)
!      call get_mo_bi_ortho_tc_two_es(hfix,p2,pfix,mo_num,hij_cache(1,2),mo_integrals_map)
      do mm = 1, mo_num
       hji_cache(mm,1) = mo_bi_ortho_tc_two_e(p1,pfix,mm,hfix)
       hji_cache(mm,2) = mo_bi_ortho_tc_two_e(p2,pfix,mm,hfix)
       hij_cache(mm,1) = mo_bi_ortho_tc_two_e(mm,hfix,p1,pfix)
       hij_cache(mm,2) = mo_bi_ortho_tc_two_e(mm,hfix,p2,pfix)
      enddo
      putj = p2
      do puti=1,mo_num
        if(lbanned(puti,ma)) cycle
        putj = p2
        if(.not. banned(puti,putj,1)) then
          hij = hij_cache(puti,1)
          if (hij /= 0.d0) then
            hij = hij * get_phase_bi(phasemask, mi, ma, hfix, pfix, puti, p1, N_int)
            !DIR$ LOOP COUNT AVG(4)
            do k=1,N_states
              tmp_row_ij(k,puti) = tmp_row_ij(k,puti) + hij * coefs(k,1)
            enddo
          endif
          hji = hji_cache(puti,1)
          if (hji /= 0.d0) then
            hji = hji * get_phase_bi(phasemask, mi, ma, hfix, pfix, puti, p1, N_int)
            !DIR$ LOOP COUNT AVG(4)
            do k=1,N_states
              tmp_row_ji(k,puti) = tmp_row_ji(k,puti) + hji * coefs(k,2)
            enddo
          endif
        end if

        putj = p1
        if(.not. banned(puti,putj,1)) then
          hij = hij_cache(puti,2)
          if (hij /= 0.d0) then
            hij = hij * get_phase_bi(phasemask, mi, ma, hfix, pfix, puti, p2, N_int)
            do k=1,N_states
              tmp_row_ij2(k,puti) = tmp_row_ij2(k,puti) + hij * coefs(k,1)
            enddo
          endif
          hji = hji_cache(puti,2)
          if (hji /= 0.d0) then
            hji = hji * get_phase_bi(phasemask, mi, ma, hfix, pfix, puti, p2, N_int)
            do k=1,N_states
              tmp_row_ji2(k,puti) = tmp_row_ji2(k,puti) + hji * coefs(k,2)
            enddo
          endif
        end if
      end do
      mat_p(:,:p2-1,p2) = mat_p(:,:p2-1,p2) + tmp_row_ij(:,:p2-1)
      mat_m(:,:p2-1,p2) = mat_m(:,:p2-1,p2) + tmp_row_ji(:,:p2-1)
      do l=p2,mo_num
        !DIR$ LOOP COUNT AVG(4)
        do k=1,N_states
          mat_p(k,p2,l) = mat_p(k,p2,l) + tmp_row_ij(k,l)
          mat_m(k,p2,l) = mat_m(k,p2,l) + tmp_row_ji(k,l)
        enddo
      enddo
      mat_p(:,:p1-1,p1) = mat_p(:,:p1-1,p1) + tmp_row_ij2(:,:p1-1)
      mat_m(:,:p1-1,p1) = mat_m(:,:p1-1,p1) + tmp_row_ji2(:,:p1-1)
      do l=p1,mo_num
        !DIR$ LOOP COUNT AVG(4)
        do k=1,N_states
          mat_p(k,p1,l) = mat_p(k,p1,l) + tmp_row_ij2(k,l)
          mat_m(k,p1,l) = mat_m(k,p1,l) + tmp_row_ji2(k,l)
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

    do i1 = 1, p(0,s1)
      ib = 1
      if(s1 == s2) ib = i1+1
      do i2 = ib, p(0,s2)
        p1 = p(i1,s1)
        p2 = p(i2,s2)
        if(bannedOrb(p1, s1) .or. bannedOrb(p2, s2) .or. banned(p1, p2, 1)) cycle
        call apply_particles(mask, s1, p1, s2, p2, det, ok, N_int)
!        call i_h_j(gen, det, N_int, hij)
        !!!! GUESS ON THE ORDER OF DETS 
        print*,'compute hij'
!        hij = 0.d0
!        hji = 0.d0
        call htilde_mu_mat_opt_bi_ortho_no_3e(gen,det,N_int, hji)
        call htilde_mu_mat_opt_bi_ortho_no_3e(det,gen,N_int, hij)
        !DIR$ LOOP COUNT AVG(4)
        do k = 1, N_states
          mat_p(k, p1, p2) = mat_p(k, p1, p2) + coefs(k,1) * hij
          mat_m(k, p1, p2) = mat_m(k, p1, p2) + coefs(k,2) * hji
        enddo
      enddo
    enddo

end subroutine get_d1

! ---

subroutine get_d0(gen, phasemask, bannedOrb, banned, mat_p, mat_m, mask, h, p, sp, coefs)

  use bitmasks
  implicit none

  integer(bit_kind), intent(in) :: gen(N_int, 2), mask(N_int, 2)
  integer(bit_kind), intent(in) :: phasemask(N_int,2)
  logical, intent(in) :: bannedOrb(mo_num, 2), banned(mo_num, mo_num,2)
  integer(bit_kind) :: det(N_int, 2)
  double precision, intent(in) :: coefs(N_states,2)
  double precision, intent(inout) :: mat_m(N_states, mo_num, mo_num)
  double precision, intent(inout) :: mat_p(N_states, mo_num, mo_num)
  integer, intent(in) :: h(0:2,2), p(0:4,2), sp

  integer :: i, j, k, s, h1, h2, p1, p2, puti, putj, mm
  double precision :: hij, phase, hji
  double precision, external :: get_phase_bi
  logical :: ok

  integer, parameter :: bant=1
  double precision, allocatable :: hij_cache1(:), hij_cache2(:)
  allocate (hij_cache1(mo_num),hij_cache2(mo_num))
  double precision, allocatable :: hji_cache1(:), hji_cache2(:)
  allocate (hji_cache1(mo_num),hji_cache2(mo_num))

  print*,'in get d0'
!  call debug_det(gen, N_int)

  if(sp == 3) then ! AB
    h1 = p(1,1)
    h2 = p(1,2)
!   print*,'in AB'
    do p1=1, mo_num
      if(bannedOrb(p1, 1)) cycle
!      call get_mo_bi_ortho_tc_two_es(p1,h2,h1,mo_num,hij_cache1,mo_integrals_map)
      do mm =1, mo_num
       hji_cache1(mm) = mo_bi_ortho_tc_two_e(mm,p1,h2,h1)
       hji_cache1(mm) = mo_bi_ortho_tc_two_e(h2,h1,mm,p1)
      enddo
      do p2=1, mo_num
        if(bannedOrb(p2,2)) cycle
        if(banned(p1, p2, bant)) cycle ! rentable?
        if(p1 == h1 .or. p2 == h2) then
!          print*,'in p1 == h1 or p2 == h2'
          call apply_particles(mask, 1,p1,2,p2, det, ok, N_int)
!          call i_h_j(gen, det, N_int, hij)
           !!! GUESS ON THE ORDER 
          call htilde_mu_mat_opt_bi_ortho_no_3e(det,gen,N_int, hji)
          call htilde_mu_mat_opt_bi_ortho_no_3e(gen,det,N_int, hij)
        else
!          print*,'ELSE '
          phase = get_phase_bi(phasemask, 1, 2, h1, p1, h2, p2, N_int)
          hij = hij_cache1(p2) * phase
          hji = hji_cache1(p2) * phase
        end if
        if (hij == 0.d0) cycle
        !DIR$ LOOP COUNT AVG(4)
        do k=1,N_states
          mat_p(k, p1, p2) = mat_p(k, p1, p2) + coefs(k,1) * hij  ! HOTSPOT
          mat_m(k, p1, p2) = mat_m(k, p1, p2) + coefs(k,2) * hji  ! HOTSPOT
        enddo
      end do
    end do

  else ! AA BB
!    print*, 'in AA BB' 
    p1 = p(1,sp)
    p2 = p(2,sp)
    do puti=1, mo_num
      if(bannedOrb(puti, sp)) cycle
      do mm = 1, mo_num
       hij_cache1(mm) = mo_bi_ortho_tc_two_e(p2,p1,mm,puti)
       hij_cache2(mm) = mo_bi_ortho_tc_two_e(p1,p2,mm,puti)
       hji_cache1(mm) = mo_bi_ortho_tc_two_e(mm,puti,p2,p1)
       hji_cache2(mm) = mo_bi_ortho_tc_two_e(mm,puti,p1,p2)
      enddo
!      call get_mo_bi_ortho_tc_two_es(puti,p2,p1,mo_num,hij_cache1,mo_integrals_map)
!      call get_mo_bi_ortho_tc_two_es(puti,p1,p2,mo_num,hij_cache2,mo_integrals_map)
      do putj=puti+1, mo_num
        if(bannedOrb(putj, sp)) cycle
        if(banned(puti, putj, bant)) cycle ! rentable?
        if(puti == p1 .or. putj == p2 .or. puti == p2 .or. putj == p1) then
          call apply_particles(mask, sp,puti,sp,putj, det, ok, N_int)
!          call i_h_j(gen, det, N_int, hij)
          !!! GUESS 
          call htilde_mu_mat_opt_bi_ortho_no_3e(gen,det,N_int, hij)
          call htilde_mu_mat_opt_bi_ortho_no_3e(det,gen,N_int, hji)
          if (hij == 0.d0.or.hji == 0.d0) cycle
        else
          hji = (mo_bi_ortho_tc_two_e(p1, p2, puti, putj) -  mo_bi_ortho_tc_two_e(p2, p1, puti, putj))
          hij = (mo_bi_ortho_tc_two_e(puti, putj, p1, p2) -  mo_bi_ortho_tc_two_e(puti, putj, p2, p1))
          if (hij == 0.d0.or.hji==0.d0) cycle
          hij = hij * get_phase_bi(phasemask, sp, sp, puti, p1 , putj, p2, N_int)
          hji = hji * get_phase_bi(phasemask, sp, sp, puti, p1 , putj, p2, N_int)
        end if
        !DIR$ LOOP COUNT AVG(4)
        do k=1,N_states
          mat_p(k, puti, putj) = mat_p(k, puti, putj) + coefs(k,1) * hij
          mat_m(k, puti, putj) = mat_m(k, puti, putj) + coefs(k,2) * hji
        enddo
      end do
    end do
  end if

!  deallocate(hij_cache1,hij_cache2)
!  deallocate(hji_cache1,hji_cache2)

end subroutine get_d0

! ---

! ___________________________________________________________________________________________________________________________________________________
! ___________________________________________________________________________________________________________________________________________________

!subroutine get_pm2(gen, phasemask, bannedOrb, banned, mat_p, mat_m, mask, h, p, sp, coefs)
!
!  use bitmasks
!
!  implicit none
!
!  integer,           intent(in)   :: h(0:2,2), p(0:4,2), sp
!  integer(bit_kind), intent(in)   :: mask(N_int, 2), gen(N_int, 2), phasemask(N_int,2)
!  logical,           intent(in)   :: bannedOrb(mo_num, 2), banned(mo_num, mo_num,2)
!  double precision,  intent(in)   :: coefs(N_states)
!  double precision, intent(inout) :: mat_p(N_states, mo_num, mo_num), mat_m(N_states, mo_num, mo_num)
!
!  integer, parameter              :: turn2d(2,3,4) = reshape((/0,0, 0,0, 0,0,  3,4, 0,0, 0,0,  2,4, 1,4, 0,0,  2,3, 1,3, 1,2 /), (/2,3,4/))
!  integer, parameter              :: turn2(2) = (/2, 1/)
!  integer, parameter              :: turn3(2,3) = reshape((/2,3,  1,3, 1,2/), (/2,3/))
!
!  integer                         :: i, j, k, tip, ma, mi, puti, putj
!  integer                         :: h1, h2, p1, p2, i1, i2
!  integer                         :: bant
!  double precision                :: hij_p, hij_m, phase
!
!  double precision, external      :: get_phase_bi
!  double precision, external      :: get_mo_bi_ortho_tc_two_e_tc_int, get_mo_bi_ortho_tc_two_e_tcdag_int 
!
!  PROVIDE mo_integrals_tc_int_map mo_integrals_tcdag_int_map
!
!  bant = 1
!
!  tip = p(0,1) * p(0,2)
!
!  ma = sp
!  if(p(0,1) > p(0,2)) ma = 1
!  if(p(0,1) < p(0,2)) ma = 2
!  mi = mod(ma, 2) + 1
!
!  if(sp == 3) then
!    if(ma == 2) bant = 2
!    if(tip == 3) then
!      puti = p(1, mi)
!      if(bannedOrb(puti, mi)) return
!      h1 = h(1, ma)
!      h2 = h(2, ma)
!
!      do i = 1, 3
!        putj = p(i, ma)
!        if(banned(putj,puti,bant)) cycle
!        i1 = turn3(1,i)
!        i2 = turn3(2,i)
!        p1 = p(i1, ma)
!        p2 = p(i2, ma)
!
!        hij_p = get_mo_bi_ortho_tc_two_e_tc_int   (p1, p2, h1, h2, mo_integrals_tc_int_map   ) &
!              - get_mo_bi_ortho_tc_two_e_tc_int   (p2, p1, h1, h2, mo_integrals_tc_int_map   )
!        hij_m = get_mo_bi_ortho_tc_two_e_tcdag_int(p1, p2, h1, h2, mo_integrals_tcdag_int_map) &
!              - get_mo_bi_ortho_tc_two_e_tcdag_int(p2, p1, h1, h2, mo_integrals_tcdag_int_map)
!
!        if( (hij_p.eq.0.d0) .and. (hij_m.eq.0.d0) ) cycle
!
!        hij_p = hij_p * get_phase_bi(phasemask, ma, ma, h1, p1, h2, p2, N_int)
!        hij_m = hij_m * get_phase_bi(phasemask, ma, ma, h1, p1, h2, p2, N_int)
!
!        if(ma == 1) then
!          !DIR$ LOOP COUNT AVG(4)
!          do k = 1, N_states
!            mat_p(k, putj, puti) = mat_p(k, putj, puti) + coefs(k) * hij_p
!            mat_m(k, putj, puti) = mat_m(k, putj, puti) + coefs(k) * hij_m
!          enddo
!        else
!          !DIR$ LOOP COUNT AVG(4)
!          do k = 1, N_states
!            mat_p(k, puti, putj) = mat_p(k, puti, putj) + coefs(k) * hij_p
!            mat_m(k, puti, putj) = mat_m(k, puti, putj) + coefs(k) * hij_m
!          enddo
!        end if
!      end do
!
!    else
!
!      h1 = h(1,1)
!      h2 = h(1,2)
!      do j = 1,2
!        putj = p(j, 2)
!        if(bannedOrb(putj, 2)) cycle
!        p2 = p(turn2(j), 2)
!        do i = 1,2
!          puti = p(i, 1)
!
!          if(banned(puti,putj,bant) .or. bannedOrb(puti,1)) cycle
!          p1 = p(turn2(i), 1)
!
!          hij_p = get_mo_bi_ortho_tc_two_e_tc_int   (p1, p2, h1, h2, mo_integrals_tc_int_map   )
!          hij_m = get_mo_bi_ortho_tc_two_e_tcdag_int(p1, p2, h1, h2, mo_integrals_tcdag_int_map)
!
!          if( (hij_p.ne.0.d0) .and. (hij_m.ne.0.d0) ) then
!            hij_p = hij_p * get_phase_bi(phasemask, 1, 2, h1, p1, h2, p2, N_int)
!            hij_m = hij_m * get_phase_bi(phasemask, 1, 2, h1, p1, h2, p2, N_int)
!            !DIR$ LOOP COUNT AVG(4)
!            do k = 1, N_states
!              mat_p(k, puti, putj) = mat_p(k, puti, putj) + coefs(k) * hij_p
!              mat_m(k, puti, putj) = mat_m(k, puti, putj) + coefs(k) * hij_m
!            enddo
!          endif
!        end do
!      end do
!    end if
!
!  else
!    if(tip == 0) then
!      h1 = h(1, ma)
!      h2 = h(2, ma)
!      do i=1,3
!        puti = p(i, ma)
!        if(bannedOrb(puti,ma)) cycle
!        do j=i+1,4
!          putj = p(j, ma)
!          if(bannedOrb(putj,ma)) cycle
!          if(banned(puti,putj,1)) cycle
!
!          i1 = turn2d(1, i, j)
!          i2 = turn2d(2, i, j)
!          p1 = p(i1, ma)
!          p2 = p(i2, ma)
!
!          hij_p = get_mo_bi_ortho_tc_two_e_tc_int   (p1, p2, h1, h2, mo_integrals_tc_int_map   ) &
!                - get_mo_bi_ortho_tc_two_e_tc_int   (p2, p1, h1, h2, mo_integrals_tc_int_map   )
!          hij_m = get_mo_bi_ortho_tc_two_e_tcdag_int(p1, p2, h1, h2, mo_integrals_tcdag_int_map) &
!                - get_mo_bi_ortho_tc_two_e_tcdag_int(p2, p1, h1, h2, mo_integrals_tcdag_int_map)
!
!          if( (hij_p.eq.0.d0) .and. (hij_m.eq.0.d0) ) cycle
!
!          hij_p = hij_p * get_phase_bi(phasemask, ma, ma, h1, p1, h2, p2, N_int)
!          hij_m = hij_m * get_phase_bi(phasemask, ma, ma, h1, p1, h2, p2, N_int)
!
!          !DIR$ LOOP COUNT AVG(4)
!          do k=1,N_states
!            mat_p(k, puti, putj) = mat_p(k, puti, putj) + coefs(k) * hij_p
!            mat_m(k, puti, putj) = mat_m(k, puti, putj) + coefs(k) * hij_m
!          enddo
!        end do
!      end do
!
!    else if(tip == 3) then
!      h1 = h(1, mi)
!      h2 = h(1, ma)
!      p1 = p(1, mi)
!      do i=1,3
!        puti = p(turn3(1,i), ma)
!        if(bannedOrb(puti,ma)) cycle
!        putj = p(turn3(2,i), ma)
!        if(bannedOrb(putj,ma)) cycle
!        if(banned(puti,putj,1)) cycle
!        p2 = p(i, ma)
!
!        hij_p = get_mo_bi_ortho_tc_two_e_tc_int   (p1, p2, h1, h2, mo_integrals_tc_int_map   )
!        hij_m = get_mo_bi_ortho_tc_two_e_tcdag_int(p1, p2, h1, h2, mo_integrals_tcdag_int_map)
!
!        if( (hij_p.eq.0.d0) .and. (hij_m.eq.0.d0) ) cycle
!
!        hij_p = hij_p * get_phase_bi(phasemask, mi, ma, h1, p1, h2, p2, N_int)
!        hij_m = hij_m * get_phase_bi(phasemask, mi, ma, h1, p1, h2, p2, N_int)
!        if (puti < putj) then
!          !DIR$ LOOP COUNT AVG(4)
!          do k=1,N_states
!            mat_p(k, puti, putj) = mat_p(k, puti, putj) + coefs(k) * hij_p
!            mat_m(k, puti, putj) = mat_m(k, puti, putj) + coefs(k) * hij_m
!          enddo
!        else
!          !DIR$ LOOP COUNT AVG(4)
!          do k=1,N_states
!            mat_p(k, putj, puti) = mat_p(k, putj, puti) + coefs(k) * hij_p
!            mat_m(k, putj, puti) = mat_m(k, putj, puti) + coefs(k) * hij_m
!          enddo
!        endif
!      end do
!    else ! tip == 4
!      puti = p(1, sp)
!      putj = p(2, sp)
!      if(.not. banned(puti,putj,1)) then
!        p1 = p(1, mi)
!        p2 = p(2, mi)
!        h1 = h(1, mi)
!        h2 = h(2, mi)
!
!        hij_p = get_mo_bi_ortho_tc_two_e_tc_int   (p1, p2, h1, h2, mo_integrals_tc_int_map   ) &
!              - get_mo_bi_ortho_tc_two_e_tc_int   (p2, p1, h1, h2, mo_integrals_tc_int_map   )
!        hij_m = get_mo_bi_ortho_tc_two_e_tcdag_int(p1, p2, h1, h2, mo_integrals_tcdag_int_map) &
!              - get_mo_bi_ortho_tc_two_e_tcdag_int(p2, p1, h1, h2, mo_integrals_tcdag_int_map)
!
!        if( (hij_p.ne.0.d0) .and. (hij_m.ne.0.d0) ) then
!          hij_p = hij_p * get_phase_bi(phasemask, mi, mi, h1, p1, h2, p2, N_int)
!          hij_m = hij_m * get_phase_bi(phasemask, mi, mi, h1, p1, h2, p2, N_int)
!          !DIR$ LOOP COUNT AVG(4)
!          do k=1,N_states
!            mat_p(k, puti, putj) = mat_p(k, puti, putj) + coefs(k) * hij_p
!            mat_m(k, puti, putj) = mat_m(k, puti, putj) + coefs(k) * hij_m
!          enddo
!        end if
!      end if
!    end if
!  end if
!
!end subroutine get_pm2
! ___________________________________________________________________________________________________________________________________________________
! ___________________________________________________________________________________________________________________________________________________
! ___________________________________________________________________________________________________________________________________________________


! ___________________________________________________________________________________________________________________________________________________
! ___________________________________________________________________________________________________________________________________________________

!subroutine get_pm1(gen, phasemask, bannedOrb, banned, mat_p, mat_m, mask, h, p, sp, coefs)
!
!  use bitmasks
!
!  implicit none
!
!  integer(bit_kind)               :: det(N_int, 2)
!  integer(bit_kind), intent(in)   :: mask(N_int, 2), gen(N_int, 2)
!  integer(bit_kind), intent(in)   :: phasemask(N_int,2)
!  logical,           intent(in)   :: bannedOrb(mo_num, 2), banned(mo_num, mo_num,2)
!  double precision,  intent(in)   :: coefs(N_states)
!  integer,           intent(in)   :: h(0:2,2), p(0:4,2), sp
!  double precision, intent(inout) :: mat_p(N_states, mo_num, mo_num), mat_m(N_states, mo_num, mo_num)
!
!  double precision, external      :: get_phase_bi
!  double precision, external      :: get_mo_bi_ortho_tc_two_e_tc_int, get_mo_bi_ortho_tc_two_e_tcdag_int
!
!  logical                         :: ok
!  logical, allocatable            :: lbanned(:,:)
!  integer                         :: bant
!  integer                         :: puti, putj, ma, mi, s1, s2, i, i1, i2, j
!  integer                         :: hfix, pfix, h1, h2, p1, p2, ib, k, l
!  double precision                :: tmp_row_ij_p (N_states, mo_num), tmp_row_ij_m (N_states, mo_num)
!  double precision                :: hij_p, hij_m, tmp_row_ij2_p(N_states, mo_num), tmp_row_ij2_m(N_states, mo_num)
!  double precision, allocatable   :: hijp_cache(:,:), hijm_cache(:,:)
!
!  integer, parameter              :: turn2(2) = (/2,1/)
!  integer, parameter              :: turn3(2,3) = reshape((/2,3,  1,3, 1,2/), (/2,3/))
!
!  PROVIDE mo_integrals_tc_int_map mo_integrals_tcdag_int_map
!
!  allocate( lbanned(mo_num, 2) )
!  allocate( hijp_cache(mo_num,2), hijm_cache(mo_num,2) )
!  lbanned = bannedOrb
!
!  do i=1, p(0,1)
!    lbanned(p(i,1), 1) = .true.
!  end do
!  do i=1, p(0,2)
!    lbanned(p(i,2), 2) = .true.
!  end do
!
!  ma = 1
!  if(p(0,2) >= 2) ma = 2
!  mi = turn2(ma)
!
!  bant = 1
!
!  if(sp == 3) then
!    !move MA
!    if(ma == 2) bant = 2
!    puti = p(1,mi)
!    hfix = h(1,ma)
!    p1 = p(1,ma)
!    p2 = p(2,ma)
!    if(.not. bannedOrb(puti, mi)) then
!
!      call get_mo_bi_ortho_tc_two_es_tc_int   (hfix, p1, p2, mo_num, hijp_cache(1,1), mo_integrals_tc_int_map   )
!      call get_mo_bi_ortho_tc_two_es_tc_int   (hfix, p2, p1, mo_num, hijp_cache(1,2), mo_integrals_tc_int_map   )
!      call get_mo_bi_ortho_tc_two_es_tcdag_int(hfix, p1, p2, mo_num, hijm_cache(1,1), mo_integrals_tcdag_int_map)
!      call get_mo_bi_ortho_tc_two_es_tcdag_int(hfix, p2, p1, mo_num, hijm_cache(1,2), mo_integrals_tcdag_int_map)
!
!      tmp_row_ij_p = 0d0
!      tmp_row_ij_m = 0d0
!      do putj=1, hfix-1
!        if(lbanned(putj, ma)) cycle
!        if(banned(putj, puti,bant)) cycle
!
!        hij_p = hijp_cache(putj,1) - hijp_cache(putj,2)
!        hij_m = hijm_cache(putj,1) - hijm_cache(putj,2)
!
!        if( (hij_p.ne.0.d0) .and. (hij_m.ne.0.d0) ) then
!          hij_p = hij_p * get_phase_bi(phasemask, ma, ma, putj, p1, hfix, p2, N_int)
!          hij_m = hij_m * get_phase_bi(phasemask, ma, ma, putj, p1, hfix, p2, N_int)
!          !DIR$ LOOP COUNT AVG(4)
!          do k=1,N_states
!            tmp_row_ij_p(k,putj) = tmp_row_ij_p(k,putj) + hij_p * coefs(k)
!            tmp_row_ij_m(k,putj) = tmp_row_ij_m(k,putj) + hij_m * coefs(k)
!          enddo
!        endif
!      end do
!      do putj=hfix+1, mo_num
!        if(lbanned(putj, ma)) cycle
!        if(banned(putj, puti,bant)) cycle
!
!        hij_p = hijp_cache(putj,2) - hijp_cache(putj,1)
!        hij_m = hijm_cache(putj,2) - hijm_cache(putj,1)
!
!        if( (hij_p.ne.0.d0) .and. (hij_m.ne.0.d0) ) then
!          hij_p = hij_p * get_phase_bi(phasemask, ma, ma, hfix, p1, putj, p2, N_int)
!          hij_m = hij_m * get_phase_bi(phasemask, ma, ma, hfix, p1, putj, p2, N_int)
!          !DIR$ LOOP COUNT AVG(4)
!          do k=1,N_states
!            tmp_row_ij_p(k,putj) = tmp_row_ij_p(k,putj) + hij_p * coefs(k)
!            tmp_row_ij_m(k,putj) = tmp_row_ij_m(k,putj) + hij_m * coefs(k)
!          enddo
!        endif
!      end do
!
!      if(ma == 1) then
!        mat_p(1:N_states,1:mo_num,puti) = mat_p(1:N_states,1:mo_num,puti) + tmp_row_ij_p(1:N_states,1:mo_num)
!        mat_m(1:N_states,1:mo_num,puti) = mat_m(1:N_states,1:mo_num,puti) + tmp_row_ij_m(1:N_states,1:mo_num)
!      else
!        do l=1,mo_num
!          !DIR$ LOOP COUNT AVG(4)
!          do k=1,N_states
!            mat_p(k,puti,l) = mat_p(k,puti,l) + tmp_row_ij_p(k,l)
!            mat_m(k,puti,l) = mat_m(k,puti,l) + tmp_row_ij_m(k,l)
!          enddo
!        enddo
!      end if
!    end if
!
!    !MOVE MI
!    pfix = p(1,mi)
!    tmp_row_ij_p  = 0d0
!    tmp_row_ij_m  = 0d0
!    tmp_row_ij2_p = 0d0
!    tmp_row_ij2_m = 0d0
!
!    call get_mo_bi_ortho_tc_two_es_tc_int   (hfix, pfix, p1, mo_num, hijp_cache(1,1), mo_integrals_tc_int_map   )
!    call get_mo_bi_ortho_tc_two_es_tc_int   (hfix, pfix, p2, mo_num, hijp_cache(1,2), mo_integrals_tc_int_map   )
!    call get_mo_bi_ortho_tc_two_es_tcdag_int(hfix, pfix, p1, mo_num, hijm_cache(1,1), mo_integrals_tcdag_int_map)
!    call get_mo_bi_ortho_tc_two_es_tcdag_int(hfix, pfix, p2, mo_num, hijm_cache(1,2), mo_integrals_tcdag_int_map)
!
!    putj = p1
!    do puti=1,mo_num !HOT
!      if(lbanned(puti,mi)) cycle
!      !p1 fixed
!      putj = p1
!      if(.not. banned(putj,puti,bant)) then
!
!        hij_p = hijp_cache(puti,2)
!        hij_m = hijm_cache(puti,2)
!
!        if( (hij_p.ne.0.d0) .and. (hij_m.ne.0.d0) ) then
!          hij_p = hij_p * get_phase_bi(phasemask, ma, mi, hfix, p2, puti, pfix, N_int)
!          hij_m = hij_m * get_phase_bi(phasemask, ma, mi, hfix, p2, puti, pfix, N_int)
!          !DIR$ LOOP COUNT AVG(4)
!          do k=1,N_states
!            tmp_row_ij_p(k,puti) = tmp_row_ij_p(k,puti) + hij_p * coefs(k)
!            tmp_row_ij_m(k,puti) = tmp_row_ij_m(k,puti) + hij_m * coefs(k)
!          enddo
!        endif
!      end if
!
!      putj = p2
!      if(.not. banned(putj,puti,bant)) then
!
!        hij_p = hijp_cache(puti,1)
!        hij_m = hijm_cache(puti,1)
!
!        if( (hij_p.ne.0.d0) .and. (hij_m.ne.0.d0) ) then
!          hij_p = hij_p * get_phase_bi(phasemask, ma, mi, hfix, p1, puti, pfix, N_int)
!          hij_m = hij_m * get_phase_bi(phasemask, ma, mi, hfix, p1, puti, pfix, N_int)
!          do k=1,N_states
!            tmp_row_ij2_p(k,puti) = tmp_row_ij2_p(k,puti) + hij_p * coefs(k)
!            tmp_row_ij2_m(k,puti) = tmp_row_ij2_m(k,puti) + hij_m * coefs(k)
!          enddo
!        endif
!      end if
!    end do
!
!    if(mi == 1) then
!      mat_p(:,:,p1) = mat_p(:,:,p1) + tmp_row_ij_p (:,:)
!      mat_p(:,:,p2) = mat_p(:,:,p2) + tmp_row_ij2_p(:,:)
!      mat_m(:,:,p1) = mat_m(:,:,p1) + tmp_row_ij_m (:,:)
!      mat_m(:,:,p2) = mat_m(:,:,p2) + tmp_row_ij2_m(:,:)
!    else
!      do l=1,mo_num
!        !DIR$ LOOP COUNT AVG(4)
!        do k=1,N_states
!          mat_p(k,p1,l) = mat_p(k,p1,l) + tmp_row_ij_p (k,l)
!          mat_p(k,p2,l) = mat_p(k,p2,l) + tmp_row_ij2_p(k,l)
!          mat_m(k,p1,l) = mat_m(k,p1,l) + tmp_row_ij_m (k,l)
!          mat_m(k,p2,l) = mat_m(k,p2,l) + tmp_row_ij2_m(k,l)
!        enddo
!      enddo
!    end if
!
!  else  ! sp /= 3
!
!    if(p(0,ma) == 3) then
!      do i=1,3
!        hfix = h(1,ma)
!        puti = p(i, ma)
!        p1 = p(turn3(1,i), ma)
!        p2 = p(turn3(2,i), ma)
!
!        call get_mo_bi_ortho_tc_two_es_tc_int   (hfix, p1, p2, mo_num, hijp_cache(1,1), mo_integrals_tc_int_map   )
!        call get_mo_bi_ortho_tc_two_es_tc_int   (hfix, p2, p1, mo_num, hijp_cache(1,2), mo_integrals_tc_int_map   )
!        call get_mo_bi_ortho_tc_two_es_tcdag_int(hfix, p1, p2, mo_num, hijm_cache(1,1), mo_integrals_tcdag_int_map)
!        call get_mo_bi_ortho_tc_two_es_tcdag_int(hfix, p2, p1, mo_num, hijm_cache(1,2), mo_integrals_tcdag_int_map)
!
!        tmp_row_ij_p = 0d0
!        tmp_row_ij_m = 0d0
!        do putj=1,hfix-1
!          if(banned(putj,puti,1)) cycle
!          if(lbanned(putj,ma)) cycle
!
!          hij_p = hijp_cache(putj,1) - hijp_cache(putj,2)
!          hij_m = hijm_cache(putj,1) - hijm_cache(putj,2)
!
!          if( (hij_p.ne.0.d0) .and. (hij_m.ne.0.d0) ) then
!            hij_p = hij_p * get_phase_bi(phasemask, ma, ma, putj, p1, hfix, p2, N_int)
!            hij_m = hij_m * get_phase_bi(phasemask, ma, ma, putj, p1, hfix, p2, N_int)
!            tmp_row_ij_p(:,putj) = tmp_row_ij_p(:,putj) + hij_p * coefs(:)
!            tmp_row_ij_m(:,putj) = tmp_row_ij_m(:,putj) + hij_m * coefs(:)
!          endif
!        end do
!        do putj=hfix+1,mo_num
!          if(banned(putj,puti,1)) cycle
!          if(lbanned(putj,ma)) cycle
!
!          hij_p = hijp_cache(putj,2) - hijp_cache(putj,1)
!          hij_m = hijm_cache(putj,2) - hijm_cache(putj,1)
!
!          if( (hij_p.ne.0.d0) .and. (hij_m.ne.0.d0) ) then
!            hij_p = hij_p * get_phase_bi(phasemask, ma, ma, hfix, p1, putj, p2, N_int)
!            hij_m = hij_m * get_phase_bi(phasemask, ma, ma, hfix, p1, putj, p2, N_int)
!            tmp_row_ij_p(:,putj) = tmp_row_ij_p(:,putj) + hij_p * coefs(:)
!            tmp_row_ij_m(:,putj) = tmp_row_ij_m(:,putj) + hij_m * coefs(:)
!          endif
!        end do
!
!        mat_p(:, :puti-1, puti) = mat_p(:, :puti-1, puti) + tmp_row_ij_p(:,:puti-1)
!        mat_m(:, :puti-1, puti) = mat_m(:, :puti-1, puti) + tmp_row_ij_m(:,:puti-1)
!        do l=puti,mo_num
!          !DIR$ LOOP COUNT AVG(4)
!          do k=1,N_states
!            mat_p(k, puti, l) = mat_p(k, puti,l) + tmp_row_ij_p(k,l)
!            mat_m(k, puti, l) = mat_m(k, puti,l) + tmp_row_ij_m(k,l)
!          enddo
!        enddo
!      end do
!    else
!      hfix = h(1,mi)
!      pfix = p(1,mi)
!      p1 = p(1,ma)
!      p2 = p(2,ma)
!      tmp_row_ij_p  = 0d0
!      tmp_row_ij_m  = 0d0
!      tmp_row_ij2_p = 0d0
!      tmp_row_ij2_m = 0d0
!
!      call get_mo_bi_ortho_tc_two_es_tc_int   (hfix, p1, pfix, mo_num, hijp_cache(1,1), mo_integrals_tc_int_map   )
!      call get_mo_bi_ortho_tc_two_es_tc_int   (hfix, p2, pfix, mo_num, hijp_cache(1,2), mo_integrals_tc_int_map   )
!      call get_mo_bi_ortho_tc_two_es_tcdag_int(hfix, p1, pfix, mo_num, hijp_cache(1,1), mo_integrals_tcdag_int_map)
!      call get_mo_bi_ortho_tc_two_es_tcdag_int(hfix, p2, pfix, mo_num, hijp_cache(1,2), mo_integrals_tcdag_int_map)
!
!      putj = p2
!      do puti=1,mo_num
!        if(lbanned(puti,ma)) cycle
!        putj = p2
!        if(.not. banned(puti,putj,1)) then
!
!          hij_p = hijp_cache(puti,1)
!          hij_m = hijm_cache(puti,1)
!
!          if( (hij_p.ne.0.d0) .and. (hij_m.ne.0.d0) ) then
!            hij_p = hij_p * get_phase_bi(phasemask, mi, ma, hfix, pfix, puti, p1, N_int)
!            hij_m = hij_m * get_phase_bi(phasemask, mi, ma, hfix, pfix, puti, p1, N_int)
!            !DIR$ LOOP COUNT AVG(4)
!            do k=1,N_states
!              tmp_row_ij_p(k,puti) = tmp_row_ij_p(k,puti) + hij_p * coefs(k)
!              tmp_row_ij_m(k,puti) = tmp_row_ij_m(k,puti) + hij_m * coefs(k)
!            enddo
!          endif
!        end if
!
!        putj = p1
!        if(.not. banned(puti,putj,1)) then
!          hij_p = hijp_cache(puti,2)
!          hij_m = hijm_cache(puti,2)
!          if( (hij_p.ne.0.d0) .and. (hij_m.ne.0.d0) ) then
!            hij_p = hij_p * get_phase_bi(phasemask, mi, ma, hfix, pfix, puti, p2, N_int)
!            hij_m = hij_m * get_phase_bi(phasemask, mi, ma, hfix, pfix, puti, p2, N_int)
!            do k=1,N_states
!              tmp_row_ij2_p(k,puti) = tmp_row_ij2_p(k,puti) + hij_p * coefs(k)
!              tmp_row_ij2_m(k,puti) = tmp_row_ij2_m(k,puti) + hij_m * coefs(k)
!            enddo
!          endif
!        end if
!      end do
!      mat_p(:,:p2-1,p2) = mat_p(:,:p2-1,p2) + tmp_row_ij_p(:,:p2-1)
!      mat_m(:,:p2-1,p2) = mat_m(:,:p2-1,p2) + tmp_row_ij_m(:,:p2-1)
!      do l=p2,mo_num
!        !DIR$ LOOP COUNT AVG(4)
!        do k=1,N_states
!          mat_p(k,p2,l) = mat_p(k,p2,l) + tmp_row_ij_p(k,l)
!          mat_m(k,p2,l) = mat_m(k,p2,l) + tmp_row_ij_m(k,l)
!        enddo
!      enddo
!      mat_p(:,:p1-1,p1) = mat_p(:,:p1-1,p1) + tmp_row_ij2_p(:,:p1-1)
!      mat_m(:,:p1-1,p1) = mat_m(:,:p1-1,p1) + tmp_row_ij2_m(:,:p1-1)
!      do l=p1,mo_num
!        !DIR$ LOOP COUNT AVG(4)
!        do k=1,N_states
!          mat_p(k,p1,l) = mat_p(k,p1,l) + tmp_row_ij2_p(k,l)
!          mat_m(k,p1,l) = mat_m(k,p1,l) + tmp_row_ij2_m(k,l)
!        enddo
!      enddo
!    end if
!  end if
!  deallocate(lbanned,hijp_cache, hijm_cache)
!
! !! MONO
!  if(sp == 3) then
!    s1 = 1
!    s2 = 2
!  else
!    s1 = sp
!    s2 = sp
!  end if
!
!  do i1 = 1, p(0,s1)
!    ib = 1
!    if(s1 == s2) ib = i1+1
!    do i2 = ib, p(0,s2)
!      p1 = p(i1,s1)
!      p2 = p(i2,s2)
!      if(bannedOrb(p1, s1) .or. bannedOrb(p2, s2) .or. banned(p1, p2, 1)) cycle
!      call apply_particles(mask, s1, p1, s2, p2, det, ok, N_int)
!
!      call htilde_mu_mat_tot   (gen, det, N_int, hij_p)
!      call htildedag_mu_mat_tot(gen, det, N_int, hij_m)
!
!      !DIR$ LOOP COUNT AVG(4)
!      do k = 1, N_states
!        mat_p(k, p1, p2) = mat_p(k, p1, p2) + coefs(k) * hij_p
!        mat_m(k, p1, p2) = mat_m(k, p1, p2) + coefs(k) * hij_m
!      enddo
!    enddo
!  enddo
!
!end subroutine get_pm1
! ___________________________________________________________________________________________________________________________________________________
! ___________________________________________________________________________________________________________________________________________________
! ___________________________________________________________________________________________________________________________________________________



! ___________________________________________________________________________________________________________________________________________________
! ___________________________________________________________________________________________________________________________________________________

!subroutine get_pm0(gen, phasemask, bannedOrb, banned, mat_p, mat_m, mask, h, p, sp, coefs)
!
!  use bitmasks
!  implicit none
!
!  integer(bit_kind)               :: det(N_int, 2)
!  integer(bit_kind), intent(in)   :: gen(N_int, 2), mask(N_int, 2)
!  integer(bit_kind), intent(in)   :: phasemask(N_int,2)
!  integer,           intent(in)   :: h(0:2,2), p(0:4,2), sp
!  logical,           intent(in)   :: bannedOrb(mo_num, 2), banned(mo_num, mo_num,2)
!  double precision,  intent(in)   :: coefs(N_states)
!  double precision, intent(inout) :: mat_p(N_states, mo_num, mo_num), mat_m(N_states, mo_num, mo_num)
!
!  double precision, external      :: get_phase_bi
!  double precision, external      :: get_mo_bi_ortho_tc_two_e_tc_int, get_mo_bi_ortho_tc_two_e_tcdag_int 
!  integer, parameter              :: bant=1
!  integer                         :: i, j, k, s, h1, h2, p1, p2, puti, putj
!  logical                         :: ok
!  double precision                :: hij_p, hij_m, phase
!  double precision, allocatable   :: hijp_cache1(:), hijp_cache2(:), hijm_cache1(:), hijm_cache2(:)
!
!  PROVIDE mo_integrals_tc_int_map mo_integrals_tcdag_int_map
!
!  allocate( hijp_cache1(mo_num) , hijp_cache2(mo_num) )
!  allocate( hijm_cache1(mo_num) , hijm_cache2(mo_num) )
!
!  if(sp == 3) then ! AB
!    h1 = p(1,1)
!    h2 = p(1,2)
!    do p1=1, mo_num
!      if(bannedOrb(p1, 1)) cycle
!
!      call get_mo_bi_ortho_tc_two_es_tc_int   (p1, h2, h1, mo_num, hijp_cache1, mo_integrals_tc_int_map   )
!      call get_mo_bi_ortho_tc_two_es_tcdag_int(p1, h2, h1, mo_num, hijm_cache1, mo_integrals_tcdag_int_map)
!
!      do p2 = 1, mo_num
!        if(bannedOrb(p2,2)) cycle
!        if(banned(p1, p2, bant)) cycle ! rentable?
!        if(p1 == h1 .or. p2 == h2) then
!          call apply_particles(mask, 1, p1, 2, p2, det, ok, N_int)
!          call htilde_mu_mat_tot   (gen, det, N_int, hij_p)
!          call htildedag_mu_mat_tot(gen, det, N_int, hij_m)
!        else
!          phase = get_phase_bi(phasemask, 1, 2, h1, p1, h2, p2, N_int)
!          hij_p = hijp_cache1(p2) * phase
!          hij_m = hijm_cache1(p2) * phase
!        end if
!        if( (hij_p.eq.0.d0).and.(hij_m.eq.0.d0) ) cycle
!        !DIR$ LOOP COUNT AVG(4)
!        do k = 1, N_states
!          mat_p(k, p1, p2) = mat_p(k, p1, p2) + coefs(k) * hij_p  ! HOTSPOT
!          mat_m(k, p1, p2) = mat_m(k, p1, p2) + coefs(k) * hij_m  ! HOTSPOT
!        enddo
!      end do
!    end do
!
!  else ! AA BB
!    p1 = p(1,sp)
!    p2 = p(2,sp)
!    do puti=1, mo_num
!      if(bannedOrb(puti, sp)) cycle
!
!      call get_mo_bi_ortho_tc_two_es_tc_int   (puti, p2, p1, mo_num, hijp_cache1, mo_integrals_tc_int_map   )
!      call get_mo_bi_ortho_tc_two_es_tc_int   (puti, p1, p2, mo_num, hijp_cache2, mo_integrals_tc_int_map   )
!      call get_mo_bi_ortho_tc_two_es_tcdag_int(puti, p2, p1, mo_num, hijm_cache1, mo_integrals_tcdag_int_map)
!      call get_mo_bi_ortho_tc_two_es_tcdag_int(puti, p1, p2, mo_num, hijm_cache2, mo_integrals_tcdag_int_map)
!
!      do putj=puti+1, mo_num
!        if(bannedOrb(putj, sp)) cycle
!        if(banned(puti, putj, bant)) cycle ! rentable?
!        if(puti == p1 .or. putj == p2 .or. puti == p2 .or. putj == p1) then
!          call apply_particles(mask, sp, puti, sp, putj, det, ok, N_int)
!          call htilde_mu_mat_tot   (gen, det, N_int, hij_p)
!          call htildedag_mu_mat_tot(gen, det, N_int, hij_m)
!          if( (hij_p.eq.0.d0).and.(hij_m.eq.0.d0) ) cycle
!        else
!
!          hij_p = get_mo_bi_ortho_tc_two_e_tc_int   (p1, p2, puti, putj, mo_integrals_tc_int_map   ) &
!                - get_mo_bi_ortho_tc_two_e_tc_int   (p2, p1, puti, putj, mo_integrals_tc_int_map   )
!          hij_m = get_mo_bi_ortho_tc_two_e_tcdag_int(p1, p2, puti, putj, mo_integrals_tcdag_int_map) &
!                - get_mo_bi_ortho_tc_two_e_tcdag_int(p2, p1, puti, putj, mo_integrals_tcdag_int_map)
!
!          if( (hij_p.eq.0.d0).and.(hij_m.eq.0.d0) ) cycle
!
!          hij_p = hij_p * get_phase_bi(phasemask, sp, sp, puti, p1 , putj, p2, N_int)
!          hij_m = hij_m * get_phase_bi(phasemask, sp, sp, puti, p1 , putj, p2, N_int)
!
!        end if
!        !DIR$ LOOP COUNT AVG(4)
!        do k=1,N_states
!          mat_p(k, puti, putj) = mat_p(k, puti, putj) + coefs(k) * hij_p
!          mat_m(k, puti, putj) = mat_m(k, puti, putj) + coefs(k) * hij_m
!        enddo
!      end do
!    end do
!  end if
!
!  deallocate( hijp_cache1 , hijp_cache2 )
!  deallocate( hijm_cache1 , hijm_cache2 )
!
!end subroutine get_pm0
! ___________________________________________________________________________________________________________________________________________________
! ___________________________________________________________________________________________________________________________________________________
! ___________________________________________________________________________________________________________________________________________________


! OLD unoptimized routines for debugging
! ======================================

subroutine get_d0_reference(gen, phasemask, bannedOrb, banned, mat, mask, h, p, sp, coefs)

  use bitmasks
  implicit none

  integer(bit_kind), intent(in) :: gen(N_int, 2), mask(N_int, 2)
  integer(bit_kind), intent(in) :: phasemask(N_int,2)
  logical, intent(in) :: bannedOrb(mo_num, 2), banned(mo_num, mo_num,2)
  integer(bit_kind) :: det(N_int, 2)
  double precision, intent(in) :: coefs(N_states)
  double precision, intent(inout) :: mat(N_states, mo_num, mo_num)
  integer, intent(in) :: h(0:2,2), p(0:4,2), sp

  integer :: i, j, s, h1, h2, p1, p2, puti, putj
  double precision :: hij, phase
  double precision, external :: get_phase_bi
  logical :: ok

  integer :: bant
  bant = 1


  if(sp == 3) then ! AB
    h1 = p(1,1)
    h2 = p(1,2)
    do p1=1, mo_num
      if(bannedOrb(p1, 1)) cycle
      do p2=1, mo_num
        if(bannedOrb(p2,2)) cycle
        if(banned(p1, p2, bant)) cycle ! rentable?
        if(p1 == h1 .or. p2 == h2) then
          call apply_particles(mask, 1,p1,2,p2, det, ok, N_int)
          call i_h_j(gen, det, N_int, hij)
        else
          phase = get_phase_bi(phasemask, 1, 2, h1, p1, h2, p2, N_int)
          hij = mo_bi_ortho_tc_two_e(p1, p2, h1, h2) * phase
        end if
        mat(:, p1, p2) = mat(:, p1, p2) + coefs(:) * hij
      end do
    end do
  else ! AA BB
    p1 = p(1,sp)
    p2 = p(2,sp)
    do puti=1, mo_num
      if(bannedOrb(puti, sp)) cycle
      do putj=puti+1, mo_num
        if(bannedOrb(putj, sp)) cycle
        if(banned(puti, putj, bant)) cycle ! rentable?
        if(puti == p1 .or. putj == p2 .or. puti == p2 .or. putj == p1) then
          call apply_particles(mask, sp,puti,sp,putj, det, ok, N_int)
          call i_h_j(gen, det, N_int, hij)
        else
          hij = (mo_bi_ortho_tc_two_e(p1, p2, puti, putj) -  mo_bi_ortho_tc_two_e(p2, p1, puti, putj))* get_phase_bi(phasemask, sp, sp, puti, p1 , putj, p2, N_int)
        end if
        mat(:, puti, putj) = mat(:, puti, putj) + coefs(:) * hij
      end do
    end do
  end if

end subroutine get_d0_reference

! ---

subroutine get_d1_reference(gen, phasemask, bannedOrb, banned, mat, mask, h, p, sp, coefs)

  use bitmasks
  implicit none

  integer(bit_kind), intent(in)  :: mask(N_int, 2), gen(N_int, 2)
  integer(bit_kind), intent(in)  :: phasemask(N_int,2)
  logical, intent(in)            :: bannedOrb(mo_num, 2), banned(mo_num, mo_num,2)
  integer(bit_kind)              :: det(N_int, 2)
  double precision, intent(in)   :: coefs(N_states)
  double precision, intent(inout) :: mat(N_states, mo_num, mo_num)
  integer, intent(in)            :: h(0:2,2), p(0:4,2), sp
  double precision               :: hij, tmp_row_ij(N_states, mo_num), tmp_row_ij2(N_states, mo_num), hji
  double precision, external     :: get_phase_bi
  logical                        :: ok

  logical, allocatable           :: lbanned(:,:)
  integer                        :: puti, putj, ma, mi, s1, s2, i, i1, i2, j
  integer                        :: hfix, pfix, h1, h2, p1, p2, ib

  integer, parameter             :: turn2(2) = (/2,1/)
  integer, parameter             :: turn3(2,3) = reshape((/2,3,  1,3, 1,2/), (/2,3/))

  integer                        :: bant


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
      tmp_row_ij = 0d0
      do putj=1, hfix-1
        if(lbanned(putj, ma) .or. banned(putj, puti,bant)) cycle
        hij = (mo_bi_ortho_tc_two_e(p1, p2, putj, hfix)-mo_bi_ortho_tc_two_e(p2,p1,putj,hfix)) * get_phase_bi(phasemask, ma, ma, putj, p1, hfix, p2, N_int)
        tmp_row_ij(1:N_states,putj) = tmp_row_ij(1:N_states,putj) + hij * coefs(1:N_states)
      end do
      do putj=hfix+1, mo_num
        if(lbanned(putj, ma) .or. banned(putj, puti,bant)) cycle
        hij = (mo_bi_ortho_tc_two_e(p1, p2, hfix, putj)-mo_bi_ortho_tc_two_e(p2,p1,hfix,putj)) * get_phase_bi(phasemask, ma, ma, hfix, p1, putj, p2, N_int)
        tmp_row_ij(1:N_states,putj) = tmp_row_ij(1:N_states,putj) + hij * coefs(1:N_states)
      end do

      if(ma == 1) then
        mat(1:N_states,1:mo_num,puti) = mat(1:N_states,1:mo_num,puti) + tmp_row_ij(1:N_states,1:mo_num)
      else
        mat(1:N_states,puti,1:mo_num) = mat(1:N_states,puti,1:mo_num) + tmp_row_ij(1:N_states,1:mo_num)
      end if
    end if

    !MOVE MI
    pfix = p(1,mi)
    tmp_row_ij = 0d0
    tmp_row_ij2 = 0d0
    do puti=1,mo_num
      if(lbanned(puti,mi)) cycle
      !p1 fixed
      putj = p1
      if(.not. banned(putj,puti,bant)) then
        hij = mo_bi_ortho_tc_two_e(p2,pfix,hfix,puti) * get_phase_bi(phasemask, ma, mi, hfix, p2, puti, pfix, N_int)
        tmp_row_ij(:,puti) = tmp_row_ij(:,puti) + hij * coefs(:)
      end if

      putj = p2
      if(.not. banned(putj,puti,bant)) then
        hij = mo_bi_ortho_tc_two_e(p1,pfix,hfix,puti) * get_phase_bi(phasemask, ma, mi, hfix, p1, puti, pfix, N_int)
        tmp_row_ij2(:,puti) = tmp_row_ij2(:,puti) + hij * coefs(:)
      end if
    end do

    if(mi == 1) then
      mat(:,:,p1) = mat(:,:,p1) + tmp_row_ij(:,:)
      mat(:,:,p2) = mat(:,:,p2) + tmp_row_ij2(:,:)
    else
      mat(:,p1,:) = mat(:,p1,:) + tmp_row_ij(:,:)
      mat(:,p2,:) = mat(:,p2,:) + tmp_row_ij2(:,:)
    end if
  else
    if(p(0,ma) == 3) then
      do i=1,3
        hfix = h(1,ma)
        puti = p(i, ma)
        p1 = p(turn3(1,i), ma)
        p2 = p(turn3(2,i), ma)
        tmp_row_ij = 0d0
        do putj=1,hfix-1
          if(lbanned(putj,ma) .or. banned(puti,putj,1)) cycle
          hij = (mo_bi_ortho_tc_two_e(p1, p2, putj, hfix)-mo_bi_ortho_tc_two_e(p2,p1,putj,hfix)) * get_phase_bi(phasemask, ma, ma, putj, p1, hfix, p2, N_int)
          tmp_row_ij(:,putj) = tmp_row_ij(:,putj) + hij * coefs(:)
        end do
        do putj=hfix+1,mo_num
          if(lbanned(putj,ma) .or. banned(puti,putj,1)) cycle
          hij = (mo_bi_ortho_tc_two_e(p1, p2, hfix, putj)-mo_bi_ortho_tc_two_e(p2,p1,hfix,putj)) * get_phase_bi(phasemask, ma, ma, hfix, p1, putj, p2, N_int)
          tmp_row_ij(:,putj) = tmp_row_ij(:,putj) + hij * coefs(:)
        end do

        mat(:, :puti-1, puti) = mat(:, :puti-1, puti) + tmp_row_ij(:,:puti-1)
        mat(:, puti, puti:) = mat(:, puti, puti:) + tmp_row_ij(:,puti:)
      end do
    else
      hfix = h(1,mi)
      pfix = p(1,mi)
      p1 = p(1,ma)
      p2 = p(2,ma)
      tmp_row_ij = 0d0
      tmp_row_ij2 = 0d0
      do puti=1,mo_num
        if(lbanned(puti,ma)) cycle
        putj = p2
        if(.not. banned(puti,putj,1)) then
          hij = mo_bi_ortho_tc_two_e(pfix, p1, hfix, puti) * get_phase_bi(phasemask, mi, ma, hfix, pfix, puti, p1, N_int)
          tmp_row_ij(:,puti) = tmp_row_ij(:,puti) + hij * coefs(:)
        end if

        putj = p1
        if(.not. banned(puti,putj,1)) then
          hij = mo_bi_ortho_tc_two_e(pfix, p2, hfix, puti) * get_phase_bi(phasemask, mi, ma, hfix, pfix, puti, p2, N_int)
          tmp_row_ij2(:,puti) = tmp_row_ij2(:,puti) + hij * coefs(:)
        end if
      end do
      mat(:,:p2-1,p2) = mat(:,:p2-1,p2) + tmp_row_ij(:,:p2-1)
      mat(:,p2,p2:) = mat(:,p2,p2:) + tmp_row_ij(:,p2:)
      mat(:,:p1-1,p1) = mat(:,:p1-1,p1) + tmp_row_ij2(:,:p1-1)
      mat(:,p1,p1:) = mat(:,p1,p1:) + tmp_row_ij2(:,p1:)
    end if
  end if
  deallocate(lbanned)

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
        call i_h_j(gen, det, N_int, hij)
        mat(:, p1, p2) = mat(:, p1, p2) + coefs(:) * hij
      end do
    end do

end subroutine get_d1_reference

! ---

subroutine get_d2_reference(gen, phasemask, bannedOrb, banned, mat, mask, h, p, sp, coefs)

  use bitmasks
  implicit none

  integer(bit_kind), intent(in) :: mask(N_int, 2), gen(N_int, 2)
  integer(bit_kind), intent(in) :: phasemask(2,N_int)
  logical, intent(in) :: bannedOrb(mo_num, 2), banned(mo_num, mo_num,2)
  double precision, intent(in) :: coefs(N_states)
  double precision, intent(inout) :: mat(N_states, mo_num, mo_num)
  integer, intent(in) :: h(0:2,2), p(0:4,2), sp

  double precision, external :: get_phase_bi

  integer :: i, j, tip, ma, mi, puti, putj
  integer :: h1, h2, p1, p2, i1, i2, mm
  double precision :: hij, phase, hji

  integer, parameter:: turn2d(2,3,4) = reshape((/0,0, 0,0, 0,0,  3,4, 0,0, 0,0,  2,4, 1,4, 0,0,  2,3, 1,3, 1,2 /), (/2,3,4/))
  integer, parameter :: turn2(2) = (/2, 1/)
  integer, parameter :: turn3(2,3) = reshape((/2,3,  1,3, 1,2/), (/2,3/))

  integer :: bant
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
        i1 = turn3(1,i)
        i2 = turn3(2,i)
        p1 = p(i1, ma)
        p2 = p(i2, ma)
        h1 = h(1, ma)
        h2 = h(2, ma)

        hij = (mo_bi_ortho_tc_two_e(p1, p2, h1, h2) - mo_bi_ortho_tc_two_e(p2,p1, h1, h2)) * get_phase_bi(phasemask, ma, ma, h1, p1, h2, p2, N_int)
        if(ma == 1) then
          mat(:, putj, puti) = mat(:, putj, puti) + coefs(:) * hij
        else
          mat(:, puti, putj) = mat(:, puti, putj) + coefs(:) * hij
        end if
      end do
    else
      h1 = h(1,1)
      h2 = h(1,2)
      do j = 1,2
        putj = p(j, 2)
        p2 = p(turn2(j), 2)
        do i = 1,2
          puti = p(i, 1)

          if(banned(puti,putj,bant)) cycle
          p1 = p(turn2(i), 1)

          hij = mo_bi_ortho_tc_two_e(p1, p2, h1, h2) * get_phase_bi(phasemask, 1, 2, h1, p1, h2, p2,N_int)
          mat(:, puti, putj) = mat(:, puti, putj) + coefs(:) * hij
        end do
      end do
    end if

  else
    if(tip == 0) then
      h1 = h(1, ma)
      h2 = h(2, ma)
      do i=1,3
      puti = p(i, ma)
      do j=i+1,4
        putj = p(j, ma)
        if(banned(puti,putj,1)) cycle

        i1 = turn2d(1, i, j)
        i2 = turn2d(2, i, j)
        p1 = p(i1, ma)
        p2 = p(i2, ma)
        hij = (mo_bi_ortho_tc_two_e(p1, p2, h1, h2) - mo_bi_ortho_tc_two_e(p2,p1, h1, h2)) * get_phase_bi(phasemask, ma, ma, h1, p1, h2, p2,N_int)
        mat(:, puti, putj) = mat(:, puti, putj) + coefs(:) * hij
      end do
      end do
    else if(tip == 3) then
      h1 = h(1, mi)
      h2 = h(1, ma)
      p1 = p(1, mi)
      do i=1,3
        puti = p(turn3(1,i), ma)
        putj = p(turn3(2,i), ma)
        if(banned(puti,putj,1)) cycle
        p2 = p(i, ma)

        hij = mo_bi_ortho_tc_two_e(p1, p2, h1, h2) * get_phase_bi(phasemask, mi, ma, h1, p1, h2, p2,N_int)
        mat(:, min(puti, putj), max(puti, putj)) = mat(:, min(puti, putj), max(puti, putj)) + coefs(:) * hij
      end do
    else ! tip == 4
      puti = p(1, sp)
      putj = p(2, sp)
      if(.not. banned(puti,putj,1)) then
        p1 = p(1, mi)
        p2 = p(2, mi)
        h1 = h(1, mi)
        h2 = h(2, mi)
        hij = (mo_bi_ortho_tc_two_e(p1, p2, h1, h2) - mo_bi_ortho_tc_two_e(p2,p1, h1, h2)) * get_phase_bi(phasemask, mi, mi, h1, p1, h2, p2,N_int)
        mat(:, puti, putj) = mat(:, puti, putj) + coefs(:) * hij
      end if
    end if
  end if

end subroutine get_d2_reference(gen, phasemask, bannedOrb, banned, mat, mask, h, p, sp, coefs)

! ---

