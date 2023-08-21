
subroutine get_d2_new(gen, phasemask, bannedOrb, banned, mat_l, mat_r, mask, h, p, sp, coefs)
  !todo: indices/conjg should be correct for complex
  use bitmasks
  implicit none

  integer(bit_kind), intent(in) :: mask(N_int, 2), gen(N_int, 2)
  integer(bit_kind), intent(in) :: phasemask(N_int,2)
  logical, intent(in) :: bannedOrb(mo_num, 2), banned(mo_num, mo_num,2)
  double precision, intent(in) :: coefs(N_states,2)
  double precision, intent(inout) :: mat_r(N_states, mo_num, mo_num)
  double precision, intent(inout) :: mat_l(N_states, mo_num, mo_num)
  integer, intent(in) :: h(0:2,2), p(0:4,2), sp

  double precision, external :: get_phase_bi

  integer :: i, j, k, tip, ma, mi, puti, putj
  integer :: h1, h2, p1, p2, i1, i2
  double precision :: phase
  double precision :: hij,hji

  integer, parameter:: turn2d(2,3,4) = reshape((/0,0, 0,0, 0,0,  3,4, 0,0, 0,0,  2,4, 1,4, 0,0,  2,3, 1,3, 1,2 /), (/2,3,4/))
  integer, parameter :: turn2(2) = (/2, 1/)
  integer, parameter :: turn3(2,3) = reshape((/2,3,  1,3, 1,2/), (/2,3/))

  integer :: bant
  bant = 1
!  print*, 'in get_d2_new'
!  call debug_det(gen,N_int)
!  print*,'coefs',coefs(1,:)

  tip = p(0,1) * p(0,2) ! number of alpha particles times number of beta particles

  ma = sp !1:(alpha,alpha); 2:(b,b); 3:(a,b)
  if(p(0,1) > p(0,2)) ma = 1 ! more alpha particles than beta particles
  if(p(0,1) < p(0,2)) ma = 2 ! fewer alpha particles than beta particles
  mi = mod(ma, 2) + 1

  if(sp == 3) then ! if one alpha and one beta xhole 
    !(where xholes refer to the ionizations from the generator, not the holes occupied in the ionized generator)
    if(ma == 2) bant = 2 ! if more beta particles than alpha particles

    if(tip == 3) then ! if 3 of one particle spin and 1 of the other particle spin
      puti = p(1, mi)
      if(bannedOrb(puti, mi)) return
      h1 = h(1, ma)
      h2 = h(2, ma)

      !! <alpha|H|psi>
      do i = 1, 3    ! loop over all 3 combinations of 2 particles with spin ma
        putj = p(i, ma)
        if(banned(putj,puti,bant)) cycle
        i1 = turn3(1,i)
        i2 = turn3(2,i)
        p1 = p(i1, ma)
        p2 = p(i2, ma)
        
     ! |G> = |psi_{gen,i}>
     ! |G'> = a_{x1} a_{x2} |G>
     ! |alpha> = a_{puti}^{\dagger} a_{putj}^{\dagger} |G'>
     ! |alpha> = t_{x1,x2}^{puti,putj} |G>
     ! hij = <psi_{selectors,i}|H|alpha>
     ! |alpha> = t_{p1,p2}^{h1,h2}|psi_{selectors,i}>
        !todo: <i|H|j>  =  (<h1,h2|p1,p2> - <h1,h2|p2,p1>) * phase
        !    <psi|H|j> +=  dconjg(c_i) * <i|H|j>
        !      <j|H|i>  =  (<p1,p2|h1,h2> - <p2,p1|h1,h2>) * phase
        !    <j|H|psi> +=  <j|H|i> * c_i
!        hij = mo_bi_ortho_tc_two_e(p1, p2, h1, h2) - mo_bi_ortho_tc_two_e(p2, p1, h1, h2)

!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!! 
        ! take the transpose of what's written above because later use the complex conjugate 
        hij = mo_bi_ortho_tc_two_e(h1, h2, p1, p2) - mo_bi_ortho_tc_two_e( h1, h2, p2, p1)
        if (hij == 0.d0) cycle

        ! take conjugate to get contribution to <alpha|H|psi> instead of <psi|H|alpha>
!        hij = dconjg(hij) * get_phase_bi(phasemask, ma, ma, h1, p1, h2, p2, N_int)
        hij = hij * get_phase_bi(phasemask, ma, ma, h1, p1, h2, p2, N_int)

        if(ma == 1) then ! if particle spins are (alpha,alpha,alpha,beta), then puti is beta and putj is alpha
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat_r(k, putj, puti) = mat_r(k, putj, puti) + coefs(k,2) * hij
          enddo
        else            ! if particle spins are (beta,beta,beta,alpha), then puti is alpha and putj is beta
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat_r(k, puti, putj) = mat_r(k, puti, putj) + coefs(k,2) * hij
          enddo
        end if
      end do
      !! <phi|H|alpha>
      do i = 1, 3    ! loop over all 3 combinations of 2 particles with spin ma
        putj = p(i, ma)
        if(banned(putj,puti,bant)) cycle
        i1 = turn3(1,i)
        i2 = turn3(2,i)
        p1 = p(i1, ma)
        p2 = p(i2, ma)
        hji = mo_bi_ortho_tc_two_e(p1, p2,h1, h2) - mo_bi_ortho_tc_two_e( p2, p1, h1, h2)
        if (hji == 0.d0) cycle
        hji = hji * get_phase_bi(phasemask, ma, ma, h1, p1, h2, p2, N_int)

        if(ma == 1) then ! if particle spins are (alpha,alpha,alpha,beta), then puti is beta and putj is alpha
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat_l(k, putj, puti) = mat_l(k, putj, puti) + coefs(k,1) * hji
          enddo
        else            ! if particle spins are (beta,beta,beta,alpha), then puti is alpha and putj is beta
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat_l(k, puti, putj) = mat_l(k, puti, putj) + coefs(k,1) * hji
          enddo
        end if
      end do
    else ! if 2 alpha and 2 beta particles
      h1 = h(1,1)
      h2 = h(1,2)
      !! <alpha|H|psi>
      do j = 1,2 ! loop over all 4 combinations of one alpha and one beta particle
        putj = p(j, 2)
        if(bannedOrb(putj, 2)) cycle
        p2 = p(turn2(j), 2)
        do i = 1,2
          puti = p(i, 1)
          if(banned(puti,putj,bant) .or. bannedOrb(puti,1)) cycle
          p1 = p(turn2(i), 1)
    ! hij = <psi_{selectors,i}|H|alpha> 
!          hij = mo_bi_ortho_tc_two_e(p1, p2, h1, h2)
!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!! 
        ! take the transpose of what's written above because later use the complex conjugate 
          hij = mo_bi_ortho_tc_two_e(h1, h2, p1, p2 )
          if (hij /= 0.d0) then
            ! take conjugate to get contribution to <alpha|H|psi> instead of <psi|H|alpha>
!            hij = dconjg(hij) * get_phase_bi(phasemask, 1, 2, h1, p1, h2, p2, N_int)
            hij = hij * get_phase_bi(phasemask, 1, 2, h1, p1, h2, p2, N_int)
            !DIR$ LOOP COUNT AVG(4)
            do k=1,N_states
              mat_r(k, puti, putj) = mat_r(k, puti, putj) + coefs(k,2) * hij
            enddo
          endif
        end do
      end do
      !! <phi|H|alpha>
      do j = 1,2 ! loop over all 4 combinations of one alpha and one beta particle
        putj = p(j, 2)
        if(bannedOrb(putj, 2)) cycle
        p2 = p(turn2(j), 2)
        do i = 1,2
          puti = p(i, 1)
          if(banned(puti,putj,bant) .or. bannedOrb(puti,1)) cycle
          p1 = p(turn2(i), 1)
          hji = mo_bi_ortho_tc_two_e( p1, p2, h1, h2)
          if (hji /= 0.d0) then
            hji = hji * get_phase_bi(phasemask, 1, 2, h1, p1, h2, p2, N_int)
            !DIR$ LOOP COUNT AVG(4)
            do k=1,N_states
              mat_l(k, puti, putj) = mat_l(k, puti, putj) + coefs(k,1) * hji
            enddo
          endif
        end do
      end do
    end if

  else ! if holes are (a,a) or (b,b)
    if(tip == 0) then ! if particles are (a,a,a,a) or (b,b,b,b)
      h1 = h(1, ma)
      h2 = h(2, ma)
      !! <alpha|H|psi>
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
!          hij = mo_bi_ortho_tc_two_e(p1, p2, h1, h2) - mo_bi_ortho_tc_two_e(p2,p1, h1, h2)
!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!! 
        ! take the transpose of what's written above because later use the complex conjugate 
          hij = mo_bi_ortho_tc_two_e(h1, h2, p1, p2) - mo_bi_ortho_tc_two_e(h1, h2, p2,p1 )
          if (hij == 0.d0) cycle

          ! take conjugate to get contribution to <alpha|H|psi> instead of <psi|H|alpha>
!          hij = dconjg(hij) * get_phase_bi(phasemask, ma, ma, h1, p1, h2, p2, N_int)
          hij = hij * get_phase_bi(phasemask, ma, ma, h1, p1, h2, p2, N_int)
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat_r(k, puti, putj) = mat_r(k, puti, putj) +coefs(k,2) * hij
          enddo
        end do
      end do
      !! <phi|H|alpha>
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
          hji = mo_bi_ortho_tc_two_e(p1, p2, h1, h2) - mo_bi_ortho_tc_two_e(p2,p1,h1, h2 )
          if (hji == 0.d0) cycle
          hji = hji * get_phase_bi(phasemask, ma, ma, h1, p1, h2, p2, N_int)
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat_l(k, puti, putj) = mat_l(k, puti, putj) +coefs(k,1) * hji
          enddo
        end do
      end do
    else if(tip == 3) then ! if particles are (a,a,a,b) (ma=1,mi=2) or (a,b,b,b) (ma=2,mi=1)
      h1 = h(1, mi)
      h2 = h(1, ma)
      p1 = p(1, mi)
      !! <alpha|H|psi>
      do i=1,3
        puti = p(turn3(1,i), ma)
        if(bannedOrb(puti,ma)) cycle
        putj = p(turn3(2,i), ma)
        if(bannedOrb(putj,ma)) cycle
        if(banned(puti,putj,1)) cycle
        p2 = p(i, ma)

!        hij = mo_bi_ortho_tc_two_e(p1, p2, h1, h2)
!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!! 
        ! take the transpose of what's written above because later use the complex conjugate 
        hij = mo_bi_ortho_tc_two_e(h1, h2,p1, p2 )
        if (hij == 0.d0) cycle

        ! take conjugate to get contribution to <alpha|H|psi> instead of <psi|H|alpha>
!        hij = dconjg(hij) * get_phase_bi(phasemask, mi, ma, h1, p1, h2, p2, N_int)
        hij = hij * get_phase_bi(phasemask, mi, ma, h1, p1, h2, p2, N_int)
        if (puti < putj) then
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat_r(k, puti, putj) = mat_r(k, puti, putj) + coefs(k,2) * hij
          enddo
        else
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat_r(k, putj, puti) = mat_r(k, putj, puti) + coefs(k,2) * hij
          enddo
        endif
      end do
      !! <phi|H|alpha>
      do i=1,3
        puti = p(turn3(1,i), ma)
        if(bannedOrb(puti,ma)) cycle
        putj = p(turn3(2,i), ma)
        if(bannedOrb(putj,ma)) cycle
        if(banned(puti,putj,1)) cycle
        p2 = p(i, ma)
        hji = mo_bi_ortho_tc_two_e(p1, p2,h1, h2)
        if (hji == 0.d0) cycle
        hji = hji * get_phase_bi(phasemask, mi, ma, h1, p1, h2, p2, N_int)
        if (puti < putj) then
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat_l(k, puti, putj) = mat_l(k, puti, putj) + coefs(k,1) * hji
          enddo
        else
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat_l(k, putj, puti) = mat_l(k, putj, puti) + coefs(k,1) * hji
          enddo
        endif
      end do
    else ! tip == 4  (a,a,b,b)
      puti = p(1, sp)
      putj = p(2, sp)
      if(.not. banned(puti,putj,1)) then
        p1 = p(1, mi)
        p2 = p(2, mi)
        h1 = h(1, mi)
        h2 = h(2, mi)
      !! <alpha|H|psi>
!        hij = (mo_bi_ortho_tc_two_e(p1, p2, h1, h2) - mo_bi_ortho_tc_two_e(p2,p1, h1, h2))
!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!! 
        ! take the transpose of what's written above because later use the complex conjugate 
        hij = (mo_bi_ortho_tc_two_e(h1, h2,p1, p2) - mo_bi_ortho_tc_two_e(h1, h2, p2,p1))
        if (hij /= 0.d0) then
          ! take conjugate to get contribution to <alpha|H|psi> instead of <psi|H|alpha>
!          hij = dconjg(hij) * get_phase_bi(phasemask, mi, mi, h1, p1, h2, p2, N_int)
          hij = hij * get_phase_bi(phasemask, mi, mi, h1, p1, h2, p2, N_int)
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat_r(k, puti, putj) = mat_r(k, puti, putj) + coefs(k,2) * hij
          enddo
        end if
      !! <phi|H|alpha>
        hji = (mo_bi_ortho_tc_two_e(p1, p2,h1, h2) - mo_bi_ortho_tc_two_e( p2,p1, h1, h2))
        if (hji /= 0.d0) then
          hji = hji * get_phase_bi(phasemask, mi, mi, h1, p1, h2, p2, N_int)
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat_l(k, puti, putj) = mat_l(k, puti, putj) + coefs(k,1) * hji
          enddo
        end if
      end if
    end if
  end if
end
