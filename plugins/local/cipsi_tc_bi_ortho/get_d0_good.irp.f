subroutine get_d0_new(gen, phasemask, bannedOrb, banned, mat_l, mat_r, mask, h, p, sp, coefs)
  !todo: indices/conjg should be okay for complex
  use bitmasks
  implicit none

  integer(bit_kind), intent(in) :: gen(N_int, 2), mask(N_int, 2)
  integer(bit_kind), intent(in) :: phasemask(N_int,2)
  logical, intent(in) :: bannedOrb(mo_num, 2), banned(mo_num, mo_num,2)
  integer(bit_kind) :: det(N_int, 2)
  double precision, intent(in) :: coefs(N_states,2)
  double precision, intent(inout) :: mat_l(N_states, mo_num, mo_num)
  double precision, intent(inout) :: mat_r(N_states, mo_num, mo_num)
  integer, intent(in) :: h(0:2,2), p(0:4,2), sp

  integer :: i, j, k, s, h1, h2, p1, p2, puti, putj, mm
  double precision :: phase
  double precision :: hij,hji
  double precision, external :: get_phase_bi
  logical :: ok

  integer, parameter :: bant=1
  double precision, allocatable :: hij_cache1(:), hij_cache2(:)
  allocate (hij_cache1(mo_num),hij_cache2(mo_num))
  double precision, allocatable :: hji_cache1(:), hji_cache2(:)
  allocate (hji_cache1(mo_num),hji_cache2(mo_num))
!  print*,'in get_d0_new'
!  call debug_det(gen,N_int)
!  print*,'coefs',coefs(1,:)

  if(sp == 3) then ! AB
    h1 = p(1,1)
    h2 = p(1,2)
    do p1=1, mo_num
      if(bannedOrb(p1, 1)) cycle
!      call get_mo_two_e_integrals_complex(p1,h2,h1,mo_num,hij_cache1,mo_integrals_map)
      do mm = 1, mo_num
       hij_cache1(mm) = mo_bi_ortho_tc_two_e(mm,p1,h2,h1)
       hji_cache1(mm) = mo_bi_ortho_tc_two_e(h2,h1,mm,p1)
      enddo
      !!!!!!!!!! <alpha|H|psi>
      do p2=1, mo_num
        if(bannedOrb(p2,2)) cycle
        if(banned(p1, p2, bant)) cycle ! rentable?
        if(p1 == h1 .or. p2 == h2) then
          call apply_particles(mask, 1,p1,2,p2, det, ok, N_int)
          ! call i_h_j_complex(gen, det, N_int, hij) ! need to take conjugate of this
!          call i_h_j_complex(det, gen, N_int, hij)
          call htilde_mu_mat_opt_bi_ortho_no_3e(det,gen,N_int, hij)
        else
          phase = get_phase_bi(phasemask, 1, 2, h1, p1, h2, p2, N_int)
          hij = hij_cache1(p2) * phase
        end if
        if (hij == (0.d0,0.d0)) cycle
        !DIR$ LOOP COUNT AVG(4)
        do k=1,N_states
          mat_r(k, p1, p2) = mat_r(k, p1, p2) + coefs(k,2) * hij  ! HOTSPOT
        enddo
      end do
      !!!!!!!!!! <phi|H|alpha>
      do p2=1, mo_num
        if(bannedOrb(p2,2)) cycle
        if(banned(p1, p2, bant)) cycle ! rentable?
        if(p1 == h1 .or. p2 == h2) then
          call apply_particles(mask, 1,p1,2,p2, det, ok, N_int)
          ! call i_h_j_complex(gen, det, N_int, hij) ! need to take conjugate of this
!          call i_h_j_complex(det, gen, N_int, hij)
          call htilde_mu_mat_opt_bi_ortho_no_3e(gen,det,N_int, hji)
        else
          phase = get_phase_bi(phasemask, 1, 2, h1, p1, h2, p2, N_int)
          hji = hji_cache1(p2) * phase
        end if
        if (hji == (0.d0,0.d0)) cycle
        !DIR$ LOOP COUNT AVG(4)
        do k=1,N_states
          mat_l(k, p1, p2) = mat_l(k, p1, p2) + coefs(k,1) * hji  ! HOTSPOT
        enddo
      end do
    end do

  else ! AA BB
    p1 = p(1,sp)
    p2 = p(2,sp)
    do puti=1, mo_num
      if(bannedOrb(puti, sp)) cycle
!      call get_mo_two_e_integrals_complex(puti,p2,p1,mo_num,hij_cache1,mo_integrals_map,mo_integrals_map_2)
!      call get_mo_two_e_integrals_complex(puti,p1,p2,mo_num,hij_cache2,mo_integrals_map,mo_integrals_map_2)
      do mm = 1, mo_num
       hij_cache1(mm) = mo_bi_ortho_tc_two_e(mm,puti,p2,p1)
       hij_cache2(mm) = mo_bi_ortho_tc_two_e(mm,puti,p1,p2)
       hji_cache1(mm) = mo_bi_ortho_tc_two_e(p2,p1,mm,puti)
       hji_cache2(mm) = mo_bi_ortho_tc_two_e(p1,p2,mm,puti)
      enddo
      !!!!!!!!!! <alpha|H|psi>
      do putj=puti+1, mo_num
        if(bannedOrb(putj, sp)) cycle
        if(banned(puti, putj, bant)) cycle ! rentable?
        if(puti == p1 .or. putj == p2 .or. puti == p2 .or. putj == p1) then
          call apply_particles(mask, sp,puti,sp,putj, det, ok, N_int)
          !call i_h_j_complex(gen, det, N_int, hij) ! need to take conjugate of this
!          call i_h_j_complex(det, gen, N_int, hij)
          call htilde_mu_mat_opt_bi_ortho_no_3e(det,gen,N_int, hij)
          if (hij == 0.d0) cycle
        else
!          hij = (mo_two_e_integral_complex(p1, p2, puti, putj) -  mo_two_e_integral_complex(p2, p1, puti, putj))
!          hij = (mo_bi_ortho_tc_two_e(p1, p2, puti, putj) -  mo_bi_ortho_tc_two_e(p2, p1, puti, putj))
          hij = (mo_bi_ortho_tc_two_e(puti, putj, p1, p2) -  mo_bi_ortho_tc_two_e(puti, putj, p2, p1))
          if (hij == 0.d0) cycle
          hij = (hij) * get_phase_bi(phasemask, sp, sp, puti, p1 , putj, p2, N_int)
        end if
        !DIR$ LOOP COUNT AVG(4)
        do k=1,N_states
          mat_r(k, puti, putj) = mat_r(k, puti, putj) + coefs(k,2) * hij
        enddo
      end do

      !!!!!!!!!! <phi|H|alpha>
      do putj=puti+1, mo_num
        if(bannedOrb(putj, sp)) cycle
        if(banned(puti, putj, bant)) cycle ! rentable?
        if(puti == p1 .or. putj == p2 .or. puti == p2 .or. putj == p1) then
          call apply_particles(mask, sp,puti,sp,putj, det, ok, N_int)
          call htilde_mu_mat_opt_bi_ortho_no_3e(gen,det,N_int, hji)
          if (hji == 0.d0) cycle
        else
          hji = (mo_bi_ortho_tc_two_e( p1, p2, puti, putj) -  mo_bi_ortho_tc_two_e( p2, p1, puti, putj))
          if (hji == 0.d0) cycle
          hji = (hji) * get_phase_bi(phasemask, sp, sp, puti, p1 , putj, p2, N_int)
        end if
        !DIR$ LOOP COUNT AVG(4)
        do k=1,N_states
          mat_l(k, puti, putj) = mat_l(k, puti, putj) + coefs(k,1) * hji
        enddo
      end do
    end do
  end if

  deallocate(hij_cache1,hij_cache2)
end

