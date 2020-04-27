subroutine get_d1_kpts_new(gen, phasemask, bannedOrb, banned, mat, mask, h, p, sp, coefs)
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

  integer :: kp1,ip1, kp2,ip2, khfix,ihfix, kputi,iputi, kputj,iputj, putj0
  integer :: kpfix, ipfix, puti0
  integer :: kputi1,kputi2,puti01,puti02
  integer :: ii0

  integer, parameter             :: turn2(2) = (/2,1/)
  integer, parameter             :: turn3(2,3) = reshape((/2,3,  1,3, 1,2/), (/2,3/))

  integer                        :: bant
  complex*16, allocatable :: hij_cache(:,:),hij_cache2(:,:)
  complex*16               :: hij, tmp_row(N_states, mo_num), tmp_row2(N_states, mo_num)
  complex*16               :: tmp_row_kpts(N_states, mo_num), tmp_row2_kpts(N_states, mo_num)
  complex*16               :: tmp_row_kpts2(N_states, mo_num_per_kpt), tmp_row2_kpts2(N_states,mo_num_per_kpt)
  complex*16 :: tmp_mat1(N_states,mo_num,mo_num), tmp_mat2(N_states,mo_num,mo_num)
  PROVIDE mo_integrals_map N_int

  allocate (lbanned(mo_num, 2))
  allocate (hij_cache(mo_num,2),hij_cache2(mo_num_per_kpt,2))
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
!    kputi = (puti-1)/mo_num_per_kpt + 1
!    khfix = (hfix-1)/mo_num_per_kpt + 1
!    kp1   =   (p1-1)/mo_num_per_kpt + 1
!    kp2   =   (p2-1)/mo_num_per_kpt + 1
!    iputi = mod(puti-1,mo_num_per_kpt) + 1
!    ihfix = mod(hfix-1,mo_num_per_kpt) + 1
!    ip1   = mod(p1-1,  mo_num_per_kpt) + 1
!    ip2   = mod(p2-1,  mo_num_per_kpt) + 1

    if(.not. bannedOrb(puti, mi)) then
      !==================
      call get_mo_two_e_integrals_complex(hfix,p1,p2,mo_num,hij_cache(1,1),mo_integrals_map,mo_integrals_map_2)
      call get_mo_two_e_integrals_complex(hfix,p2,p1,mo_num,hij_cache(1,2),mo_integrals_map,mo_integrals_map_2)
      !==================
!      call get_mo_two_e_integrals_kpts(hfix,ihfix,khfix,p1,ip1,kp1,p2,ip2,kp2,mo_num_per_kpt,hij_cache2(1,1),mo_integrals_map,mo_integrals_map_2)
!      call get_mo_two_e_integrals_kpts(hfix,ihfix,khfix,p2,ip2,kp2,p1,ip1,kp1,mo_num_per_kpt,hij_cache2(1,2),mo_integrals_map,mo_integrals_map_2)
      tmp_row = (0.d0,0.d0)
!      tmp_row_kpts2 = (0.d0,0.d0)
!      kputj = kconserv(kp1,kp2,khfix)
!      putj0 = (kputj-1)*mo_num_per_kpt
      !==================
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
      !===========================
      ! begin kpts testing
!      do putj = putj0+1, hfix-1
!        iputj = putj-putj0
!        if(lbanned(putj, ma)) cycle
!        if(banned(putj, puti,bant)) cycle
!        hij = hij_cache2(iputj,1) - hij_cache2(iputj,2)
!        if (hij /= (0.d0,0.d0)) then
!          hij = hij * get_phase_bi(phasemask, ma, ma, putj, p1, hfix, p2, N_int)
!          !DIR$ LOOP COUNT AVG(4)
!          do k=1,N_states
!            !tmp_row_kpts(k,putj) = tmp_row_kpts(k,putj) + hij * coefs(k)
!            tmp_row_kpts2(k,iputj) = tmp_row_kpts2(k,iputj) + hij * coefs(k)
!          enddo
!        endif
!      end do
!      do putj = hfix+1,putj0+mo_num_per_kpt
!        iputj = putj - putj0
!        if(lbanned(putj, ma)) cycle
!        if(banned(putj, puti,bant)) cycle
!        hij = hij_cache2(iputj,2) - hij_cache2(iputj,1)
!        if (hij /= (0.d0,0.d0)) then
!          hij = hij * get_phase_bi(phasemask, ma, ma, hfix, p1, putj, p2, N_int)
!          !DIR$ LOOP COUNT AVG(4)
!          do k=1,N_states
!            !tmp_row_kpts(k,putj) = tmp_row_kpts(k,putj) + hij * coefs(k)
!            tmp_row_kpts2(k,iputj) = tmp_row_kpts2(k,iputj) + hij * coefs(k)
!          enddo
!        endif
!      end do
!      ! end kpts testing
!      !===========================================================
!      !print*,'tmp_row_k,tmp_row'
!      !do ii0=1,mo_num
!      !  if (cdabs(tmp_row_kpts(1,ii0)-tmp_row(1,ii0)).gt.1.d-12) then
!      !    print'((A),4(I5),2(2(E25.15),2X))','WarNInG, ',ii0,hfix,p1,p2,tmp_row_kpts(1,ii0),tmp_row(1,ii0)
!      !  endif
!      !enddo
!      !===========================================================
!      tmp_mat1 = (0.d0,0.d0)
!      tmp_mat2 = (0.d0,0.d0)
      !===========================================================
      if(ma == 1) then
      !===========================================================
!        tmp_mat1(1:N_states,1:mo_num,puti) = tmp_mat1(1:N_states,1:mo_num,puti) + tmp_row(1:N_states,1:mo_num)
!        tmp_mat2(1:N_states,putj0+1:putj0+mo_num_per_kpt,puti) = tmp_mat2(1:N_states,putj0+1:putj0+mo_num_per_kpt,puti) + &
!                                         tmp_row_kpts2(1:N_states,1:mo_num_per_kpt)
!      !===========================================================
        mat(1:N_states,1:mo_num,puti) = mat(1:N_states,1:mo_num,puti) + tmp_row(1:N_states,1:mo_num)
!        mat(1:N_states,putj0+1:putj0+mo_num_per_kpt,puti) = mat(1:N_states,putj0+1:putj0+mo_num_per_kpt,puti) + &
!                                         tmp_row_kpts2(1:N_states,1:mo_num_per_kpt)
      else
      !===========================================================
!        do l=1,mo_num
!          !DIR$ LOOP COUNT AVG(4)
!          do k=1,N_states
!            tmp_mat1(k,puti,l) = tmp_mat1(k,puti,l) + tmp_row(k,l)
!          enddo
!        enddo
!        do l=1,mo_num_per_kpt
!          !DIR$ LOOP COUNT AVG(4)
!          do k=1,N_states
!            tmp_mat2(k,puti,l+putj0) = tmp_mat2(k,puti,l+putj0) + tmp_row_kpts2(k,l)
!          enddo
!        enddo
      !===========================================================
        do l=1,mo_num
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat(k,puti,l) = mat(k,puti,l) + tmp_row(k,l)
          enddo
        enddo
        !do l=1,mo_num_per_kpt
        !  !DIR$ LOOP COUNT AVG(4)
        !  do k=1,N_states
        !    mat(k,puti,l+putj0) = mat(k,puti,l+putj0) + tmp_row_kpts2(k,l)
        !  enddo
        !enddo
      end if
      !===========================================================
        !do k=1,N_states
        !  do l=1,mo_num
        !    do ii0=1,mo_num
        !      if (cdabs(tmp_mat2(k,l,ii0)-tmp_mat1(k,l,ii0)).gt.1.d-12) then
        !        print'((A),6(I5),2(2(E25.15),2X))','WarNInG 4a, ',k,l,ii0,hfix,p1,p2,tmp_mat2(k,l,ii0),tmp_mat1(k,l,ii0)
        !    !  else if ((cdabs(tmp_row_kpts(1,ii0))+cdabs(tmp_row(1,ii0))).gt.1.d-12) then
        !    !    print'((A),4(I5),2(2(E25.15),2X))','WarNInG 1b, ',ii0,hfix,pfix,p2,tmp_row_kpts(1,ii0),tmp_row(1,ii0)
        !      endif
        !    enddo
        !  enddo
        !enddo
      !===========================================================
    end if

    !MOVE MI
    pfix = p(1,mi)
!    kpfix = (pfix-1)/mo_num_per_kpt + 1
!    ipfix = mod(pfix-1,mo_num_per_kpt) + 1
    tmp_row = (0.d0,0.d0)
    tmp_row2 = (0.d0,0.d0)
!    !tmp_row_kpts = (0.d0,0.d0)
!    !tmp_row2_kpts = (0.d0,0.d0)
!    tmp_row_kpts2 = (0.d0,0.d0)
!    tmp_row2_kpts2 = (0.d0,0.d0)
    !===========================================================
    call get_mo_two_e_integrals_complex(hfix,pfix,p1,mo_num,hij_cache(1,1),mo_integrals_map,mo_integrals_map_2)
    call get_mo_two_e_integrals_complex(hfix,pfix,p2,mo_num,hij_cache(1,2),mo_integrals_map,mo_integrals_map_2)
    !===========================================================
!    call get_mo_two_e_integrals_kpts(hfix,ihfix,khfix,pfix,ipfix,kpfix,p1,ip1,kp1,mo_num_per_kpt,hij_cache2(1,1),mo_integrals_map,mo_integrals_map_2)
!    call get_mo_two_e_integrals_kpts(hfix,ihfix,khfix,pfix,ipfix,kpfix,p2,ip2,kp2,mo_num_per_kpt,hij_cache2(1,2),mo_integrals_map,mo_integrals_map_2)
    putj = p1
    !============
    !begin ref
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
!    !end ref
!    !===================
!    !begin kpts
!    if (kp1.eq.kp2) then
!    !if (.False.) then
!      kputi1 = kconserv(kpfix,kp1,khfix)
!      kputi2 = kputi1
!      puti01 = (kputi1-1)*mo_num_per_kpt
!      puti02 = puti01
!      do iputi=1,mo_num_per_kpt !HOT
!        puti = puti01 + iputi
!        if(lbanned(puti,mi)) cycle
!        !p1 fixed
!        putj = p1
!        if(.not. banned(putj,puti,bant)) then
!          hij = hij_cache2(iputi,2)
!          if (hij /= (0.d0,0.d0)) then
!            hij = hij * get_phase_bi(phasemask, ma, mi, hfix, p2, puti, pfix, N_int)
!            !DIR$ LOOP COUNT AVG(4)
!            do k=1,N_states
!              tmp_row_kpts2(k,iputi) = tmp_row_kpts2(k,iputi) + hij * coefs(k)
!              !tmp_row_kpts(k,puti) = tmp_row_kpts(k,puti) + hij * coefs(k)
!            enddo
!          endif
!        end if
!!      enddo
!!        
!        putj = p2
!!      do puti=1,mo_num !HOT
!        if(.not. banned(putj,puti,bant)) then
!          hij = hij_cache2(iputi,1)
!          if (hij /= (0.d0,0.d0)) then
!            hij = hij * get_phase_bi(phasemask, ma, mi, hfix, p1, puti, pfix, N_int)
!            do k=1,N_states
!              tmp_row2_kpts2(k,iputi) = tmp_row2_kpts2(k,iputi) + hij * coefs(k)
!              !tmp_row2_kpts(k,puti) = tmp_row2_kpts(k,puti) + hij * coefs(k)
!            enddo
!          endif
!        end if
!      end do
!    else !kp1.ne.kp2
!      kputi2 = kconserv(kpfix,kp2,khfix)
!      puti02 = (kputi2-1)*mo_num_per_kpt
!      putj = p1
!      do iputi=1,mo_num_per_kpt !HOT
!        puti = puti02 + iputi
!        if(lbanned(puti,mi)) cycle
!        !p1 fixed
!        if(.not. banned(putj,puti,bant)) then
!          hij = hij_cache2(iputi,2)
!          if (hij /= (0.d0,0.d0)) then
!            hij = hij * get_phase_bi(phasemask, ma, mi, hfix, p2, puti, pfix, N_int)
!            !DIR$ LOOP COUNT AVG(4)
!            do k=1,N_states
!              tmp_row_kpts2(k,iputi) = tmp_row_kpts2(k,iputi) + hij * coefs(k)
!              !tmp_row_kpts(k,puti) = tmp_row_kpts(k,puti) + hij * coefs(k)
!            enddo
!          endif
!        end if
!      enddo
!!        
!      putj = p2
!      kputi1 = kconserv(kpfix,kp1,khfix)
!      puti01 = (kputi1-1)*mo_num_per_kpt
!      do iputi=1,mo_num_per_kpt !HOT
!        puti = puti01 + iputi
!        if(lbanned(puti,mi)) cycle
!        if(.not. banned(putj,puti,bant)) then
!          hij = hij_cache2(iputi,1)
!          if (hij /= (0.d0,0.d0)) then
!            hij = hij * get_phase_bi(phasemask, ma, mi, hfix, p1, puti, pfix, N_int)
!            do k=1,N_states
!              tmp_row2_kpts2(k,iputi) = tmp_row2_kpts2(k,iputi) + hij * coefs(k)
!              !tmp_row2_kpts(k,puti) = tmp_row2_kpts(k,puti) + hij * coefs(k)
!            enddo
!          endif
!        end if
!      end do
!    endif
!    !end kpts
!    !===================
!    !test printing
!    !print'((A),5(I5))','kpt info1: ',kconserv(kpfix,kp2,khfix),khfix,kpfix,kp2,kputi2
!    !print'((A),5(I5))','kpt info2: ',kconserv(kpfix,kp1,khfix),khfix,kpfix,kp1,kputi1
!    !do ii0=1,mo_num
!    !  if (cdabs(tmp_row_kpts(1,ii0)-tmp_row(1,ii0)).gt.1.d-12) then
!    !    print'((A),4(I5),2(2(E25.15),2X))','WarNInG 1a, ',ii0,hfix,pfix,p2,tmp_row_kpts(1,ii0),tmp_row(1,ii0)
!    !!  else if ((cdabs(tmp_row_kpts(1,ii0))+cdabs(tmp_row(1,ii0))).gt.1.d-12) then
!    !!    print'((A),4(I5),2(2(E25.15),2X))','WarNInG 1b, ',ii0,hfix,pfix,p2,tmp_row_kpts(1,ii0),tmp_row(1,ii0)
!    !  endif
!    !  if (cdabs(tmp_row2_kpts(1,ii0)-tmp_row2(1,ii0)).gt.1.d-12) then
!    !    print'((A),4(I5),2(2(E25.15),2X))','WarNInG 2a, ',ii0,hfix,pfix,p1,tmp_row2_kpts(1,ii0),tmp_row2(1,ii0)
!    !!  else if ((cdabs(tmp_row2_kpts(1,ii0))+cdabs(tmp_row2(1,ii0))).gt.1.d-12) then
!    !!    print'((A),4(I5),2(2(E25.15),2X))','WarNInG 2b, ',ii0,hfix,pfix,p1,tmp_row2_kpts(1,ii0),tmp_row2(1,ii0)
!    !  endif
!    !enddo
!    !===================
!
!    tmp_mat1 = (0.d0,0.d0)
!    tmp_mat2 = (0.d0,0.d0)
    if(mi == 1) then
!    !===================
!      tmp_mat1(:,:,p1) = tmp_mat1(:,:,p1) + tmp_row(:,:)
!      tmp_mat1(:,:,p2) = tmp_mat1(:,:,p2) + tmp_row2(:,:)
!      tmp_mat2(:,puti02+1:puti02+mo_num_per_kpt,p1) = tmp_mat2(:,puti02+1:puti02+mo_num_per_kpt,p1) + tmp_row_kpts2(:,:)
!      tmp_mat2(:,puti01+1:puti01+mo_num_per_kpt,p2) = tmp_mat2(:,puti01+1:puti01+mo_num_per_kpt,p2) + tmp_row2_kpts2(:,:)
!    !===================
      mat(:,:,p1) = mat(:,:,p1) + tmp_row(:,:)
      mat(:,:,p2) = mat(:,:,p2) + tmp_row2(:,:)
!      mat(:,puti02+1:puti02+mo_num_per_kpt,p1) = mat(:,puti02+1:puti02+mo_num_per_kpt,p1) + tmp_row_kpts2(:,:)
!      mat(:,puti01+1:puti01+mo_num_per_kpt,p2) = mat(:,puti01+1:puti01+mo_num_per_kpt,p2) + tmp_row2_kpts2(:,:)
    else
    !===================
!      do l=1,mo_num
!        !DIR$ LOOP COUNT AVG(4)
!        do k=1,N_states
!          tmp_mat1(k,p1,l) = tmp_mat1(k,p1,l) + tmp_row(k,l)
!          tmp_mat1(k,p2,l) = tmp_mat1(k,p2,l) + tmp_row2(k,l)
!        enddo
!      enddo
!      do l=1,mo_num_per_kpt
!        !DIR$ LOOP COUNT AVG(4)
!        do k=1,N_states
!          tmp_mat2(k,p1,l+puti02) = tmp_mat2(k,p1,l+puti02) + tmp_row_kpts2(k,l)
!          tmp_mat2(k,p2,l+puti01) = tmp_mat2(k,p2,l+puti01) + tmp_row2_kpts2(k,l)
!        enddo
!      enddo
    !===================
      do l=1,mo_num
        !DIR$ LOOP COUNT AVG(4)
        do k=1,N_states
          mat(k,p1,l) = mat(k,p1,l) + tmp_row(k,l)
          mat(k,p2,l) = mat(k,p2,l) + tmp_row2(k,l)
        enddo
      enddo
!      do l=1,mo_num_per_kpt
!        !DIR$ LOOP COUNT AVG(4)
!        do k=1,N_states
!          mat(k,p1,l+puti02) = mat(k,p1,l+puti02) + tmp_row_kpts2(k,l)
!          mat(k,p2,l+puti01) = mat(k,p2,l+puti01) + tmp_row2_kpts2(k,l)
!        enddo
!      enddo
    end if
      !===========================================================
!        do k=1,N_states
!          do l=1,mo_num
!            do ii0=1,mo_num
!              if (cdabs(tmp_mat2(k,l,ii0)-tmp_mat1(k,l,ii0)).gt.1.d-12) then
!                print'((A),7(I5),2(2(E25.15),2X))','WarNInG 5a, ',k,l,ii0,hfix,pfix,p1,p2,tmp_mat2(k,l,ii0),tmp_mat1(k,l,ii0)
!            !  else if ((cdabs(tmp_row_kpts(1,ii0))+cdabs(tmp_row(1,ii0))).gt.1.d-12) then
!            !    print'((A),4(I5),2(2(E25.15),2X))','WarNInG 1b, ',ii0,hfix,pfix,p2,tmp_row_kpts(1,ii0),tmp_row(1,ii0)
!              endif
!            enddo
!          enddo
!        enddo
      !===========================================================
    !todo: kpts okay up to this point in get_d1_complex

  else  ! sp /= 3

    if(p(0,ma) == 3) then
      do i=1,3
        hfix = h(1,ma)
        puti = p(i, ma)
        p1 = p(turn3(1,i), ma)
        p2 = p(turn3(2,i), ma)
!        kputi = (puti-1)/mo_num_per_kpt + 1
!        khfix = (hfix-1)/mo_num_per_kpt + 1
!        kp1   =   (p1-1)/mo_num_per_kpt + 1
!        kp2   =   (p2-1)/mo_num_per_kpt + 1
!        iputi = mod(puti-1,mo_num_per_kpt) + 1
!        ihfix = mod(hfix-1,mo_num_per_kpt) + 1
!        ip1   = mod(p1-1,  mo_num_per_kpt) + 1
!        ip2   = mod(p2-1,  mo_num_per_kpt) + 1
        call get_mo_two_e_integrals_complex(hfix,p1,p2,mo_num,hij_cache(1,1),mo_integrals_map,mo_integrals_map_2)
        call get_mo_two_e_integrals_complex(hfix,p2,p1,mo_num,hij_cache(1,2),mo_integrals_map,mo_integrals_map_2)
!        call get_mo_two_e_integrals_kpts(hfix,ihfix,khfix,p1,ip1,kp1,p2,ip2,kp2,mo_num_per_kpt,hij_cache2(1,1),mo_integrals_map,mo_integrals_map_2)
!        call get_mo_two_e_integrals_kpts(hfix,ihfix,khfix,p2,ip2,kp2,p1,ip1,kp1,mo_num_per_kpt,hij_cache2(1,2),mo_integrals_map,mo_integrals_map_2)
        tmp_row = (0.d0,0.d0)
        !tmp_row_kpts = (0.d0,0.d0)
!        tmp_row_kpts2 = (0.d0,0.d0)
        !===================
        !begin ref
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
        !end ref
        !=================
        !begin kpts
!        kputj = kconserv(kp1,kp2,khfix)
!        putj0 = (kputj-1)*mo_num_per_kpt
!        do putj = putj0+1,hfix-1
!          iputj = putj - putj0
!          if(banned(putj,puti,1)) cycle
!          if(lbanned(putj,ma)) cycle
!          hij = hij_cache2(iputj,1) - hij_cache2(iputj,2)
!          if (hij /= (0.d0,0.d0)) then
!            hij = hij * get_phase_bi(phasemask, ma, ma, putj, p1, hfix, p2, N_int)
!            !tmp_row_kpts(:,putj) = tmp_row_kpts(:,putj) + hij * coefs(:)
!            tmp_row_kpts2(:,iputj) = tmp_row_kpts2(:,iputj) + hij * coefs(:)
!          endif
!        end do
!        do putj=hfix+1,putj0+mo_num_per_kpt
!          iputj = putj - putj0
!          if(banned(putj,puti,1)) cycle
!          if(lbanned(putj,ma)) cycle
!          hij = hij_cache2(iputj,2) - hij_cache2(iputj,1)
!          if (hij /= (0.d0,0.d0)) then
!            hij = hij * get_phase_bi(phasemask, ma, ma, hfix, p1, putj, p2, N_int)
!            !tmp_row_kpts(:,putj) = tmp_row_kpts(:,putj) + hij * coefs(:)
!            tmp_row_kpts2(:,iputj) = tmp_row_kpts2(:,iputj) + hij * coefs(:)
!          endif
!        end do
!
!        !end kpts
!    !do ii0=1,mo_num
!    !  if (cdabs(tmp_row_kpts(1,ii0)-tmp_row(1,ii0)).gt.1.d-12) then
!    !    print'((A),4(I5),2(2(E25.15),2X))','WarNInG 1a, ',ii0,hfix,pfix,p2,tmp_row_kpts(1,ii0),tmp_row(1,ii0)
!    !!  else if ((cdabs(tmp_row_kpts(1,ii0))+cdabs(tmp_row(1,ii0))).gt.1.d-12) then
!    !!    print'((A),4(I5),2(2(E25.15),2X))','WarNInG 1b, ',ii0,hfix,pfix,p2,tmp_row_kpts(1,ii0),tmp_row(1,ii0)
!    !  endif
!    !enddo
!        !=================
!        tmp_mat1 = (0.d0,0.d0)
!        tmp_mat2 = (0.d0,0.d0)
!        tmp_mat1(:, :puti-1, puti) = tmp_mat1(:, :puti-1, puti) + tmp_row(:,:puti-1)
!        do l=puti,mo_num
!          !DIR$ LOOP COUNT AVG(4)
!          do k=1,N_states
!            tmp_mat1(k, puti, l) = tmp_mat1(k, puti,l) + tmp_row(k,l)
!          enddo
!        enddo
!        !=================
!        if (kputj.lt.kputi) then
!          tmp_mat2(1:N_states,putj0+1:putj0+mo_num_per_kpt,puti) =  &
!                  tmp_mat2(1:N_states,putj0+1:putj0+mo_num_per_kpt,puti) + &
!                  tmp_row_kpts2(1:N_states,1:mo_num_per_kpt)
!        else if (kputj.gt.kputi) then
!          do l=1,mo_num_per_kpt
!            !DIR$ LOOP COUNT AVG(4)
!            do k=1,N_states
!              tmp_mat2(k, puti, l+putj0) = tmp_mat2(k, puti,l+putj0) + tmp_row_kpts2(k,l)
!            enddo
!          enddo
!        else !kputj == kputi
!          tmp_mat2(1:N_states,putj0+1:puti-1,puti) =  &
!                  tmp_mat2(1:N_states,putj0+1:puti-1,puti) + &
!                  tmp_row_kpts2(1:N_states,1:iputi-1)
!          do l=iputi,mo_num_per_kpt
!            !DIR$ LOOP COUNT AVG(4)
!            do k=1,N_states
!              tmp_mat2(k, puti, l+putj0) = tmp_mat2(k, puti,l+putj0) + tmp_row_kpts2(k,l)
!            enddo
!          enddo
!        endif
!        !=================
!        do k=1,N_states
!          do l=1,mo_num
!            do ii0=1,mo_num
!              if (cdabs(tmp_mat2(k,l,ii0)-tmp_mat1(k,l,ii0)).gt.1.d-12) then
!                print'((A),6(I5),2(2(E25.15),2X))','WarNInG 3a, ',k,l,ii0,hfix,p1,p2,tmp_mat2(k,l,ii0),tmp_mat1(k,l,ii0)
!            !  else if ((cdabs(tmp_row_kpts(1,ii0))+cdabs(tmp_row(1,ii0))).gt.1.d-12) then
!            !    print'((A),4(I5),2(2(E25.15),2X))','WarNInG 1b, ',ii0,hfix,pfix,p2,tmp_row_kpts(1,ii0),tmp_row(1,ii0)
!              endif
!            enddo
!          enddo
!        enddo

        !=================
        mat(:, :puti-1, puti) = mat(:, :puti-1, puti) + tmp_row(:,:puti-1)
        do l=puti,mo_num
          !DIR$ LOOP COUNT AVG(4)
          do k=1,N_states
            mat(k, puti, l) = mat(k, puti,l) + tmp_row(k,l)
          enddo
        enddo
        !!=================
        !!todo: check for iputi=1,2
        !if (kputj.lt.kputi) then
        !  mat(1:N_states,putj0+1:putj0+mo_num_per_kpt,puti) =  &
        !          mat(1:N_states,putj0+1:putj0+mo_num_per_kpt,puti) + &
        !          tmp_row_kpts2(1:N_states,1:mo_num_per_kpt)
        !else if (kputj.gt.kputi) then
        !  do l=1,mo_num_per_kpt
        !    !DIR$ LOOP COUNT AVG(4)
        !    do k=1,N_states
        !      mat(k, puti, l+putj0) = mat(k, puti,l+putj0) + tmp_row_kpts2(k,l)
        !    enddo
        !  enddo
        !else !kputj == kputi
        !  mat(1:N_states,putj0+1:puti-1,puti) =  &
        !          mat(1:N_states,putj0+1:puti-1,puti) + &
        !          tmp_row_kpts2(1:N_states,1:iputi-1)
        !  do l=iputi,mo_num_per_kpt
        !    !DIR$ LOOP COUNT AVG(4)
        !    do k=1,N_states
        !      mat(k, puti, l+putj0) = mat(k, puti,l+putj0) + tmp_row_kpts2(k,l)
        !    enddo
        !  enddo
        !endif
      end do
    else
      hfix = h(1,mi)
      pfix = p(1,mi)
      p1 = p(1,ma)
      p2 = p(2,ma)
!      kpfix = (pfix-1)/mo_num_per_kpt + 1
!      khfix = (hfix-1)/mo_num_per_kpt + 1
!      kp1   =   (p1-1)/mo_num_per_kpt + 1
!      kp2   =   (p2-1)/mo_num_per_kpt + 1
!      ipfix = mod(pfix-1,mo_num_per_kpt) + 1
!      ihfix = mod(hfix-1,mo_num_per_kpt) + 1
!      ip1   = mod(p1-1,  mo_num_per_kpt) + 1
!      ip2   = mod(p2-1,  mo_num_per_kpt) + 1
      tmp_row = (0.d0,0.d0)
      tmp_row2 = (0.d0,0.d0)
      !tmp_row_kpts = (0.d0,0.d0)
      !tmp_row2_kpts = (0.d0,0.d0)
      call get_mo_two_e_integrals_complex(hfix,p1,pfix,mo_num,hij_cache(1,1),mo_integrals_map,mo_integrals_map_2)
      call get_mo_two_e_integrals_complex(hfix,p2,pfix,mo_num,hij_cache(1,2),mo_integrals_map,mo_integrals_map_2)
      !call get_mo_two_e_integrals_kpts(hfix,ihfix,khfix,p1,ip1,kp1,pfix,ipfix,kpfix,mo_num_per_kpt,hij_cache2(1,1),mo_integrals_map,mo_integrals_map_2)
      !call get_mo_two_e_integrals_kpts(hfix,ihfix,khfix,p2,ip2,kp2,pfix,ipfix,kpfix,mo_num_per_kpt,hij_cache2(1,2),mo_integrals_map,mo_integrals_map_2)
      !===============
      !begin ref
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
      !end ref
      !===============
      !begin kpts
      !todo: combine if kp1==kp2
 !     putj = p2
 !     kputi1 = kconserv(kp1,kpfix,khfix)
 !     puti01 = (kputi1-1)*mo_num_per_kpt
 !     do iputi=1,mo_num_per_kpt
 !       puti = puti01 + iputi
 !       if(lbanned(puti,ma)) cycle
 !       if(.not. banned(puti,putj,1)) then
 !         hij = hij_cache2(iputi,1)
 !         if (hij /= (0.d0,0.d0)) then
 !           hij = hij * get_phase_bi(phasemask, mi, ma, hfix, pfix, puti, p1, N_int)
 !           !DIR$ LOOP COUNT AVG(4)
 !           do k=1,N_states
 !             tmp_row_kpts(k,puti) = tmp_row_kpts(k,puti) + hij * coefs(k)
 !           enddo
 !         endif
 !       end if
 !     enddo
 !     putj = p1
 !     kputi2 = kconserv(kp2,kpfix,khfix)
 !     puti02 = (kputi2-1)*mo_num_per_kpt
 !     do iputi=1,mo_num_per_kpt
 !       puti = puti02 + iputi
 !       if(lbanned(puti,ma)) cycle
 !       if(.not. banned(puti,putj,1)) then
 !         hij = hij_cache2(iputi,2)
 !         if (hij /= (0.d0,0.d0)) then
 !           hij = hij * get_phase_bi(phasemask, mi, ma, hfix, pfix, puti, p2, N_int)
 !           do k=1,N_states
 !             tmp_row2_kpts(k,puti) = tmp_row2_kpts(k,puti) + hij * coefs(k)
 !           enddo
 !         endif
 !       end if
 !     end do
 !     !end kpts
 !     !===============
 !   !test printing
 !   !print'((A),5(I5))','kpt info1: ',kconserv(kpfix,kp2,khfix),khfix,kpfix,kp2,kputi2
 !   !print'((A),5(I5))','kpt info2: ',kconserv(kpfix,kp1,khfix),khfix,kpfix,kp1,kputi1
 !   do ii0=1,mo_num
 !     if (cdabs(tmp_row_kpts(1,ii0)-tmp_row(1,ii0)).gt.1.d-12) then
 !       print'((A),4(I5),2(2(E25.15),2X))','WarNInG 1a, ',ii0,hfix,p1,pfix,tmp_row_kpts(1,ii0),tmp_row(1,ii0)
 !   !  else if ((cdabs(tmp_row_kpts(1,ii0))+cdabs(tmp_row(1,ii0))).gt.1.d-12) then
 !   !    print'((A),4(I5),2(2(E25.15),2X))','WarNInG 1b, ',ii0,hfix,pfix,p2,tmp_row_kpts(1,ii0),tmp_row(1,ii0)
 !     endif
 !     if (cdabs(tmp_row2_kpts(1,ii0)-tmp_row2(1,ii0)).gt.1.d-12) then
 !       print'((A),4(I5),2(2(E25.15),2X))','WarNInG 2a, ',ii0,hfix,p2,pfix,tmp_row2_kpts(1,ii0),tmp_row2(1,ii0)
 !   !  else if ((cdabs(tmp_row2_kpts(1,ii0))+cdabs(tmp_row2(1,ii0))).gt.1.d-12) then
 !   !    print'((A),4(I5),2(2(E25.15),2X))','WarNInG 2b, ',ii0,hfix,pfix,p1,tmp_row2_kpts(1,ii0),tmp_row2(1,ii0)
 !     endif
 !   enddo
    !===================
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


