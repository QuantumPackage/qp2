subroutine get_excitation_degree(key1,key2,degree,Nint)
  use bitmasks
  include 'utils/constants.include.F'
  implicit none
  BEGIN_DOC
  ! This function calculates the excitation degree between two
  ! determinants, which is half the number of bits that are different between the two
  ! determinants. The function takes four arguments: 
  !
  !  * key1: An integer array of length Nint*2, representing the first determinant.
  !
  !  * key2: An integer array of length Nint*2, representing the second determinant.
  !
  !  * degree: An integer, passed by reference, that will store the calculated excitation degree.
  !
  !  * Nint: An integer representing the number of integers in each of the key1 and key2 arrays.
  ! 
  ! It starts a select case block that depends on the value of Nint. 
  ! In each case, the function first calculates the bitwise XOR of each
  ! corresponding pair of elements in key1 and key2, storing the results in the
  ! xorvec array. It then calculates the number of bits set (using the popcnt
  ! function) for each element in xorvec, and sums these counts up. This sum is
  ! stored in the degree variable.  
  ! Finally, the degree variable is right-shifted by 1 bit to divide the result by 2.
  END_DOC

  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key1(Nint*2)
  integer(bit_kind), intent(in)  :: key2(Nint*2)
  integer, intent(out)           :: degree

  integer(bit_kind)              :: xorvec(2*N_int_max)
  integer                        :: l

  ASSERT (Nint > 0)

  select case (Nint)

    case (1)
      xorvec(1) = xor( key1(1), key2(1))
      xorvec(2) = xor( key1(2), key2(2))
      degree = popcnt(xorvec(1))+popcnt(xorvec(2))

    case (2)
      xorvec(1) = xor( key1(1), key2(1))
      xorvec(2) = xor( key1(2), key2(2))
      xorvec(3) = xor( key1(3), key2(3))
      xorvec(4) = xor( key1(4), key2(4))
      degree = sum(popcnt(xorvec(1:4)))

    case (3)
      do l=1,6
        xorvec(l) = xor( key1(l), key2(l))
      enddo
      degree = sum(popcnt(xorvec(1:6)))

    case (4)
      do l=1,8
        xorvec(l) = xor( key1(l), key2(l))
      enddo
      degree = sum(popcnt(xorvec(1:8)))

    case default
      integer :: lmax
      lmax = shiftl(Nint,1)
      do l=1,lmax
        xorvec(l) = xor( key1(l), key2(l))
      enddo
      degree = sum(popcnt(xorvec(1:lmax)))

  end select

  degree = shiftr(degree,1)

end


subroutine get_excitation(det1,det2,exc,degree,phase,Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns the excitation operators between two determinants and the phase.
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: det1(Nint,2)
  integer(bit_kind), intent(in)  :: det2(Nint,2)
  integer, intent(out)           :: exc(0:2,2,2)
  integer, intent(out)           :: degree
  double precision, intent(out)  :: phase
  ! exc(number,hole/particle,spin)
  ! ex :
  ! exc(0,1,1) = number of holes alpha
  ! exc(0,2,1) = number of particle alpha
  ! exc(0,2,2) = number of particle beta
  ! exc(1,2,1) = first particle alpha
  ! exc(1,1,1) = first hole     alpha
  ! exc(1,2,2) = first particle beta
  ! exc(1,1,2) = first hole     beta
  ! E_pq : T^alpha_pq + T^beta_pq
  ! T^alpha_pq  : exc(0,1,1) = 1
  !               exc(0,2,1) = 1
  !               exc(1,1,1) = q
  !               exc(1,2,1) = p

  ! T^beta_pq   : exc(0,1,2) = 1
  !               exc(0,2,2) = 1
  !               exc(1,1,2) = q
  !               exc(1,2,2) = p

  ASSERT (Nint > 0)

  !DIR$ FORCEINLINE
  call get_excitation_degree(det1,det2,degree,Nint)
  select case (degree)

    case (3:)
      degree = -1
      return

    case (2)
      call get_double_excitation(det1,det2,exc,phase,Nint)
      return

    case (1)
      call get_single_excitation(det1,det2,exc,phase,Nint)
      return

    case(0)
      ! Avoid uninitialized phase
      phase = 1d0 
      return

  end select
end

subroutine decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Decodes the exc arrays returned by get_excitation.
  ! h1,h2 : Holes
  ! p1,p2 : Particles
  ! s1,s2 : Spins (1:alpha, 2:beta)
  ! degree : Degree of excitation
  END_DOC
  integer, intent(in)            :: exc(0:2,2,2),degree
  integer, intent(out)           :: h1,h2,p1,p2,s1,s2
  ASSERT (degree > 0)
  ASSERT (degree < 3)

  select case(degree)
    case(2)
      if (exc(0,1,1) == 2) then
        h1 = exc(1,1,1)
        h2 = exc(2,1,1)
        p1 = exc(1,2,1)
        p2 = exc(2,2,1)
        s1 = 1
        s2 = 1
      else if (exc(0,1,2) == 2) then
        h1 = exc(1,1,2)
        h2 = exc(2,1,2)
        p1 = exc(1,2,2)
        p2 = exc(2,2,2)
        s1 = 2
        s2 = 2
      else
        h1 = exc(1,1,1)
        h2 = exc(1,1,2)
        p1 = exc(1,2,1)
        p2 = exc(1,2,2)
        s1 = 1
        s2 = 2
      endif
    case(1)
      if (exc(0,1,1) == 1) then
        h1 = exc(1,1,1)
        h2 = 0
        p1 = exc(1,2,1)
        p2 = 0
        s1 = 1
        s2 = 0
      else
        h1 = exc(1,1,2)
        h2 = 0
        p1 = exc(1,2,2)
        p2 = 0
        s1 = 2
        s2 = 0
      endif
    case(0)
      h1 = 0
      p1 = 0
      h2 = 0
      p2 = 0
      s1 = 0
      s2 = 0
  end select
end

subroutine get_double_excitation(det1,det2,exc,phase,Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns the two excitation operators between two doubly excited determinants and the phase.
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: det1(Nint,2)
  integer(bit_kind), intent(in)  :: det2(Nint,2)
  integer, intent(out)           :: exc(0:2,2,2)
  double precision, intent(out)  :: phase
  integer                        :: tz
  integer                        :: l, ispin, idx_hole, idx_particle, ishift
  integer                        :: nperm
  integer                        :: i,j,k,m,n
  integer                        :: high, low
  integer                        :: a,b,c,d
  integer(bit_kind)              :: hole, particle, tmp
  double precision, parameter    :: phase_dble(0:1) = (/ 1.d0, -1.d0 /)

  ASSERT (Nint > 0)
  nperm = 0
  exc(0,1,1) = 0
  exc(0,2,1) = 0
  exc(0,1,2) = 0
  exc(0,2,2) = 0
  do ispin = 1,2
    idx_particle = 0
    idx_hole = 0
    ishift = 1-bit_kind_size
    do l=1,Nint
      ishift = ishift + bit_kind_size
      if (det1(l,ispin) == det2(l,ispin)) then
        cycle
      endif
      tmp = xor( det1(l,ispin), det2(l,ispin) )
      particle = iand(tmp, det2(l,ispin))
      hole     = iand(tmp, det1(l,ispin))
      do while (particle /= 0_bit_kind)
        tz = trailz(particle)
        idx_particle = idx_particle + 1
        exc(0,2,ispin) = exc(0,2,ispin) + 1
        exc(idx_particle,2,ispin) = tz+ishift
        particle = iand(particle,particle-1_bit_kind)
      enddo
      if (iand(exc(0,1,ispin),exc(0,2,ispin))==2) then  ! exc(0,1,ispin)==2 or exc(0,2,ispin)==2
        exit
      endif
      do while (hole /= 0_bit_kind)
        tz = trailz(hole)
        idx_hole = idx_hole + 1
        exc(0,1,ispin) = exc(0,1,ispin) + 1
        exc(idx_hole,1,ispin) = tz+ishift
        hole = iand(hole,hole-1_bit_kind)
      enddo
      if (iand(exc(0,1,ispin),exc(0,2,ispin))==2) then ! exc(0,1,ispin)==2 or exc(0,2,ispin)
        exit
      endif
    enddo

    select case (exc(0,1,ispin))
      case(0)
        cycle

      case(1)

        high = max(exc(1,1,ispin), exc(1,2,ispin))-1
        low  = min(exc(1,1,ispin), exc(1,2,ispin))

        ASSERT (low >= 0)
        ASSERT (high > 0)

        k = shiftr(high,bit_kind_shift)+1
        j = shiftr(low,bit_kind_shift)+1
        m = iand(high,bit_kind_size-1)
        n = iand(low,bit_kind_size-1)

        if (j==k) then
          nperm = nperm + popcnt(iand(det1(j,ispin),           &
              iand( shiftl(1_bit_kind,m)-1_bit_kind,            &
                    not(shiftl(1_bit_kind,n))+1_bit_kind)) )
        else
          nperm = nperm + popcnt(                                    &
               iand(det1(j,ispin),                                   &
                    iand(not(0_bit_kind),                            &
                         (not(shiftl(1_bit_kind,n)) + 1_bit_kind) ))) &
               + popcnt(iand(det1(k,ispin),                          &
                             (shiftl(1_bit_kind,m) - 1_bit_kind ) ))

          do i=j+1,k-1
            nperm = nperm + popcnt(det1(i,ispin))
          end do

        endif

      case (2)

        do l=1,2
          high = max(exc(l,1,ispin), exc(l,2,ispin))-1
          low  = min(exc(l,1,ispin), exc(l,2,ispin))

          ASSERT (low > 0)
          ASSERT (high > 0)

          k = shiftr(high,bit_kind_shift)+1
          j = shiftr(low,bit_kind_shift)+1
          m = iand(high,bit_kind_size-1)
          n = iand(low,bit_kind_size-1)

          if (j==k) then
            nperm = nperm + popcnt(iand(det1(j,ispin),           &
                iand( shiftl(1_bit_kind,m)-1_bit_kind,            &
                    not(shiftl(1_bit_kind,n))+1_bit_kind)) )
          else
            nperm = nperm + popcnt(                                    &
                 iand(det1(j,ispin),                                   &
                      iand(not(0_bit_kind),                            &
                           (not(shiftl(1_bit_kind,n)) + 1_bit_kind) ))) &
                 + popcnt(iand(det1(k,ispin),                          &
                               (shiftl(1_bit_kind,m) - 1_bit_kind ) ))

            do i=j+1,k-1
              nperm = nperm + popcnt(det1(i,ispin))
            end do

          endif

        enddo

        a = min(exc(1,1,ispin), exc(1,2,ispin))
        b = max(exc(1,1,ispin), exc(1,2,ispin))
        c = min(exc(2,1,ispin), exc(2,2,ispin))
        d = max(exc(2,1,ispin), exc(2,2,ispin))
        if ((a<c) .and. (c<b) .and. (b<d)) then
          nperm = nperm + 1
        endif
        exit
    end select

  enddo
  phase = phase_dble(iand(nperm,1))

end

subroutine get_phasemask_bit(det1, pm, Nint)
  use bitmasks
  implicit none
  integer, intent(in) :: Nint
  integer(bit_kind), intent(in) :: det1(Nint,2)
  integer(bit_kind), intent(out) :: pm(Nint,2)
  integer(bit_kind) :: tmp
 integer :: ispin, i
  do ispin=1,2
  tmp = 0_8
  do i=1,Nint
    pm(i,ispin) = xor(det1(i,ispin), shiftl(det1(i,ispin), 1))
    pm(i,ispin) = xor(pm(i,ispin), shiftl(pm(i,ispin), 2))
    pm(i,ispin) = xor(pm(i,ispin), shiftl(pm(i,ispin), 4))
    pm(i,ispin) = xor(pm(i,ispin), shiftl(pm(i,ispin), 8))
    pm(i,ispin) = xor(pm(i,ispin), shiftl(pm(i,ispin), 16))
    pm(i,ispin) = xor(pm(i,ispin), shiftl(pm(i,ispin), 32))
    pm(i,ispin) = xor(pm(i,ispin), tmp)
    if(iand(popcnt(det1(i,ispin)), 1) == 1) tmp = not(tmp)
  end do
  end do
end



subroutine get_single_excitation(det1,det2,exc,phase,Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns the excitation operator between two singly excited determinants and the phase.
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: det1(Nint,2)
  integer(bit_kind), intent(in)  :: det2(Nint,2)
  integer, intent(out)           :: exc(0:2,2,2)
  double precision, intent(out)  :: phase
  integer                        :: tz
  integer                        :: l, ispin, idx_hole, idx_particle, ishift
  integer                        :: nperm
  integer                        :: i,j,k,m,n
  integer                        :: high, low
  integer                        :: a,b,c,d
  integer(bit_kind)              :: hole, particle, tmp
  double precision, parameter    :: phase_dble(0:1) = (/ 1.d0, -1.d0 /)

  ASSERT (Nint > 0)
  nperm = 0
  exc(0,1,1) = 0
  exc(0,2,1) = 0
  exc(0,1,2) = 0
  exc(0,2,2) = 0
  do ispin = 1,2
    ishift = 1-bit_kind_size
    do l=1,Nint
      ishift = ishift + bit_kind_size
      if (det1(l,ispin) == det2(l,ispin)) then
        cycle
      endif
      tmp = xor( det1(l,ispin), det2(l,ispin) )
      particle = iand(tmp, det2(l,ispin))
      hole     = iand(tmp, det1(l,ispin))
      if (particle /= 0_bit_kind) then
        tz = trailz(particle)
        exc(0,2,ispin) = 1
        exc(1,2,ispin) = tz+ishift
      endif
      if (hole /= 0_bit_kind) then
        tz = trailz(hole)
        exc(0,1,ispin) = 1
        exc(1,1,ispin) = tz+ishift
      endif

      if ( iand(exc(0,1,ispin),exc(0,2,ispin)) /= 1) then  ! exc(0,1,ispin)/=1 and exc(0,2,ispin) /= 1
        cycle
      endif

      high = max(exc(1,1,ispin), exc(1,2,ispin))-1
      low  = min(exc(1,1,ispin), exc(1,2,ispin))

      ASSERT (low >= 0)
      ASSERT (high > 0)

      k = shiftr(high,bit_kind_shift)+1
      j = shiftr(low,bit_kind_shift)+1
      m = iand(high,bit_kind_size-1)
      n = iand(low,bit_kind_size-1)

      if (j==k) then
        nperm = nperm + popcnt(iand(det1(j,ispin),           &
            iand( shiftl(1_bit_kind,m)-1_bit_kind,            &
                  not(shiftl(1_bit_kind,n))+1_bit_kind)) )
      else
        nperm = nperm + popcnt(                                    &
             iand(det1(j,ispin),                                   &
                  iand(not(0_bit_kind),                            &
                       (not(shiftl(1_bit_kind,n)) + 1_bit_kind) ))) &
             + popcnt(iand(det1(k,ispin),                          &
                           (shiftl(1_bit_kind,m) - 1_bit_kind ) ))

        do i=j+1,k-1
          nperm = nperm + popcnt(det1(i,ispin))
        end do

      endif

      phase = phase_dble(iand(nperm,1))
      return

    enddo
  enddo

end

subroutine get_single_excitation_cfg(cfg1,cfg2,p,q,Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns the excitation operator between two singly excited configurations.
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: cfg1(Nint,2)
  integer(bit_kind), intent(in)  :: cfg2(Nint,2)
  integer, intent(out)           :: p, q
  integer                        :: tz
  integer                        :: l, ispin, idx_hole, idx_particle, ishift
  integer                        :: nperm
  integer                        :: i,j,k,m,n
  integer                        :: high, low
  integer                        :: a,b,c,d
  integer(bit_kind)              :: hole, particle, tmp
  integer                        :: exc(0:2,2,2)

  ASSERT (Nint > 0)
  nperm = 0
  p = 0
  q = 0
  exc(0,1,1) = 0
  exc(0,2,1) = 0
  exc(0,1,2) = 0
  exc(0,2,2) = 0
  do ispin = 1,2
    ishift = 1-bit_kind_size
    do l=1,Nint
      ishift = ishift + bit_kind_size
      if (cfg1(l,ispin) == cfg2(l,ispin)) then
        cycle
      endif
      tmp = xor( cfg1(l,ispin), cfg2(l,ispin) )
      particle = iand(tmp, cfg2(l,ispin))
      hole     = iand(tmp, cfg1(l,ispin))
      if (particle /= 0_bit_kind) then
        tz = trailz(particle)
        exc(0,2,ispin) = 1
        exc(1,2,ispin) = tz+ishift
        !print *,"part ",tz+ishift, " ispin=",ispin
      endif
      if (hole /= 0_bit_kind) then
        tz = trailz(hole)
        exc(0,1,ispin) = 1
        exc(1,1,ispin) = tz+ishift
        !print *,"hole ",tz+ishift, " ispin=",ispin
      endif

      if ( iand(exc(0,1,ispin),exc(0,2,ispin)) /= 1) then  ! exc(0,1,ispin)/=1 and exc(0,2,ispin) /= 1
        cycle
      endif

      high = max(exc(1,1,ispin), exc(1,2,ispin))-1
      low  = min(exc(1,1,ispin), exc(1,2,ispin))

      ASSERT (low >= 0)
      ASSERT (high > 0)

      k = shiftr(high,bit_kind_shift)+1
      j = shiftr(low,bit_kind_shift)+1
      m = iand(high,bit_kind_size-1)
      n = iand(low,bit_kind_size-1)

      if (j==k) then
        nperm = nperm + popcnt(iand(cfg1(j,ispin),           &
            iand( shiftl(1_bit_kind,m)-1_bit_kind,            &
                  not(shiftl(1_bit_kind,n))+1_bit_kind)) )
      else
        nperm = nperm + popcnt(                                    &
             iand(cfg1(j,ispin),                                   &
                  iand(not(0_bit_kind),                            &
                       (not(shiftl(1_bit_kind,n)) + 1_bit_kind) ))) &
             + popcnt(iand(cfg1(k,ispin),                          &
                           (shiftl(1_bit_kind,m) - 1_bit_kind ) ))

        do i=j+1,k-1
          nperm = nperm + popcnt(cfg1(i,ispin))
        end do

      endif

      ! Set p and q
      q = max(exc(1,1,1),exc(1,1,2))
      p = max(exc(1,2,1),exc(1,2,2))
      return

    enddo
  enddo
end

subroutine bitstring_to_list_ab( string, list, n_elements, Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Gives the indices(+1) of the bits set to 1 in the bit string
  ! For alpha/beta determinants.
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: string(Nint,2)
  integer, intent(out)           :: list(Nint*bit_kind_size,2)
  integer, intent(out)           :: n_elements(2)

  integer                        :: i, j, ishift
  integer(bit_kind)              :: l

  n_elements(1) = 0
  n_elements(2) = 0
  ishift = 1
  do i=1,Nint
    l = string(i,1)
    do while (l /= 0_bit_kind)
      j = trailz(l)
      n_elements(1) = n_elements(1)+1
      l = ibclr(l,j)
      list(n_elements(1),1) = ishift+j
    enddo
    l = string(i,2)
    do while (l /= 0_bit_kind)
      j = trailz(l)
      n_elements(2) = n_elements(2)+1
      l = ibclr(l,j)
      list(n_elements(2),2) = ishift+j
    enddo
    ishift = ishift + bit_kind_size
  enddo

end


subroutine i_H_j_s2(key_i,key_j,Nint,hij,s2)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns $\langle i|H|j \rangle$ and $\langle i|S^2|j \rangle$
  ! where $i$ and $j$ are determinants.
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(out)  :: hij, s2

  integer                        :: exc(0:2,2,2)
  integer                        :: degree
  double precision               :: get_two_e_integral
  integer                        :: m,n,p,q
  integer                        :: i,j,k
  integer                        :: occ(Nint*bit_kind_size,2)
  double precision               :: diag_H_mat_elem, phase
  integer                        :: n_occ_ab(2)
  PROVIDE mo_two_e_integrals_in_map mo_integrals_map big_array_exchange_integrals

  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  ASSERT (sum(popcnt(key_i(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(key_i(:,2))) == elec_beta_num)
  ASSERT (sum(popcnt(key_j(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(key_j(:,2))) == elec_beta_num)

  hij = 0.d0
  s2 = 0d0
  !DIR$ FORCEINLINE
  call get_excitation_degree(key_i,key_j,degree,Nint)
  integer :: spin
  select case (degree)
    case (2)
      call get_double_excitation(key_i,key_j,exc,phase,Nint)
      ! Single alpha, single beta
      if (exc(0,1,1) == 1) then
        if ( (exc(1,1,1) == exc(1,2,2)).and.(exc(1,1,2) == exc(1,2,1)) ) then
          s2 =  -phase
        endif
        if(exc(1,1,1) == exc(1,2,2) )then
          hij = phase * big_array_exchange_integrals(exc(1,1,1),exc(1,1,2),exc(1,2,1))
        else if (exc(1,2,1) ==exc(1,1,2))then
          hij = phase * big_array_exchange_integrals(exc(1,2,1),exc(1,1,1),exc(1,2,2))
        else
          hij = phase*get_two_e_integral(                          &
              exc(1,1,1),                                              &
              exc(1,1,2),                                              &
              exc(1,2,1),                                              &
              exc(1,2,2) ,mo_integrals_map)
        endif
      ! Double alpha
      else if (exc(0,1,1) == 2) then
        hij = phase*(get_two_e_integral(                         &
            exc(1,1,1),                                              &
            exc(2,1,1),                                              &
            exc(1,2,1),                                              &
            exc(2,2,1) ,mo_integrals_map) -                          &
            get_two_e_integral(                                  &
            exc(1,1,1),                                              &
            exc(2,1,1),                                              &
            exc(2,2,1),                                              &
            exc(1,2,1) ,mo_integrals_map) )
      ! Double beta
      else if (exc(0,1,2) == 2) then
        hij = phase*(get_two_e_integral(                         &
            exc(1,1,2),                                              &
            exc(2,1,2),                                              &
            exc(1,2,2),                                              &
            exc(2,2,2) ,mo_integrals_map) -                          &
            get_two_e_integral(                                  &
            exc(1,1,2),                                              &
            exc(2,1,2),                                              &
            exc(2,2,2),                                              &
            exc(1,2,2) ,mo_integrals_map) )
      endif
    case (1)
      call get_single_excitation(key_i,key_j,exc,phase,Nint)
      !DIR$ FORCEINLINE
      call bitstring_to_list_ab(key_i, occ, n_occ_ab, Nint)
      ! Single alpha
      if (exc(0,1,1) == 1) then
        m = exc(1,1,1)
        p = exc(1,2,1)
        spin = 1
      ! Single beta
      else
        m = exc(1,1,2)
        p = exc(1,2,2)
        spin = 2
      endif
      call get_single_excitation_from_fock(key_i,key_j,p,m,spin,phase,hij)

    case (0)
      double precision, external :: diag_S_mat_elem
      s2 = diag_S_mat_elem(key_i,Nint)
      hij = diag_H_mat_elem(key_i,Nint)
  end select
end



subroutine i_H_j(key_i,key_j,Nint,hij)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns $\langle i|H|j \rangle$ where $i$ and $j$ are determinants.
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(out)  :: hij

  integer                        :: exc(0:2,2,2)
  integer                        :: degree
  double precision               :: get_two_e_integral
  integer                        :: m,n,p,q
  integer                        :: i,j,k
  integer                        :: occ(Nint*bit_kind_size,2)
  double precision               :: diag_H_mat_elem, phase
  integer                        :: n_occ_ab(2)
  PROVIDE mo_two_e_integrals_in_map mo_integrals_map big_array_exchange_integrals

  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  ASSERT (sum(popcnt(key_i(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(key_i(:,2))) == elec_beta_num)
  ASSERT (sum(popcnt(key_j(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(key_j(:,2))) == elec_beta_num)


  hij = 0.d0
  !DIR$ FORCEINLINE
  call get_excitation_degree(key_i,key_j,degree,Nint)
  integer :: spin
  select case (degree)
    case (2)
      call get_double_excitation(key_i,key_j,exc,phase,Nint)
      if (exc(0,1,1) == 1) then
        ! Single alpha, single beta
        if(exc(1,1,1) == exc(1,2,2) )then
          hij = phase * big_array_exchange_integrals(exc(1,1,1),exc(1,1,2),exc(1,2,1))
        else if (exc(1,2,1) ==exc(1,1,2))then
          hij = phase * big_array_exchange_integrals(exc(1,2,1),exc(1,1,1),exc(1,2,2))
        else
          hij = phase*get_two_e_integral(                          &
              exc(1,1,1),                                              &
              exc(1,1,2),                                              &
              exc(1,2,1),                                              &
              exc(1,2,2) ,mo_integrals_map)
        endif
      else if (exc(0,1,1) == 2) then
        ! Double alpha
        hij = phase*(get_two_e_integral(                         &
            exc(1,1,1),                                              &
            exc(2,1,1),                                              &
            exc(1,2,1),                                              &
            exc(2,2,1) ,mo_integrals_map) -                          &
            get_two_e_integral(                                  &
            exc(1,1,1),                                              &
            exc(2,1,1),                                              &
            exc(2,2,1),                                              &
            exc(1,2,1) ,mo_integrals_map) )
      else if (exc(0,1,2) == 2) then
        ! Double beta
        hij = phase*(get_two_e_integral(                         &
            exc(1,1,2),                                              &
            exc(2,1,2),                                              &
            exc(1,2,2),                                              &
            exc(2,2,2) ,mo_integrals_map) -                          &
            get_two_e_integral(                                  &
            exc(1,1,2),                                              &
            exc(2,1,2),                                              &
            exc(2,2,2),                                              &
            exc(1,2,2) ,mo_integrals_map) )
      endif
    case (1)
      call get_single_excitation(key_i,key_j,exc,phase,Nint)
      !DIR$ FORCEINLINE
      if (exc(0,1,1) == 1) then
        ! Single alpha
        m = exc(1,1,1)
        p = exc(1,2,1)
        spin = 1
      else
        ! Single beta
        m = exc(1,1,2)
        p = exc(1,2,2)
        spin = 2
      endif
      call get_single_excitation_from_fock(key_i,key_j,p,m,spin,phase,hij)

    case (0)
      hij = diag_H_mat_elem(key_i,Nint)
  end select
end

subroutine i_H_j_verbose(key_i,key_j,Nint,hij,hmono,hdouble,phase)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns $\langle i|H|j \rangle$ where $i$ and $j$ are determinants.
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(out)  :: hij,hmono,hdouble,phase

  integer                        :: exc(0:2,2,2)
  integer                        :: degree
  double precision               :: get_two_e_integral
  integer                        :: m,n,p,q
  integer                        :: i,j,k
  integer                        :: occ(Nint*bit_kind_size,2)
  double precision               :: diag_H_mat_elem
  integer                        :: n_occ_ab(2)
  logical                        :: has_mipi(Nint*bit_kind_size)
  double precision               :: mipi(Nint*bit_kind_size), miip(Nint*bit_kind_size)
  PROVIDE mo_two_e_integrals_in_map mo_integrals_map

  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  ASSERT (sum(popcnt(key_i(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(key_i(:,2))) == elec_beta_num)
  ASSERT (sum(popcnt(key_j(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(key_j(:,2))) == elec_beta_num)

  hij = 0.d0
  hmono = 0.d0
  hdouble = 0.d0
  !DIR$ FORCEINLINE
  call get_excitation_degree(key_i,key_j,degree,Nint)
  select case (degree)
    case (2)
      call get_double_excitation(key_i,key_j,exc,phase,Nint)
      if (exc(0,1,1) == 1) then
        ! Single alpha, single beta
        hij = phase*get_two_e_integral(                          &
            exc(1,1,1),                                              &
            exc(1,1,2),                                              &
            exc(1,2,1),                                              &
            exc(1,2,2) ,mo_integrals_map)
      else if (exc(0,1,1) == 2) then
        ! Double alpha
        hij = phase*(get_two_e_integral(                         &
            exc(1,1,1),                                              &
            exc(2,1,1),                                              &
            exc(1,2,1),                                              &
            exc(2,2,1) ,mo_integrals_map) -                          &
            get_two_e_integral(                                  &
            exc(1,1,1),                                              &
            exc(2,1,1),                                              &
            exc(2,2,1),                                              &
            exc(1,2,1) ,mo_integrals_map) )

      else if (exc(0,1,2) == 2) then
        ! Double beta
        hij = phase*(get_two_e_integral(                         &
            exc(1,1,2),                                              &
            exc(2,1,2),                                              &
            exc(1,2,2),                                              &
            exc(2,2,2) ,mo_integrals_map) -                          &
            get_two_e_integral(                                  &
            exc(1,1,2),                                              &
            exc(2,1,2),                                              &
            exc(2,2,2),                                              &
            exc(1,2,2) ,mo_integrals_map) )
      endif
    case (1)
      call get_single_excitation(key_i,key_j,exc,phase,Nint)
      !DIR$ FORCEINLINE
      call bitstring_to_list_ab(key_i, occ, n_occ_ab, Nint)
      has_mipi = .False.
      if (exc(0,1,1) == 1) then
        ! Single alpha
        m = exc(1,1,1)
        p = exc(1,2,1)
        do k = 1, elec_alpha_num
          i = occ(k,1)
          if (.not.has_mipi(i)) then
            mipi(i) = get_two_e_integral(m,i,p,i,mo_integrals_map)
            miip(i) = get_two_e_integral(m,i,i,p,mo_integrals_map)
            has_mipi(i) = .True.
          endif
        enddo
        do k = 1, elec_beta_num
          i = occ(k,2)
          if (.not.has_mipi(i)) then
            mipi(i) = get_two_e_integral(m,i,p,i,mo_integrals_map)
            has_mipi(i) = .True.
          endif
        enddo

        do k = 1, elec_alpha_num
          hdouble = hdouble + mipi(occ(k,1)) - miip(occ(k,1))
        enddo
        do k = 1, elec_beta_num
          hdouble = hdouble + mipi(occ(k,2))
        enddo

      else
        ! Single beta
        m = exc(1,1,2)
        p = exc(1,2,2)
        do k = 1, elec_beta_num
          i = occ(k,2)
          if (.not.has_mipi(i)) then
            mipi(i) = get_two_e_integral(m,i,p,i,mo_integrals_map)
            miip(i) = get_two_e_integral(m,i,i,p,mo_integrals_map)
            has_mipi(i) = .True.
          endif
        enddo
        do k = 1, elec_alpha_num
          i = occ(k,1)
          if (.not.has_mipi(i)) then
            mipi(i) = get_two_e_integral(m,i,p,i,mo_integrals_map)
            has_mipi(i) = .True.
          endif
        enddo

        do k = 1, elec_alpha_num
          hdouble = hdouble + mipi(occ(k,1))
        enddo
        do k = 1, elec_beta_num
          hdouble = hdouble + mipi(occ(k,2)) - miip(occ(k,2))
        enddo

      endif
      hmono = mo_one_e_integrals(m,p)
      hij = phase*(hdouble + hmono)

    case (0)
      phase = 1.d0
      hij = diag_H_mat_elem(key_i,Nint)
  end select
end

subroutine create_minilist(key_mask, fullList, miniList, idx_miniList, N_fullList, N_miniList, Nint)
  use bitmasks
  implicit none

  integer, intent(in)                      :: N_fullList
  integer, intent(in)                      :: Nint
  integer(bit_kind), intent(in)            :: fullList(Nint, 2, N_fullList)
  integer(bit_kind),intent(out)            :: miniList(Nint, 2, N_fullList)
  integer,intent(out)                      :: idx_miniList(N_fullList), N_miniList
  integer(bit_kind)                        :: key_mask(Nint, 2)
  integer                                  :: ni, k, i, n_a, n_b, e_a, e_b


  n_a = popcnt(key_mask(1,1))
  n_b = popcnt(key_mask(1,2))
  do ni=2,nint
    n_a = n_a + popcnt(key_mask(ni,1))
    n_b = n_b + popcnt(key_mask(ni,2))
  end do

  if(n_a == 0) then
    N_miniList = N_fullList
    do k=1,N_fullList
      do ni=1,Nint
        miniList(ni,1,k) = fullList(ni,1,k)
        miniList(ni,2,k) = fullList(ni,2,k)
      enddo
    enddo
    do i=1,N_fullList
      idx_miniList(i) = i
    end do
    return
  end if

  N_miniList = 0

  integer :: e_ab
  e_ab = n_a+n_b
  do i=1,N_fullList
    e_a = e_ab - popcnt(iand(fullList(1, 1, i), key_mask(1, 1))) &
               - popcnt(iand(fullList(1, 2, i), key_mask(1, 2)))
    do ni=2,nint
      e_a = e_a - popcnt(iand(fullList(ni, 1, i), key_mask(ni, 1))) &
                - popcnt(iand(fullList(ni, 2, i), key_mask(ni, 2)))
    end do

    if(e_a > 2) then
      cycle
    endif

    N_miniList = N_miniList + 1
    miniList(1,1,N_miniList) = fullList(1,1,i)
    miniList(1,2,N_miniList) = fullList(1,2,i)
    do ni=2,Nint
      miniList(ni,1,N_miniList) = fullList(ni,1,i)
      miniList(ni,2,N_miniList) = fullList(ni,2,i)
    enddo
    idx_miniList(N_miniList) = i

  end do
end subroutine

subroutine create_minilist_find_previous(key_mask, fullList, miniList, N_fullList, N_miniList, fullMatch, Nint)
  use bitmasks
  implicit none

  integer, intent(in)                      :: N_fullList
  integer, intent(in)                      :: Nint
  integer(bit_kind), intent(in)            :: fullList(Nint, 2, N_fullList)
  integer(bit_kind),intent(out)            :: miniList(Nint, 2, N_fullList)
  integer(bit_kind), allocatable           :: subList(:,:,:)
  logical,intent(out)                      :: fullMatch
  integer,intent(out)                      :: N_miniList
  integer(bit_kind)                        :: key_mask(Nint, 2)
  integer                                  :: ni, i, k, l, N_subList

  allocate (subList(Nint, 2, N_fullList))

  fullMatch = .false.
  N_miniList = 0
  N_subList = 0
  l = popcnt(key_mask(1,1)) + popcnt(key_mask(1,2))
  do ni = 2,Nint
    l = l + popcnt(key_mask(ni,1)) + popcnt(key_mask(ni,2))
  end do

  if(l == 0) then
    N_miniList = N_fullList
    do k=1,N_fullList
      do ni=1,Nint
        miniList(ni,1,k) = fullList(ni,1,k)
        miniList(ni,2,k) = fullList(ni,2,k)
      enddo
    enddo
  else
    do i=N_fullList,1,-1
      k = l
      do ni=1,nint
        k -= popcnt(iand(key_mask(ni,1), fullList(ni,1,i))) + popcnt(iand(key_mask(ni,2), fullList(ni,2,i)))
      end do
      if(k == 2) then
        N_subList += 1
        do ni=1,Nint
          subList(ni,1,N_subList) = fullList(ni,1,i)
          subList(ni,2,N_subList) = fullList(ni,2,i)
        enddo
      else if(k == 1) then
        N_minilist += 1
        do ni=1,Nint
          miniList(ni,1,N_minilist) = fullList(ni,1,i)
          miniList(ni,2,N_minilist) = fullList(ni,2,i)
        enddo
      else if(k == 0) then
        N_minilist += 1
        do ni=1,Nint
          miniList(ni,1,N_minilist) = fullList(ni,1,i)
          miniList(ni,2,N_minilist) = fullList(ni,2,i)
        enddo
!         fullMatch = .true.
!         return
      end if
    end do
  end if

  if(N_subList > 0) then
    do k=1,N_subList
      do ni=1,Nint
        miniList(ni,1,N_minilist+k) = sublist(ni,1,k)
        miniList(ni,2,N_minilist+k) = sublist(ni,2,k)
      enddo
    enddo
    N_minilist = N_minilist + N_subList
  end if

  deallocate(sublist)
end subroutine


subroutine i_H_psi(key,keys,coef,Nint,Ndet,Ndet_max,Nstate,i_H_psi_array)
  use bitmasks
  implicit none
  BEGIN_DOC
! Computes $\langle i|H|Psi \rangle  = \sum_J c_J \langle i | H | J \rangle$.
!
! Uses filter_connected_i_H_psi0 to get all the $|J \rangle$ to which $|i \rangle$
! is connected.
! The i_H_psi_minilist is much faster but requires to build the
! minilists.
  END_DOC
  integer, intent(in)            :: Nint, Ndet,Ndet_max,Nstate
  integer(bit_kind), intent(in)  :: keys(Nint,2,Ndet)
  integer(bit_kind), intent(in)  :: key(Nint,2)
  double precision, intent(in)   :: coef(Ndet_max,Nstate)
  double precision, intent(out)  :: i_H_psi_array(Nstate)

  integer                        :: i, ii,j
  double precision               :: phase
  integer                        :: exc(0:2,2,2)
  double precision               :: hij
  integer, allocatable           :: idx(:)

  ASSERT (Nint > 0)
  ASSERT (N_int == Nint)
  ASSERT (Nstate > 0)
  ASSERT (Ndet > 0)
  ASSERT (Ndet_max >= Ndet)
  allocate(idx(0:Ndet))

  i_H_psi_array = 0.d0

  call filter_connected_i_H_psi0(keys,key,Nint,Ndet,idx)
  if (Nstate == 1) then

    do ii=1,idx(0)
      i = idx(ii)
      !DIR$ FORCEINLINE
      call i_H_j(keys(1,1,i),key,Nint,hij)
      i_H_psi_array(1) = i_H_psi_array(1) + coef(i,1)*hij
    enddo

  else

    do ii=1,idx(0)
      i = idx(ii)
      !DIR$ FORCEINLINE
      call i_H_j(keys(1,1,i),key,Nint,hij)
      do j = 1, Nstate
        i_H_psi_array(j) = i_H_psi_array(j) + coef(i,j)*hij
      enddo
    enddo

  endif

end

subroutine i_H_psi_minilist(key,keys,idx_key,N_minilist,coef,Nint,Ndet,Ndet_max,Nstate,i_H_psi_array)
  use bitmasks
  implicit none
  integer, intent(in)            :: Nint, Ndet,Ndet_max,Nstate,idx_key(Ndet), N_minilist
  integer(bit_kind), intent(in)  :: keys(Nint,2,Ndet)
  integer(bit_kind), intent(in)  :: key(Nint,2)
  double precision, intent(in)   :: coef(Ndet_max,Nstate)
  double precision, intent(out)  :: i_H_psi_array(Nstate)

  integer                        :: i, ii,j, i_in_key, i_in_coef
  double precision               :: phase
  integer                        :: exc(0:2,2,2)
  double precision               :: hij
  integer, allocatable           :: idx(:)
  BEGIN_DOC
! Computes $\langle i|H|\Psi \rangle = \sum_J c_J \langle i|H|J\rangle$.
!
! Uses filter_connected_i_H_psi0 to get all the $|J \rangle$ to which $|i \rangle$
! is connected. The $|J\rangle$ are searched in short pre-computed lists.
  END_DOC

  ASSERT (Nint > 0)
  ASSERT (N_int == Nint)
  ASSERT (Nstate > 0)
  ASSERT (Ndet > 0)
  ASSERT (Ndet_max >= Ndet)
  allocate(idx(0:Ndet))
  i_H_psi_array = 0.d0

  call filter_connected_i_H_psi0(keys,key,Nint,N_minilist,idx)
  if (Nstate == 1) then

    do ii=1,idx(0)
      i_in_key = idx(ii)
      i_in_coef = idx_key(idx(ii))
      !DIR$ FORCEINLINE
      call i_H_j(keys(1,1,i_in_key),key,Nint,hij)
      ! TODO : Cache misses
      i_H_psi_array(1) = i_H_psi_array(1) + coef(i_in_coef,1)*hij
    enddo

  else

    do ii=1,idx(0)
      i_in_key = idx(ii)
      i_in_coef = idx_key(idx(ii))
      !DIR$ FORCEINLINE
      call i_H_j(keys(1,1,i_in_key),key,Nint,hij)
      do j = 1, Nstate
        i_H_psi_array(j) = i_H_psi_array(j) + coef(i_in_coef,j)*hij
      enddo
    enddo

  endif

end




subroutine get_excitation_degree_vector_single(key1,key2,degree,Nint,sze,idx)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Applies get_excitation_degree to an array of determinants and returns only
  ! the single excitations.
  END_DOC
  integer, intent(in)            :: Nint, sze
  integer(bit_kind), intent(in)  :: key1(Nint,2,sze)
  integer(bit_kind), intent(in)  :: key2(Nint,2)
  integer, intent(out)           :: degree(sze)
  integer, intent(out)           :: idx(0:sze)

  integer                        :: i,l,d,m

  ASSERT (Nint > 0)
  ASSERT (sze > 0)

  l=1
  if (Nint==1) then

    !DIR$ LOOP COUNT (1000)
    do i=1,sze
      d = popcnt(xor( key1(1,1,i), key2(1,1))) +       &
          popcnt(xor( key1(1,2,i), key2(1,2)))
      if (d > 2) then
        cycle
      else
        degree(l) = shiftr(d,1)
        idx(l) = i
        l = l+1
      endif
    enddo

  else if (Nint==2) then

    !DIR$ LOOP COUNT (1000)
    do i=1,sze
      d = popcnt(xor( key1(1,1,i), key2(1,1))) +                     &
          popcnt(xor( key1(1,2,i), key2(1,2))) +                     &
          popcnt(xor( key1(2,1,i), key2(2,1))) +                     &
          popcnt(xor( key1(2,2,i), key2(2,2)))
      if (d > 2) then
        cycle
      else
        degree(l) = shiftr(d,1)
        idx(l)    = i
        l         = l+1
      endif
    enddo

  else if (Nint==3) then

    !DIR$ LOOP COUNT (1000)
    do i=1,sze
      d = popcnt(xor( key1(1,1,i), key2(1,1))) +                     &
          popcnt(xor( key1(1,2,i), key2(1,2))) +                     &
          popcnt(xor( key1(2,1,i), key2(2,1))) +                     &
          popcnt(xor( key1(2,2,i), key2(2,2))) +                     &
          popcnt(xor( key1(3,1,i), key2(3,1))) +                     &
          popcnt(xor( key1(3,2,i), key2(3,2)))
      if (d > 2) then
        cycle
      else
        degree(l) = shiftr(d,1)
        idx(l)    = i
        l         = l+1
      endif
    enddo

  else

    !DIR$ LOOP COUNT (1000)
    do i=1,sze
      d = 0
      !DIR$ LOOP COUNT MIN(4)
      do m=1,Nint
        d = d + popcnt(xor( key1(m,1,i), key2(m,1)))                 &
              + popcnt(xor( key1(m,2,i), key2(m,2)))
      enddo
      if (d > 2) then
        cycle
      else
        degree(l) = shiftr(d,1)
        idx(l)    = i
        l         = l+1
      endif
    enddo

  endif
  idx(0) = l-1
end


subroutine get_excitation_degree_vector_single_or_exchange(key1,key2,degree,Nint,sze,idx)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Applies get_excitation_degree to an array of determinants and return only the
  ! single excitations and the connections through exchange integrals.
  END_DOC
  integer, intent(in)            :: Nint, sze
  integer(bit_kind), intent(in)  :: key1(Nint,2,sze)
  integer(bit_kind), intent(in)  :: key2(Nint,2)
  integer, intent(out)           :: degree(sze)
  integer, intent(out)           :: idx(0:sze)
  integer(bit_kind)              :: key_tmp(Nint,2)

  integer                        :: i,l,d,m
  integer                        :: exchange_1,exchange_2

  ASSERT (Nint > 0)
  ASSERT (sze > 0)

  l=1
  if (Nint==1) then

    !DIR$ LOOP COUNT (1000)
    do i=1,sze
      d = popcnt(xor( key1(1,1,i), key2(1,1))) +       &
          popcnt(xor( key1(1,2,i), key2(1,2)))
      key_tmp(1,1) = xor(key1(1,1,i),key2(1,1))
      key_tmp(1,2) = xor(key1(1,2,i),key2(1,2))
      if(popcnt(key_tmp(1,1)) .ge.3 .or. popcnt(key_tmp(1,2)) .ge.3 )cycle !! no double excitations of same spin
      if (d > 4)cycle
      if (d ==4)then
       if(popcnt(xor(key_tmp(1,1),key_tmp(1,2))) == 0)then
        degree(l) = shiftr(d,1)
        idx(l) = i
        l = l+1
       else
        cycle
       endif
!      pause
      else
        degree(l) = shiftr(d,1)
        idx(l) = i
        l = l+1
      endif
    enddo
  else

    print*, 'get_excitation_degree_vector_single_or_exchange not yet implemented for N_int > 1 ...'
    stop

  endif
  idx(0) = l-1
end




subroutine get_excitation_degree_vector_double_alpha_beta(key1,key2,degree,Nint,sze,idx)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Applies get_excitation_degree to an array of determinants and return only the
  ! single excitations and the connections through exchange integrals.
  END_DOC
  integer, intent(in)            :: Nint, sze
  integer(bit_kind), intent(in)  :: key1(Nint,2,sze)
  integer(bit_kind), intent(in)  :: key2(Nint,2)
  integer, intent(out)           :: degree(sze)
  integer, intent(out)           :: idx(0:sze)
  integer(bit_kind)              :: key_tmp(Nint,2)

  integer                        :: i,l,d,m
  integer :: degree_alpha, degree_beta

  ASSERT (Nint > 0)
  ASSERT (sze > 0)

  l=1
  if (Nint==1) then

    !DIR$ LOOP COUNT (1000)
    do i=1,sze
      d = popcnt(xor( key1(1,1,i), key2(1,1))) +       &
          popcnt(xor( key1(1,2,i), key2(1,2)))
      if (d .ne.4)cycle
      key_tmp(1,1) = xor(key1(1,1,i),key2(1,1))
      key_tmp(1,2) = xor(key1(1,2,i),key2(1,2))
      degree_alpha = popcnt(key_tmp(1,1))
      degree_beta  = popcnt(key_tmp(1,2))
      if(degree_alpha .ge.3 .or. degree_beta .ge.3 )cycle !! no double excitations of same spin
      degree(l) = shiftr(d,1)
      idx(l) = i
      l = l+1
    enddo
  else if (Nint==2) then

    !DIR$ LOOP COUNT (1000)
    do i=1,sze
      d = popcnt(xor( key1(1,1,i), key2(1,1))) +                     &
          popcnt(xor( key1(1,2,i), key2(1,2))) +                     &
          popcnt(xor( key1(2,1,i), key2(2,1))) +                     &
          popcnt(xor( key1(2,2,i), key2(2,2)))
      if (d .ne.4)cycle
      key_tmp(1,1) = xor(key1(1,1,i),key2(1,1))
      key_tmp(1,2) = xor(key1(1,2,i),key2(1,2))
      key_tmp(2,1) = xor(key1(2,1,i),key2(2,1))
      key_tmp(2,2) = xor(key1(2,2,i),key2(2,2))
      degree_alpha = popcnt(key_tmp(1,1)) + popcnt(key_tmp(2,1))
      degree_beta  = popcnt(key_tmp(1,2)) + popcnt(key_tmp(2,2))
      if(degree_alpha .ge.3 .or. degree_beta .ge.3 )cycle !! no double excitations of same spin
      degree(l) = shiftr(d,1)
      idx(l) = i
        l = l+1
    enddo

  else if (Nint==3) then

    !DIR$ LOOP COUNT (1000)
    do i=1,sze
      d = popcnt(xor( key1(1,1,i), key2(1,1))) +                     &
          popcnt(xor( key1(1,2,i), key2(1,2))) +                     &
          popcnt(xor( key1(2,1,i), key2(2,1))) +                     &
          popcnt(xor( key1(2,2,i), key2(2,2))) +                     &
          popcnt(xor( key1(3,1,i), key2(3,1))) +                     &
          popcnt(xor( key1(3,2,i), key2(3,2)))
      if (d .ne.4)cycle
      key_tmp(1,1) = xor(key1(1,1,i),key2(1,1))
      key_tmp(1,2) = xor(key1(1,2,i),key2(1,2))
      key_tmp(2,1) = xor(key1(2,1,i),key2(2,1))
      key_tmp(2,2) = xor(key1(2,2,i),key2(2,2))
      key_tmp(3,1) = xor(key1(3,1,i),key2(3,1))
      key_tmp(3,2) = xor(key1(3,2,i),key2(3,2))
      degree_alpha = popcnt(key_tmp(1,1)) + popcnt(key_tmp(2,1)) + popcnt(key_tmp(3,1))
      degree_beta  = popcnt(key_tmp(1,2)) + popcnt(key_tmp(2,2)) + popcnt(key_tmp(3,2))
      if(degree_alpha .ge.3 .or. degree_beta .ge.3 )cycle !! no double excitations of same spin
      degree(l) = shiftr(d,1)
      idx(l) = i
      l = l+1
    enddo

  else

    !DIR$ LOOP COUNT (1000)
    do i=1,sze
      d = 0
      degree_alpha = 0
      degree_beta  = 0
      !DIR$ LOOP COUNT MIN(4)
      do m=1,Nint
        d = d + popcnt(xor( key1(m,1,i), key2(m,1)))                 &
              + popcnt(xor( key1(m,2,i), key2(m,2)))
        key_tmp(m,1) = xor(key1(m,1,i),key2(m,1))
        key_tmp(m,2) = xor(key1(m,2,i),key2(m,2))
        degree_alpha += popcnt(key_tmp(m,1))
        degree_beta  += popcnt(key_tmp(m,2))
      enddo
      if(degree_alpha .ge.3 .or. degree_beta .ge.3 )cycle !! no double excitations of same spin
      degree(l) = shiftr(d,1)
      idx(l) = i
      l = l+1
    enddo

  endif
  idx(0) = l-1
end


subroutine get_excitation_degree_vector_single_or_exchange_verbose(key1,key2,degree,Nint,sze,idx)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Applies get_excitation_degree to an array of determinants and return only the single
  ! excitations and the connections through exchange integrals.
  END_DOC
  integer, intent(in)            :: Nint, sze
  integer(bit_kind), intent(in)  :: key1(Nint,2,sze)
  integer(bit_kind), intent(in)  :: key2(Nint,2)
  integer, intent(out)           :: degree(sze)
  integer, intent(out)           :: idx(0:sze)

  integer                        :: i,l,d,m
  integer                        :: exchange_1,exchange_2

  ASSERT (Nint > 0)
  ASSERT (sze > 0)

  l=1
  if (Nint==1) then

    !DIR$ LOOP COUNT (1000)
    do i=1,sze
      d = popcnt(xor( key1(1,1,i), key2(1,1))) +       &
          popcnt(xor( key1(1,2,i), key2(1,2)))
      exchange_1 = popcnt(xor(ior(key1(1,1,i),key1(1,2,i)),ior(key2(1,1),key2(1,2))))
      exchange_2 = popcnt(ior(xor(key1(1,1,i),key2(1,1)),xor(key1(1,2,i),key2(1,2))))
      if(i==99)then
      integer(bit_kind) :: key_test(N_int,2)
      key_test(1,2) = 0_bit_kind
       call debug_det(key2,N_int)
      key_test(1,1) = ior(key2(1,1),key2(1,2))
       call debug_det(key_test,N_int)
      key_test(1,1) = ior(key1(1,1,i),key1(1,2,i))
       call debug_det(key1(1,1,i),N_int)
       call debug_det(key_test,N_int)
       key_test(1,1) = xor(ior(key1(1,1,i),key1(1,2,i)),ior(key2(1,1),key2(1,2)))
       call debug_det(key_test,N_int)
       print*, exchange_1 , exchange_2
       stop
      endif
      if (d > 4)cycle
      if (d ==4)then
       if(exchange_1 .eq. 0 ) then
        degree(l) = shiftr(d,1)
        idx(l) = i
        l = l+1
       else if (exchange_1 .eq. 2 .and. exchange_2.eq.2)then
        degree(l) = shiftr(d,1)
        idx(l) = i
        l = l+1
       else
        cycle
       endif
!      pause
      else
        degree(l) = shiftr(d,1)
        idx(l) = i
        l = l+1
      endif
    enddo
  else if (Nint==2) then

    !DIR$ LOOP COUNT (1000)
    do i=1,sze
      d = popcnt(xor( key1(1,1,i), key2(1,1))) +                     &
          popcnt(xor( key1(1,2,i), key2(1,2))) +                     &
          popcnt(xor( key1(2,1,i), key2(2,1))) +                     &
          popcnt(xor( key1(2,2,i), key2(2,2)))
      exchange_1 = popcnt(xor(iand(key1(1,1,i),key1(1,2,i)),iand(key2(1,2),key2(1,2))))   +  &
                   popcnt(xor(iand(key1(2,1,i),key1(2,2,i)),iand(key2(2,2),key2(2,2))))
      exchange_2 = popcnt(iand(xor(key1(1,1,i),key2(1,1)),xor(key1(1,2,i),key2(1,2))))    +  &
                   popcnt(iand(xor(key1(2,1,i),key2(2,1)),xor(key1(2,2,i),key2(2,2))))
      if (d > 4)cycle
      if (d ==4)then
       if(exchange_1 .eq. 0 ) then
        degree(l) = shiftr(d,1)
        idx(l) = i
        l = l+1
       else if (exchange_1 .eq. 2 .and. exchange_2.eq.2)then
        degree(l) = shiftr(d,1)
        idx(l) = i
        l = l+1
       else
        cycle
       endif
!      pause
      else
        degree(l) = shiftr(d,1)
        idx(l) = i
        l = l+1
      endif
    enddo

  else if (Nint==3) then

    !DIR$ LOOP COUNT (1000)
    do i=1,sze
      d = popcnt(xor( key1(1,1,i), key2(1,1))) +                     &
          popcnt(xor( key1(1,2,i), key2(1,2))) +                     &
          popcnt(xor( key1(2,1,i), key2(2,1))) +                     &
          popcnt(xor( key1(2,2,i), key2(2,2))) +                     &
          popcnt(xor( key1(3,1,i), key2(3,1))) +                     &
          popcnt(xor( key1(3,2,i), key2(3,2)))
      exchange_1 = popcnt(xor(iand(key1(1,1,i),key1(1,2,i)),iand(key2(1,1),key2(1,2))))   +  &
                   popcnt(xor(iand(key1(2,1,i),key1(2,2,i)),iand(key2(2,1),key2(2,2))))   +  &
                   popcnt(xor(iand(key1(3,1,i),key1(3,2,i)),iand(key2(3,1),key2(3,2))))
      exchange_2 = popcnt(iand(xor(key1(1,1,i),key2(1,1)),xor(key1(1,2,i),key2(1,2))))    +  &
                   popcnt(iand(xor(key1(2,1,i),key2(2,1)),xor(key1(2,2,i),key2(2,2))))    +  &
                   popcnt(iand(xor(key1(3,1,i),key2(3,1)),xor(key1(3,2,i),key2(3,2))))
      if (d > 4)cycle
      if (d ==4)then
       if(exchange_1 .eq. 0 ) then
        degree(l) = shiftr(d,1)
        idx(l) = i
        l = l+1
       else if (exchange_1 .eq. 2 .and. exchange_2.eq.2)then
        degree(l) = shiftr(d,1)
        idx(l) = i
        l = l+1
       else
        cycle
       endif
!      pause
      else
        degree(l) = shiftr(d,1)
        idx(l) = i
        l = l+1
      endif
    enddo

  else

    !DIR$ LOOP COUNT (1000)
    do i=1,sze
      d = 0
      exchange_1 = 0
      exchange_2 = 0
      !DIR$ LOOP COUNT MIN(4)
      do m=1,Nint
        d = d + popcnt(xor( key1(m,1,i), key2(m,1)))                 &
              + popcnt(xor( key1(m,2,i), key2(m,2)))
        exchange_1 += popcnt(xor(iand(key1(m,1,i),key1(m,2,i)),iand(key2(m,1),key2(m,2))))
        exchange_2 += popcnt(iand(xor(key1(m,1,i),key2(m,1)),xor(key1(m,2,i),key2(m,2))))
      enddo
      if (d > 4)cycle
      if (d ==4)then
       if(exchange_1 .eq. 0 ) then
        degree(l) = shiftr(d,1)
        idx(l) = i
        l = l+1
       else if (exchange_1 .eq. 2 .and. exchange_2.eq.2)then
        degree(l) = shiftr(d,1)
        idx(l) = i
        l = l+1
       else
        cycle
       endif
      else
        degree(l) = shiftr(d,1)
        idx(l) = i
        l = l+1
      endif
    enddo

  endif
  idx(0) = l-1
end


subroutine get_excitation_degree_vector(key1,key2,degree,Nint,sze,idx)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Applies get_excitation_degree to an array of determinants.
  END_DOC
  integer, intent(in)            :: Nint, sze
  integer(bit_kind), intent(in)  :: key1(Nint,2,sze)
  integer(bit_kind), intent(in)  :: key2(Nint,2)
  integer, intent(out)           :: degree(sze)
  integer, intent(out)           :: idx(0:sze)

  integer                        :: i,l,d,m

  ASSERT (Nint > 0)
  ASSERT (sze > 0)

  l=1
  if (Nint==1) then

    do i=1,sze
      d = popcnt(xor( key1(1,1,i), key2(1,1))) +       &
          popcnt(xor( key1(1,2,i), key2(1,2)))
      if (d > 4) then
        cycle
      else
        degree(l) = shiftr(d,1)
        idx(l) = i
        l = l+1
      endif
    enddo

  else if (Nint==2) then

    do i=1,sze
      d = popcnt(xor( key1(1,1,i), key2(1,1))) +                     &
          popcnt(xor( key1(1,2,i), key2(1,2))) +                     &
          popcnt(xor( key1(2,1,i), key2(2,1))) +                     &
          popcnt(xor( key1(2,2,i), key2(2,2)))
      if (d > 4) then
        cycle
      else
        degree(l) = shiftr(d,1)
        idx(l)    = i
        l         = l+1
      endif
    enddo

  else if (Nint==3) then

    do i=1,sze
      d = popcnt(xor( key1(1,1,i), key2(1,1))) +                     &
          popcnt(xor( key1(1,2,i), key2(1,2))) +                     &
          popcnt(xor( key1(2,1,i), key2(2,1))) +                     &
          popcnt(xor( key1(2,2,i), key2(2,2))) +                     &
          popcnt(xor( key1(3,1,i), key2(3,1))) +                     &
          popcnt(xor( key1(3,2,i), key2(3,2)))
      if (d > 4) then
        cycle
      else
        degree(l) = shiftr(d,1)
        idx(l)    = i
        l         = l+1
      endif
    enddo

  else

    do i=1,sze
      d = 0
      do m=1,Nint
        d = d + popcnt(xor( key1(m,1,i), key2(m,1)))                 &
              + popcnt(xor( key1(m,2,i), key2(m,2)))
      enddo
      if (d > 4) then
        cycle
      else
        degree(l) = shiftr(d,1)
        idx(l)    = i
        l         = l+1
      endif
    enddo

  endif
  idx(0) = l-1
end




double precision function diag_H_mat_elem_fock(det_ref,det_pert,fock_diag_tmp,Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes $\langle i|H|i \rangle$ when $i$ is at most a double excitation from
  ! a reference.
  END_DOC
  integer,intent(in)             :: Nint
  integer(bit_kind),intent(in)   :: det_ref(Nint,2), det_pert(Nint,2)
  double precision, intent(in)   :: fock_diag_tmp(2,mo_num+1)

  integer                        :: degree
  double precision               :: phase, E0
  integer                        :: exc(0:2,2,2)
  integer                        :: h1, p1, h2, p2, s1, s2

  call get_excitation_degree(det_ref,det_pert,degree,Nint)
  E0 = fock_diag_tmp(1,mo_num+1)
  if (degree == 2) then
    call get_double_excitation(det_ref,det_pert,exc,phase,Nint)
    call decode_exc(exc,2,h1,p1,h2,p2,s1,s2)

    if ( (s1 == 1).and.(s2 == 1) ) then      ! alpha/alpha
      diag_H_mat_elem_fock = E0 &
        - fock_diag_tmp(1,h1) &
        + ( fock_diag_tmp(1,p1) - mo_two_e_integrals_jj_anti(h1,p1) ) &
        - ( fock_diag_tmp(1,h2) - mo_two_e_integrals_jj_anti(h1,h2)   &
            + mo_two_e_integrals_jj_anti(p1,h2) )                     &
        + ( fock_diag_tmp(1,p2) - mo_two_e_integrals_jj_anti(h1,p2)   &
            + mo_two_e_integrals_jj_anti(p1,p2) - mo_two_e_integrals_jj_anti(h2,p2) )

    else if ( (s1 == 2).and.(s2 == 2) ) then ! beta/beta
      diag_H_mat_elem_fock = E0 &
        - fock_diag_tmp(2,h1) &
        + ( fock_diag_tmp(2,p1) - mo_two_e_integrals_jj_anti(h1,p1) ) &
        - ( fock_diag_tmp(2,h2) - mo_two_e_integrals_jj_anti(h1,h2)   &
            + mo_two_e_integrals_jj_anti(p1,h2) )                     &
        + ( fock_diag_tmp(2,p2) - mo_two_e_integrals_jj_anti(h1,p2)   &
            + mo_two_e_integrals_jj_anti(p1,p2) - mo_two_e_integrals_jj_anti(h2,p2) )

    else                                    ! alpha/beta
      diag_H_mat_elem_fock = E0 &
        - fock_diag_tmp(1,h1) &
        + ( fock_diag_tmp(1,p1) - mo_two_e_integrals_jj_anti(h1,p1) ) &
        - ( fock_diag_tmp(2,h2) - mo_two_e_integrals_jj(h1,h2)        &
            + mo_two_e_integrals_jj(p1,h2) )                          &
        + ( fock_diag_tmp(2,p2) - mo_two_e_integrals_jj(h1,p2)        &
            + mo_two_e_integrals_jj(p1,p2) - mo_two_e_integrals_jj_anti(h2,p2) )

    endif

  else if (degree == 1) then
    call get_single_excitation(det_ref,det_pert,exc,phase,Nint)
    call decode_exc(exc,1,h1,p1,h2,p2,s1,s2)
    if (s1 == 1) then
      diag_H_mat_elem_fock = E0 - fock_diag_tmp(1,h1) &
        + ( fock_diag_tmp(1,p1) - mo_two_e_integrals_jj_anti(h1,p1) )
    else
      diag_H_mat_elem_fock = E0 - fock_diag_tmp(2,h1) &
        + ( fock_diag_tmp(2,p1) - mo_two_e_integrals_jj_anti(h1,p1) )
    endif

  else if (degree == 0) then
    diag_H_mat_elem_fock = E0
  else
    STOP 'Bug in diag_H_mat_elem_fock'
  endif
end

double precision function diag_H_mat_elem(det_in,Nint)
  implicit none
  BEGIN_DOC
  ! Computes $\langle i|H|i \rangle$.
  END_DOC
  integer,intent(in)             :: Nint
  integer(bit_kind),intent(in)   :: det_in(Nint,2)

  integer(bit_kind)              :: hole(Nint,2)
  integer(bit_kind)              :: particle(Nint,2)
  integer                        :: i, nexc(2), ispin
  integer                        :: occ_particle(Nint*bit_kind_size,2)
  integer                        :: occ_hole(Nint*bit_kind_size,2)
  integer(bit_kind)              :: det_tmp(Nint,2)
  integer                        :: na, nb

  ASSERT (Nint > 0)
  ASSERT (sum(popcnt(det_in(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(det_in(:,2))) == elec_beta_num)

  nexc(1) = 0
  nexc(2) = 0
  do i=1,Nint
    hole(i,1)     = xor(det_in(i,1),ref_bitmask(i,1))
    hole(i,2)     = xor(det_in(i,2),ref_bitmask(i,2))
    particle(i,1) = iand(hole(i,1),det_in(i,1))
    particle(i,2) = iand(hole(i,2),det_in(i,2))
    hole(i,1)     = iand(hole(i,1),ref_bitmask(i,1))
    hole(i,2)     = iand(hole(i,2),ref_bitmask(i,2))
    nexc(1)       = nexc(1) + popcnt(hole(i,1))
    nexc(2)       = nexc(2) + popcnt(hole(i,2))
  enddo

  diag_H_mat_elem = ref_bitmask_energy
  if (nexc(1)+nexc(2) == 0) then
    return
  endif

  !call debug_det(det_in,Nint)
  integer                        :: tmp(2)
  !DIR$ FORCEINLINE
  call bitstring_to_list_ab(particle, occ_particle, tmp, Nint)
  ASSERT (tmp(1) == nexc(1)) ! Number of particles alpha
  ASSERT (tmp(2) == nexc(2)) ! Number of particle beta 
  !DIR$ FORCEINLINE
  call bitstring_to_list_ab(hole, occ_hole, tmp, Nint)
  ASSERT (tmp(1) == nexc(1)) ! Number of holes alpha
  ASSERT (tmp(2) == nexc(2)) ! Number of holes beta 

  det_tmp = ref_bitmask
  do ispin=1,2
    na = elec_num_tab(ispin)
    nb = elec_num_tab(iand(ispin,1)+1)
    do i=1,nexc(ispin)
      !DIR$ FORCEINLINE
      call ac_operator( occ_particle(i,ispin), ispin, det_tmp, diag_H_mat_elem, Nint,na,nb)
      !DIR$ FORCEINLINE
      call a_operator ( occ_hole    (i,ispin), ispin, det_tmp, diag_H_mat_elem, Nint,na,nb)
    enddo
  enddo
end

subroutine a_operator(iorb,ispin,key,hjj,Nint,na,nb)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Needed for :c:func:`diag_H_mat_elem`.
  END_DOC
  integer, intent(in)            :: iorb, ispin, Nint
  integer, intent(inout)         :: na, nb
  integer(bit_kind), intent(inout) :: key(Nint,2)
  double precision, intent(inout) :: hjj

  integer                        :: occ(Nint*bit_kind_size,2)
  integer                        :: other_spin
  integer                        :: k,l,i
  integer                        :: tmp(2)

  ASSERT (iorb > 0)
  ASSERT (ispin > 0)
  ASSERT (ispin < 3)
  ASSERT (Nint > 0)

  k = shiftr(iorb-1,bit_kind_shift)+1
  ASSERT (k>0)
  l = iorb - shiftl(k-1,bit_kind_shift)-1
  key(k,ispin) = ibclr(key(k,ispin),l)
  other_spin = iand(ispin,1)+1

  !DIR$ FORCEINLINE
  call bitstring_to_list_ab(key, occ, tmp, Nint)
  na = na-1

  hjj = hjj - mo_one_e_integrals(iorb,iorb)

  ! Same spin
  do i=1,na
    hjj = hjj - mo_two_e_integrals_jj_anti(occ(i,ispin),iorb)
  enddo

  ! Opposite spin
  do i=1,nb
    hjj = hjj - mo_two_e_integrals_jj(occ(i,other_spin),iorb)
  enddo

end


subroutine ac_operator(iorb,ispin,key,hjj,Nint,na,nb)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Needed for :c:func:`diag_H_mat_elem`.
  END_DOC
  integer, intent(in)            :: iorb, ispin, Nint
  integer, intent(inout)         :: na, nb
  integer(bit_kind), intent(inout) :: key(Nint,2)
  double precision, intent(inout) :: hjj

  integer                        :: occ(Nint*bit_kind_size,2)
  integer                        :: other_spin
  integer                        :: k,l,i

  if (iorb < 1) then
    print *,  irp_here, ': iorb < 1'
    print *,  iorb, mo_num
    stop -1
  endif
  if (iorb > mo_num) then
    print *,  irp_here, ': iorb > mo_num'
    print *,  iorb, mo_num
    stop -1
  endif

  ASSERT (ispin > 0)
  ASSERT (ispin < 3)
  ASSERT (Nint > 0)

  integer                        :: tmp(2)
  !DIR$ FORCEINLINE
  call bitstring_to_list_ab(key, occ, tmp, Nint)
  ASSERT (tmp(1) == elec_alpha_num)
  ASSERT (tmp(2) == elec_beta_num)

  k = shiftr(iorb-1,bit_kind_shift)+1
  ASSERT (k >0)
  l = iorb - shiftl(k-1,bit_kind_shift)-1
  ASSERT (l >= 0)
  key(k,ispin) = ibset(key(k,ispin),l)
  other_spin = iand(ispin,1)+1

  hjj = hjj + mo_one_e_integrals(iorb,iorb)

  ! Same spin
  do i=1,na
    hjj = hjj + mo_two_e_integrals_jj_anti(occ(i,ispin),iorb)
  enddo

  ! Opposite spin
  do i=1,nb
    hjj = hjj + mo_two_e_integrals_jj(occ(i,other_spin),iorb)
  enddo
  na = na+1
end


subroutine get_phase(key1,key2,phase,Nint)
  use bitmasks
  implicit none
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key1(Nint,2), key2(Nint,2)
  double precision, intent(out)  :: phase
  BEGIN_DOC
! Returns the phase between key1 and key2.
  END_DOC
  integer                        :: exc(0:2, 2, 2), degree

  !DIR$ FORCEINLINE
  call get_excitation(key1, key2, exc, degree, phase, Nint)
end



! Spin-determinant routines
! -------------------------

subroutine get_excitation_degree_spin(key1,key2,degree,Nint)
  use bitmasks
  include 'utils/constants.include.F'
  implicit none
  BEGIN_DOC
  ! Returns the excitation degree between two determinants.
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key1(Nint)
  integer(bit_kind), intent(in)  :: key2(Nint)
  integer, intent(out)           :: degree

  integer(bit_kind)              :: xorvec(N_int_max)
  integer                        :: l

  ASSERT (Nint > 0)

  select case (Nint)

    case (1)
      xorvec(1) = xor( key1(1), key2(1))
      degree = popcnt(xorvec(1))

    case (2)
      xorvec(1) = xor( key1(1), key2(1))
      xorvec(2) = xor( key1(2), key2(2))
      degree = popcnt(xorvec(1))+popcnt(xorvec(2))

    case (3)
      xorvec(1) = xor( key1(1), key2(1))
      xorvec(2) = xor( key1(2), key2(2))
      xorvec(3) = xor( key1(3), key2(3))
      degree = sum(popcnt(xorvec(1:3)))

    case (4)
      xorvec(1) = xor( key1(1), key2(1))
      xorvec(2) = xor( key1(2), key2(2))
      xorvec(3) = xor( key1(3), key2(3))
      xorvec(4) = xor( key1(4), key2(4))
      degree = sum(popcnt(xorvec(1:4)))

    case default
      do l=1,Nint
        xorvec(l) = xor( key1(l), key2(l))
      enddo
      degree = sum(popcnt(xorvec(1:Nint)))

  end select

  degree = shiftr(degree,1)

end


subroutine get_excitation_spin(det1,det2,exc,degree,phase,Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns the excitation operators between two determinants and the phase.
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: det1(Nint)
  integer(bit_kind), intent(in)  :: det2(Nint)
  integer, intent(out)           :: exc(0:2,2)
  integer, intent(out)           :: degree
  double precision, intent(out)  :: phase
  ! exc(number,hole/particle)
  ! ex :
  ! exc(0,1) = number of holes
  ! exc(0,2) = number of particles
  ! exc(1,2) = first particle
  ! exc(1,1) = first hole

  ASSERT (Nint > 0)

  !DIR$ FORCEINLINE
  call get_excitation_degree_spin(det1,det2,degree,Nint)
  select case (degree)

    case (3:)
      degree = -1
      return

    case (2)
      call get_double_excitation_spin(det1,det2,exc,phase,Nint)
      return

    case (1)
      call get_single_excitation_spin(det1,det2,exc,phase,Nint)
      return

    case(0)
      return

  end select
end

subroutine decode_exc_spin(exc,h1,p1,h2,p2)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Decodes the exc arrays returned by get_excitation.
  !
  ! h1,h2 : Holes
  !
  ! p1,p2 : Particles
  END_DOC
  integer, intent(in)            :: exc(0:2,2)
  integer, intent(out)           :: h1,h2,p1,p2

  select case (exc(0,1))
    case(2)
      h1 = exc(1,1)
      h2 = exc(2,1)
      p1 = exc(1,2)
      p2 = exc(2,2)
    case(1)
      h1 = exc(1,1)
      h2 = 0
      p1 = exc(1,2)
      p2 = 0
    case default
      h1 = 0
      p1 = 0
      h2 = 0
      p2 = 0
  end select
end


subroutine get_double_excitation_spin(det1,det2,exc,phase,Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns the two excitation operators between two doubly excited spin-determinants
  ! and the phase.
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: det1(Nint)
  integer(bit_kind), intent(in)  :: det2(Nint)
  integer, intent(out)           :: exc(0:2,2)
  double precision, intent(out)  :: phase
  integer                        :: tz
  integer                        :: l, idx_hole, idx_particle, ishift
  integer                        :: nperm
  integer                        :: i,j,k,m,n
  integer                        :: high, low
  integer                        :: a,b,c,d
  integer(bit_kind)              :: hole, particle, tmp
  double precision, parameter    :: phase_dble(0:1) = (/ 1.d0, -1.d0 /)

  ASSERT (Nint > 0)
  nperm = 0
  exc(0,1) = 0
  exc(0,2) = 0

  idx_particle = 0
  idx_hole = 0
  ishift = 1-bit_kind_size
  do l=1,Nint
    ishift = ishift + bit_kind_size
    if (det1(l) == det2(l)) then
      cycle
    endif
    tmp = xor( det1(l), det2(l) )
    particle = iand(tmp, det2(l))
    hole     = iand(tmp, det1(l))
    do while (particle /= 0_bit_kind)
      tz = trailz(particle)
      idx_particle = idx_particle + 1
      exc(0,2) = exc(0,2) + 1
      exc(idx_particle,2) = tz+ishift
      particle = iand(particle,particle-1_bit_kind)
    enddo
    if (iand(exc(0,1),exc(0,2))==2) then  ! exc(0,1)==2 or exc(0,2)==2
      exit
    endif
    do while (hole /= 0_bit_kind)
      tz = trailz(hole)
      idx_hole = idx_hole + 1
      exc(0,1) = exc(0,1) + 1
      exc(idx_hole,1) = tz+ishift
      hole = iand(hole,hole-1_bit_kind)
    enddo
    if (iand(exc(0,1),exc(0,2))==2) then ! exc(0,1)==2 or exc(0,2)==2
      exit
    endif
  enddo

  select case (exc(0,1))

    case(1)
      low  = min(exc(1,1), exc(1,2))
      high = max(exc(1,1), exc(1,2))

      ASSERT (low > 0)
      j = shiftr(low-1,bit_kind_shift)+1   ! Find integer in array(Nint)
      n = iand(low-1,bit_kind_size-1)+1        ! mod(low,bit_kind_size)
      ASSERT (high > 0)
      k = shiftr(high-1,bit_kind_shift)+1
      m = iand(high-1,bit_kind_size-1)+1

      if (j==k) then
        nperm = nperm + popcnt(iand(det1(j),                         &
            iand( ibset(0_bit_kind,m-1)-1_bit_kind,                  &
            ibclr(-1_bit_kind,n)+1_bit_kind ) ))
      else
        nperm = nperm + popcnt(iand(det1(k),                         &
            ibset(0_bit_kind,m-1)-1_bit_kind))
        if (n < bit_kind_size) then
          nperm = nperm + popcnt(iand(det1(j), ibclr(-1_bit_kind,n) +1_bit_kind))
        endif
        do i=j+1,k-1
          nperm = nperm + popcnt(det1(i))
        end do
      endif

    case (2)

      do i=1,2
        low  = min(exc(i,1), exc(i,2))
        high = max(exc(i,1), exc(i,2))

        ASSERT (low > 0)
        j = shiftr(low-1,bit_kind_shift)+1   ! Find integer in array(Nint)
        n = iand(low-1,bit_kind_size-1)+1        ! mod(low,bit_kind_size)
        ASSERT (high > 0)
        k = shiftr(high-1,bit_kind_shift)+1
        m = iand(high-1,bit_kind_size-1)+1

        if (j==k) then
          nperm = nperm + popcnt(iand(det1(j),                       &
              iand( ibset(0_bit_kind,m-1)-1_bit_kind,                &
              ibclr(-1_bit_kind,n)+1_bit_kind ) ))
        else
          nperm = nperm + popcnt(iand(det1(k),                       &
              ibset(0_bit_kind,m-1)-1_bit_kind))
          if (n < bit_kind_size) then
            nperm = nperm + popcnt(iand(det1(j), ibclr(-1_bit_kind,n) +1_bit_kind))
          endif
          do l=j+1,k-1
            nperm = nperm + popcnt(det1(l))
          end do
        endif

      enddo

      a = min(exc(1,1), exc(1,2))
      b = max(exc(1,1), exc(1,2))
      c = min(exc(2,1), exc(2,2))
      d = max(exc(2,1), exc(2,2))
      if (c>a .and. c<b .and. d>b) then
        nperm = nperm + 1
      endif
  end select

  phase = phase_dble(iand(nperm,1))

end

subroutine get_single_excitation_spin(det1,det2,exc,phase,Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns the excitation operator between two singly excited determinants and the phase.
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: det1(Nint)
  integer(bit_kind), intent(in)  :: det2(Nint)
  integer, intent(out)           :: exc(0:2,2)
  double precision, intent(out)  :: phase
  integer                        :: tz
  integer                        :: l, idx_hole, idx_particle, ishift
  integer                        :: nperm
  integer                        :: i,j,k,m,n
  integer                        :: high, low
  integer                        :: a,b,c,d
  integer(bit_kind)              :: hole, particle, tmp
  double precision, parameter    :: phase_dble(0:1) = (/ 1.d0, -1.d0 /)

  ASSERT (Nint > 0)
  nperm = 0
  exc(0,1) = 0
  exc(0,2) = 0

  ishift = 1-bit_kind_size
  do l=1,Nint
    ishift = ishift + bit_kind_size
    if (det1(l) == det2(l)) then
      cycle
    endif
    tmp = xor( det1(l), det2(l) )
    particle = iand(tmp, det2(l))
    hole     = iand(tmp, det1(l))
    if (particle /= 0_bit_kind) then
      tz = trailz(particle)
      exc(0,2) = 1
      exc(1,2) = tz+ishift
    endif
    if (hole /= 0_bit_kind) then
      tz = trailz(hole)
      exc(0,1) = 1
      exc(1,1) = tz+ishift
    endif

    if ( iand(exc(0,1),exc(0,2)) /= 1) then  ! exc(0,1)/=1 and exc(0,2) /= 1
      cycle
    endif

    low = min(exc(1,1),exc(1,2))
    high = max(exc(1,1),exc(1,2))

    ASSERT (low > 0)
    j = shiftr(low-1,bit_kind_shift)+1   ! Find integer in array(Nint)
    n = iand(low-1,bit_kind_size-1)+1      ! mod(low,bit_kind_size)
    ASSERT (high > 0)
    k = shiftr(high-1,bit_kind_shift)+1
    m = iand(high-1,bit_kind_size-1)+1
    if (j==k) then
      nperm = popcnt(iand(det1(j),                                   &
          iand(ibset(0_bit_kind,m-1)-1_bit_kind,ibclr(-1_bit_kind,n)+1_bit_kind)))
    else
      nperm = nperm + popcnt(iand(det1(k),ibset(0_bit_kind,m-1)-1_bit_kind))
      if (n < bit_kind_size) then
        nperm = nperm + popcnt(iand(det1(j),ibclr(-1_bit_kind,n)+1_bit_kind))
      endif
      do i=j+1,k-1
        nperm = nperm + popcnt(det1(i))
      end do
    endif
    phase = phase_dble(iand(nperm,1))
    return

  enddo
end

subroutine i_H_j_single_spin(key_i,key_j,Nint,spin,hij)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns $\langle i|H|j \rangle$ where $i$ and $j$ are determinants differing by
  ! a single excitation.
  END_DOC
  integer, intent(in)            :: Nint, spin
  integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(out)  :: hij

  integer                        :: exc(0:2,2)
  double precision               :: phase

  PROVIDE big_array_exchange_integrals mo_two_e_integrals_in_map

  call get_single_excitation_spin(key_i(1,spin),key_j(1,spin),exc,phase,Nint)
  call get_single_excitation_from_fock(key_i,key_j,exc(1,1),exc(1,2),spin,phase,hij)
end

subroutine i_H_j_double_spin(key_i,key_j,Nint,hij)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns $\langle i|H|j \rangle$ where $i$ and $j$ are determinants differing by
  ! a same-spin double excitation.
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key_i(Nint), key_j(Nint)
  double precision, intent(out)  :: hij

  integer                        :: exc(0:2,2)
  double precision               :: phase
  double precision, external     :: get_two_e_integral

  PROVIDE big_array_exchange_integrals mo_two_e_integrals_in_map
  call get_double_excitation_spin(key_i,key_j,exc,phase,Nint)
  hij = phase*(get_two_e_integral(                             &
      exc(1,1),                                                    &
      exc(2,1),                                                    &
      exc(1,2),                                                    &
      exc(2,2), mo_integrals_map) -                                &
      get_two_e_integral(                                      &
      exc(1,1),                                                    &
      exc(2,1),                                                    &
      exc(2,2),                                                    &
      exc(1,2), mo_integrals_map) )
end

subroutine i_H_j_double_alpha_beta(key_i,key_j,Nint,hij)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Returns $\langle i|H|j \rangle$ where $i$ and $j$ are determinants differing by
  ! an opposite-spin double excitation.
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(out)  :: hij

  integer                        :: exc(0:2,2,2)
  double precision               :: phase, phase2
  double precision, external     :: get_two_e_integral

  PROVIDE big_array_exchange_integrals mo_two_e_integrals_in_map

  call get_single_excitation_spin(key_i(1,1),key_j(1,1),exc(0,1,1),phase,Nint)
  call get_single_excitation_spin(key_i(1,2),key_j(1,2),exc(0,1,2),phase2,Nint)
  phase = phase*phase2
  if (exc(1,1,1) == exc(1,2,2)) then
    hij = phase * big_array_exchange_integrals(exc(1,1,1),exc(1,1,2),exc(1,2,1))
  else if (exc(1,2,1) == exc(1,1,2)) then
    hij = phase * big_array_exchange_integrals(exc(1,2,1),exc(1,1,1),exc(1,2,2))
  else
    hij = phase*get_two_e_integral(                              &
        exc(1,1,1),                                                  &
        exc(1,1,2),                                                  &
        exc(1,2,1),                                                  &
        exc(1,2,2) ,mo_integrals_map)
  endif
end


subroutine connected_to_hf(key_i,yes_no)
 implicit none 
 use bitmasks
 integer(bit_kind), intent(in)  :: key_i(N_int,2)
 logical , intent(out) :: yes_no
 double precision :: hij,phase
 integer          :: exc(0:2,2,2)
 integer          :: degree
 integer          :: m,p
 yes_no = .True.
 call get_excitation_degree(ref_bitmask,key_i,degree,N_int)
 if(degree == 2)then
  call i_H_j(ref_bitmask,key_i,N_int,hij)
  if(dabs(hij) .lt. thresh_sym)then
   yes_no = .False. 
  endif
 else if(degree == 1)then  
  call get_single_excitation(ref_bitmask,key_i,exc,phase,N_int)
  ! Single alpha
  if (exc(0,1,1) == 1) then
    m = exc(1,1,1)
    p = exc(1,2,1)
  ! Single beta
  else
    m = exc(1,1,2)
    p = exc(1,2,2)
  endif
  hij = mo_one_e_integrals(m,p)
  if(dabs(hij) .lt. thresh_sym)then
   yes_no = .False. 
  endif
 else if(degree == 0)then
  yes_no = .True.
 endif
end
