use bitmasks

BEGIN_PROVIDER [ integer, N_det_selectors]
  implicit none
  BEGIN_DOC
  ! For Single reference wave functions, the number of selectors is 1 : the
  ! Hartree-Fock determinant
  END_DOC
  N_det_selectors = N_det
END_PROVIDER

 BEGIN_PROVIDER [ integer(bit_kind), psi_selectors, (N_int,2,psi_selectors_size) ]
&BEGIN_PROVIDER [ double precision, psi_selectors_coef, (psi_selectors_size,N_states) ]
  implicit none
  BEGIN_DOC
  ! Determinants on which we apply <i|H|psi> for perturbation.
  ! The selectors are equivalent to Selectors_full, but in a different
  ! order. The Generators_CAS determinants appear first, then all the
  ! others.
  END_DOC
  integer                        :: i, k, l, m
  logical                        :: good
  integer, external :: number_of_holes,number_of_particles

  do i=1,N_det_generators
    do k=1,N_int
      psi_selectors(k,1,i) = psi_det_generators(k,1,i)
      psi_selectors(k,2,i) = psi_det_generators(k,2,i)
    enddo
  enddo
  do k=1,N_states
    do i=1,N_det_generators
      psi_selectors_coef(i,k) = psi_coef_generators(i,k)
    enddo
  enddo

  m=N_det_generators

  do i=1,N_det
    good = ( number_of_holes(psi_det_sorted(1,1,i)) ==0).and.(number_of_particles(psi_det_sorted(1,1,i))==0 )
    if (.not.good) then
      m = m+1
      do k=1,N_int
        psi_selectors(k,1,m) = psi_det_sorted(k,1,i)
        psi_selectors(k,2,m) = psi_det_sorted(k,2,i)
      enddo
      psi_selectors_coef(m,:) = psi_coef_sorted(i,:)
    endif
  enddo
  if (N_det /= m) then
    print *,  N_det, m
    stop 'Selectors_CASSD : N_det /= m'
  endif
END_PROVIDER


