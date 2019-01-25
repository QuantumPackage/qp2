BEGIN_TEMPLATE

subroutine pt2_epstein_nesbet ($arguments)
  use bitmasks
  implicit none
  $declarations

  BEGIN_DOC
  ! Compute the standard Epstein-Nesbet perturbative first order coefficient and
  ! second order energetic contribution for the various N_st states.
  !
  ! `c_pert(i)` = $\\frac{\langle i|H|\\alpha \\rangle}{ E_n - \\langle \\alpha|H|\\alpha \\rangle }$.
  !
  ! `e_2_pert(i)` = $\\frac{\\langle i|H|\\alpha \\rangle^2}{ E_n - \\langle \\alpha|H|\\alpha \\rangle }$.
  !
  END_DOC

  integer                        :: i,j
  double precision               :: diag_H_mat_elem_fock, h
  double precision               :: i_H_psi_array(N_st)
  PROVIDE  selection_criterion

  ASSERT (Nint == N_int)
  ASSERT (Nint > 0)
  !call i_H_psi(det_pert,psi_selectors,psi_selectors_coef,Nint,N_det_selectors,psi_selectors_size,N_st,i_H_psi_array)
  call i_H_psi_minilist(det_pert,minilist,idx_minilist,N_minilist,psi_selectors_coef,Nint,N_minilist,psi_selectors_size,N_st,i_H_psi_array)


  h = diag_H_mat_elem_fock(det_ref,det_pert,fock_diag_tmp,Nint)
  do i =1,N_st
    if(electronic_energy(i)>h.and.electronic_energy(i).ne.0.d0)then
      c_pert(i) = -1.d0
      e_2_pert(i) = selection_criterion*selection_criterion_factor*2.d0
    else if  (dabs(electronic_energy(i) - h) > 1.d-6) then
        c_pert(i) = i_H_psi_array(i) / (electronic_energy(i) - h)
        H_pert_diag(i) = h*c_pert(i)*c_pert(i)
        e_2_pert(i) = c_pert(i) * i_H_psi_array(i)
    else
      c_pert(i) = -1.d0
      e_2_pert(i) = -dabs(i_H_psi_array(i))
      H_pert_diag(i) = h
    endif
  enddo

end

subroutine pt2_qdpt ($arguments)
  use bitmasks
  implicit none
  $declarations

  BEGIN_DOC
  ! Computes the QDPT first order coefficient and second order energetic contribution
  ! for the various N_st states.
  !
  ! `c_pert(i)` = $\\frac{\\langle i|H|\\alpha \\rangle}{\\langle i|H|i \\rangle - \\langle \\alpha|H|\\alpha \\rangle}$.
  !
  END_DOC

  integer                        :: i,j
  double precision               :: diag_H_mat_elem_fock, h, E, diag_H_mat_elem, hij
  double precision               :: i_H_psi_array(N_st)
  integer :: degree
  double precision :: delta_E
  PROVIDE  selection_criterion

  ASSERT (Nint == N_int)
  ASSERT (Nint > 0)
  !call i_H_psi(det_pert,psi_selectors,psi_selectors_coef,Nint,N_det_selectors,psi_selectors_size,N_st,i_H_psi_array)
  call i_H_psi_minilist(det_pert,minilist,idx_minilist,N_minilist,psi_selectors_coef,Nint,N_minilist,psi_selectors_size,N_st,i_H_psi_array)


  h = diag_H_mat_elem_fock(det_ref,det_pert,fock_diag_tmp,Nint)
  c_pert = 0.d0
  do j=1,N_det_selectors
    call get_excitation_degree(det_ref, psi_selectors(1,1,j), degree, Nint)
    if (degree > 2) then
      E = diag_H_mat_elem(psi_selectors(1,1,j),Nint)
    else
      E = diag_H_mat_elem_fock(det_ref,det_ref,fock_diag_tmp,Nint)
    endif
    delta_E = E-h
!    delta_E = electronic_energy(1) - h
    call i_H_j(psi_selectors(1,1,j),det_pert,Nint,hij)
    if (dabs(delta_e) > 1.d-3) then
        do i =1,N_st
          c_pert(i) += psi_selectors_coef(j,i) * hij / delta_e
        enddo
    endif
  enddo
  do i =1,N_st
    e_2_pert(i) = c_pert(i)*i_H_psi_array(i)
    H_pert_diag(i) = h*c_pert(i)*c_pert(i)
  enddo

end


subroutine pt2_epstein_nesbet_2x2 ($arguments)
  use bitmasks
  implicit none
  $declarations

  BEGIN_DOC
  ! Computes the Epstein-Nesbet 2x2 diagonalization coefficient and energetic contribution
  ! for the various N_st states.
  !
  ! `e_2_pert(i)` = $\\frac{1}{2} ( \\langle \\alpha|H|\\alpha \\rangle -  E_n) - \\sqrt{ (\\langle \\alpha|H|\\alpha \\rangle -  E_n)^2 + 4 \\langle i|H|\\alpha \\rangle^2 }$.
  !
  ! `c_pert(i)` = `e_2_pert(i)` $\\times \\frac{1}{ \\langle i|H|\\alpha \\rangle}$.
  !
  END_DOC

  integer                        :: i,j
  double precision               :: diag_H_mat_elem_fock,delta_e, h
  double precision               :: i_H_psi_array(N_st)
  ASSERT (Nint == N_int)
  ASSERT (Nint > 0)

   call i_H_psi(det_pert,psi_selectors,psi_selectors_coef,Nint,N_det_selectors,psi_selectors_size,N_st,i_H_psi_array)
  !call i_H_psi_minilist(det_pert,minilist,idx_minilist,N_minilist,psi_selectors_coef,Nint,N_minilist,psi_selectors_size,N_st,i_H_psi_array)

  h = diag_H_mat_elem_fock(det_ref,det_pert,fock_diag_tmp,Nint)
  do i =1,N_st
    if (i_H_psi_array(i) /= 0.d0) then
      delta_e = h - electronic_energy(i)
      if (delta_e > 0.d0) then
        e_2_pert(i) = 0.5d0 * (delta_e - dsqrt(delta_e * delta_e + 4.d0 * i_H_psi_array(i) * i_H_psi_array(i)))
      else
        e_2_pert(i) = 0.5d0 * (delta_e + dsqrt(delta_e * delta_e + 4.d0 * i_H_psi_array(i) * i_H_psi_array(i)))
      endif
      if (dabs(i_H_psi_array(i)) > 1.d-6) then
        c_pert(i) = e_2_pert(i)/i_H_psi_array(i)
      else
        c_pert(i) = 0.d0
      endif
      H_pert_diag(i) = h*c_pert(i)*c_pert(i)
    else
      e_2_pert(i) = 0.d0
      c_pert(i) = 0.d0
      H_pert_diag(i) = 0.d0
    endif
  enddo

end



subroutine pt2_epstein_nesbet_2x2_no_ci_diag($arguments)
  use bitmasks
  implicit none
  $declarations

  BEGIN_DOC
  ! compute the Epstein-Nesbet 2x2 diagonalization coefficient and energetic contribution
  !
  ! for the various N_st states.
  !
  ! e_2_pert(i) = 0.5 * (( <det_pert|H|det_pert> -  E(i) )  - sqrt( ( <det_pert|H|det_pert> -  E(i)) ^2 + 4 <psi(i)|H|det_pert>^2  )
  !
  ! c_pert(i) = e_2_pert(i)/ <psi(i)|H|det_pert>
  !
  END_DOC

  integer                        :: i,j
  double precision               :: diag_H_mat_elem_fock,delta_e, h
  double precision               :: i_H_psi_array(N_st)
  ASSERT (Nint == N_int)
  ASSERT (Nint > 0)
  PROVIDE psi_energy

   call i_H_psi(det_pert,psi_selectors,psi_selectors_coef,Nint,N_det_selectors,psi_selectors_size,N_st,i_H_psi_array)

  h = diag_H_mat_elem_fock(det_ref,det_pert,fock_diag_tmp,Nint)
  do i =1,N_st
    if (i_H_psi_array(i) /= 0.d0) then
      delta_e = h - psi_energy(i)
      if (delta_e > 0.d0) then
        e_2_pert(i) = 0.5d0 * (delta_e - dsqrt(delta_e * delta_e + 4.d0 * i_H_psi_array(i) * i_H_psi_array(i)))
      else
        e_2_pert(i) = 0.5d0 * (delta_e + dsqrt(delta_e * delta_e + 4.d0 * i_H_psi_array(i) * i_H_psi_array(i)))
      endif
      if (dabs(i_H_psi_array(i)) > 1.d-6) then
        c_pert(i) = e_2_pert(i)/i_H_psi_array(i)
      else
        c_pert(i) = 0.d0
      endif
      H_pert_diag(i) = h*c_pert(i)*c_pert(i)
    else
      e_2_pert(i) = 0.d0
      c_pert(i) = 0.d0
      H_pert_diag(i) = 0.d0
    endif
  enddo

end



subroutine pt2_moller_plesset ($arguments)
  use bitmasks
  implicit none
  $declarations

  BEGIN_DOC
  ! Computes the standard Moller-Plesset perturbative first order coefficient and second
  ! order energetic contribution for the various N_st states.
  !
  ! `c_pert(i)` = $\\frac{\\langle i|H|\\alpha \\rangle}{\\text{difference of orbital energies}}$.
  !
  ! `e_2_pert(i)` = $\\frac{\\langle i|H|\\alpha \\rangle^2}{\\text{difference of orbital energies}}$.
  !
  END_DOC

  integer                        :: i,j
  double precision               :: diag_H_mat_elem_fock
  integer                        :: exc(0:2,2,2)
  integer                        :: degree
  double precision               :: phase,delta_e,h
  double precision               :: i_H_psi_array(N_st)
  integer                        :: h1,h2,p1,p2,s1,s2
  ASSERT (Nint == N_int)
  ASSERT (Nint > 0)
  call get_excitation(ref_bitmask,det_pert,exc,degree,phase,Nint)
  if (degree == 2) then
    call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
    delta_e = (Fock_matrix_diag_mo(h1) - Fock_matrix_diag_mo(p1)) + &
              (Fock_matrix_diag_mo(h2) - Fock_matrix_diag_mo(p2))
  else if (degree == 1) then
    call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
    delta_e = Fock_matrix_diag_mo(h1) - Fock_matrix_diag_mo(p1)
  else
    delta_e = 0.d0
  endif

  if (dabs(delta_e) > 1.d-10) then
        delta_e = 1.d0/delta_e
    call i_H_psi_minilist(det_pert,minilist,idx_minilist,N_minilist,psi_selectors_coef,Nint,N_minilist,psi_selectors_size,N_st,i_H_psi_array)
    h = diag_H_mat_elem_fock(det_ref,det_pert,fock_diag_tmp,Nint)
  else
    i_H_psi_array(:) = 0.d0
    h = 0.d0
  endif
  do i =1,N_st
    H_pert_diag(i) = h
    c_pert(i) = i_H_psi_array(i) *delta_e
    e_2_pert(i) = c_pert(i) * i_H_psi_array(i)
  enddo

end

subroutine pt2_dummy ($arguments)
  use bitmasks
  implicit none
  $declarations

  BEGIN_DOC
  ! Dummy perturbation to add all connected determinants.
  END_DOC

  integer                        :: i,j
  double precision               :: diag_H_mat_elem_fock, h
  double precision               :: i_H_psi_array(N_st)
  PROVIDE  selection_criterion

  call i_H_psi_minilist(det_pert,minilist,idx_minilist,N_minilist,psi_selectors_coef,Nint,N_minilist,psi_selectors_size,N_st,i_H_psi_array)

  h = diag_H_mat_elem_fock(det_ref,det_pert,fock_diag_tmp,Nint)
  do i =1,N_st
    if (i_H_psi_array(i) /= 0.d0) then
     c_pert(i) = i_H_psi_array(i) / (electronic_energy(i) - h)
     H_pert_diag(i) = h*c_pert(i)*c_pert(i)
     e_2_pert(i) = 1.d0
    else
      c_pert(i) = 0.d0
      e_2_pert(i) = 0.d0
      H_pert_diag(i) = 0.d0
    endif
  enddo

end



SUBST [ arguments, declarations ]

electronic_energy,det_ref,det_pert,fock_diag_tmp,c_pert,e_2_pert,H_pert_diag,Nint,ndet,N_st,minilist,idx_minilist,N_minilist ;

    integer, intent(in)             :: Nint
    integer, intent(in)             :: ndet
    integer, intent(in)             :: N_st
    integer, intent(in)             :: N_minilist
    integer(bit_kind), intent(in)   :: det_ref (Nint,2)
    integer(bit_kind), intent(in)   :: det_pert(Nint,2)
    double precision , intent(in)   :: fock_diag_tmp(2,mo_num+1)
    double precision , intent(in)    :: electronic_energy(N_st)
    double precision , intent(out)  :: c_pert(N_st)
    double precision , intent(out)  :: e_2_pert(N_st)
    double precision, intent(out)   :: H_pert_diag(N_st)
    integer, intent(in)             :: idx_minilist(0:N_det_selectors)
    integer(bit_kind), intent(in)   :: minilist(Nint,2,N_det_selectors)
;;


END_TEMPLATE

! Note : If the arguments are changed here, they should also be changed accordingly in
! the perturbation.template.f file.

