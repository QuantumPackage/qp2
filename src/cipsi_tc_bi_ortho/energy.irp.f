BEGIN_PROVIDER [ logical, initialize_pt2_E0_denominator ]
 implicit none
 BEGIN_DOC
 ! If true, initialize pt2_E0_denominator
 END_DOC
 initialize_pt2_E0_denominator = .True.
END_PROVIDER

BEGIN_PROVIDER [ double precision, pt2_E0_denominator, (N_states) ]
 implicit none
 BEGIN_DOC
 ! E0 in the denominator of the PT2
 END_DOC
 integer :: i,j

  pt2_E0_denominator = eigval_right_tc_bi_orth

! if (initialize_pt2_E0_denominator) then
!   if (h0_type == "EN") then
!     pt2_E0_denominator(1:N_states) = psi_energy(1:N_states)
!   else if (h0_type == "HF") then
!     do i=1,N_states
!       j = maxloc(abs(psi_coef(:,i)),1)
!       pt2_E0_denominator(i) = psi_det_hii(j)
!     enddo
!   else if (h0_type == "Barycentric") then
!     pt2_E0_denominator(1:N_states) = barycentric_electronic_energy(1:N_states)
!   else if (h0_type == "CFG") then
!     pt2_E0_denominator(1:N_states) = psi_energy(1:N_states)
!   else
!     print *,  h0_type, ' not implemented'
!     stop
!   endif
!  do i=1,N_states
!    call write_double(6,pt2_E0_denominator(i)+nuclear_repulsion, 'PT2 Energy denominator')
!  enddo
! else
!   pt2_E0_denominator = -huge(1.d0)
! endif

END_PROVIDER


BEGIN_PROVIDER [ double precision, pt2_overlap, (N_states, N_states) ]
 implicit none
 BEGIN_DOC
 ! Overlap between the perturbed wave functions
 END_DOC
 pt2_overlap(1:N_states,1:N_states) = 0.d0
END_PROVIDER

