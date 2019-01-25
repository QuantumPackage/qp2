BEGIN_PROVIDER [ double precision, mo_energy_expval, (N_states,mo_num,2,2)]
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Third index is spin.
  ! Fourth index is 1:creation, 2:annihilation
  END_DOC
  integer                        :: i,j,k
  integer                        :: ispin, istate
  integer                        :: hp
  double precision               :: norm_out(N_states)

  integer, parameter             :: hole_particle(2) = (/ -1, 1 /)
  double precision               :: energies(n_states)

  integer(bit_kind), allocatable :: psi_in_out(:,:,:)
  double precision, allocatable  :: psi_in_out_coef(:,:)
  double precision               :: E0(N_states), norm
  double precision, parameter    :: t=1.d-3

  allocate (psi_in_out(N_int,2,N_det),psi_in_out_coef(N_det,N_states))
  mo_energy_expval = 0.d0

  psi_in_out_coef(1:N_det,1:N_states) = psi_coef(1:N_det,1:N_states)
  psi_in_out(1:N_int,1:2,1:N_det) = psi_det(1:N_int,1:2,1:N_det)

  ! Truncate the wave function
  do istate=1,N_states
    norm = 0.d0
    do k=1,N_det
      if (dabs(psi_in_out_coef(k,istate)) < t) then
        psi_in_out_coef(k,istate) = 0.d0
      endif
      norm = norm + psi_in_out_coef(k,istate)*psi_in_out_coef(k,istate)
    enddo
    ASSERT (norm  > 0.d0)
    norm = 1.d0/dsqrt(norm)
    psi_in_out_coef(1:N_det,istate) = psi_in_out_coef(1:N_det,istate) * norm
    call au0_h_au0(E0,psi_in_out,psi_in_out_coef,N_det,size(psi_in_out_coef,1))
  enddo


  do hp=1,2
    do ispin=1,2
      do i=1,mo_num
        psi_in_out_coef(1:N_det,1:N_states) = psi_coef(1:N_det,1:N_states)
        psi_in_out(1:N_int,1:2,1:N_det) = psi_det(1:N_int,1:2,1:N_det)
        call apply_exc_to_psi(i,hole_particle(hp),ispin,             &
            norm_out,psi_in_out,psi_in_out_coef, N_det,N_det,N_det,N_states)

        ! Truncate the wave function
        do istate=1,N_states
          norm = 0.d0
          do k=1,N_det
            if (dabs(psi_in_out_coef(k,istate)) < t) then
              psi_in_out_coef(k,istate) = 0.d0
            endif
            norm = norm + psi_in_out_coef(k,istate)*psi_in_out_coef(k,istate)
          enddo
          if (norm == 0.d0) then
            cycle
          endif
          norm = 1.d0/dsqrt(norm)
          psi_in_out_coef(1:N_det,istate) = psi_in_out_coef(1:N_det,istate) * norm
        enddo
        call au0_h_au0(energies,psi_in_out,psi_in_out_coef,N_det,size(psi_in_out_coef,1))
        mo_energy_expval(1:N_states,i,ispin,hp) = energies(1:N_states) - E0(1:N_states)
        print *,  i, ispin, real(energies(1)), real(E0(1))
      enddo
    enddo

  enddo
  mo_energy_expval(1:N_states,1:mo_num,1:2,1) = -mo_energy_expval(1:N_states,1:mo_num,1:2,1)

END_PROVIDER


subroutine au0_h_au0(energies,psi_in,psi_in_coef,ndet,dim_psi_coef)
  use bitmasks
  implicit none
  integer, intent(in)            :: ndet,dim_psi_coef
  integer(bit_kind), intent(in)  :: psi_in(N_int,2,Ndet)
  double precision,  intent(in)  :: psi_in_coef(dim_psi_coef,N_states)
  double precision,  intent(out) :: energies(N_states)

  integer                        :: i,j, istate
  double precision               :: hij,accu
  double precision, allocatable  :: psi_coef_tmp(:)

  energies(1:N_states) = 0.d0
  do i = 1, Ndet
    if(sum(dabs(psi_in_coef(i,1:N_states)))==0.d0) then
      cycle
    endif
    call diag_H_mat_elem_au0_h_au0(psi_in(1,1,i),N_int,hij)
    do istate=1,N_states
      energies(istate) += psi_in_coef(i,istate) * psi_in_coef(i,istate) * hij
    enddo
    do j = i+1, Ndet
      if(sum(dabs(psi_in_coef(j,1:N_states)))==0.d0) then
        cycle
      endif
      call i_H_j(psi_in(1,1,i),psi_in(1,1,j),N_int,hij)
      hij = hij+hij
      do istate=1,N_states
        energies(istate) = energies(istate) + psi_in_coef(i,istate) * psi_in_coef(j,istate) * hij
      enddo
    enddo
  enddo
end

subroutine diag_H_mat_elem_au0_h_au0(det_in,Nint,hii)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes $\langle i|H|i \rangle$ for any determinant $|i\rangle$.
  ! Used for wave functions with an additional electron.
  END_DOC
  integer,intent(in)             :: Nint
  integer(bit_kind),intent(in)   :: det_in(Nint,2)
  double precision, intent(out)  :: hii

  integer                        :: i, j, iorb, jorb
  integer                        :: occ(Nint*bit_kind_size,2)
  integer                        :: elec_num_tab_local(2)

  hii = 0.d0
  call bitstring_to_list(det_in(1,1), occ(1,1), elec_num_tab_local(1), Nint)
  call bitstring_to_list(det_in(1,2), occ(1,2), elec_num_tab_local(2), Nint)

  ! alpha - alpha
  do i = 1, elec_num_tab_local(1)
    iorb =  occ(i,1)
    hii += mo_one_e_integrals(iorb,iorb)
    do j = i+1, elec_num_tab_local(1)
      jorb = occ(j,1)
      hii +=  mo_two_e_integrals_jj_anti(jorb,iorb)
    enddo
  enddo

  ! beta - beta
  do i = 1, elec_num_tab_local(2)
    iorb =  occ(i,2)
    hii += mo_one_e_integrals(iorb,iorb)
    do j = i+1, elec_num_tab_local(2)
      jorb = occ(j,2)
      hii +=  mo_two_e_integrals_jj_anti(jorb,iorb)
    enddo
  enddo

  ! alpha - beta
  do i = 1, elec_num_tab_local(2)
    iorb =  occ(i,2)
    do j = 1, elec_num_tab_local(1)
      jorb = occ(j,1)
      hii +=  mo_two_e_integrals_jj(jorb,iorb)
    enddo
  enddo

end
