program simple_cis
 implicit none
 BEGIN_DOC
! Program that extracts the :option:`determinants n_states` lowest
! states of the Hamiltonian within the set of Slater determinants stored
! in the |EZFIO| directory.
!
! If :option:`determinants s2_eig` = |true|, it will retain only states
! which correspond to the desired value of
! :option:`determinants expected_s2`.
!
 END_DOC
 !read_wf = .True.
 !touch read_wf
 call add_singles
 call routine
end

subroutine add_singles
  implicit none
  !truncate to 1 det and add singles?
  integer(bit_kind) :: refdet(N_int,2), refdet_neg(N_int,2), hmask(N_int,2), pmask(N_int,2)
  integer(bit_kind), allocatable :: newdets(:,:,:)
  integer :: i, ki, ispin, ih, ip, idet
  integer :: refdetlist(N_int*bit_kind_size,2), n_elements(2)
  integer :: hlist(N_int*bit_kind_size,2), n_el_h(2)
  integer :: plist(N_int*bit_kind_size,2), n_el_p(2)
  integer :: nexc(2), nh(2), np(2), ndet_cis
  integer :: exc(0:2,2,2)
  logical :: ok
  complex*16, allocatable :: newcoef(:,:)


  refdet = psi_det_sorted(:,:,1)
  do i=1,N_int
    refdet_neg(i,1) = not(refdet(i,1))
    refdet_neg(i,2) = not(refdet(i,2))
  enddo
  call bitstring_to_list_ab(refdet, refdetlist, n_elements, N_int)
  
  ! count ndet for CIS
  nexc = 0
  do ki = 1, kpt_num
    nh = 0
    np = 0
    do i = 1, N_int
      do ispin = 1, 2
        hmask(i,ispin) = iand(reunion_of_inact_act_bitmask_kpts(i,ispin,ki),refdet(i,ispin))
        pmask(i,ispin) = iand(reunion_of_act_virt_bitmask_kpts(i,ispin,ki),refdet_neg(i,ispin))
        nh(ispin) += popcnt(hmask(i,ispin))
        np(ispin) += popcnt(pmask(i,ispin))
      enddo
    enddo
    do ispin = 1, 2
      nexc(ispin) += nh(ispin)*np(ispin)
    enddo
  enddo
  ndet_cis = 1 + nexc(1) + nexc(2)

  allocate(newdets(N_int,2,ndet_cis), newcoef(ndet_cis,n_states))
  idet = 1
  newdets(:,:,1) = refdet
  do ki = 1, kpt_num
    nh = 0
    np = 0
    do i = 1, N_int
      do ispin = 1, 2
        hmask(i,ispin) = iand(reunion_of_inact_act_bitmask_kpts(i,ispin,ki),refdet(i,ispin))
        pmask(i,ispin) = iand(reunion_of_act_virt_bitmask_kpts(i,ispin,ki),refdet_neg(i,ispin))
        nh(ispin) += popcnt(hmask(i,ispin))
        np(ispin) += popcnt(pmask(i,ispin))
      enddo
    enddo
    call bitstring_to_list_ab(hmask, hlist, n_el_h, N_int)
    call bitstring_to_list_ab(pmask, plist, n_el_p, N_int)
    do ispin = 1, 2
      do ih = 1, n_el_h(ispin)
        do ip = 1, n_el_p(ispin)
          exc = 0
          exc(0,1,ispin) = 1
          exc(0,2,ispin) = 1
          exc(1,2,ispin) = plist(ip,ispin)
          exc(1,1,ispin) = hlist(ih,ispin)
          idet += 1
          call apply_excitation(refdet, exc, newdets(:,:,idet), ok, N_int)
          if (.not.ok) then
            print *, irp_here, 'exc not ok'
            STOP -1
          endif
        enddo
      enddo
    enddo
  enddo
  newcoef = 0.d0
  do i = 1, n_states
    newcoef(i,i) = 1.d0
  enddo

  call save_wavefunction_general_complex(ndet_cis, n_states, newdets, ndet_cis, newcoef)

end

subroutine routine
 implicit none
 read_wf = .True.
 touch read_wf
 call diagonalize_ci
 print*,'N_det = ',N_det
 if (is_complex) then
  call save_wavefunction_general_complex(N_det,N_states,psi_det_sorted,size(psi_coef_sorted_complex,1),psi_coef_sorted_complex)
 else
  call save_wavefunction_general(N_det,N_states,psi_det_sorted,size(psi_coef_sorted,1),psi_coef_sorted)
 endif
end
