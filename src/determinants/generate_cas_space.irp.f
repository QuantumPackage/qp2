subroutine generate_cas_space
  use bitmasks
  implicit none
  BEGIN_DOC
! Generates the CAS space
  END_DOC
  integer :: i, n_alpha_act, n_beta_act, N_det_tot
  integer(bit_kind) :: o(N_int)
  integer :: mo_list(elec_alpha_num)

  integer :: k,n,m
  integer(bit_kind) :: u, t, t1, t2
  integer(bit_kind), parameter :: one = 1_bit_kind
  integer(bit_kind), allocatable :: det_a(:,:), det_b(:,:)
  integer(bit_kind), allocatable :: psi_det_local(:,:,:)
  double precision, allocatable  :: psi_coef_local(:,:)

  call list_to_bitstring(o, list_core_inact, n_core_inact_orb, N_int)

  ! Count number of active electrons
  n_alpha_act = 0
  n_beta_act = 0
  do i=1, n_act_orb
    if (list_act(i) <= elec_alpha_num) then
      n_alpha_act += 1
    endif
    if (list_act(i) <= elec_beta_num) then
      n_beta_act += 1
    endif
  enddo
  if (n_act_orb > 64) then
     stop 'More than 64 active MOs'
  endif

  print *, ''
  print *, 'CAS(', n_alpha_act+n_beta_act, ', ', n_act_orb, ')'
  print *, ''

  n_det_alpha_unique = int(binom_int(n_act_orb, n_alpha_act),4)
  n_det_beta_unique  = int(binom_int(n_act_orb, n_beta_act),4)
  if (binom_int(n_act_orb, n_alpha_act) > huge(1) .or. &
      binom_int(n_act_orb, n_beta_act) > huge(1)) then
    print *, n_act_orb, n_alpha_act, n_beta_act
    print *,  irp_here, ': Too many determinants to generate'
    stop 1
  endif
  N_det_tot = n_det_alpha_unique * n_det_beta_unique

  allocate(det_a(N_int,n_det_alpha_unique), det_b(N_int,n_det_beta_unique))

  ! All alpha determinants: weight n_alpha_act over the n_act_orb active MOs
  n = n_alpha_act
  u = shiftl(one,n) - one
  k = 0
  do while (u < shiftl(one,n_act_orb))
    k += 1
    call bitstring_to_list(u, mo_list, m, 1)
    do i=1,m
      mo_list(i) = list_act( mo_list(i) )
    enddo
    call list_to_bitstring(det_a(1,k), mo_list, m, N_int)
    do i=1,N_int
      det_a(i,k) = ior(det_a(i,k), o(i))
    enddo
    t = ior(u,u-1)
    t1 = t+1
    t2 = shiftr((iand(not(t),t1)-1), trailz(u)+1)
    u = ior(t1,t2)
  enddo
  if (k /= n_det_alpha_unique) stop 'CAS alpha determinant count mismatch'

  ! All beta determinants (identical pattern with n_beta_act / det_b)
  n = n_beta_act
  u = shiftl(one,n) - one
  k = 0
  do while (u < shiftl(one,n_act_orb))
    k += 1
    call bitstring_to_list(u, mo_list, m, 1)
    do i=1,m
      mo_list(i) = list_act( mo_list(i) )
    enddo
    call list_to_bitstring(det_b(1,k), mo_list, m, N_int)
    do i=1,N_int
      det_b(i,k) = ior(det_b(i,k), o(i))
    enddo
    t = ior(u,u-1)
    t1 = t+1
    t2 = shiftr((iand(not(t),t1)-1), trailz(u)+1)
    u = ior(t1,t2)
  enddo
  if (k /= n_det_beta_unique) stop 'CAS beta determinant count mismatch'

  ! Assemble the complete wave function (alpha-major x beta)
  allocate(psi_det_local(N_int,2,N_det_tot))
  allocate(psi_coef_local(N_det_tot,N_states))
  psi_coef_local = 1.d0
  k = 0
  do m=1,n_det_beta_unique
    do n=1,n_det_alpha_unique
      k += 1
      psi_det_local(:,1,k) = det_a(:,n)
      psi_det_local(:,2,k) = det_b(:,m)
    enddo
  enddo

  ! Persist to the EZFIO determinants group, then update the in-memory WF
  call save_wavefunction_general(N_det_tot,N_states,psi_det_local,&
        N_det_tot,psi_coef_local)
  read_wf = .True.
  N_det = N_det_tot
  psi_det_size = N_det_tot
  psi_det = psi_det_local
  psi_coef = psi_coef_local
  TOUCH N_det psi_det_size psi_det psi_coef

  print *, 'Ndet = ', N_det

end
