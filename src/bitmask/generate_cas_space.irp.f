subroutine generate_cas_space
  use bitmasks
  implicit none
  BEGIN_DOC
! Generates the CAS space
  END_DOC
  integer :: i, sze, ncore, n_alpha_act, n_beta_act
  integer(bit_kind) :: o(N_int)
  integer(bit_kind) :: u
  integer :: mo_list(elec_alpha_num)

  integer :: k,n,m
  integer(bit_kind) :: t, t1, t2

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

  n_det_alpha_unique = binom_int(n_act_orb, n_alpha_act)
  TOUCH n_det_alpha_unique

  n = n_alpha_act
  u = shiftl(1_bit_kind,n) - 1_bit_kind

  k=0
  do while (u < shiftl(1_bit_kind,n_act_orb))
    k = k+1
    call bitstring_to_list(u, mo_list, m, 1)
    do i=1,m
      mo_list(i) = list_act( mo_list(i) )
    enddo
    call list_to_bitstring(psi_det_alpha_unique(1,k), mo_list, m, N_int)
    do i=1,N_int
      psi_det_alpha_unique(i,k) = ior(psi_det_alpha_unique(i,k), o(i))
    enddo
    t = ior(u,u-1)
    t1 = t+1
    t2 = shiftr((iand(not(t),t1)-1), trailz(u)+1)
    u = ior(t1,t2)
  enddo

  n_det_beta_unique = binom_int(n_act_orb, n_beta_act)
  TOUCH n_det_beta_unique

  n = n_beta_act
  u = shiftl(1_bit_kind,n) -1_bit_kind

  k=0
  do while (u < shiftl(1_bit_kind,n_act_orb))
    k = k+1
    call bitstring_to_list(u, mo_list, m, 1)
    do i=1,m
      mo_list(i) = list_act( mo_list(i) )
    enddo
    call list_to_bitstring(psi_det_beta_unique(1,k), mo_list, m, N_int)
    do i=1,N_int
      psi_det_beta_unique(i,k) = ior(psi_det_beta_unique(i,k), o(i))
    enddo
    t = ior(u,u-1)
    t1 = t+1
    t2 = shiftr((iand(not(t),t1)-1), trailz(u)+1)
    u = ior(t1,t2)
  enddo

  call generate_all_alpha_beta_det_products

  print *, 'Ndet = ', N_det

end

