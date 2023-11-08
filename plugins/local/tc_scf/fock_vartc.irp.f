
! ---

 BEGIN_PROVIDER [ double precision, two_e_vartc_integral_alpha, (ao_num, ao_num)]
&BEGIN_PROVIDER [ double precision, two_e_vartc_integral_beta , (ao_num, ao_num)]

  implicit none
  integer                       :: i, j, k, l
  double precision              :: density, density_a, density_b, I_coul, I_kjli
  double precision              :: t0, t1
  double precision, allocatable :: tmp_a(:,:), tmp_b(:,:)

  two_e_vartc_integral_alpha = 0.d0
  two_e_vartc_integral_beta  = 0.d0

 !$OMP PARALLEL DEFAULT (NONE)                                                                           &
 !$OMP PRIVATE (i, j, k, l, density_a, density_b, density, tmp_a, tmp_b, I_coul, I_kjli)                 &
 !$OMP SHARED  (ao_num, TCSCF_density_matrix_ao_alpha, TCSCF_density_matrix_ao_beta, ao_two_e_vartc_tot, &
 !$OMP         two_e_vartc_integral_alpha, two_e_vartc_integral_beta)

  allocate(tmp_a(ao_num,ao_num), tmp_b(ao_num,ao_num))
  tmp_a = 0.d0
  tmp_b = 0.d0

 !$OMP DO
  do j = 1, ao_num
    do l = 1, ao_num
      density_a = TCSCF_density_matrix_ao_alpha(l,j)
      density_b = TCSCF_density_matrix_ao_beta (l,j)
      density   = density_a + density_b                      
      do i = 1, ao_num
        do k = 1, ao_num

          I_coul = density * ao_two_e_vartc_tot(k,i,l,j)
          I_kjli = ao_two_e_vartc_tot(k,j,l,i)

          tmp_a(k,i) += I_coul - density_a * I_kjli
          tmp_b(k,i) += I_coul - density_b * I_kjli
        enddo
      enddo
    enddo
  enddo
 !$OMP END DO NOWAIT

 !$OMP CRITICAL
  do i = 1, ao_num
    do j = 1, ao_num
      two_e_vartc_integral_alpha(j,i) += tmp_a(j,i)
      two_e_vartc_integral_beta (j,i) += tmp_b(j,i)
    enddo
  enddo
 !$OMP END CRITICAL

  deallocate(tmp_a, tmp_b)
 !$OMP END PARALLEL

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, Fock_matrix_vartc_ao_alpha, (ao_num, ao_num)]

  implicit none

  Fock_matrix_vartc_ao_alpha = ao_one_e_integrals_tc_tot + two_e_vartc_integral_alpha 

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, Fock_matrix_vartc_ao_beta, (ao_num, ao_num)]

  implicit none

  Fock_matrix_vartc_ao_beta = ao_one_e_integrals_tc_tot + two_e_vartc_integral_beta 

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, Fock_matrix_vartc_mo_alpha, (mo_num, mo_num) ]

  implicit none

  call ao_to_mo_bi_ortho( Fock_matrix_vartc_ao_alpha, size(Fock_matrix_vartc_ao_alpha, 1) &
                        , Fock_matrix_vartc_mo_alpha, size(Fock_matrix_vartc_mo_alpha, 1) )
  if(three_body_h_tc) then
    Fock_matrix_vartc_mo_alpha += fock_3e_uhf_mo_a
  endif

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, Fock_matrix_vartc_mo_beta, (mo_num,mo_num) ]

  implicit none

  call ao_to_mo_bi_ortho( Fock_matrix_vartc_ao_beta, size(Fock_matrix_vartc_ao_beta, 1) &
                        , Fock_matrix_vartc_mo_beta, size(Fock_matrix_vartc_mo_beta, 1) )
  if(three_body_h_tc) then
    Fock_matrix_vartc_mo_beta += fock_3e_uhf_mo_b
  endif

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, grad_vartc]

  implicit none
  integer          :: i, k
  double precision :: grad_left, grad_right

  grad_left  = 0.d0
  grad_right = 0.d0

  do i = 1, elec_beta_num ! doc --> SOMO
    do k = elec_beta_num+1, elec_alpha_num
      grad_left  = max(grad_left , dabs(Fock_matrix_vartc_mo_tot(k,i)))
      grad_right = max(grad_right, dabs(Fock_matrix_vartc_mo_tot(i,k)))
    enddo
  enddo

  do i = 1, elec_beta_num ! doc --> virt 
    do k = elec_alpha_num+1, mo_num
      grad_left  = max(grad_left , dabs(Fock_matrix_vartc_mo_tot(k,i)))
      grad_right = max(grad_right, dabs(Fock_matrix_vartc_mo_tot(i,k)))
    enddo
  enddo

  do i = elec_beta_num+1, elec_alpha_num ! SOMO --> virt 
    do k = elec_alpha_num+1, mo_num
      grad_left  = max(grad_left , dabs(Fock_matrix_vartc_mo_tot(k,i)))
      grad_right = max(grad_right, dabs(Fock_matrix_vartc_mo_tot(i,k)))
    enddo
  enddo

  grad_vartc = grad_left + grad_right

END_PROVIDER 

! ---

BEGIN_PROVIDER [ double precision, Fock_matrix_vartc_ao_tot, (ao_num, ao_num) ]

  implicit none

  call mo_to_ao_bi_ortho( Fock_matrix_vartc_mo_tot, size(Fock_matrix_vartc_mo_tot, 1) &
                        , Fock_matrix_vartc_ao_tot, size(Fock_matrix_vartc_ao_tot, 1) )

END_PROVIDER

! ---

 BEGIN_PROVIDER [ double precision, Fock_matrix_vartc_mo_tot, (mo_num,mo_num) ]
&BEGIN_PROVIDER [ double precision, Fock_matrix_vartc_diag_mo_tot, (mo_num)]

  implicit none
  integer :: i, j, n

  if(elec_alpha_num == elec_beta_num) then
    Fock_matrix_vartc_mo_tot = Fock_matrix_vartc_mo_alpha
  else

    do j = 1, elec_beta_num
      ! F-K
      do i = 1, elec_beta_num !CC
         Fock_matrix_vartc_mo_tot(i,j) = 0.5d0*(Fock_matrix_vartc_mo_alpha(i,j)+Fock_matrix_vartc_mo_beta(i,j))&
            - (Fock_matrix_vartc_mo_beta(i,j) - Fock_matrix_vartc_mo_alpha(i,j))
      enddo
      ! F+K/2
      do i = elec_beta_num+1, elec_alpha_num  !CA
         Fock_matrix_vartc_mo_tot(i,j) = 0.5d0*(Fock_matrix_vartc_mo_alpha(i,j)+Fock_matrix_vartc_mo_beta(i,j))&
             + 0.5d0*(Fock_matrix_vartc_mo_beta(i,j) - Fock_matrix_vartc_mo_alpha(i,j))
      enddo
      ! F
      do i = elec_alpha_num+1, mo_num !CV
         Fock_matrix_vartc_mo_tot(i,j) = 0.5d0*(Fock_matrix_vartc_mo_alpha(i,j)+Fock_matrix_vartc_mo_beta(i,j))
      enddo
    enddo

    do j = elec_beta_num+1, elec_alpha_num
      ! F+K/2
      do i = 1, elec_beta_num !AC
         Fock_matrix_vartc_mo_tot(i,j) = 0.5d0*(Fock_matrix_vartc_mo_alpha(i,j)+Fock_matrix_vartc_mo_beta(i,j))&
             + 0.5d0*(Fock_matrix_vartc_mo_beta(i,j) - Fock_matrix_vartc_mo_alpha(i,j))
      enddo
      ! F
      do i = elec_beta_num+1, elec_alpha_num !AA
         Fock_matrix_vartc_mo_tot(i,j) = 0.5d0*(Fock_matrix_vartc_mo_alpha(i,j)+Fock_matrix_vartc_mo_beta(i,j))
      enddo
      ! F-K/2
      do i = elec_alpha_num+1, mo_num !AV
         Fock_matrix_vartc_mo_tot(i,j) = 0.5d0*(Fock_matrix_vartc_mo_alpha(i,j)+Fock_matrix_vartc_mo_beta(i,j))&
             - 0.5d0*(Fock_matrix_vartc_mo_beta(i,j) - Fock_matrix_vartc_mo_alpha(i,j))
      enddo
    enddo
    
    do j = elec_alpha_num+1, mo_num
      ! F
      do i = 1, elec_beta_num !VC
         Fock_matrix_vartc_mo_tot(i,j) = 0.5d0*(Fock_matrix_vartc_mo_alpha(i,j)+Fock_matrix_vartc_mo_beta(i,j))
      enddo
      ! F-K/2
      do i = elec_beta_num+1, elec_alpha_num !VA
         Fock_matrix_vartc_mo_tot(i,j) = 0.5d0*(Fock_matrix_vartc_mo_alpha(i,j)+Fock_matrix_vartc_mo_beta(i,j))&
             - 0.5d0*(Fock_matrix_vartc_mo_beta(i,j) - Fock_matrix_vartc_mo_alpha(i,j))
      enddo
      ! F+K
      do i = elec_alpha_num+1, mo_num !VV
         Fock_matrix_vartc_mo_tot(i,j) = 0.5d0*(Fock_matrix_vartc_mo_alpha(i,j)+Fock_matrix_vartc_mo_beta(i,j)) &
             + (Fock_matrix_vartc_mo_beta(i,j) - Fock_matrix_vartc_mo_alpha(i,j))
      enddo
    enddo
    if(three_body_h_tc)then
      ! C-O
      do j = 1, elec_beta_num
        do i = elec_beta_num+1, elec_alpha_num
          Fock_matrix_vartc_mo_tot(i,j) += 0.5d0*(fock_a_tot_3e_bi_orth(i,j) + fock_b_tot_3e_bi_orth(i,j))
          Fock_matrix_vartc_mo_tot(j,i) += 0.5d0*(fock_a_tot_3e_bi_orth(j,i) + fock_b_tot_3e_bi_orth(j,i))
        enddo
      enddo
      ! C-V
      do j = 1, elec_beta_num
        do i = elec_alpha_num+1, mo_num
          Fock_matrix_vartc_mo_tot(i,j) += 0.5d0*(fock_a_tot_3e_bi_orth(i,j) + fock_b_tot_3e_bi_orth(i,j))
          Fock_matrix_vartc_mo_tot(j,i) += 0.5d0*(fock_a_tot_3e_bi_orth(j,i) + fock_b_tot_3e_bi_orth(j,i))
        enddo
      enddo
      ! O-V
      do j = elec_beta_num+1, elec_alpha_num
        do i = elec_alpha_num+1, mo_num
          Fock_matrix_vartc_mo_tot(i,j) += 0.5d0*(fock_a_tot_3e_bi_orth(i,j) + fock_b_tot_3e_bi_orth(i,j))
          Fock_matrix_vartc_mo_tot(j,i) += 0.5d0*(fock_a_tot_3e_bi_orth(j,i) + fock_b_tot_3e_bi_orth(j,i))
        enddo
      enddo
    endif

  endif

  do i = 1, mo_num
    Fock_matrix_vartc_diag_mo_tot(i) = Fock_matrix_vartc_mo_tot(i,i)
  enddo

  if(frozen_orb_scf)then
    integer :: iorb, jorb
    do i = 1, n_core_orb
     iorb = list_core(i)
     do j = 1, n_act_orb
      jorb = list_act(j)
      Fock_matrix_vartc_mo_tot(iorb,jorb) = 0.d0
      Fock_matrix_vartc_mo_tot(jorb,iorb) = 0.d0
     enddo
    enddo
  endif

  if(no_oa_or_av_opt)then
    do i = 1, n_act_orb
      iorb = list_act(i)
      do j = 1, n_inact_orb
        jorb = list_inact(j)
        Fock_matrix_vartc_mo_tot(iorb,jorb) = 0.d0
        Fock_matrix_vartc_mo_tot(jorb,iorb) = 0.d0
      enddo
      do j = 1, n_virt_orb
        jorb = list_virt(j)
        Fock_matrix_vartc_mo_tot(iorb,jorb) = 0.d0
        Fock_matrix_vartc_mo_tot(jorb,iorb) = 0.d0
      enddo
      do j = 1, n_core_orb
        jorb = list_core(j)
        Fock_matrix_vartc_mo_tot(iorb,jorb) = 0.d0
        Fock_matrix_vartc_mo_tot(jorb,iorb) = 0.d0                                                                                                                 
      enddo
    enddo
  endif

  !call check_sym(Fock_matrix_vartc_mo_tot, mo_num)
  !do i = 1, mo_num
  !  write(*,'(100(F15.8, I4))') Fock_matrix_vartc_mo_tot(i,:)
  !enddo

END_PROVIDER

! ---

