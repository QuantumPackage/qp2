
! ---

  use bitmasks

! ---

BEGIN_PROVIDER [integer, index_HF_psi_det]                                                                                                            

  implicit none
  integer :: i, degree

  do i = 1, N_det
    call get_excitation_degree(HF_bitmask, psi_det(1,1,i), degree, N_int)
    if(degree == 0) then
      index_HF_psi_det = i
      exit
    endif
  enddo

END_PROVIDER

! ---

subroutine diagonalize_CI_tc()

  BEGIN_DOC
  !  Replace the coefficients of the |CI| states by the coefficients of the
  !  eigenstates of the |CI| matrix.
  END_DOC

  implicit none
  integer :: i, j

  do j = 1, N_states
    do i = 1, N_det
      psi_l_coef_bi_ortho(i,j) = leigvec_tc_bi_orth(i,j)
      psi_r_coef_bi_ortho(i,j) = reigvec_tc_bi_orth(i,j)
    enddo
  enddo

  SOFT_TOUCH psi_l_coef_bi_ortho psi_r_coef_bi_ortho

end

! ---

 BEGIN_PROVIDER [double precision, eigval_right_tc_bi_orth, (N_states)      ]
&BEGIN_PROVIDER [double precision, eigval_left_tc_bi_orth , (N_states)      ]
&BEGIN_PROVIDER [double precision, reigvec_tc_bi_orth     , (N_det,N_states)]
&BEGIN_PROVIDER [double precision, leigvec_tc_bi_orth     , (N_det,N_states)]
&BEGIN_PROVIDER [double precision, s2_eigvec_tc_bi_orth   , (N_states)      ]
&BEGIN_PROVIDER [double precision, norm_ground_left_right_bi_orth           ]

  BEGIN_DOC
  ! eigenvalues, right and left eigenvectors of the transcorrelated Hamiltonian on the BI-ORTHO basis 
  END_DOC

  implicit none
  integer                       :: i, idx_dress, j, istate, k
  integer                       :: i_good_state, i_other_state, i_state
  integer                       :: n_real_tc_bi_orth_eigval_right, igood_r, igood_l
  logical                       :: converged, dagger
  double precision, parameter   :: alpha = 0.1d0
  integer,          allocatable :: index_good_state_array(:)
  integer,          allocatable :: iorder(:)
  logical,          allocatable :: good_state_array(:)
  double precision, allocatable :: reigvec_tc_bi_orth_tmp(:,:), leigvec_tc_bi_orth_tmp(:,:),eigval_right_tmp(:)
  double precision, allocatable :: s2_values_tmp(:), H_prime(:,:), expect_e(:)
  double precision, allocatable :: coef_hf_r(:),coef_hf_l(:)
  double precision, allocatable :: Stmp(:,:)

  PROVIDE N_det N_int

  if(N_det .le. N_det_max_full) then

    allocate(reigvec_tc_bi_orth_tmp(N_det,N_det), leigvec_tc_bi_orth_tmp(N_det,N_det), eigval_right_tmp(N_det), expect_e(N_det))
    allocate(H_prime(N_det,N_det), s2_values_tmp(N_det))

    H_prime(1:N_det,1:N_det) = htilde_matrix_elmt_bi_ortho(1:N_det,1:N_det)
    if(s2_eig) then
      H_prime(1:N_det,1:N_det) += alpha * S2_matrix_all_dets(1:N_det,1:N_det)
      do j = 1, N_det
        H_prime(j,j) = H_prime(j,j) - alpha*expected_s2
      enddo
    endif

    call non_hrmt_real_diag(N_det, H_prime, leigvec_tc_bi_orth_tmp, reigvec_tc_bi_orth_tmp, n_real_tc_bi_orth_eigval_right, eigval_right_tmp)
!    do i = 1, N_det
!     call get_H_tc_s2_l0_r0(leigvec_tc_bi_orth_tmp(1,i),reigvec_tc_bi_orth_tmp(1,i),1,N_det,expect_e(i), s2_values_tmp(i))
!    enddo
    call get_H_tc_s2_l0_r0(leigvec_tc_bi_orth_tmp,reigvec_tc_bi_orth_tmp,N_det,N_det,expect_e, s2_values_tmp)

    allocate(index_good_state_array(N_det),good_state_array(N_det))
    i_state = 0
    good_state_array = .False.

    if(s2_eig) then

      if(only_expected_s2) then
        do j = 1, N_det
         ! Select at least n_states states with S^2 values closed to "expected_s2"
!         print*,'s2_values_tmp(j) = ',s2_values_tmp(j),eigval_right_tmp(j),expect_e(j)
          if(dabs(s2_values_tmp(j) - expected_s2).le.0.5d0)then
            i_state +=1
            index_good_state_array(i_state) = j
            good_state_array(j) = .True.
          endif
          if(i_state.eq.N_states) then
            exit
          endif
        enddo
      else
        do j = 1, N_det
          index_good_state_array(j) = j
          good_state_array(j) = .True.
        enddo
      endif

      if(i_state .ne. 0) then
        ! Fill the first "i_state" states that have a correct S^2 value
        do j = 1, i_state
          do i = 1, N_det
            reigvec_tc_bi_orth(i,j) = reigvec_tc_bi_orth_tmp(i,index_good_state_array(j))
            leigvec_tc_bi_orth(i,j) = leigvec_tc_bi_orth_tmp(i,index_good_state_array(j))
          enddo
          eigval_right_tc_bi_orth(j) = expect_e(index_good_state_array(j))
          eigval_left_tc_bi_orth(j)  = expect_e(index_good_state_array(j))
          s2_eigvec_tc_bi_orth(j)    = s2_values_tmp(index_good_state_array(j))
        enddo
        i_other_state = 0
        do j = 1, N_det
          if(good_state_array(j))cycle
          i_other_state +=1
          if(i_state+i_other_state.gt.n_states)then
            exit
          endif
          do i = 1, N_det
            reigvec_tc_bi_orth(i,i_state+i_other_state) = reigvec_tc_bi_orth_tmp(i,j)
            leigvec_tc_bi_orth(i,i_state+i_other_state) = leigvec_tc_bi_orth_tmp(i,j)
          enddo
          eigval_right_tc_bi_orth(i_state+i_other_state) = eigval_right_tmp(j)
          eigval_left_tc_bi_orth (i_state+i_other_state) = eigval_right_tmp(j)
          s2_eigvec_tc_bi_orth(i_state+i_other_state)    = s2_values_tmp(i_state+i_other_state)
        enddo
      else ! istate == 0
        print*,''
        print*,'!!!!!!!!   WARNING  !!!!!!!!!'
        print*,'  Within the ',N_det,'determinants selected'
        print*,'  and the ',N_states_diag,'states requested'
        print*,'  We did not find only states with S^2 values close to ',expected_s2
        print*,'  We will then set the first N_states eigenvectors of the H matrix'
        print*,'  as the CI_eigenvectors'
        print*,'  You should consider more states and maybe ask for s2_eig to be .True. or just enlarge the CI space'
        print*,''
        do j = 1, min(N_states_diag, N_det)
          do i = 1, N_det
            leigvec_tc_bi_orth(i,j) = leigvec_tc_bi_orth_tmp(i,j)
            reigvec_tc_bi_orth(i,j) = reigvec_tc_bi_orth_tmp(i,j)
          enddo
          eigval_right_tc_bi_orth(j) = eigval_right_tmp(j)
          eigval_left_tc_bi_orth (j) = eigval_right_tmp(j)
          s2_eigvec_tc_bi_orth(j)    = s2_values_tmp(j)
        enddo
      endif ! istate .ne. 0

    else ! s2_eig

      allocate(coef_hf_r(N_det),coef_hf_l(N_det),iorder(N_det))
      do i = 1,N_det
        iorder(i) = i
        coef_hf_r(i) = -dabs(reigvec_tc_bi_orth_tmp(index_HF_psi_det,i))
      enddo
      call dsort(coef_hf_r,iorder,N_det)
      igood_r = iorder(1)
      print*,'igood_r, coef_hf_r = ',igood_r,coef_hf_r(1)
      do i = 1,N_det
        iorder(i) = i
        coef_hf_l(i) = -dabs(leigvec_tc_bi_orth_tmp(index_HF_psi_det,i))
      enddo
      call dsort(coef_hf_l,iorder,N_det)
      igood_l = iorder(1)
      print*,'igood_l, coef_hf_l = ',igood_l,coef_hf_l(1)
       
      if(igood_r.ne.igood_l .and. igood_r.ne.1) then
        print *,''
        print *,'Warning, the left and right eigenvectors are "not the same" '
        print *,'Warning, the ground state is not dominated by HF...'
        print *,'State with largest RIGHT coefficient of HF ',igood_r
        print *,'coef of HF in RIGHT eigenvector = ',reigvec_tc_bi_orth_tmp(index_HF_psi_det,igood_r)
        print *,'State with largest LEFT  coefficient of HF ',igood_l
        print *,'coef of HF in LEFT  eigenvector = ',leigvec_tc_bi_orth_tmp(index_HF_psi_det,igood_l)
      endif

      if(state_following_tc) then
        print *,'Following the states with the largest coef on HF'
        print *,'igood_r,igood_l',igood_r,igood_l
        i = igood_r
        eigval_right_tc_bi_orth(1) = eigval_right_tmp(i)
        do j = 1, N_det
          reigvec_tc_bi_orth(j,1) = reigvec_tc_bi_orth_tmp(j,i)
        enddo
        i = igood_l
        eigval_left_tc_bi_orth(1)  = eigval_right_tmp(i)
        do j = 1, N_det
          leigvec_tc_bi_orth(j,1) = leigvec_tc_bi_orth_tmp(j,i)
        enddo
      else 
        do i = 1, N_states
          eigval_right_tc_bi_orth(i) = eigval_right_tmp(i)
          eigval_left_tc_bi_orth(i)  = eigval_right_tmp(i)
          do j = 1, N_det
            reigvec_tc_bi_orth(j,i) = reigvec_tc_bi_orth_tmp(j,i)
            leigvec_tc_bi_orth(j,i) = leigvec_tc_bi_orth_tmp(j,i)
          enddo
        enddo
      endif

    endif

  else ! n_det > N_det_max_full

    double precision, allocatable :: H_jj(:),vec_tmp(:,:)
    external                         H_tc_u_0_opt
    external                         H_tc_dagger_u_0_opt
    external                         H_tc_s2_dagger_u_0_opt
    external                         H_tc_s2_u_0_opt
    external                         H_tc_s2_dagger_u_0_with_pure_three_omp
    external                         H_tc_s2_u_0_with_pure_three_omp

    allocate(H_jj(N_det),vec_tmp(N_det,n_states_diag))

    do i = 1, N_det
      call htilde_mu_mat_opt_bi_ortho_tot(psi_det(1,1,i), psi_det(1,1,i), N_int, H_jj(i))
    enddo

    print*,'---------------------------------'
    print*,'---------------------------------'
    print*,'Computing the left-eigenvector '
    print*,'---------------------------------'
    print*,'---------------------------------'
    !!!! Preparing the left-eigenvector
    vec_tmp = 0.d0
    do istate = 1, N_states
      vec_tmp(1:N_det,istate) = psi_l_coef_bi_ortho(1:N_det,istate)
    enddo
    do istate = N_states+1, n_states_diag
      vec_tmp(istate,istate) = 1.d0
    enddo
    integer :: n_it_max,i_it
    n_it_max = 1
    converged = .False.
    i_it = 0
    do while (.not.converged)
      if(.not.pure_three_body_h_tc)then
       call davidson_hs2_nonsym_b1space(vec_tmp, H_jj, s2_eigvec_tc_bi_orth, eigval_left_tc_bi_orth, N_det, n_states, n_states_diag, n_it_max, converged, H_tc_s2_dagger_u_0_opt)
      else 
       call davidson_hs2_nonsym_b1space(vec_tmp, H_jj, s2_eigvec_tc_bi_orth, eigval_left_tc_bi_orth, N_det, n_states, n_states_diag, n_it_max, converged, H_tc_s2_dagger_u_0_with_pure_three_omp)
      endif
      i_it += 1
      if(i_it .gt. 5) exit
    enddo
    do istate = 1, N_states
      leigvec_tc_bi_orth(1:N_det,istate) = vec_tmp(1:N_det,istate)
    enddo

    print*,'---------------------------------'
    print*,'---------------------------------'
    print*,'Computing the right-eigenvector '
    print*,'---------------------------------'
    print*,'---------------------------------'
    !!!! Preparing the right-eigenvector
    vec_tmp = 0.d0
    do istate = 1, N_states
      vec_tmp(1:N_det,istate) = psi_r_coef_bi_ortho(1:N_det,istate)
    enddo
    do istate = N_states+1, n_states_diag
      vec_tmp(istate,istate) = 1.d0
    enddo
    !call davidson_general_ext_rout_nonsym_b1space(vec_tmp, H_jj, eigval_right_tc_bi_orth, N_det, n_states, n_states_diag, converged, H_tc_u_0_opt)
    converged = .False.
    i_it = 0
    do while (.not. converged)
      if(.not.pure_three_body_h_tc)then
       call davidson_hs2_nonsym_b1space(vec_tmp, H_jj, s2_eigvec_tc_bi_orth, eigval_right_tc_bi_orth, N_det, n_states, n_states_diag, n_it_max, converged, H_tc_s2_u_0_opt)
      else
       call davidson_hs2_nonsym_b1space(vec_tmp, H_jj, s2_eigvec_tc_bi_orth, eigval_right_tc_bi_orth, N_det, n_states, n_states_diag, n_it_max, converged, H_tc_s2_u_0_with_pure_three_omp)
      endif
      i_it += 1
      if(i_it .gt. 5) exit
    enddo
    do istate = 1, N_states
      reigvec_tc_bi_orth(1:N_det,istate) = vec_tmp(1:N_det,istate)
    enddo

    deallocate(H_jj)
  endif

  call bi_normalize(leigvec_tc_bi_orth, reigvec_tc_bi_orth, size(reigvec_tc_bi_orth, 1), N_det, N_states)
  ! check bi-orthogonality
  allocate(Stmp(N_states,N_states))
  call dgemm( 'T', 'N', N_states, N_states, N_det, 1.d0                                                              &
        , leigvec_tc_bi_orth(1,1), size(leigvec_tc_bi_orth, 1), reigvec_tc_bi_orth(1,1), size(reigvec_tc_bi_orth, 1) &
        , 0.d0, Stmp(1,1), size(Stmp, 1) )
  print *, ' overlap matrix between states:'
  do i = 1, N_states
    write(*,'(1000(F16.10,X))') Stmp(i,:)
  enddo
  deallocate(Stmp)

  print*,'leigvec_tc_bi_orth(1,1),reigvec_tc_bi_orth(1,1) = ', leigvec_tc_bi_orth(1,1), reigvec_tc_bi_orth(1,1)
  do i = 1, N_states
    norm_ground_left_right_bi_orth = 0.d0
    do j = 1, N_det
      norm_ground_left_right_bi_orth += leigvec_tc_bi_orth(j,i) * reigvec_tc_bi_orth(j,i)
    enddo
    print*,' state      ', i
    print*,' norm l/r = ', norm_ground_left_right_bi_orth
    print*,' <S2>     = ', s2_eigvec_tc_bi_orth(i)
  enddo

  double precision, allocatable :: buffer(:,:)
  allocate(buffer(N_det,N_states))
  do k = 1, N_states
    do i = 1, N_det
      psi_l_coef_bi_ortho(i,k) = leigvec_tc_bi_orth(i,k)
      buffer(i,k) = leigvec_tc_bi_orth(i,k)
    enddo
  enddo
  TOUCH psi_l_coef_bi_ortho
  call ezfio_set_tc_bi_ortho_psi_l_coef_bi_ortho(buffer)
  do k = 1, N_states
    do i = 1, N_det
      psi_r_coef_bi_ortho(i,k) = reigvec_tc_bi_orth(i,k)
      buffer(i,k) = reigvec_tc_bi_orth(i,k)
    enddo
  enddo
  TOUCH psi_r_coef_bi_ortho
  call ezfio_set_tc_bi_ortho_psi_r_coef_bi_ortho(buffer)
  deallocate(buffer)
!  print*,'After diag'
!  do i = 1, N_det! old version
!   print*,'i',i,psi_l_coef_bi_ortho(i,1),psi_r_coef_bi_ortho(i,1)
!   call debug_det(psi_det(1,1,i),N_int) 
!  enddo 

END_PROVIDER 



subroutine bi_normalize(u_l, u_r, n, ld, nstates)

  BEGIN_DOC
  !!!! Normalization of the scalar product of the left/right eigenvectors
  END_DOC

  implicit none
  integer,          intent(in)    :: n, ld, nstates
  double precision, intent(inout) :: u_l(ld,nstates), u_r(ld,nstates)
  integer                         :: i, j
  double precision                :: accu, tmp

  do i = 1, nstates

    !!!! Normalization of right eigenvectors |Phi>
    accu = 0.d0
    do j = 1, n
      accu += u_r(j,i) * u_r(j,i)
    enddo
    accu = 1.d0/dsqrt(accu)
    print*,'accu_r = ',accu
    do j = 1, n
      u_r(j,i) *= accu
    enddo
    tmp = u_r(1,i) / dabs(u_r(1,i))
    do j = 1, n
      u_r(j,i) *= tmp
    enddo

    !!!! Adaptation of the norm of the left eigenvector such that <chi|Phi> = 1
    accu = 0.d0
    do j = 1, n
      accu += u_l(j,i) * u_r(j,i)
      !print*,j, u_l(j,i) , u_r(j,i)
    enddo
    print*,'accu_lr = ', accu
    if(accu.gt.0.d0)then
      accu = 1.d0/dsqrt(accu)
    else
      accu = 1.d0/dsqrt(-accu)
    endif
    tmp = (u_l(1,i) * u_r(1,i) )/dabs(u_l(1,i) * u_r(1,i))
    do j = 1, n
      u_l(j,i) *= accu * tmp
      u_r(j,i) *= accu
    enddo

  enddo

end

