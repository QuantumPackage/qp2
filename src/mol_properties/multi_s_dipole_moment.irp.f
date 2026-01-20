! Providers for the dipole moments along x,y,z and the total dipole
! moments.

! The dipole moment along the x axis is:
! \begin{align*}
! \mu_x = < \Psi_m | \sum_i x_i + \sum_A Z_A R_A | \Psi_n >
! \end{align*}
! where $i$ is used for the electrons and $A$ for the nuclei.
! $Z_A$ the charge of the nucleus $A$ and $R_A$ its position in the
! space.

! And it can be computed using the (transition, if n /= m) density
! matrix as a expectation value
! \begin{align*}
! <\Psi_n|x| \Psi_m > = \sum_p \gamma_{pp}^{nm} < \phi_p | x | \phi_p >
!      + \sum_{pq, p \neq q} \gamma_{pq}^{nm} < \phi_p |x | \phi_q > +  < \Psi_m | \sum_A Z_A R_A | \Psi_n >
! \end{align*}



 BEGIN_PROVIDER [double precision, multi_s_dipole_moment  , (N_states, N_states)]
&BEGIN_PROVIDER [double precision, multi_s_x_dipole_moment, (N_states, N_states)]
&BEGIN_PROVIDER [double precision, multi_s_y_dipole_moment, (N_states, N_states)]
&BEGIN_PROVIDER [double precision, multi_s_z_dipole_moment, (N_states, N_states)]

  implicit none

  BEGIN_DOC
  ! Providers for :
  ! <\Psi_m|\mu_x|\Psi_n>
  ! <\Psi_m|\mu_y|\Psi_n>
  ! <\Psi_m|\mu_z|\Psi_n>
  ! ||\mu|| = \sqrt{\mu_x^2 + \mu_y^2 + \mu_z^2}
  !
  ! <\Psi_n|x| \Psi_m > = \sum_p \gamma_{pp}^{nm} \bra{\phi_p} x \ket{\phi_p} 
  !   + \sum_{pq, p \neq q} \gamma_{pq}^{nm} \bra{\phi_p} x \ket{\phi_q}
  ! \Psi: wf
  ! n,m indexes for the states
  ! p,q: general spatial MOs 
  ! gamma^{nm}: density matrix \bra{\Psi^n} a^{\dagger}_a a_i \ket{\Psi^m}
  END_DOC
  USE OMP_LIB
  integer          :: istate, jstate ! States
  integer          :: i, j           ! general spatial MOs
  double precision :: nuclei_part_x, nuclei_part_y, nuclei_part_z
  double precision :: mem_tot_tr_dm
  integer :: nthreads
  nthreads = OMP_GET_MAX_THREADS()
  mem_tot_tr_dm = 8.d0*mo_num*mo_num*n_states*n_states*1d-9*(nthreads+1)
 
  multi_s_x_dipole_moment = 0.d0
  multi_s_y_dipole_moment = 0.d0
  multi_s_z_dipole_moment = 0.d0

  if(mem_tot_tr_dm .lt. 0.9d0 * qp_max_mem) then
 
    do jstate = 1, N_states
      do istate = 1, N_states
        do i = 1, mo_num  
          do j = 1, mo_num  
            multi_s_x_dipole_moment(istate,jstate) -= one_e_tr_dm_mo(j,i,istate,jstate) * mo_prop_dipole_x(j,i)  
            multi_s_y_dipole_moment(istate,jstate) -= one_e_tr_dm_mo(j,i,istate,jstate) * mo_prop_dipole_y(j,i) 
            multi_s_z_dipole_moment(istate,jstate) -= one_e_tr_dm_mo(j,i,istate,jstate) * mo_prop_dipole_z(j,i) 
          enddo
        enddo 
      enddo
    enddo

  else

    ! no enouph memory
    ! on the fly scheme
    print*,'Computing on the fly the various dipole matrix elements '

    PROVIDE psi_det_alpha_unique psi_det_beta_unique

    integer           :: l, k_a, k_b
    integer           :: occ(N_int*bit_kind_size,2)
    integer           :: h1, h2, p1, p2, degree
    integer           :: exc(0:2,2), n_occ(2)
    integer           :: krow, kcol, lrow, lcol
    integer(bit_kind) :: tmp_det(N_int,2), tmp_det2(N_int)
    double precision  :: ck, ckl, phase

    !$OMP PARALLEL DEFAULT(NONE)                                                      &
    !$OMP PRIVATE(j, l, k_a, k_b, istate, jstate, occ, ck, ckl, h1, h2, p1, p2, exc,  & 
    !$OMP         phase, degree, n_occ, krow, kcol, lrow, lcol, tmp_det, tmp_det2)    &
    !$OMP SHARED(N_int, N_states, elec_alpha_num, elec_beta_num, N_det,               &
    !$OMP        psi_bilinear_matrix_rows, psi_bilinear_matrix_columns,               &
    !$OMP        psi_bilinear_matrix_transp_rows, psi_bilinear_matrix_transp_columns, &
    !$OMP        psi_det_alpha_unique, psi_det_beta_unique,                           &
    !$OMP        psi_bilinear_matrix_values, psi_bilinear_matrix_transp_values,       &
    !$OMP        mo_prop_dipole_x, mo_prop_dipole_y, mo_prop_dipole_z,                               &
    !$OMP        multi_s_x_dipole_moment, multi_s_y_dipole_moment, multi_s_z_dipole_moment)
    !$OMP DO COLLAPSE(2)
    do istate = 1, N_states
      do jstate = 1, N_states

        do k_a = 1, N_det
          krow = psi_bilinear_matrix_rows   (k_a)
          kcol = psi_bilinear_matrix_columns(k_a)
  
          tmp_det(1:N_int,1) = psi_det_alpha_unique(1:N_int,krow)
          tmp_det(1:N_int,2) = psi_det_beta_unique (1:N_int,kcol)
  
          ! Diagonal part
          call bitstring_to_list_ab(tmp_det, occ, n_occ, N_int)
          ck = psi_bilinear_matrix_values(k_a,istate)*psi_bilinear_matrix_values(k_a,jstate)
          do l = 1, elec_alpha_num
            j = occ(l,1)
            multi_s_x_dipole_moment(istate,jstate) -= ck * mo_prop_dipole_x(j,j) 
            multi_s_y_dipole_moment(istate,jstate) -= ck * mo_prop_dipole_y(j,j) 
            multi_s_z_dipole_moment(istate,jstate) -= ck * mo_prop_dipole_z(j,j) 
          enddo
  
          if (k_a == N_det) cycle
          l = k_a + 1
          lrow = psi_bilinear_matrix_rows   (l)
          lcol = psi_bilinear_matrix_columns(l)
          ! Fix beta determinant, loop over alphas
          do while (lcol == kcol)
            tmp_det2(:) = psi_det_alpha_unique(:,lrow)
            call get_excitation_degree_spin(tmp_det(1,1), tmp_det2, degree, N_int)
            if (degree == 1) then
              exc = 0
              call get_single_excitation_spin(tmp_det(1,1), tmp_det2, exc, phase, N_int)
              call decode_exc_spin(exc, h1, p1, h2, p2)
              ckl = psi_bilinear_matrix_values(k_a,istate)*psi_bilinear_matrix_values(l,jstate) * phase
              multi_s_x_dipole_moment(istate,jstate) -= ckl * mo_prop_dipole_x(h1,p1) 
              multi_s_y_dipole_moment(istate,jstate) -= ckl * mo_prop_dipole_y(h1,p1) 
              multi_s_z_dipole_moment(istate,jstate) -= ckl * mo_prop_dipole_z(h1,p1) 
              ckl = psi_bilinear_matrix_values(k_a,jstate)*psi_bilinear_matrix_values(l,istate) * phase
              multi_s_x_dipole_moment(istate,jstate) -= ckl * mo_prop_dipole_x(p1,h1) 
              multi_s_y_dipole_moment(istate,jstate) -= ckl * mo_prop_dipole_y(p1,h1) 
              multi_s_z_dipole_moment(istate,jstate) -= ckl * mo_prop_dipole_z(p1,h1) 
            endif
            l = l+1
            if (l > N_det) exit
            lrow = psi_bilinear_matrix_rows   (l)
            lcol = psi_bilinear_matrix_columns(l)
          enddo
        enddo ! k_a
  
        do k_b = 1, N_det
          krow = psi_bilinear_matrix_transp_rows   (k_b)
          kcol = psi_bilinear_matrix_transp_columns(k_b)
      
          tmp_det(1:N_int,1) = psi_det_alpha_unique(1:N_int,krow)
          tmp_det(1:N_int,2) = psi_det_beta_unique (1:N_int,kcol)
      
          ! Diagonal part
          call bitstring_to_list_ab(tmp_det, occ, n_occ, N_int)
          ck = psi_bilinear_matrix_transp_values(k_b,istate)*psi_bilinear_matrix_transp_values(k_b,jstate)
          do l = 1, elec_beta_num
            j = occ(l,2)
            multi_s_x_dipole_moment(istate,jstate) -= ck * mo_prop_dipole_x(j,j) 
            multi_s_y_dipole_moment(istate,jstate) -= ck * mo_prop_dipole_y(j,j) 
            multi_s_z_dipole_moment(istate,jstate) -= ck * mo_prop_dipole_z(j,j) 
          enddo
      
          if (k_b == N_det) cycle
          l = k_b+1
          lrow = psi_bilinear_matrix_transp_rows   (l)
          lcol = psi_bilinear_matrix_transp_columns(l)
          ! Fix beta determinant, loop over alphas
          do while (lrow == krow)
            tmp_det2(:) = psi_det_beta_unique(:,lcol)
            call get_excitation_degree_spin(tmp_det(1,2), tmp_det2, degree, N_int)
            if (degree == 1) then
              exc = 0
              call get_single_excitation_spin(tmp_det(1,2), tmp_det2, exc, phase, N_int)
              call decode_exc_spin(exc, h1, p1, h2, p2)
              ckl = psi_bilinear_matrix_transp_values(k_b,istate)*psi_bilinear_matrix_transp_values(l,jstate) * phase
              multi_s_x_dipole_moment(istate,jstate) -= ckl * mo_prop_dipole_x(h1,p1) 
              multi_s_y_dipole_moment(istate,jstate) -= ckl * mo_prop_dipole_y(h1,p1) 
              multi_s_z_dipole_moment(istate,jstate) -= ckl * mo_prop_dipole_z(h1,p1) 
              ckl = psi_bilinear_matrix_transp_values(k_b,jstate)*psi_bilinear_matrix_transp_values(l,istate) * phase
              multi_s_x_dipole_moment(istate,jstate) -= ckl * mo_prop_dipole_x(p1,h1) 
              multi_s_y_dipole_moment(istate,jstate) -= ckl * mo_prop_dipole_y(p1,h1) 
              multi_s_z_dipole_moment(istate,jstate) -= ckl * mo_prop_dipole_z(p1,h1) 
            endif
            l = l+1
            if (l > N_det) exit
            lrow = psi_bilinear_matrix_transp_rows   (l)
            lcol = psi_bilinear_matrix_transp_columns(l)
          enddo
        enddo ! k_b

      enddo ! istate
    enddo ! jstate
    !$OMP END DO
    !$OMP END PARALLEL

  endif ! memory condition
 
  ! Nuclei part
  nuclei_part_x = 0.d0
  nuclei_part_y = 0.d0
  nuclei_part_z = 0.d0
 
  do i = 1,nucl_num 
    nuclei_part_x += nucl_charge(i) * nucl_coord(i,1) 
    nuclei_part_y += nucl_charge(i) * nucl_coord(i,2) 
    nuclei_part_z += nucl_charge(i) * nucl_coord(i,3) 
  enddo
 
  ! Only if istate = jstate, otherwise 0 by the orthogonality of the states
  do istate = 1, N_states
    multi_s_x_dipole_moment(istate,istate) += nuclei_part_x
    multi_s_y_dipole_moment(istate,istate) += nuclei_part_y
    multi_s_z_dipole_moment(istate,istate) += nuclei_part_z
  enddo
 
  ! d = <Psi|r|Psi>
  do jstate = 1, N_states
    do istate = 1, N_states
      multi_s_dipole_moment(istate,jstate) = &
        dsqrt(multi_s_x_dipole_moment(istate,jstate)**2 & 
            + multi_s_y_dipole_moment(istate,jstate)**2 &
            + multi_s_z_dipole_moment(istate,jstate)**2) 
    enddo
  enddo

END_PROVIDER

! ---

 BEGIN_PROVIDER [double precision, multi_s_x_dipole_moment_eigenvec, (N_states, N_states)]
&BEGIN_PROVIDER [double precision, multi_s_y_dipole_moment_eigenvec, (N_states, N_states)]
&BEGIN_PROVIDER [double precision, multi_s_z_dipole_moment_eigenvec, (N_states, N_states)]
&BEGIN_PROVIDER [double precision, multi_s_x_dipole_moment_eigenval,           (N_states)]
&BEGIN_PROVIDER [double precision, multi_s_y_dipole_moment_eigenval,           (N_states)]
&BEGIN_PROVIDER [double precision, multi_s_z_dipole_moment_eigenval,           (N_states)]

  implicit none
  double precision, allocatable :: eigval(:), eigvec(:,:), A(:,:)

  PROVIDE multi_s_x_dipole_moment multi_s_y_dipole_moment multi_s_z_dipole_moment

  allocate(A(N_states,N_states), eigvec(N_states,N_states), eigval(N_states))

  A = multi_s_x_dipole_moment
  call lapack_diag(eigval(1), eigvec(1,1), A(1,1), N_states, N_states)
  multi_s_x_dipole_moment_eigenval = eigval
  multi_s_x_dipole_moment_eigenvec = eigvec

  A = multi_s_y_dipole_moment
  call lapack_diag(eigval(1), eigvec(1,1), A(1,1), N_states, N_states)
  multi_s_y_dipole_moment_eigenval = eigval
  multi_s_y_dipole_moment_eigenvec = eigvec

  A = multi_s_z_dipole_moment
  call lapack_diag(eigval(1), eigvec(1,1), A(1,1), N_states, N_states)
  multi_s_z_dipole_moment_eigenval = eigval
  multi_s_z_dipole_moment_eigenvec = eigvec

  deallocate(A, eigvec, eigval)

END_PROVIDER

! ---


