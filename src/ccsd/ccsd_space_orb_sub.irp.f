subroutine run_ccsd_space_orb
  use gpu

  implicit none

  integer :: i,j,k,l,a,b,c,d,tmp_a,tmp_b,tmp_c,tmp_d
  integer :: u,v,gam,beta,tmp_gam,tmp_beta
  integer :: nb_iter
  double precision :: get_two_e_integral
  double precision :: uncorr_energy,energy, max_elem, max_r, max_r1, max_r2,ta,tb
  logical :: not_converged

  type(gpu_double4) :: t2, r2, tau, tau_x
  type(gpu_double2) :: t1, r1
  type(gpu_double2) :: H_oo, H_vv, H_vo

  type(gpu_double2) :: d_cc_space_f_oo, d_cc_space_f_vo
  type(gpu_double2) :: d_cc_space_f_ov, d_cc_space_f_vv

  type(gpu_double3) :: d_cc_space_v_oo_chol, d_cc_space_v_vo_chol
  type(gpu_double3) :: d_cc_space_v_ov_chol, d_cc_space_v_vv_chol

  type(gpu_double4) :: d_cc_space_v_oovv, d_cc_space_v_voov, d_cc_space_v_ovov
  type(gpu_double4) :: d_cc_space_v_oovo, d_cc_space_v_vooo, d_cc_space_v_oooo
  type(gpu_double4) :: d_cc_space_v_vvoo, d_cc_space_v_ovvo, d_cc_space_v_ovoo

  double precision, allocatable :: all_err(:,:), all_t(:,:)
  integer, allocatable          :: list_occ(:), list_vir(:)
  integer(bit_kind)             :: det(N_int,2)
  integer                       :: nO, nV, nOa, nVa

  call set_multiple_levels_omp(.False.)

  if (do_mo_cholesky) then
    PROVIDE cholesky_mo_transp
    FREE cholesky_ao
  else
    PROVIDE mo_two_e_integrals_in_map
  endif

  det = psi_det(:,:,cc_ref)
  print*,'Reference determinant:'
  call print_det(det,N_int)

  nOa = cc_nOa
  nVa = cc_nVa

  ! Check that the reference is a closed shell determinant
  if (cc_ref_is_open_shell) then
    call abort
  endif

  ! Number of occ/vir spatial orb
  nO = nOa
  nV = nVa

  allocate(list_occ(nO),list_vir(nV))
  list_occ = cc_list_occ
  list_vir = cc_list_vir
  ! Debug
  !call extract_list_orb_space(det,nO,nV,list_occ,list_vir)
  !print*,'occ',list_occ
  !print*,'vir',list_vir

  ! GPU arrays
  call gpu_allocate(d_cc_space_f_oo, nO, nO)
  call gpu_allocate(d_cc_space_f_vo, nV, nO)
  call gpu_allocate(d_cc_space_f_ov, nO, nV)
  call gpu_allocate(d_cc_space_f_vv, nV, nV)

  call gpu_upload(cc_space_f_oo, d_cc_space_f_oo)
  call gpu_upload(cc_space_f_vo, d_cc_space_f_vo)
  call gpu_upload(cc_space_f_ov, d_cc_space_f_ov)
  call gpu_upload(cc_space_f_vv, d_cc_space_f_vv)

!  FREE cc_space_f_oo
!  FREE cc_space_f_vo
!  FREE cc_space_f_vv

  if (do_mo_cholesky) then
    call gpu_allocate(d_cc_space_v_oo_chol, cholesky_mo_num, nO, nO)
    call gpu_allocate(d_cc_space_v_ov_chol, cholesky_mo_num, nO, nV)
    call gpu_allocate(d_cc_space_v_vo_chol, cholesky_mo_num, nV, nO)
    call gpu_allocate(d_cc_space_v_vv_chol, cholesky_mo_num, nV, nV)

    call gpu_upload(cc_space_v_oo_chol, d_cc_space_v_oo_chol)
    call gpu_upload(cc_space_v_ov_chol, d_cc_space_v_ov_chol)
    call gpu_upload(cc_space_v_vo_chol, d_cc_space_v_vo_chol)
    call gpu_upload(cc_space_v_vv_chol, d_cc_space_v_vv_chol)

!    FREE cc_space_v_oo_chol
!    FREE cc_space_v_ov_chol
!    FREE cc_space_v_vo_chol
!    FREE cc_space_v_vv_chol
  endif

  call gpu_allocate(d_cc_space_v_oovv, nO, nO, nV, nV)
  call gpu_allocate(d_cc_space_v_voov, nV, nO, nO, nV)
  call gpu_allocate(d_cc_space_v_ovov, nO, nV, nO, nV)
  call gpu_allocate(d_cc_space_v_oovo, nO, nO, nV, nO)
  call gpu_allocate(d_cc_space_v_ovvo, nO, nV, nV, nO)
  call gpu_allocate(d_cc_space_v_vooo, nV, nO, nO, nO)
  call gpu_allocate(d_cc_space_v_oooo, nO, nO, nO, nO)
  call gpu_allocate(d_cc_space_v_vvoo, nV, nV, nO, nO)
  call gpu_allocate(d_cc_space_v_ovoo, nO, nV, nO, nO)

  call gpu_upload(cc_space_v_oovv, d_cc_space_v_oovv)
  call gpu_upload(cc_space_v_voov, d_cc_space_v_voov)
  call gpu_upload(cc_space_v_ovov, d_cc_space_v_ovov)
  call gpu_upload(cc_space_v_oovo, d_cc_space_v_oovo)
  call gpu_upload(cc_space_v_ovvo, d_cc_space_v_ovvo)
  call gpu_upload(cc_space_v_vooo, d_cc_space_v_vooo)
  call gpu_upload(cc_space_v_oooo, d_cc_space_v_oooo)
  call gpu_upload(cc_space_v_vvoo, d_cc_space_v_vvoo)
  call gpu_upload(cc_space_v_ovoo, d_cc_space_v_ovoo)

!  FREE cc_space_v_voov
!  FREE cc_space_v_ovov
!  FREE cc_space_v_oovo
!  FREE cc_space_v_oovv
!  FREE cc_space_v_vooo
!  FREE cc_space_v_oooo
!  FREE cc_space_v_vvoo
!  FREE cc_space_v_ovvo
!  FREE cc_space_v_ovoo

  call gpu_allocate(t2, nO,nO,nV,nV)
  call gpu_allocate(r2, nO,nO,nV,nV)
  call gpu_allocate(tau, nO,nO,nV,nV)
  call gpu_allocate(tau_x, nO,nO,nV,nV)
  call gpu_allocate(t1, nO,nV)
  call gpu_allocate(r1, nO,nV)
  call gpu_allocate(H_oo, nO, nO)
  call gpu_allocate(H_vo, nV, nO)
  call gpu_allocate(H_vv, nV, nV)

  if (cc_update_method == 'diis') then
    double precision :: rss, diis_mem, extra_mem
    double precision, external :: memory_of_double
    diis_mem = 2.d0*memory_of_double(nO*nV)*(1.d0+nO*nV)
    call resident_memory(rss)
    do while (cc_diis_depth > 1)
      if (rss + diis_mem * cc_diis_depth > qp_max_mem) then
        cc_diis_depth = cc_diis_depth - 1
      else
        exit
      endif
    end do
    if (cc_diis_depth <= 1) then
      print *,  'Not enough memory for DIIS'
      stop -1
    endif
    print *,  'DIIS size  ', cc_diis_depth

    allocate(all_err(nO*nV+nO*nO*nV*(nV*1_8),cc_diis_depth), all_t(nO*nV+nO*nO*nV*(nV*1_8),cc_diis_depth))
    !$OMP PARALLEL PRIVATE(i,j) DEFAULT(SHARED)
    !$OMP DO COLLAPSE(2)
    do j=1,cc_diis_depth
      do i=1, size(all_err,1)
        all_err(i,j) = 0d0
        all_t(i,j)   = 0d0
      enddo
    enddo
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL
  endif

  if (elec_alpha_num /= elec_beta_num) then
    print*, 'Only for closed shell systems'
    print*, 'elec_alpha_num=',elec_alpha_num
    print*, 'elec_beta_num =',elec_beta_num
    print*, 'abort'
    call abort
  endif

  ! Init
  double precision, allocatable :: h_t1(:,:), h_t2(:,:,:,:)
  allocate(h_t1(nO,nV), h_t2(nO,nO,nV,nV))

  call guess_t1(nO,nV,cc_space_f_o,cc_space_f_v,cc_space_f_ov,h_t1)
  call gpu_upload(h_t1, t1)

  call guess_t2(nO,nV,cc_space_f_o,cc_space_f_v,cc_space_v_oovv,h_t2)
  call gpu_upload(h_t2, t2)


  call update_tau_space(nO,nV,h_t1,t1,t2,tau)
  call update_tau_x_space(nO,nV,tau,tau_x)
  call det_energy(det,uncorr_energy)
  print*,'Det energy', uncorr_energy

  call ccsd_energy_space_x(nO,nV,d_cc_space_v_oovv,d_cc_space_f_vo,tau_x,t1,energy)
  print*,'Guess energy', uncorr_energy+energy, energy

  nb_iter = 0
  not_converged = .True.
  max_r1 = 0d0
  max_r2 = 0d0

  write(*,'(A77)') ' -----------------------------------------------------------------------------'
  write(*,'(A77)') ' |   It.  |       E(CCSD) (Ha) | Correlation (Ha) |  Conv. T1  |  Conv. T2  |'
  write(*,'(A77)') ' -----------------------------------------------------------------------------'
  call wall_time(ta)

  do while (not_converged)

    ! Residue
    if (do_mo_cholesky) then
      call compute_H_oo_chol(nO,nV,tau_x,d_cc_space_f_oo, d_cc_space_v_ov_chol,d_cc_space_v_vo_chol,H_oo)
      call compute_H_vv_chol(nO,nV,tau_x,d_cc_space_f_vv, d_cc_space_v_ov_chol,H_vv)
      call compute_H_vo_chol(nO,nV,t1,d_cc_space_f_vo, d_cc_space_v_ov_chol,d_cc_space_v_vo_chol, H_vo)

      call compute_r1_space_chol(nO,nV,t1,t2,tau,H_oo,H_vv,H_vo,r1,max_r1,d_cc_space_f_ov,d_cc_space_f_vo, &
           d_cc_space_v_voov, d_cc_space_v_ovov, d_cc_space_v_oovo, d_cc_space_v_vo_chol, d_cc_space_v_vv_chol)
      call compute_r2_space_chol(nO,nV,t1,t2,tau,H_oo,H_vv, &
           d_cc_space_v_oovv, d_cc_space_v_vooo, d_cc_space_v_oooo, d_cc_space_v_oovo, d_cc_space_v_ovvo, d_cc_space_v_ovoo, &
           d_cc_space_v_ovov, d_cc_space_v_vvoo, d_cc_space_v_oo_chol, d_cc_space_v_ov_chol, d_cc_space_v_vo_chol, d_cc_space_v_vv_chol, &
           d_cc_space_f_vo, &
           r2, max_r2)
    else
      call compute_H_oo(nO,nV,t1%f,t2%f,tau%f,H_oo%f)
      call compute_H_vv(nO,nV,t1%f,t2%f,tau%f,H_vv%f)
      call compute_H_vo(nO,nV,t1%f,t2%f,H_vo%f)

      call compute_r1_space(nO,nV,t1%f,t2%f,tau%f,H_oo%f,H_vv%f,H_vo%f,r1%f,max_r1)
      call compute_r2_space(nO,nV,t1%f,t2%f,tau%f,H_oo%f,H_vv%f,H_vo%f,r2%f,max_r2)
    endif
    max_r = max(max_r1,max_r2)

    ! Update
    if (cc_update_method == 'diis') then
      call update_t_ccsd_diis_v3(nO,nV,nb_iter,cc_space_f_o,cc_space_f_v,r1%f,r2%f,t1%f,t2%f,all_err,all_t)

    ! Standard update as T = T - Delta
    elseif (cc_update_method == 'none') then
      call update_t1(nO,nV,cc_space_f_o,cc_space_f_v,r1%f,t1%f)
      call update_t2(nO,nV,cc_space_f_o,cc_space_f_v,r2%f,t2%f)
    else
      print*,'Unkown cc_method_method: '//cc_update_method
      call abort
    endif

    call update_tau_space(nO,nV,t1%f,t1,t2,tau)
    call update_tau_x_space(nO,nV,tau,tau_x)

    ! Energy
    call ccsd_energy_space_x(nO,nV,d_cc_space_v_oovv,d_cc_space_f_vo,tau_x,t1,energy)
    write(*,'(A3,I6,A3,F18.12,A3,F16.12,A3,ES10.2,A3,ES10.2,A2)') ' | ',nb_iter,' | ', uncorr_energy+energy,' | ', energy,' | ', max_r1,' | ', max_r2,' |'

    nb_iter = nb_iter + 1
    if (max_r < cc_thresh_conv .or. nb_iter > cc_max_iter) then
      not_converged = .False.
    endif

  enddo
  write(*,'(A77)') ' -----------------------------------------------------------------------------'
  call wall_time(tb)
  print*,'Time: ',tb-ta, ' s'
  print*,''
  if (max_r < cc_thresh_conv) then
    write(*,'(A30,I6,A11)') ' Successful convergence after ', nb_iter, ' iterations'
  else
    write(*,'(A26,I6,A11)') ' Failed convergence after ', nb_iter, ' iterations'
  endif
  print*,''
  write(*,'(A15,F18.12,A3)') ' E(CCSD)     = ', uncorr_energy+energy, ' Ha'
  write(*,'(A15,F18.12,A3)') ' Correlation = ', energy, ' Ha'
  write(*,'(A15,ES10.2,A3)')' Conv        = ', max_r
  print*,''

  if (write_amplitudes) then
    call write_t1(nO,nV,t1%f)
    call write_t2(nO,nV,t2%f)
    call ezfio_set_utils_cc_io_amplitudes('Read')
  endif

  ! Deallocation
  if (cc_update_method == 'diis') then
    deallocate(all_err,all_t)
  endif

  call gpu_deallocate(H_oo)
  call gpu_deallocate(H_vv)
  call gpu_deallocate(H_vo)

  call gpu_deallocate(r1)
  call gpu_deallocate(r2)
  call gpu_deallocate(tau)
  call gpu_deallocate(tau_x)

  ! CCSD(T)
  double precision :: e_t, e_t_err
  e_t = 0.d0

  if (cc_par_t .and. elec_alpha_num + elec_beta_num > 2) then

    ! New
    e_t = uncorr_energy + energy ! For print in (T) call
    e_t_err = 0.d0

    print*,'Computing (T) correction...'
    call wall_time(ta)

    call ccsd_par_t_space_stoch(nO,nV,t1%f,t2%f,cc_space_f_o,cc_space_f_v &
         ,cc_space_v_vvvo,cc_space_v_vvoo,cc_space_v_vooo,e_t, e_t_err)

    call wall_time(tb)
    print*,'Time: ',tb-ta, ' s'

    print*,''
    write(*,'(A15,F18.12,A7,F18.12)') ' E(CCSD(T))  = ', uncorr_energy + energy + e_t, ' Ha +/- ', e_t_err
    write(*,'(A15,F18.12,A7,F18.12)') ' E(T)        = ', e_t, ' Ha +/- ', e_t_err
    write(*,'(A15,F18.12,A7,F18.12)') ' Correlation = ', energy + e_t, ' Ha +/- ', e_t_err
    print*,''
  endif

  call save_energy(uncorr_energy + energy, e_t)

  deallocate(h_t1, h_t2)
  if (do_mo_cholesky) then
    call gpu_deallocate(d_cc_space_v_oo_chol)
    call gpu_deallocate(d_cc_space_v_ov_chol)
    call gpu_deallocate(d_cc_space_v_vo_chol)
    call gpu_deallocate(d_cc_space_v_vv_chol)
  endif

  call gpu_deallocate(d_cc_space_v_oovv)
  call gpu_deallocate(d_cc_space_v_voov)
  call gpu_deallocate(d_cc_space_v_ovov)
  call gpu_deallocate(d_cc_space_v_oovo)
  call gpu_deallocate(d_cc_space_v_ovvo)
  call gpu_deallocate(d_cc_space_v_vooo)
  call gpu_deallocate(d_cc_space_v_oooo)
  call gpu_deallocate(d_cc_space_v_vvoo)
  call gpu_deallocate(d_cc_space_v_ovoo)

  call gpu_deallocate(d_cc_space_f_oo)
  call gpu_deallocate(d_cc_space_f_vo)
  call gpu_deallocate(d_cc_space_f_ov)
  call gpu_deallocate(d_cc_space_f_vv)

  call gpu_deallocate(t1)
  call gpu_deallocate(t2)

end

! Energy

subroutine ccsd_energy_space_x(nO,nV,d_cc_space_v_oovv,d_cc_space_f_vo,tau_x,t1,energy)
  use gpu
  implicit none

  integer, intent(in)            :: nO, nV
  type(gpu_double4), intent(in)  :: tau_x, d_cc_space_v_oovv
  type(gpu_double2), intent(in)  :: t1, d_cc_space_f_vo
  double precision, intent(out)  :: energy

  ! internal
  integer :: i,j,a,b
  double precision :: e

  type(gpu_stream) :: s1, s2
  call gpu_stream_create(s1)
  call gpu_stream_create(s2)

  call gpu_set_stream(blas_handle,s1)
  call gpu_ddot(blas_handle, nO*nV, d_cc_space_f_vo%f(1,1), 1, t1%f(1,1), 1, e)

  call gpu_set_stream(blas_handle,s2)
  call gpu_ddot_64(blas_handle, nO*nO*nV*nV*1_8, tau_x%f(1,1,1,1), 1_8, d_cc_space_v_oovv%f(1,1,1,1), 1_8, energy)
  call gpu_set_stream(blas_handle,gpu_default_stream)

  call gpu_synchronize()
  call gpu_stream_destroy(s1)
  call gpu_stream_destroy(s2)

   energy = energy + 2.d0*e

end

! Tau

subroutine update_tau_space(nO,nV,h_t1,t1,t2,tau)
  use gpu
  implicit none

  ! in
  integer, intent(in)           :: nO, nV
  double precision, intent(in)  :: h_t1(nO,nV)
  type(gpu_double2), intent(in) :: t1
  type(gpu_double4), intent(in) :: t2

  ! out
  type(gpu_double4) :: tau

  ! internal
  integer                       :: i,j,a,b

  type(gpu_stream) :: stream(nV)

  !$OMP PARALLEL &
  !$OMP SHARED(nO,nV,tau,t2,t1,h_t1,stream,blas_handle) &
  !$OMP PRIVATE(i,j,a,b) &
  !$OMP DEFAULT(NONE)
  !$OMP DO
  do b=1,nV
    call gpu_stream_create(stream(b))
    call gpu_set_stream(blas_handle,stream(b))
    do j=1,nO
      call gpu_dgeam(blas_handle, 'N', 'N', nO, nV, &
         1.d0, t2%f(1,j,1,b), nO*nO, &
         h_t1(j,b), t1%f(1,1), nO, &
         tau%f(1,j,1,b), nO*nO)
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  call gpu_synchronize()

  do b=1,nV
    call gpu_stream_destroy(stream(b))
  enddo
  call gpu_set_stream(blas_handle,gpu_default_stream)


end

subroutine update_tau_x_space(nO,nV,tau,tau_x)
  use gpu
  implicit none

  ! in
  integer, intent(in)         :: nO, nV
  type(gpu_double4), intent(in)  :: tau

  ! out
  type(gpu_double4) :: tau_x

  ! internal
  integer                       :: i,j,a,b

  type(gpu_stream) :: stream(nV)

  do a=1,nV
    call gpu_stream_create(stream(a))
  enddo

  !$OMP PARALLEL &
  !$OMP SHARED(nO,nV,tau,tau_x,stream,blas_handle) &
  !$OMP PRIVATE(a,b) &
  !$OMP DEFAULT(NONE)
  !$OMP DO
  do b=1,nV
    do a=1,nV
      call gpu_set_stream(blas_handle,stream(a))
      call gpu_dgeam(blas_handle, 'N', 'N', nO, nO, &
          2.d0, tau%f(1,1,a,b), nO, &
         -1.d0, tau%f(1,1,b,a), nO, &
         tau_x%f(1,1,a,b), nO)
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  call gpu_set_stream(blas_handle,gpu_default_stream)
  call gpu_synchronize()

  do b=1,nV
    call gpu_stream_destroy(stream(b))
  enddo


end

! R1

subroutine compute_r1_space(nO,nV,t1,t2,tau,H_oo,H_vv,H_vo,r1,max_r1)

  implicit none

  ! in
  integer, intent(in)           :: nO, nV
  double precision, intent(in)  :: t1(nO,nV), t2(nO,nO,nV,nV), tau(nO,nO,nV,nV)
  double precision, intent(in)  :: H_oo(nO,nO), H_vv(nV,nV), H_vo(nV,nO)

  ! out
  double precision, intent(out) :: r1(nO,nV), max_r1

  ! internal
  integer                       :: u,i,j,beta,a,b

  !$omp parallel &
  !$omp shared(nO,nV,r1,cc_space_f_ov) &
  !$omp private(u,beta) &
  !$omp default(none)
  !$omp do
  do beta = 1, nV
    do u = 1, nO
      r1(u,beta) = cc_space_f_ov(u,beta)
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  ! r1(u,beta) = r1(u,beta) - 2d0 * cc_space_f_vo(a,i) * t1(i,beta) * t1(u,a)
  ! cc_space_f_vo(a,i) * t1(i,beta) -> X1(nV,nV), O(nV*nV*nO)
  ! X1(a,beta) * t1(u,a) -> O(nO*nV*nV)
  ! cc_space_f_vo(a,i) * t1(u,a)    -> X1(nO,nO), O(nO*nO*nV)
  ! X1(i,u) * t1(i,beta) -> O(nO*nO*nV)
  !do beta = 1, nV
  !  do u = 1, nO
  !    do i = 1, nO
  !      do a = 1, nV
  !        r1(u,beta) = r1(u,beta) - 2d0 * cc_space_f_vo(a,i) * t1(i,beta) * t1(u,a)
  !      enddo
  !    enddo
  !  enddo
  !enddo
  double precision, allocatable :: X_oo(:,:)
  allocate(X_oo(nO,nO))
  call dgemm('N','N', nO, nO, nV, &
             -2d0, t1    , size(t1,1), &
                   cc_space_f_vo, size(cc_space_f_vo,1), &
              0d0, X_oo  , size(X_oo,1))

  call dgemm('T','N', nO, nV, nO, &
             1d0, X_oo, size(X_oo,2), &
                  t1  , size(t1,1), &
             1d0, r1  , size(r1,1))
  deallocate(X_oo)

  ! r1(u,beta) = r1(u,beta) + H_vv(a,beta) * t1(u,a)
  !do beta = 1, nV
  !  do u = 1, nO
  !    do a = 1, nV
  !      r1(u,beta) = r1(u,beta) + H_vv(a,beta) * t1(u,a)
  !    enddo
  !  enddo
  !enddo
  call dgemm('N','N', nO, nV, nV, &
             1d0, t1  , size(t1,1), &
                  H_vv, size(H_vv,1), &
             1d0, r1  , size(r1,1))

  ! r1(u,beta) = r1(u,beta) - H_oo(u,i) * t1(i,beta)
  !do beta = 1, nV
  !  do u = 1, nO
  !    do i = 1, nO
  !      r1(u,beta) = r1(u,beta) - H_oo(u,i) * t1(i,beta)
  !    enddo
  !  enddo
  !enddo
  call dgemm('N','N', nO, nV, nO, &
             -1d0, H_oo, size(H_oo,1), &
                   t1  , size(t1,1), &
              1d0, r1, size(r1,1))

  !r1(u,beta) = r1(u,beta) + H_vo(a,i) * (2d0 * t2(i,u,a,beta) - t2(u,i,a,beta) + t1(u,a) * t1(i,beta))
  ! <=>
  ! r1(u,beta) = r1(u,beta) + H_vo(a,i) * X(a,i,u,beta)
  !do beta = 1, nV
  !  do u = 1, nO
  !    do i = 1, nO
  !      do a = 1, nV
  !        r1(u,beta) = r1(u,beta) + H_vo(a,i) * &
  !        (2d0 * t2(i,u,a,beta) - t2(u,i,a,beta) + t1(u,a) * t1(i,beta))
  !      enddo
  !    enddo
  !  enddo
  !enddo
  double precision, allocatable :: X_voov(:,:,:,:)
  allocate(X_voov(nV, nO, nO, nV))

  !$omp parallel &
  !$omp shared(nO,nV,X_voov,t2,t1) &
  !$omp private(u,beta,i,a) &
  !$omp default(none)
  !$omp do
  do beta = 1, nV
    do u = 1, nO
      do i = 1, nO
        do a = 1, nV
          X_voov(a,i,u,beta) = 2d0 * t2(i,u,a,beta) - t2(u,i,a,beta) + t1(u,a) * t1(i,beta)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  call dgemv('T', nV*nO, nO*nV, &
             1d0, X_voov, size(X_voov,1) * size(X_voov,2), &
                  H_vo  , 1, &
             1d0, r1    , 1)

  deallocate(X_voov)

  ! r1(u,beta) = r1(u,beta) + (2d0 * cc_space_v_voov(a,u,i,beta) - cc_space_v_ovov(u,a,i,beta)) * t1(i,a)
  ! <=>
  ! r1(u,beta) = r1(u,beta) + X(i,a,u,beta)
  !do beta = 1, nV
  !  do u = 1, nO
  !    do i = 1, nO
  !      do a = 1, nV
  !        r1(u,beta) = r1(u,beta) + (2d0 * cc_space_v_voov(a,u,i,beta) - cc_space_v_ovov(u,a,i,beta)) * t1(i,a)
  !      enddo
  !    enddo
  !  enddo
  !enddo
  double precision, allocatable :: X_ovov(:,:,:,:)
  allocate(X_ovov(nO, nV, nO, nV))

  !$omp parallel &
  !$omp shared(nO,nV,cc_space_v_ovov,cc_space_v_voov,X_ovov) &
  !$omp private(u,beta,i,a) &
  !$omp default(none)
  !$omp do
  do beta = 1, nV
    do u = 1, nO
      do a = 1, nv
        do i = 1, nO
          X_ovov(i,a,u,beta) = 2d0 * cc_space_v_voov(a,u,i,beta) - cc_space_v_ovov(u,a,i,beta)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  call dgemv('T', nO*nV, nO*nV, &
             1d0, X_ovov, size(X_ovov,1) * size(X_ovov,2), &
                  t1     , 1, &
             1d0, r1     , 1)

  deallocate(X_ovov)

  ! r1(u,beta) = r1(u,beta) + (2d0 * cc_space_v_vvov(a,b,i,beta) - cc_space_v_vvov(b,a,i,beta)) * tau(i,u,a,b)
  ! r1(u,beta) = r1(u,beta) + W(a,b,i,beta) * T(u,a,b,i)
  !do beta = 1, nV
  !  do u = 1, nO
  !    do i = 1, nO
  !      do a = 1, nV
  !        do b = 1, nV
  !          r1(u,beta) = r1(u,beta) + (2d0 * cc_space_v_vvov(a,b,i,beta) - cc_space_v_vvov(b,a,i,beta)) * tau(i,u,a,b)
  !        enddo
  !      enddo
  !    enddo
  !  enddo
  !enddo

  integer :: iblock, block_size, nVmax
  double precision, allocatable :: W_vvov(:,:,:,:), T_vvoo(:,:,:,:)
  block_size = 8
  allocate(W_vvov(nV,nV,nO,block_size), T_vvoo(nV,nV,nO,nO))

  !$omp parallel &
  !$omp shared(nO,nV,cc_space_v_vvov,W_vvov,T_vvoo,tau) &
  !$omp private(b,beta,i,a) &
  !$omp default(none)
  !$omp do
  do u = 1, nO
    do i = 1, nO
      do b = 1, nV
        do a = 1, nV
          T_vvoo(a,b,i,u) = tau(i,u,a,b)
        enddo
      enddo
    enddo
  enddo
  !$omp end do nowait
  !$omp end parallel

  do iblock = 1, nV, block_size
    nVmax = min(block_size,nV-iblock+1)
    !$omp parallel &
    !$omp shared(nO,nV,cc_space_v_vvov,W_vvov,T_vvoo,tau,nVmax,iblock) &
    !$omp private(b,i,a,beta) &
    !$omp default(none)
    !$omp do collapse(2)
    do beta = iblock, iblock + nVmax - 1
      do i = 1, nO
        do b = 1, nV
          do a = 1, nV
            W_vvov(a,b,i,beta-iblock+1) = 2d0 * cc_space_v_vvov(a,b,i,beta) - cc_space_v_vvov(b,a,i,beta)
          enddo
        enddo
      enddo
    enddo
    !$omp end do nowait
    !$omp end parallel

    call dgemm('T','N',nO,nVmax,nO*nV*nV, &
             1d0, T_vvoo, nV*nV*nO, &
                  W_vvov, nO*nV*nV, &
             1d0, r1(1,iblock), nO)
  enddo

  deallocate(W_vvov,T_vvoo)

  ! r1(u,beta) = r1(u,beta) - (2d0 * cc_space_v_vooo(a,u,i,j) - cc_space_v_vooo(a,u,j,i)) * tau(i,j,a,beta)
  ! r1(u,beta) = r1(u,beta) - W(i,j,a,u) * tau(i,j,a,beta)
  !do beta = 1, nV
  !  do u = 1, nO
  !    do i = 1, nO
  !      do j = 1, nO
  !        do a = 1, nV
  !          r1(u,beta) = r1(u,beta) - (2d0 * cc_space_v_vooo(a,u,i,j) - cc_space_v_vooo(a,u,j,i)) * tau(i,j,a,beta)
  !        enddo
  !      enddo
  !    enddo
  !  enddo
  !enddo
  double precision, allocatable :: W_oovo(:,:,:,:)
  allocate(W_oovo(nO,nO,nV,nO))

  !$omp parallel &
  !$omp shared(nO,nV,cc_space_v_vooo,W_oovo) &
  !$omp private(u,a,i,j) &
  !$omp default(none)
  do u = 1, nO
    !$omp do
    do a = 1, nV
      do j = 1, nO
        do i = 1, nO
          W_oovo(i,j,a,u) = 2d0 * cc_space_v_vooo(a,u,i,j) - cc_space_v_vooo(a,u,j,i)
        enddo
      enddo
    enddo
    !$omp end do nowait
  enddo
  !$omp end parallel

  call dgemm('T','N', nO, nV, nO*nO*nV, &
             -1d0, W_oovo, size(W_oovo,1) * size(W_oovo,2) * size(W_oovo,3), &
                   tau   , size(tau,1) * size(tau,2) * size(tau,3), &
              1d0, r1    , size(r1,1))

  deallocate(W_oovo)

  max_r1 = 0d0
  do a = 1, nV
    do i = 1, nO
      max_r1 = max(dabs(r1(i,a)), max_r1)
    enddo
  enddo

  ! Change the sign for consistency with the code in spin orbitals
  !$omp parallel &
  !$omp shared(nO,nV,r1) &
  !$omp private(a,i) &
  !$omp default(none)
  !$omp do
  do a = 1, nV
    do i = 1, nO
      r1(i,a) = -r1(i,a)
    enddo
  enddo
  !$omp end do
  !$omp end parallel

end

! H_oo

subroutine compute_H_oo(nO,nV,t1,t2,tau,H_oo)

  implicit none

  integer, intent(in)           :: nO,nV
  double precision, intent(in)  :: t1(nO, nV)
  double precision, intent(in)  :: t2(nO, nO, nV, nV)
  double precision, intent(in)  :: tau(nO, nO, nV, nV)
  double precision, intent(out) :: H_oo(nO, nO)

  integer :: a,tmp_a,k,b,l,c,d,tmp_c,tmp_d,i,j,u

  !H_oo = 0d0

  !do i = 1, nO
  !  do u = 1, nO
  !    H_oo(u,i) = cc_space_f_oo(u,i)

  !    do j = 1, nO
  !      do a = 1, nV
  !        do b = 1, nV
  !          !H_oo(u,i) = H_oo(u,i) + (2d0 * cc_space_v_vvoo(a,b,i,j) - cc_space_v_vvoo(b,a,i,j)) * tau(u,j,a,b)
  !          !H_oo(u,i) = H_oo(u,i) + cc_space_w_vvoo(a,b,i,j) * tau(u,j,a,b)
  !          H_oo(u,i) = H_oo(u,i) + cc_space_w_oovv(i,j,a,b) * tau(u,j,a,b)
  !        enddo
  !      enddo
  !    enddo
  !
  !  enddo
  !enddo

  ! H_oo(u,i) = cc_space_f_oo(u,i)
  !$omp parallel &
  !$omp shared(nO,H_oo,cc_space_f_oo) &
  !$omp private(i,u) &
  !$omp default(none)
  !$omp do
  do i = 1, nO
    do u = 1, nO
      H_oo(u,i) = cc_space_f_oo(u,i)
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  ! H_oo(u,i) += cc_space_w_oovv(i,j,a,b) * tau(u,j,a,b)
  ! H_oo(u,i) += tau(u,j,a,b) * cc_space_w_oovv(i,j,a,b)
  call dgemm('N','T', nO, nO, nO*nV*nV,       &
             1d0, tau     , size(tau,1),      &
                  cc_space_w_oovv, size(cc_space_w_oovv,1), &
             1d0, H_oo    , size(H_oo,1))

end

! H_vv

subroutine compute_H_vv(nO,nV,t1,t2,tau,H_vv)

  implicit none

  integer, intent(in)           :: nO,nV
  double precision, intent(in)  :: t1(nO, nV)
  double precision, intent(in)  :: t2(nO, nO, nV, nV)
  double precision, intent(in)  :: tau(nO, nO, nV, nV)
  double precision, intent(out) :: H_vv(nV, nV)

  integer :: a,tmp_a,b,k,l,c,d,tmp_c,tmp_d,i,j,u, beta

  !H_vv = 0d0

  !do beta = 1, nV
  !  do a = 1, nV
  !    H_vv(a,beta) = cc_space_f_vv(a,beta)

  !    do j = 1, nO
  !      do i = 1, nO
  !        do b = 1, nV
  !          !H_vv(a,beta) = H_vv(a,beta) - (2d0 * cc_space_v_vvoo(a,b,i,j) - cc_space_v_vvoo(a,b,j,i)) * tau(i,j,beta,b)
  !          H_vv(a,beta) = H_vv(a,beta) - cc_space_w_vvoo(a,b,i,j) * tau(i,j,beta,b)
  !        enddo
  !      enddo
  !    enddo
  !
  !  enddo
  !enddo

  double precision, allocatable :: tmp_tau(:,:,:,:)

  allocate(tmp_tau(nV,nO,nO,nV))

  ! H_vv(a,beta) = cc_space_f_vv(a,beta)
  !$omp parallel &
  !$omp shared(nV,nO,H_vv,cc_space_f_vv,tmp_tau,tau) &
  !$omp private(a,beta,i,j,b) &
  !$omp default(none)
  !$omp do
  do beta = 1, nV
    do a = 1, nV
      H_vv(a,beta) = cc_space_f_vv(a,beta)
    enddo
  enddo
  !$omp end do nowait

  ! H_vv(a,beta) = H_vv(a,beta) - cc_space_w_vvoo(a,b,i,j) * tau(i,j,beta,b)
  ! H_vv(a,beta) = H_vv(a,beta) - cc_space_w_vvoo(a,b,i,j) * tmp_tau(b,i,j,beta)

  !$omp do
  do beta = 1, nV
    do j = 1, nO
      do i = 1, nO
        do b = 1, nV
          tmp_tau(b,i,j,beta) = tau(i,j,beta,b)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  call dgemm('N','N',nV,nV,nO*nO*nV,           &
             -1d0, cc_space_w_vvoo, size(cc_space_w_vvoo,1), &
                   tmp_tau , size(tmp_tau,1) * size(tmp_tau,2) * size(tmp_tau,3), &
              1d0, H_vv    , size(H_vv,1))

  deallocate(tmp_tau)

end

! H_vo

subroutine compute_H_vo(nO,nV,t1,t2,H_vo)

  implicit none

  integer, intent(in)           :: nO,nV
  double precision, intent(in)  :: t1(nO, nV)
  double precision, intent(in)  :: t2(nO, nO, nV, nV)
  double precision, intent(out) :: H_vo(nV, nO)

  integer :: a,tmp_a,b,k,l,c,d,tmp_c,tmp_d,i,j,u, beta

  !H_vo = 0d0

  !do i = 1, nO
  !  do a = 1, nV
  !    H_vo(a,i) = cc_space_f_vo(a,i)

  !    do j = 1, nO
  !      do b = 1, nV
  !        !H_vo(a,i) = H_vo(a,i) + (2d0 * cc_space_v_vvoo(a,b,i,j) - cc_space_v_vvoo(b,a,i,j)) * t1(j,b)
  !        H_vo(a,i) = H_vo(a,i) + cc_space_w_vvoo(a,b,i,j) * t1(j,b)
  !      enddo
  !    enddo
  !
  !  enddo
  !enddo

  double precision, allocatable :: w(:,:,:,:)

  allocate(w(nV,nO,nO,nV))

  !$omp parallel &
  !$omp shared(nV,nO,H_vo,cc_space_f_vo,w,cc_space_w_vvoo,t1) &
  !$omp private(a,beta,i,j,b) &
  !$omp default(none)
  !$omp do
  do i = 1, nO
    do a = 1, nV
      H_vo(a,i) = cc_space_f_vo(a,i)
    enddo
  enddo
  !$omp end do nowait

  ! H_vo(a,i) = H_vo(a,i) + cc_space_w_vvoo(a,b,i,j) * t1(j,b)
  ! H_vo(a,i) = H_vo(a,i) + w(a,i,j,b) * t1(j,b)

  !$omp do
  do b = 1, nV
    do j = 1, nO
      do i = 1, nO
        do a = 1, nV
          w(a,i,j,b) = cc_space_w_vvoo(a,b,i,j)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  call dgemv('N',nV*nO, nO*nV, &
             1d0, w   , size(w,1) * size(w,2), &
                  t1  , 1, &
             1d0, H_vo, 1)

  deallocate(w)

end

! R2

subroutine compute_r2_space(nO,nV,t1,t2,tau,H_oo,H_vv,H_vo,r2,max_r2)

  implicit none

  ! in
  integer, intent(in)           :: nO, nV
  double precision, intent(in)  :: t1(nO,nV), t2(nO,nO,nV,nV), tau(nO,nO,nV,nV)
  double precision, intent(in)  :: H_oo(nO,nO), H_vv(nV,nV), H_vo(nV,nO)

  ! out
  double precision, intent(out) :: r2(nO,nO,nV,nV), max_r2

  ! internal
  double precision, allocatable :: g_occ(:,:), g_vir(:,:), J1(:,:,:,:), K1(:,:,:,:)
  double precision, allocatable :: A1(:,:,:,:), B1_gam(:,:,:)
  integer                       :: u,v,i,j,beta,gam,a,b

  allocate(g_occ(nO,nO), g_vir(nV,nV))
  allocate(J1(nO,nV,nV,nO), K1(nO,nV,nO,nV))
  allocate(A1(nO,nO,nO,nO))

  call compute_g_occ(nO,nV,t1,t2,H_oo,g_occ)
  call compute_g_vir(nO,nV,t1,t2,H_vv,g_vir)
  call compute_A1(nO,nV,t1,t2,tau,A1)
  call compute_J1(nO,nV,t1,t2,cc_space_v_ovvo,cc_space_v_ovoo, &
       cc_space_v_vvvo,cc_space_v_vvoo,J1)
  call compute_K1(nO,nV,t1,t2,cc_space_v_ovoo,cc_space_v_vvoo, &
       cc_space_v_ovov,cc_space_v_vvov,K1)

  ! Residual
  !r2 = 0d0

  !$omp parallel &
  !$omp shared(nO,nV,r2,cc_space_v_oovv) &
  !$omp private(u,v,gam,beta) &
  !$omp default(none)
  !$omp do
  do gam = 1, nV
    do beta = 1, nV
      do v = 1, nO
        do u = 1, nO
         r2(u,v,beta,gam) = cc_space_v_oovv(u,v,beta,gam)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  !do gam = 1, nV
  !  do beta = 1, nV
  !    do v = 1, nO
  !      do u = 1, nO
  !       do j = 1, nO
  !         do i = 1, nO
  !           r2(u,v,beta,gam) = r2(u,v,beta,gam) &
  !           + A1(u,v,i,j) * tau(i,j,beta,gam)
  !         enddo
  !       enddo
  !      enddo
  !    enddo
  !  enddo
  !enddo

  call dgemm('N','N',nO*nO,nV*nV,nO*nO, &
             1d0, A1, size(A1,1) * size(A1,2), &
                  tau, size(tau,1) * size(tau,2), &
             1d0, r2, size(r2,1) * size(r2,2))

  !do gam = 1, nV
  !  do beta = 1, nV
  !    do v = 1, nO
  !      do u = 1, nO
  !       do a = 1, nV
  !         do b = 1, nv
  !           r2(u,v,beta,gam) = r2(u,v,beta,gam) &
  !           + B1(a,b,beta,gam) * tau(u,v,a,b)
  !         enddo
  !       enddo
  !      enddo
  !    enddo
  !  enddo
  !enddo

!  allocate(B1(nV,nV,nV,nV))
!  call compute_B1(nO,nV,t1,t2,B1)
!  call dgemm('N','N',nO*nO,nV*nV,nV*nV, &
!              1d0, tau, size(tau,1) * size(tau,2), &
!                   B1 , size(B1_gam,1) * size(B1_gam,2), &
!              1d0, r2, size(r2,1) * size(r2,2))
  allocate(B1_gam(nV,nV,nV))
  do gam=1,nV
    call compute_B1_gam(nO,nV,t1,t2,B1_gam,gam)
    call dgemm('N','N',nO*nO,nV,nV*nV, &
                1d0, tau, size(tau,1) * size(tau,2), &
                     B1_gam        , size(B1_gam,1) * size(B1_gam,2), &
                1d0, r2(1,1,1,gam), size(r2,1) * size(r2,2))
  enddo
  deallocate(B1_gam)


  !do gam = 1, nV
  !  do beta = 1, nV
  !    do v = 1, nO
  !      do u = 1, nO
  !       do a = 1, nV
  !         r2(u,v,beta,gam) = r2(u,v,beta,gam) &
  !         + g_vir(a,beta) * t2(u,v,a,gam) &
  !         + g_vir(a,gam)  * t2(v,u,a,beta) ! P
  !       enddo
  !      enddo
  !    enddo
  !  enddo
  !enddo
  double precision, allocatable :: X_oovv(:,:,:,:),Y_oovv(:,:,:,:)
  allocate(X_oovv(nO,nO,nV,nV),Y_oovv(nO,nO,nV,nV))

  !$omp parallel &
  !$omp shared(nO,nV,t2,X_oovv) &
  !$omp private(u,v,gam,a) &
  !$omp default(none)
  !$omp do
  do a = 1, nV
    do gam = 1, nV
      do v = 1, nO
        do u = 1, nO
          X_oovv(u,v,gam,a) = t2(u,v,gam,a)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  call dgemm('N','N',nO*nO*nV,nV,nV, &
             1d0, X_oovv, size(X_oovv,1) * size(X_oovv,2) * size(X_oovv,3), &
                  g_vir, size(g_vir,1), &
             0d0, Y_oovv, size(Y_oovv,1) * size(Y_oovv,2) * size(Y_oovv,3))

  !$omp parallel &
  !$omp shared(nO,nV,r2,Y_oovv) &
  !$omp private(u,v,gam,beta) &
  !$omp default(none)
  !$omp do
  do gam = 1, nV
    do beta = 1, nV
      do v = 1, nO
        do u = 1, nO
          r2(u,v,beta,gam) = r2(u,v,beta,gam) + Y_oovv(u,v,beta,gam) + Y_oovv(v,u,gam,beta)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  !do gam = 1, nV
  !  do beta = 1, nV
  !    do v = 1, nO
  !      do u = 1, nO
  !       do i = 1, nO
  !         r2(u,v,beta,gam) = r2(u,v,beta,gam) &
  !         - g_occ(u,i) * t2(i,v,beta,gam) &
  !         - g_occ(v,i) * t2(i,u,gam,beta) ! P
  !       enddo
  !      enddo
  !    enddo
  !  enddo
  !enddo

  call dgemm('N','N',nO,nO*nV*nV,nO, &
             1d0, g_occ , size(g_occ,1), &
                  t2    , size(t2,1), &
             0d0, X_oovv, size(X_oovv,1))

  !$omp parallel &
  !$omp shared(nO,nV,r2,X_oovv) &
  !$omp private(u,v,gam,beta) &
  !$omp default(none)
  !$omp do
  do gam = 1, nV
    do beta = 1, nV
      do v = 1, nO
        do u = 1, nO
          r2(u,v,beta,gam) = r2(u,v,beta,gam) - X_oovv(u,v,beta,gam) - X_oovv(v,u,gam,beta)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  deallocate(X_oovv)

  !do gam = 1, nV
  !  do beta = 1, nV
  !    do v = 1, nO
  !      do u = 1, nO
  !        do a = 1, nV
  !          r2(u,v,beta,gam) = r2(u,v,beta,gam) &
  !          + cc_space_v_ovvv(u,a,beta,gam) * t1(v,a) &
  !          + cc_space_v_ovvv(v,a,gam,beta) * t1(u,a) ! P
  !        enddo
  !      enddo
  !    enddo
  !  enddo
  !enddo

  double precision, allocatable :: X_vovv(:,:,:,:)
  allocate(X_vovv(nV,nO,nV,nV))

  !$omp parallel &
  !$omp shared(nO,nV,X_vovv,cc_space_v_ovvv) &
  !$omp private(u,a,gam,beta) &
  !$omp default(none)
  !$omp do
  do gam = 1, nV
    do beta = 1, nV
      do u = 1, nO
        do a = 1, nV
          X_vovv(a,u,beta,gam) = cc_space_v_ovvv(u,a,beta,gam)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  call dgemm('N','N',nO,nO*nV*nV,nV, &
             1d0, t1    , size(t1,1), &
                  X_vovv, size(X_vovv,1), &
             0d0, Y_oovv, size(Y_oovv,1))

  !$omp parallel &
  !$omp shared(nO,nV,r2,Y_oovv) &
  !$omp private(u,v,gam,beta) &
  !$omp default(none)
  !$omp do
  do gam = 1, nV
    do beta = 1, nV
      do v = 1, nO
        do u = 1, nO
          r2(u,v,beta,gam) = r2(u,v,beta,gam) + Y_oovv(v,u,beta,gam) + Y_oovv(u,v,gam,beta)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  !do gam = 1, nV
  !  do beta = 1, nV
  !    do v = 1, nO
  !      do u = 1, nO
  !        do a = 1, nV
  !          do i = 1, nO
  !           r2(u,v,beta,gam) = r2(u,v,beta,gam) &
  !           - cc_space_v_ovov(u,a,i,gam)  * t1(i,beta) * t1(v,a) &
  !           - cc_space_v_ovov(v,a,i,beta) * t1(i,gam)  * t1(u,a) ! P
  !          enddo
  !        enddo
  !      enddo
  !    enddo
  !  enddo
  !enddo
  double precision, allocatable :: X_vovo(:,:,:,:), Y_vovv(:,:,:,:)
  allocate(X_vovo(nV,nO,nV,nO), Y_vovv(nV,nO,nV,nV),X_oovv(nO,nO,nV,nV))

  !$omp parallel &
  !$omp shared(nO,nV,X_vovo,cc_space_v_ovov) &
  !$omp private(u,v,gam,i) &
  !$omp default(none)
  do i = 1, nO
    !$omp do
    do gam = 1, nV
      do u = 1, nO
        do a = 1, nV
          X_vovo(a,u,gam,i) = cc_space_v_ovov(u,a,i,gam)
        enddo
      enddo
    enddo
    !$omp end do nowait
  enddo
  !$omp end parallel

  call dgemm('N','N',nV*nO*nV,nV,nO, &
              1d0, X_vovo, size(X_vovo,1) * size(X_vovo,2) * size(X_vovo,3), &
                   t1    , size(t1,1), &
              0d0, Y_vovv, size(Y_vovv,1) * size(Y_vovv,2) * size(Y_vovv,3))

  call dgemm('N','N',nO,nO*nV*nV,nV, &
             1d0, t1, size(t1,1), &
                  Y_vovv, size(Y_vovv,1), &
             0d0, X_oovv, size(X_oovv,1))

  !$omp parallel &
  !$omp shared(nO,nV,r2,X_oovv) &
  !$omp private(u,v,gam,beta) &
  !$omp default(none)
  !$omp do
  do gam = 1, nV
    do beta = 1, nV
      do v = 1, nO
        do u = 1, nO
          r2(u,v,beta,gam) = r2(u,v,beta,gam) - X_oovv(v,u,gam,beta) - X_oovv(u,v,beta,gam)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  deallocate(X_vovo,Y_vovv)

  !do gam = 1, nV
  !  do beta = 1, nV
  !    do v = 1, nO
  !      do u = 1, nO
  !       do i = 1, nO
  !         r2(u,v,beta,gam) = r2(u,v,beta,gam) &
  !         - cc_space_v_oovo(u,v,beta,i) * t1(i,gam) &
  !         - cc_space_v_oovo(v,u,gam,i)  * t1(i,beta) ! P
  !        enddo
  !      enddo
  !    enddo
  !  enddo
  !enddo

  call dgemm('N','N',nO*nO*nV,nV,nO, &
             1d0, cc_space_v_oovo, size(cc_space_v_oovo,1) * size(cc_space_v_oovo,2) * size(cc_space_v_oovo,3), &
                  t1 , size(t1,1), &
             0d0, X_oovv, size(X_oovv,1) * size(X_oovv,2) * size(X_oovv,3))

  !$omp parallel &
  !$omp shared(nO,nV,r2,X_oovv) &
  !$omp private(u,v,gam,beta) &
  !$omp default(none)
  !$omp do
  do gam = 1, nV
    do beta = 1, nV
      do v = 1, nO
        do u = 1, nO
           r2(u,v,beta,gam) = r2(u,v,beta,gam) - X_oovv(u,v,beta,gam) - X_oovv(v,u,gam,beta)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel


  !do gam = 1, nV
  !  do beta = 1, nV
  !    do v = 1, nO
  !      do u = 1, nO
  !       do i = 1, nO
  !         do a = 1, nV
  !           r2(u,v,beta,gam) = r2(u,v,beta,gam) &
  !           - cc_space_v_ovvo(u,a,beta,i) * t1(v,a) * t1(i,gam) &
  !           - cc_space_v_ovvo(v,a,gam,i)  * t1(u,a) * t1(i,beta) ! P
  !         enddo
  !       enddo
  !      enddo
  !    enddo
  !  enddo
  !enddo

  double precision, allocatable :: Y_oovo(:,:,:,:)
  allocate(X_vovo(nV,nO,nV,nO), Y_oovo(nO,nO,nV,nO))

  !$omp parallel &
  !$omp shared(nO,nV,X_vovo,cc_space_v_ovvo) &
  !$omp private(a,v,gam,i) &
  !$omp default(none)
  do i = 1, nO
    !$omp do
    do gam = 1, nV
      do v = 1, nO
        do a = 1, nV
          X_vovo(a,v,gam,i) = cc_space_v_ovvo(v,a,gam,i)
        enddo
      enddo
    enddo
    !$omp end do nowait
  enddo
  !$omp end parallel

  call dgemm('N','N',nO,nO*nV*nO,nV, &
             1d0, t1, size(t1,1), &
                  X_vovo, size(X_vovo,1), &
             0d0, Y_oovo, size(Y_oovo,1))

  call dgemm('N','N',nO*nO*nV, nV, nO, &
             1d0, Y_oovo, size(Y_oovo,1) * size(Y_oovo,2) * size(Y_oovo,3), &
                  t1    , size(t1,1), &
             0d0, X_oovv, size(X_oovv,1) * size(X_oovv,2) * size(X_oovv,3))

  !$omp parallel &
  !$omp shared(nO,nV,r2,X_oovv) &
  !$omp private(u,v,gam,beta) &
  !$omp default(none)
  !$omp do
  do gam = 1, nV
    do beta = 1, nV
      do v = 1, nO
        do u = 1, nO
           r2(u,v,beta,gam) = r2(u,v,beta,gam) - X_oovv(u,v,gam,beta) - X_oovv(v,u,beta,gam)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  deallocate(X_vovo,Y_oovo)

  !do gam = 1, nV
  !  do beta = 1, nV
  !    do v = 1, nO
  !      do u = 1, nO
  !       do a = 1, nV
  !         do i = 1, nO
  !           r2(u,v,beta,gam) = r2(u,v,beta,gam) &
  !           + 0.5d0 * (2d0 * J1(u,a,beta,i) - K1(u,a,i,beta)) * &
  !             (2d0 * t2(i,v,a,gam) - t2(i,v,gam,a)) &
  !           + 0.5d0 * (2d0 * J1(v,a,gam,i)  - K1(v,a,i,gam)) * &
  !             (2d0 * t2(i,u,a,beta) - t2(i,u,beta,a)) ! P
  !         enddo
  !       enddo
  !      enddo
  !    enddo
  !  enddo
  !enddo

  double precision, allocatable :: X_ovvo(:,:,:,:), Y_voov(:,:,:,:), Z_ovov(:,:,:,:)
  allocate(X_ovvo(nO,nV,nV,nO), Y_voov(nV,nO,nO,nV),Z_ovov(nO,nV,nO,nV))
  !$omp parallel &
  !$omp shared(nO,nV,X_ovvo,Y_voov,K1,J1,t2) &
  !$omp private(u,v,gam,beta,i,a) &
  !$omp default(none)
  do i = 1, nO
    !$omp do
    do a = 1, nV
      do beta = 1, nV
        do u = 1, nO
          X_ovvo(u,beta,a,i) = (J1(u,a,beta,i) - 0.5d0 * K1(u,a,i,beta))
        enddo
      enddo
    enddo
    !$omp end do nowait
  enddo

  !$omp do
  do gam = 1, nV
    do v = 1, nO
      do i = 1, nO
        do a = 1, nV
          Y_voov(a,i,v,gam) = 2d0 * t2(i,v,a,gam) - t2(i,v,gam,a)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  call dgemm('N','N', nO*nV,nO*nV,nV*nO, &
             1d0, X_ovvo, size(X_ovvo,1) * size(X_ovvo,2), &
                  Y_voov, size(Y_voov,1) * size(Y_voov,2), &
             0d0, Z_ovov, size(Z_ovov,1) * size(Z_ovov,2))

  !$omp parallel &
  !$omp shared(nO,nV,r2,Z_ovov) &
  !$omp private(u,v,gam,beta) &
  !$omp default(none)
  !$omp do
  do gam = 1, nV
    do beta = 1, nV
      do v = 1, nO
        do u = 1, nO
           r2(u,v,beta,gam) = r2(u,v,beta,gam) + Z_ovov(u,beta,v,gam) + Z_ovov(v,gam,u,beta)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  deallocate(X_ovvo,Y_voov)

  !do gam = 1, nV
  !  do beta = 1, nV
  !    do v = 1, nO
  !      do u = 1, nO
  !       do a = 1, nV
  !         do i = 1, nO
  !           r2(u,v,beta,gam) = r2(u,v,beta,gam) &
  !           - 0.5d0 * K1(u,a,i,beta) * t2(i,v,gam,a) &
  !           - 0.5d0 * K1(v,a,i,gam)  * t2(i,u,beta,a) !P
  !         enddo
  !       enddo
  !      enddo
  !    enddo
  !  enddo
  !enddo
  double precision, allocatable :: X_ovov(:,:,:,:),Y_ovov(:,:,:,:)
  allocate(X_ovov(nO,nV,nO,nV),Y_ovov(nO,nV,nO,nV))
  !$omp parallel &
  !$omp shared(nO,nV,r2,K1,X_ovov,Y_ovov,t2) &
  !$omp private(u,a,i,beta,gam) &
  !$omp default(none)
  !$omp do
  do beta = 1, nV
    do u = 1, nO
      do a = 1, nV
        do i = 1, nO
          X_ovov(i,a,u,beta) = 0.5d0 * K1(u,a,i,beta)
        enddo
      enddo
    enddo
  enddo
  !$omp end do nowait

  !$omp do
  do gam = 1, nV
    do v = 1, nO
      do a = 1, nV
        do i = 1, nO
          Y_ovov(i,a,v,gam) = t2(i,v,gam,a)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  call dgemm('T','N',nO*nV,nO*nV,nO*nV, &
             1d0, X_ovov, size(X_ovov,1) * size(X_ovov,2), &
                  Y_ovov, size(Y_ovov,1) * size(Y_ovov,2), &
             0d0, Z_ovov, size(Y_ovov,1) * size(Y_ovov,2))

  !$omp parallel &
  !$omp shared(nO,nV,r2,Z_ovov) &
  !$omp private(u,v,gam,beta) &
  !$omp default(none)
  !$omp do
  do gam = 1, nV
    do beta = 1, nV
      do v = 1, nO
        do u = 1, nO
           r2(u,v,beta,gam) = r2(u,v,beta,gam) - Z_ovov(u,beta,v,gam) - Z_ovov(v,gam,u,beta)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  !do gam = 1, nV
  !  do beta = 1, nV
  !    do v = 1, nO
  !      do u = 1, nO
  !       do a = 1, nV
  !         do i = 1, nO
  !           r2(u,v,beta,gam) = r2(u,v,beta,gam) &
  !           - K1(u,a,i,gam)  * t2(i,v,beta,a) &
  !           - K1(v,a,i,beta) * t2(i,u,gam,a) ! P
  !         enddo
  !       enddo
  !      enddo
  !    enddo
  !  enddo
  !enddo

  !$omp parallel &
  !$omp shared(nO,nV,K1,X_ovov,Y_ovov,t2) &
  !$omp private(u,v,gam,beta,i,a) &
  !$omp default(none)
  !$omp do
  do a = 1, nV
    do i = 1, nO
      do gam = 1, nV
        do u = 1, nO
          X_ovov(u,gam,i,a) = K1(u,a,i,gam)
        enddo
      enddo
    enddo
  enddo
  !$omp end do nowait

  !$omp do
  do beta = 1, nV
    do v = 1, nO
      do a = 1, nV
        do i = 1, nO
          Y_ovov(i,a,v,beta) = t2(i,v,beta,a)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  call dgemm('N','N',nO*nV,nO*nV,nO*nV, &
             1d0, X_ovov, size(X_ovov,1) * size(X_ovov,2), &
                  Y_ovov, size(Y_ovov,1) * size(Y_ovov,2), &
             0d0, Z_ovov, size(Y_ovov,1) * size(Y_ovov,2))

  !$omp parallel &
  !$omp shared(nO,nV,r2,Z_ovov) &
  !$omp private(u,v,gam,beta) &
  !$omp default(none)
  !$omp do
  do gam = 1, nV
    do beta = 1, nV
      do v = 1, nO
        do u = 1, nO
           r2(u,v,beta,gam) = r2(u,v,beta,gam) - Z_ovov(u,gam,v,beta) - Z_ovov(v,beta,u,gam)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  deallocate(X_ovov,Y_ovov,Z_ovov)

  ! Change the sign for consistency with the code in spin orbitals
  !$omp parallel &
  !$omp shared(nO,nV,r2) &
  !$omp private(i,j,a,b) &
  !$omp default(none)
  !$omp do
  do b = 1, nV
    do a = 1, nV
      do j = 1, nO
        do i = 1, nO
          r2(i,j,a,b) = -r2(i,j,a,b)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  max_r2 = 0d0
  do b = 1, nV
    do a = 1, nV
      do j = 1, nO
        do i = 1, nO
          max_r2 = max(r2(i,j,a,b), max_r2)
        enddo
      enddo
    enddo
  enddo

  deallocate(g_occ,g_vir,J1,K1,A1)

end

! A1

subroutine compute_A1(nO,nV,t1,t2,tau,A1)

  implicit none

  integer, intent(in)           :: nO,nV
  double precision, intent(in)  :: t1(nO, nV)
  double precision, intent(in)  :: t2(nO, nO, nV, nV)
  double precision, intent(in)  :: tau(nO, nO, nV, nV)
  double precision, intent(out) :: A1(nO, nO, nO, nO)

  integer :: a,tmp_a,b,k,l,c,d,tmp_c,tmp_d,i,j,u,v, beta

  !A1 = 0d0

  !do j = 1, nO
  !  do i = 1, nO
  !    do v = 1, nO
  !      do u = 1, nO
  !        A1(u,v,i,j) = cc_space_v_oooo(u,v,i,j)

  !        do a = 1, nV
  !          A1(u,v,i,j) = A1(u,v,i,j) &
  !          + cc_space_v_ovoo(u,a,i,j) * t1(v,a) &
  !          + cc_space_v_vooo(a,v,i,j) * t1(u,a)
  !
  !          do b = 1, nV
  !            A1(u,v,i,j) = A1(u,v,i,j) + cc_space_v_vvoo(a,b,i,j) * tau(u,v,a,b)
  !          enddo
  !        enddo
  !
  !      enddo
  !    enddo
  !  enddo
  !enddo

  double precision, allocatable :: X_vooo(:,:,:,:), Y_oooo(:,:,:,:)
  allocate(X_vooo(nV,nO,nO,nO), Y_oooo(nO,nO,nO,nO))

  ! A1(u,v,i,j) = cc_space_v_oooo(u,v,i,j)
  !$omp parallel &
  !$omp shared(nO,nV,A1,cc_space_v_oooo,cc_space_v_ovoo,X_vooo) &
  !$omp private(u,v,i,j) &
  !$omp default(none)
  !$omp do collapse(2)
  do j = 1, nO
    do i = 1, nO
      do v = 1, nO
        do u = 1, nO
          A1(u,v,i,j) = cc_space_v_oooo(u,v,i,j)
        enddo
      enddo
    enddo
  enddo
  !$omp end do nowait

  ! A1(u,v,i,j) += cc_space_v_ovoo(u,a,i,j) * t1(v,a) &

  !$omp do collapse(2)
  do j = 1, nO
    do i = 1, nO
      do u = 1, nO
        do a = 1, nV
          X_vooo(a,u,i,j) = cc_space_v_ovoo(u,a,i,j)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  call dgemm('N','N', nO, nO*nO*nO, nV, &
             1d0, t1    , size(t1,1), &
                  X_vooo, size(X_vooo,1), &
             0d0, Y_oooo, size(Y_oooo,1))

  !$omp parallel &
  !$omp shared(nO,nV,A1,Y_oooo) &
  !$omp private(u,v,i,j) &
  !$omp default(none)
  !$omp do collapse(2)
  do j = 1, nO
    do i = 1, nO
      do v = 1, nO
        do u = 1, nO
          A1(u,v,i,j) = A1(u,v,i,j) + Y_oooo(v,u,i,j)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  deallocate(X_vooo,Y_oooo)

  ! A1(u,v,i,j) += cc_space_v_vooo(a,v,i,j) * t1(u,a)
  call dgemm('N','N', nO, nO*nO*nO, nV, &
             1d0, t1      , size(t1,1), &
                  cc_space_v_vooo, size(cc_space_v_vooo,1), &
             1d0, A1      , size(A1,1))

  ! A1(u,v,i,j) += cc_space_v_vvoo(a,b,i,j) * tau(u,v,a,b)
  call dgemm('N','N', nO*nO, nO*nO, nV*nV, &
             1d0, tau     , size(tau,1) * size(tau,2), &
                  cc_space_v_vvoo, size(cc_space_v_vvoo,1) * size(cc_space_v_vvoo,2), &
             1d0, A1      , size(A1,1) * size(A1,2))

end

! B1

subroutine compute_B1_gam(nO,nV,t1,t2,B1,gam)

  implicit none

  integer, intent(in)           :: nO,nV,gam
  double precision, intent(in)  :: t1(nO, nV)
  double precision, intent(in)  :: t2(nO, nO, nV, nV)
  double precision, intent(out) :: B1(nV, nV, nV)

  integer :: a,tmp_a,b,k,l,c,d,tmp_c,tmp_d,i,j,u,v, beta

!  do beta = 1, nV
!    do b = 1, nV
!      do a = 1, nV
!        B1(a,b,beta) = cc_space_v_vvvv(a,b,beta,gam)
!
!        do i = 1, nO
!          B1(a,b,beta) = B1(a,b,beta) &
!          - cc_space_v_vvvo(a,b,beta,i) * t1(i,gam) &
!          - cc_space_v_vvov(a,b,i,gam) * t1(i,beta)
!        enddo
!
!      enddo
!    enddo
!  enddo

  double precision, allocatable :: X_vvvo(:,:,:), Y_vvvv(:,:,:)

  allocate(X_vvvo(nV,nV,nO), Y_vvvv(nV,nV,nV))
!  ! B1(a,b,beta,gam) = cc_space_v_vvvv(a,b,beta,gam)

  call gen_v_space(cc_nVa,cc_nVa,cc_nVa,1, &
     cc_list_vir,cc_list_vir,cc_list_vir,cc_list_vir(gam), B1)


  !$omp parallel &
  !$omp shared(nO,nV,B1,cc_space_v_vvvv,cc_space_v_vvov,X_vvvo,gam) &
  !$omp private(a,b,beta) &
  !$omp default(none)

!  !$omp do
!    do beta = 1, nV
!      do b = 1, nV
!        do a = 1, nV
!          B1(a,b,beta) = cc_space_v_vvvv(a,b,beta,gam)
!        enddo
!      enddo
!    enddo
!  !$omp end do nowait

  do i = 1, nO
    !$omp do
      do b = 1, nV
        do a = 1, nV
          X_vvvo(a,b,i) = cc_space_v_vvov(a,b,i,gam)
        enddo
      enddo
    !$omp end do
  enddo
  !$omp end parallel

!  ! B1(a,b,beta) -= cc_space_v_vvvo(a,b,beta,i) * t1(i,gam) &
  call dgemm('N','N', nV*nV*nV, 1, nO, &
             -1d0, cc_space_v_vvvo, size(cc_space_v_vvvo,1) * size(cc_space_v_vvvo,2) * size(cc_space_v_vvvo,3), &
                   t1(1,gam), size(t1,1), &
              1d0, B1      , size(B1,1) * size(B1,2) * size(B1,3))

  ! B1(a,b,beta,gam) -= cc_space_v_vvov(a,b,i,gam) * t1(i,beta)
  call dgemm('N','N', nV*nV, nV, nO, &
             -1d0, X_vvvo, size(X_vvvo,1) * size(X_vvvo,2), &
                   t1    , size(t1,1), &
              0d0, Y_vvvv, size(Y_vvvv,1) * size(Y_vvvv,2))

  !$omp parallel &
  !$omp shared(nV,B1,Y_vvvv,gam) &
  !$omp private(a,b,beta) &
  !$omp default(none)
  !$omp do
  do beta = 1, nV
    do b = 1, nV
      do a = 1, nV
        B1(a,b,beta) = B1(a,b,beta) + Y_vvvv(a,b,beta)
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  deallocate(X_vvvo,Y_vvvv)

end

subroutine compute_B1(nO,nV,t1,t2,B1)

  implicit none

  integer, intent(in)           :: nO,nV
  double precision, intent(in)  :: t1(nO, nV)
  double precision, intent(in)  :: t2(nO, nO, nV, nV)
  double precision, intent(out) :: B1(nV, nV, nV, nV)

  integer :: a,tmp_a,b,k,l,c,d,tmp_c,tmp_d,i,j,u,v, beta, gam

  !B1 = 0d0

  !do gam = 1, nV
  !  do beta = 1, nV
  !    do b = 1, nV
  !      do a = 1, nV
  !        B1(a,b,beta,gam) = cc_space_v_vvvv(a,b,beta,gam)

  !        do i = 1, nO
  !          B1(a,b,beta,gam) = B1(a,b,beta,gam) &
  !          - cc_space_v_vvvo(a,b,beta,i) * t1(i,gam) &
  !          - cc_space_v_vvov(a,b,i,gam) * t1(i,beta)
  !        enddo
  !
  !      enddo
  !    enddo
  !  enddo
  !enddo

  double precision, allocatable :: X_vvvo(:,:,:,:), Y_vvvv(:,:,:,:)
  allocate(X_vvvo(nV,nV,nV,nO), Y_vvvv(nV,nV,nV,nV))

  ! B1(a,b,beta,gam) = cc_space_v_vvvv(a,b,beta,gam)
  !$omp parallel &
  !$omp shared(nO,nV,B1,cc_space_v_vvvv,cc_space_v_vvov,X_vvvo) &
  !$omp private(a,b,beta,gam) &
  !$omp default(none)
  !$omp do
  do gam = 1, nV
    do beta = 1, nV
      do b = 1, nV
        do a = 1, nV
          B1(a,b,beta,gam) = cc_space_v_vvvv(a,b,beta,gam)
        enddo
      enddo
    enddo
  enddo
  !$omp end do nowait
  do i = 1, nO
    !$omp do
    do gam = 1, nV
      do b = 1, nV
        do a = 1, nV
          X_vvvo(a,b,gam,i) = cc_space_v_vvov(a,b,i,gam)
        enddo
      enddo
    enddo
    !$omp end do nowait
  enddo
  !$omp end parallel

  ! B1(a,b,beta,gam) -= cc_space_v_vvvo(a,b,beta,i) * t1(i,gam) &
  call dgemm('N','N', nV*nV*nV, nV, nO, &
             -1d0, cc_space_v_vvvo, size(cc_space_v_vvvo,1) * size(cc_space_v_vvvo,2) * size(cc_space_v_vvvo,3), &
                   t1      , size(t1,1), &
              1d0, B1      , size(B1,1) * size(B1,2) * size(B1,3))


  ! B1(a,b,beta,gam) -= cc_space_v_vvov(a,b,i,gam) * t1(i,beta)
  call dgemm('N','N', nV*nV*nV, nV, nO, &
             -1d0, X_vvvo, size(X_vvvo,1) * size(X_vvvo,2) * size(X_vvvo,3), &
                   t1    , size(t1,1), &
              0d0, Y_vvvv, size(Y_vvvv,1) * size(Y_vvvv,2) * size(Y_vvvv,3))

  !$omp parallel &
  !$omp shared(nV,B1,Y_vvvv) &
  !$omp private(a,b,beta,gam) &
  !$omp default(none)
  !$omp do
  do gam = 1, nV
    do beta = 1, nV
      do b = 1, nV
        do a = 1, nV
          B1(a,b,beta,gam) = B1(a,b,beta,gam) + Y_vvvv(a,b,gam,beta)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  deallocate(X_vvvo,Y_vvvv)

end

! g_occ

subroutine compute_g_occ(nO,nV,t1,t2,H_oo,g_occ)

  implicit none

  integer, intent(in)           :: nO,nV
  double precision, intent(in)  :: t1(nO, nV), H_oo(nO, nO)
  double precision, intent(in)  :: t2(nO, nO, nV, nV)
  double precision, intent(out) :: g_occ(nO, nO)

  integer :: a,tmp_a,b,k,l,c,d,tmp_c,tmp_d,i,j,u,v, beta, gam

  !g_occ = 0d0

  !do i = 1, nO
  !  do u = 1, nO
  !    g_occ(u,i) = H_oo(u,i)
  !
  !    do a = 1, nV
  !      g_occ(u,i) = g_occ(u,i) + cc_space_f_vo(a,i) * t1(u,a)
  !
  !      do j = 1, nO
  !        g_occ(u,i) = g_occ(u,i) + (2d0 * cc_space_v_ovoo(u,a,i,j) - cc_space_v_ovoo(u,a,j,i)) * t1(j,a)
  !      enddo
  !
  !    enddo
  !  enddo
  !enddo

  call dgemm('N','N',nO,nO,nV, &
             1d0, t1, size(t1,1), &
                  cc_space_f_vo, size(cc_space_f_vo,1), &
             0d0, g_occ, size(g_occ,1))

  !$omp parallel &
  !$omp shared(nO,nV,g_occ,H_oo, cc_space_v_ovoo,t1) &
  !$omp private(i,j,a,u) &
  !$omp default(none)
  !$omp do
  do i = 1, nO
    do u = 1, nO
      g_occ(u,i) = g_occ(u,i) + H_oo(u,i)
    enddo
  enddo
  !$omp end do

  !$omp do
  do i = 1, nO
    do j = 1, nO
      do a = 1, nV
        do u = 1, nO
          g_occ(u,i) = g_occ(u,i) + (2d0 * cc_space_v_ovoo(u,a,i,j) - cc_space_v_ovoo(u,a,j,i)) * t1(j,a)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

end

! g_vir

subroutine compute_g_vir(nO,nV,t1,t2,H_vv,g_vir)

  implicit none

  integer, intent(in)           :: nO,nV
  double precision, intent(in)  :: t1(nO, nV), H_vv(nV, nV)
  double precision, intent(in)  :: t2(nO, nO, nV, nV)
  double precision, intent(out) :: g_vir(nV, nV)

  integer :: a,tmp_a,b,k,l,c,d,tmp_c,tmp_d,i,j,u,v, beta, gam

  !g_vir = 0d0

  !do beta = 1, nV
  !  do a = 1, nV
  !    g_vir(a,beta) = H_vv(a,beta)
  !
  !    do i = 1, nO
  !      g_vir(a,beta) = g_vir(a,beta) - cc_space_f_vo(a,i) * t1(i,beta)
  !
  !      do b = 1, nV
  !        g_vir(a,beta) = g_vir(a,beta) + (2d0 * cc_space_v_vvvo(a,b,beta,i) - cc_space_v_vvvo(b,a,beta,i)) * t1(i,b)
  !      enddo
  !
  !    enddo
  !  enddo
  !enddo

  call dgemm('N','N',nV,nV,nO, &
             -1d0, cc_space_f_vo , size(cc_space_f_vo,1), &
                   t1   , size(t1,1), &
              0d0, g_vir, size(g_vir,1))

  !$omp parallel &
  !$omp shared(nO,nV,g_vir,H_vv, cc_space_v_vvvo,t1) &
  !$omp private(i,b,a,beta) &
  !$omp default(none)
  !$omp do
  do beta = 1, nV
    do a = 1, nV
      g_vir(a,beta) = g_vir(a,beta) + H_vv(a,beta)
    enddo
  enddo
  !$omp end do

  !$omp do
  do beta = 1, nV
    do i = 1, nO
      do b = 1, nV
        do a = 1, nV
          g_vir(a,beta) = g_vir(a,beta) + (2d0 * cc_space_v_vvvo(a,b,beta,i) - cc_space_v_vvvo(b,a,beta,i)) * t1(i,b)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

end

! J1

subroutine compute_J1(nO,nV,t1,t2,v_ovvo,v_ovoo,v_vvvo,v_vvoo,J1)

  implicit none

  integer, intent(in)           :: nO,nV
  double precision, intent(in)  :: t1(nO, nV)
  double precision, intent(in)  :: t2(nO, nO, nV, nV)
  double precision, intent(in)  :: v_ovvo(nO,nV,nV,nO), v_ovoo(nO,nV,nO,nO)
  double precision, intent(in)  :: v_vvvo(nV,nV,nV,nO), v_vvoo(nV,nV,nO,nO)
  double precision, intent(out) :: J1(nO, nV, nV, nO)

  integer :: a,tmp_a,b,k,l,c,d,tmp_c,tmp_d,i,j,u,v, beta, gam

  !J1 = 0d0

  !do i = 1, nO
  !  do beta = 1, nV
  !    do a = 1, nV
  !      do u = 1, nO
  !       J1(u,a,beta,i) = cc_space_v_ovvo(u,a,beta,i)

  !        do j = 1, nO
  !          J1(u,a,beta,i) = J1(u,a,beta,i) &
  !          - cc_space_v_ovoo(u,a,j,i) * t1(j,beta)
  !        enddo

  !        do b = 1, nV
  !          J1(u,a,beta,i) = J1(u,a,beta,i) &
  !          + cc_space_v_vvvo(b,a,beta,i) * t1(u,b)
  !        enddo

  !        do j = 1, nO
  !          do b = 1, nV
  !           J1(u,a,beta,i) = J1(u,a,beta,i) &
  !           - cc_space_v_vvoo(a,b,i,j) * (0.5d0 * t2(u,j,b,beta) + t1(u,b) * t1(j,beta)) &
  !           + 0.5d0 * (2d0 * cc_space_v_vvoo(a,b,i,j) - cc_space_v_vvoo(b,a,i,j)) * t2(u,j,beta,b)
  !          enddo
  !        enddo
  !
  !      enddo
  !    enddo
  !  enddo
  !enddo

  double precision, allocatable :: X_ovoo(:,:,:,:), Y_ovov(:,:,:,:)
  allocate(X_ovoo(nO,nV,nO,nO),Y_ovov(nO,nV,nO,nV))

  !$omp parallel &
  !$omp shared(nO,nV,J1,v_ovvo,v_ovoo,X_ovoo) &
  !$omp private(i,j,a,u,beta) &
  !$omp default(none)
  do i = 1, nO
    !$omp do
    do beta = 1, nV
      do a = 1, nV
        do u = 1, nO
          J1(u,a,beta,i) = v_ovvo(u,a,beta,i)
        enddo
      enddo
    enddo
    !$omp end do nowait
  enddo

  !$omp do collapse(2)
  do j = 1, nO
    do i = 1, nO
      do a = 1, nV
        do u = 1, nO
          X_ovoo(u,a,i,j) = v_ovoo(u,a,j,i)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  call dgemm('N','N',nO*nV*nO,nV,nO, &
            -1d0, X_ovoo, size(X_ovoo,1) * size(X_ovoo,2) * size(X_ovoo,3), &
                  t1    , size(t1,1), &
             0d0, Y_ovov, size(Y_ovov,1) * size(Y_ovov,2) * size(Y_ovov,3))

  !$omp parallel &
  !$omp shared(nO,nV,J1,Y_ovov) &
  !$omp private(i,beta,a,u) &
  !$omp default(none)
  do i = 1, nO
    !$omp do
    do beta = 1, nV
      do a = 1, nV
        do u = 1, nO
          J1(u,a,beta,i) = J1(u,a,beta,i) + Y_ovov(u,a,i,beta)
        enddo
      enddo
    enddo
    !$omp end do nowait
  enddo
  !$omp end parallel
  deallocate(X_ovoo)

  ! v_vvvo(b,a,beta,i) * t1(u,b)
  call dgemm('N','N',nO,nV*nV*nO,nV, &
             1d0, t1    , size(t1,1), &
                  v_vvvo, size(v_vvvo,1), &
             1d0, J1    , size(J1,1))

  !- cc_space_v_vvoo(a,b,i,j) * (0.5d0 * t2(u,j,b,beta) + t1(u,b) * t1(j,beta)) &
  double precision, allocatable :: X_voov(:,:,:,:), Z_ovvo(:,:,:,:)
  allocate(X_voov(nV,nO,nO,nV), Z_ovvo(nO,nV,nV,nO))
  !$omp parallel &
  !$omp shared(nO,nV,t2,t1,Y_ovov,X_voov,v_vvoo) &
  !$omp private(i,beta,a,u,b,j) &
  !$omp default(none)
  !$omp do
  do b = 1, nV
    do j = 1, nO
      do beta = 1, nV
        do u = 1, nO
          Y_ovov(u,beta,j,b) = 0.5d0 * t2(u,j,b,beta) + t1(u,b) * t1(j,beta)
        enddo
      enddo
    enddo
  enddo
  !$omp end do nowait

  !$omp do
  do b = 1, nV
    do j = 1, nO
      do i = 1, nO
        do a = 1, nV
          X_voov(a,i,j,b) = v_vvoo(a,b,i,j)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  call dgemm('N','T',nO*nV,nV*nO,nO*nV, &
             -1d0, Y_ovov, size(Y_ovov,1) * size(Y_ovov,2), &
                   X_voov, size(X_voov,1) * size(X_voov,2), &
              0d0, Z_ovvo, size(Z_ovvo,1) * size(Z_ovvo,2))
  deallocate(X_voov)

  double precision, allocatable :: X_ovvo(:,:,:,:), Y_vovo(:,:,:,:)
  allocate(X_ovvo(nO,nV,nV,nO),Y_vovo(nV,nO,nV,nO))
  !$omp parallel &
  !$omp shared(nO,nV,J1,Z_ovvo,t2,Y_vovo,v_vvoo,X_ovvo) &
  !$omp private(i,beta,a,u,j,b) &
  !$omp default(none)
  do i = 1, nO
    !$omp do
    do beta = 1, nV
      do a = 1, nV
        do u = 1, nO
          J1(u,a,beta,i) = J1(u,a,beta,i) + Z_ovvo(u,beta,a,i)
        enddo
      enddo
    enddo
    !$omp end do nowait
  enddo

  !+ 0.5d0 * (2d0 * cc_space_v_vvoo(a,b,i,j) - cc_space_v_vvoo(b,a,i,j)) * t2(u,j,beta,b)
  do j = 1, nO
    !$omp do
    do b = 1, nV
      do i = 1, nO
        do a = 1, nV
          Y_vovo(a,i,b,j) = 0.5d0 * (2d0 * v_vvoo(a,b,i,j) - v_vvoo(b,a,i,j))
        enddo
      enddo
    enddo
    !$omp end do nowait
  enddo

  do j = 1, nO
    !$omp do
    do b = 1, nV
      do beta = 1, nV
        do u = 1, nO
          X_ovvo(u,beta,b,j) = t2(u,j,beta,b)
        enddo
      enddo
    enddo
    !$omp end do nowait
  enddo
  !$omp end parallel

  call dgemm('N','T',nO*nV,nV*nO,nV*nO, &
             1d0, X_ovvo, size(X_ovvo,1) * size(X_ovvo,2), &
                  Y_vovo, size(Y_vovo,1) * size(Y_vovo,2), &
             0d0, Z_ovvo, size(Z_ovvo,1) * size(Z_ovvo,2))

  !$omp parallel &
  !$omp shared(nO,nV,J1,Z_ovvo) &
  !$omp private(i,beta,a,u) &
  !$omp default(none)
  do i = 1, nO
    !$omp do
    do beta = 1, nV
      do a = 1, nV
        do u = 1, nO
          J1(u,a,beta,i) = J1(u,a,beta,i) + Z_ovvo(u,beta,a,i)
        enddo
      enddo
    enddo
    !$omp end do nowait
  enddo
  !$omp end parallel

  deallocate(X_ovvo,Z_ovvo,Y_ovov)

end

! K1

subroutine compute_K1(nO,nV,t1,t2,v_ovoo,v_vvoo,v_ovov,v_vvov,K1)

  implicit none

  integer, intent(in)           :: nO,nV
  double precision, intent(in)  :: t1(nO, nV)
  double precision, intent(in)  :: t2(nO, nO, nV, nV)
  double precision, intent(in)  :: v_vvoo(nV,nV,nO,nO), v_ovov(nO,nV,nO,nV)
  double precision, intent(in)  :: v_vvov(nV,nV,nO,nV), v_ovoo(nO,nV,nO,nO)
  double precision, intent(out) :: K1(nO, nV, nO, nV)

  double precision, allocatable :: X(:,:,:,:), Y(:,:,:,:), Z(:,:,:,:)

  integer :: a,tmp_a,b,k,l,c,d,tmp_c,tmp_d,i,j,u,v, beta, gam

  !K1 = 0d0

  !do beta = 1, nV
  !  do i = 1, nO
  !    do a = 1, nV
  !      do u = 1, nO
  !        K1(u,a,i,beta) = cc_space_v_ovov(u,a,i,beta)

  !        do j = 1, nO
  !          K1(u,a,i,beta) = K1(u,a,i,beta) &
  !          - cc_space_v_ovoo(u,a,i,j) * t1(j,beta)
  !        enddo

  !        do b = 1, nV
  !          K1(u,a,i,beta) = K1(u,a,i,beta) &
  !          + cc_space_v_vvov(b,a,i,beta) * t1(u,b)
  !        enddo

  !        do j = 1, nO
  !          do b = 1, nV
  !           K1(u,a,i,beta) = K1(u,a,i,beta) &
  !           - cc_space_v_vvoo(b,a,i,j) * (0.5d0 * t2(u,j,b,beta) + t1(u,b) * t1(j,beta))
  !          enddo
  !        enddo
  !
  !      enddo
  !    enddo
  !  enddo
  !enddo

  allocate(X(nV,nO,nV,nO),Y(nO,nV,nV,nO),Z(nO,nV,nV,nO))

  !$omp parallel &
  !$omp shared(nO,nV,K1,X,Y,v_vvoo,v_ovov,t1,t2) &
  !$omp private(i,beta,a,u,j,b) &
  !$omp default(none)
  !$omp do
  do beta = 1, nV
    do i = 1, nO
      do a = 1, nV
        do u = 1, nO
          K1(u,a,i,beta) = v_ovov(u,a,i,beta)
        enddo
      enddo
    enddo
  enddo
  !$omp end do nowait

  do i = 1, nO
    !$omp do
    do a = 1, nV
      do j = 1, nO
        do b = 1, nV
          X(b,j,a,i) = - v_vvoo(b,a,i,j)
        enddo
      enddo
    enddo
    !$omp end do nowait
  enddo

  do j = 1, nO
    !$omp do
    do b = 1, nV
      do beta = 1, nV
        do u = 1, nO
          Y(u,beta,b,j) = 0.5d0 * t2(u,j,b,beta) + t1(u,b) * t1(j,beta)
        enddo
      enddo
    enddo
    !$omp end do
  enddo
  !$omp end parallel

  call dgemm('N','N',nO*nV*nO,nV,nO, &
            -1d0, v_ovoo, size(v_ovoo,1) * size(v_ovoo,2) * size(v_ovoo,3), &
                  t1    , size(t1,1), &
            1d0, K1    , size(K1,1) * size(K1,2) * size(K1,3))

  call dgemm('N','N',nO,nV*nO*nV,nV, &
             1d0, t1    , size(t1,1), &
                  v_vvov, size(v_vvov,1), &
             1d0, K1    , size(K1,1))

  ! Y(u,beta,b,j) * X(b,j,a,i) = Z(u,beta,a,i)
  call dgemm('N','N',nV*nO,nO*nV,nV*nO, &
             1d0, Y, size(Y,1) * size(Y,2), &
                  X, size(X,1) * size(X,2), &
             0d0, Z, size(Z,1) * size(Z,2))

  !$omp parallel &
  !$omp shared(nO,nV,K1,Z) &
  !$omp private(i,beta,a,u) &
  !$omp default(none)
  !$omp do
   do beta = 1, nV
    do i = 1, nO
      do a = 1, nV
        do u = 1, nO
          K1(u,a,i,beta) = K1(u,a,i,beta) + Z(u,beta,a,i)
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  deallocate(X,Y,Z)

end
