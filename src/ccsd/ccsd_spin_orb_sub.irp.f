! Code

subroutine run_ccsd_spin_orb

  implicit none

  BEGIN_DOC
  ! CCSD in spin orbitals
  END_DOC

  double precision, allocatable :: t1(:,:), t2(:,:,:,:), tau(:,:,:,:), tau_t(:,:,:,:)
  double precision, allocatable :: r1(:,:), r2(:,:,:,:)
  double precision, allocatable :: cF_oo(:,:), cF_ov(:,:), cF_vv(:,:)
  double precision, allocatable :: cW_oooo(:,:,:,:), cW_ovvo(:,:,:,:), cW_vvvv(:,:,:,:)
  
  double precision, allocatable :: f_oo(:,:), f_ov(:,:), f_vv(:,:), f_o(:), f_v(:)
  double precision, allocatable :: v_oooo(:,:,:,:), v_vooo(:,:,:,:), v_ovoo(:,:,:,:)
  double precision, allocatable :: v_oovo(:,:,:,:), v_ooov(:,:,:,:), v_vvoo(:,:,:,:)
  double precision, allocatable :: v_vovo(:,:,:,:), v_voov(:,:,:,:), v_ovvo(:,:,:,:)
  double precision, allocatable :: v_ovov(:,:,:,:), v_oovv(:,:,:,:), v_vvvo(:,:,:,:)
  double precision, allocatable :: v_vvov(:,:,:,:), v_vovv(:,:,:,:), v_ovvv(:,:,:,:)
  double precision, allocatable :: v_vvvv(:,:,:,:)

  double precision, allocatable :: all_err(:,:), all_t(:,:)

  logical                       :: not_converged
  integer, allocatable          :: list_occ(:,:), list_vir(:,:)
  integer                       :: nO,nV,nOa,nOb,nVa,nVb,nO_m,nV_m,nO_S(2),nV_S(2),n_spin(4)
  integer                       :: nb_iter, i,j,a,b
  double precision              :: uncorr_energy, energy, max_r, max_r1, max_r2, cc, ta, tb,ti,tf,tbi,tfi
  integer(bit_kind)             :: det(N_int,2)

  det = psi_det(:,:,cc_ref)
  print*,'Reference determinant:'
  call print_det(det,N_int)
 
  ! Extract number of occ/vir alpha/beta spin orbitals
  !call extract_n_spin(det,n_spin)
  nOa = cc_nOa !n_spin(1)
  nOb = cc_nOb !n_spin(2)
  nVa = cc_nVa !n_spin(3)
  nVb = cc_nVb !n_spin(4)

  ! Total number of occ/vir spin orb
  nO = cc_nOab !nOa + nOb
  nV = cc_nVab !nVa + nVb
  ! Debug
  !print*,nO,nV

  ! Number of occ/vir spin orb per spin
  nO_S = cc_nO_S !(/nOa,nOb/)
  nV_S = cc_nV_S !(/nVa,nVb/)
  ! Debug
  !print*,nO_S,nV_S

  ! Maximal number of occ/vir 
  nO_m = cc_nO_m !max(nOa, nOb)
  nV_m = cc_nV_m !max(nVa, nVb)
  ! Debug
  !print*,nO_m,nV_m
  
  allocate(list_occ(nO_m,2), list_vir(nV_m,2))
  list_occ = cc_list_occ_spin
  list_vir = cc_list_vir_spin
  ! Debug
  !call extract_list_orb_spin(det,nO_m,nV_m,list_occ,list_vir)
  !print*,list_occ(:,1)
  !print*,list_occ(:,2)
  !print*,list_vir(:,1)
  !print*,list_vir(:,2)

  ! Allocation
  allocate(t1(nO,nV), t2(nO,nO,nV,nV), tau(nO,nO,nV,nV), tau_t(nO,nO,nV,nV))
  allocate(r1(nO,nV), r2(nO,nO,nV,nV))
  allocate(cF_oo(nO,nO), cF_ov(nO,nV), cF_vv(nV,nV))
  allocate(cW_oooo(nO,nO,nO,nO), cW_ovvo(nO,nV,nV,nO))!, cW_vvvv(nV,nV,nV,nV))
  allocate(v_oooo(nO,nO,nO,nO))
  !allocate(v_vooo(nV,nO,nO,nO))
  allocate(v_ovoo(nO,nV,nO,nO))
  allocate(v_oovo(nO,nO,nV,nO))
  allocate(v_ooov(nO,nO,nO,nV))
  allocate(v_vvoo(nV,nV,nO,nO))
  !allocate(v_vovo(nV,nO,nV,nO))
  !allocate(v_voov(nV,nO,nO,nV))
  allocate(v_ovvo(nO,nV,nV,nO))
  allocate(v_ovov(nO,nV,nO,nV))
  allocate(v_oovv(nO,nO,nV,nV))
  !allocate(v_vvvo(nV,nV,nV,nO))
  !allocate(v_vvov(nV,nV,nO,nV))
  !allocate(v_vovv(nV,nO,nV,nV))
  !allocate(v_ovvv(nO,nV,nV,nV))
  !allocate(v_vvvv(nV,nV,nV,nV))
  allocate(f_o(nO), f_v(nV))
  allocate(f_oo(nO, nO))
  allocate(f_ov(nO, nV))
  allocate(f_vv(nV, nV))
  
  ! Allocation for the diis
  if (cc_update_method == 'diis') then
    allocate(all_err(nO*nV+nO*nO*nV*nV,cc_diis_depth), all_t(nO*nV+nO*nO*nV*nV,cc_diis_depth))
    all_err = 0d0
    all_t   = 0d0
  endif

  ! Fock elements
  call gen_f_spin(det, nO_m,nO_m, nO_S,nO_S, list_occ,list_occ, nO,nO, f_oo)
  call gen_f_spin(det, nO_m,nV_m, nO_S,nV_S, list_occ,list_vir, nO,nV, f_ov)
  call gen_f_spin(det, nV_m,nV_m, nV_S,nV_S, list_vir,list_vir, nV,nV, f_vv)

  ! Diag elements
  do i = 1, nO
    f_o(i) = f_oo(i,i)
  enddo
  do i = 1, nV
    f_v(i) = f_vv(i,i)
  enddo

  ! Bi electronic integrals from list
  call wall_time(ti)
  ! OOOO
  call gen_v_spin(nO_m,nO_m,nO_m,nO_m, nO_S,nO_S,nO_S,nO_S, list_occ,list_occ,list_occ,list_occ, nO,nO,nO,nO, v_oooo)

  ! OOO V
  !call gen_v_spin(nV_m,nO_m,nO_m,nO_m, nV_S,nO_S,nO_S,nO_S, list_vir,list_occ,list_occ,list_occ, nV,nO,nO,nO, v_vooo)
  call gen_v_spin(nO_m,nV_m,nO_m,nO_m, nO_S,nV_S,nO_S,nO_S, list_occ,list_vir,list_occ,list_occ, nO,nV,nO,nO, v_ovoo)
  call gen_v_spin(nO_m,nO_m,nV_m,nO_m, nO_S,nO_S,nV_S,nO_S, list_occ,list_occ,list_vir,list_occ, nO,nO,nV,nO, v_oovo)
  call gen_v_spin(nO_m,nO_m,nO_m,nV_m, nO_S,nO_S,nO_S,nV_S, list_occ,list_occ,list_occ,list_vir, nO,nO,nO,nV, v_ooov)

  ! OO VV
  call gen_v_spin(nV_m,nV_m,nO_m,nO_m, nV_S,nV_S,nO_S,nO_S, list_vir,list_vir,list_occ,list_occ, nV,nV,nO,nO, v_vvoo)
  !call gen_v_spin(nV_m,nO_m,nV_m,nO_m, nV_S,nO_S,nV_S,nO_S, list_vir,list_occ,list_vir,list_occ, nV,nO,nV,nO, v_vovo)
  !call gen_v_spin(nV_m,nO_m,nO_m,nV_m, nV_S,nO_S,nO_S,nV_S, list_vir,list_occ,list_occ,list_vir, nV,nO,nO,nV, v_voov)
  call gen_v_spin(nO_m,nV_m,nV_m,nO_m, nO_S,nV_S,nV_S,nO_S, list_occ,list_vir,list_vir,list_occ, nO,nV,nV,nO, v_ovvo)
  call gen_v_spin(nO_m,nV_m,nO_m,nV_m, nO_S,nV_S,nO_S,nV_S, list_occ,list_vir,list_occ,list_vir, nO,nV,nO,nV, v_ovov)
  call gen_v_spin(nO_m,nO_m,nV_m,nV_m, nO_S,nO_S,nV_S,nV_S, list_occ,list_occ,list_vir,list_vir, nO,nO,nV,nV, v_oovv)

  ! O VVV
  !call gen_v_spin(nV_m,nV_m,nV_m,nO_m, nV_S,nV_S,nV_S,nO_S, list_vir,list_vir,list_vir,list_occ, nV,nV,nV,nO, v_vvvo)
  !call gen_v_spin(nV_m,nV_m,nO_m,nV_m, nV_S,nV_S,nO_S,nV_S, list_vir,list_vir,list_occ,list_vir, nV,nV,nO,nV, v_vvov)
  !call gen_v_spin(nV_m,nO_m,nV_m,nV_m, nV_S,nO_S,nV_S,nV_S, list_vir,list_occ,list_vir,list_vir, nV,nO,nV,nV, v_vovv)
  !call gen_v_spin(nO_m,nV_m,nV_m,nV_m, nO_S,nV_S,nV_S,nV_S, list_occ,list_vir,list_vir,list_vir, nO,nV,nV,nV, v_ovvv)

  ! VVVV
  !call gen_v_spin(nV_m,nV_m,nV_m,nV_m, nV_S,nV_S,nV_S,nV_S, list_vir,list_vir,list_vir,list_vir, nV,nV,nV,nV, v_vvvv)
  call wall_time(tf)
  if (cc_dev) then
    print*,'Load bi elec int:',tf-ti,'s'
  endif

  ! Init of T
  t1 = 0d0
  call guess_t1(nO,nV,f_o,f_v,f_ov,t1)
  call guess_t2(nO,nV,f_o,f_v,v_oovv,t2)
  call compute_tau_spin(nO,nV,t1,t2,tau)
  call compute_tau_t_spin(nO,nV,t1,t2,tau_t)
  
  ! Loop init
  nb_iter = 0
  not_converged = .True.
  r1 = 0d0
  r2 = 0d0
  max_r1 = 0d0
  max_r2 = 0d0

  call det_energy(det,uncorr_energy)
  print*,'Det energy', uncorr_energy
  call ccsd_energy_spin(nO,nV,t1,t2,F_ov,v_oovv,energy)
  print*,'guess energy', uncorr_energy+energy, energy
  
  write(*,'(A77)') ' -----------------------------------------------------------------------------'
  write(*,'(A77)') ' |   It.  |       E(CCSD) (Ha) | Correlation (Ha) |  Conv. T1  |  Conv. T2  |'
  write(*,'(A77)') ' -----------------------------------------------------------------------------'

  call wall_time(ta)

  ! Loop
  do while (not_converged)

    ! Intermediates
    call wall_time(tbi)
    call wall_time(ti)
    call compute_cF_oo(nO,nV,t1,tau_t,F_oo,F_ov,v_ooov,v_oovv,cF_oo)
    call compute_cF_ov(nO,nV,t1,F_ov,v_oovv,cF_ov)
    call compute_cF_vv(nO,nV,t1,tau_t,F_ov,F_vv,v_oovv,cF_vv)
    call wall_time(tf)
    if (cc_dev) then
      print*,'Compute cFs:',tf-ti,'s'
    endif
    
    call wall_time(ti)
    call compute_cW_oooo(nO,nV,t1,t2,tau,v_oooo,v_ooov,v_oovv,cW_oooo)
    call compute_cW_ovvo(nO,nV,t1,t2,tau,v_ovvo,v_oovo,v_oovv,cW_ovvo)
    !call compute_cW_vvvv(nO,nV,t1,t2,tau,v_vvvv,v_vovv,v_oovv,cW_vvvv)
    call wall_time(tf)
    if (cc_dev) then
      print*,'Compute cFs:',tf-ti,'s'
    endif

    ! Residuals
    call wall_time(ti)
    call compute_r1_spin(nO,nV,t1,t2,f_o,f_v,F_ov,cF_oo,cF_ov,cF_vv,v_oovo,v_ovov,r1)
    call wall_time(tf)
    if (cc_dev) then
      print*,'Compute r1:',tf-ti,'s'
    endif
    call wall_time(ti)
    call compute_r2_spin(nO,nV,t1,t2,tau,f_o,f_v,cF_oo,cF_ov,cF_vv,cW_oooo,cW_ovvo,v_ovoo,v_oovv,v_ovvo,r2)
    call wall_time(tf)
    if (cc_dev) then
      print*,'Compute r2:',tf-ti,'s'
    endif

    ! Max elements in the residuals
    max_r1 = maxval(abs(r1(:,:)))
    max_r2 = maxval(abs(r2(:,:,:,:)))
    max_r  = max(max_r1,max_r2)

    call wall_time(ti)
    ! Update
    if (cc_update_method == 'diis') then
      !call update_t_ccsd(nO,nV,nb_iter,f_o,f_v,r1,r2,t1,t2,all_err1,all_err2,all_t1,all_t2)
      !call update_t_ccsd_diis(nO,nV,nb_iter,f_o,f_v,r1,r2,t1,t2,all_err1,all_err2,all_t1,all_t2)
      call update_t_ccsd_diis_v3(nO,nV,nb_iter,f_o,f_v,r1,r2,t1,t2,all_err,all_t)

    ! Standard update as T = T - Delta
    elseif (cc_update_method == 'none') then
      call update_t1(nO,nV,f_o,f_v,r1,t1)
      call update_t2(nO,nV,f_o,f_v,r2,t2)
    else
      print*,'Unkonw cc_method_method: '//cc_update_method
    endif

    call compute_tau_spin(nO,nV,t1,t2,tau)
    call compute_tau_t_spin(nO,nV,t1,t2,tau_t)
    call wall_time(tf)
    if (cc_dev) then
      print*,'Update:',tf-ti,'s'
    endif

    ! Print
    call ccsd_energy_spin(nO,nV,t1,t2,F_ov,v_oovv,energy)
    call wall_time(tfi)
    
    write(*,'(A3,I6,A3,F18.12,A3,F16.12,A3,1pE10.2,A3,1pE10.2,A2)') ' | ',nb_iter,' | ', &
         uncorr_energy+energy,' | ', energy,' | ', max_r1,' | ', max_r2,' |'
    if (cc_dev) then
      print*,'Total:',tfi-tbi,'s'
    endif

    ! Convergence
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
  write(*,'(A15,1pE10.2,A3)')' Conv        = ', max_r
  print*,''

  if (write_amplitudes) then
    call write_t1(nO,nV,t1)
    call write_t2(nO,nV,t2)
    call ezfio_set_utils_cc_io_amplitudes('Read')
  endif

  ! Deallocate
  if (cc_update_method == 'diis') then
     deallocate(all_err,all_t)
  endif
  deallocate(tau,tau_t)
  deallocate(r1,r2)
  deallocate(cF_oo,cF_ov,cF_vv)
  deallocate(cW_oooo,cW_ovvo)!,cW_vvvv)
  deallocate(v_oooo)
  deallocate(v_ovoo,v_oovo)
  deallocate(v_ovvo,v_ovov,v_oovv)
  
  double precision :: t_corr
  t_corr = 0.d0
  if (cc_par_t .and. elec_alpha_num  +elec_beta_num > 2) then
    print*,'CCSD(T) calculation...'
    call wall_time(ta)
    !allocate(v_vvvo(nV,nV,nV,nO))
    !call gen_v_spin(cc_nV_m,cc_nV_m,cc_nV_m,cc_nO_m, &
    !   cc_nV_S,cc_nV_S,cc_nV_S,cc_nO_S, &
    !   cc_list_vir_spin,cc_list_vir_spin,cc_list_vir_spin,cc_list_occ_spin, &
    !   nV,nV,nV,nO, v_vvvo)

    !call ccsd_par_t_spin(nO,nV,t1,t2,f_o,f_v,f_ov,v_ooov,v_vvoo,v_vvvo,t_corr)
    call ccsd_par_t_spin_v2(nO,nV,t1,t2,f_o,f_v,f_ov,v_ooov,v_vvoo,t_corr)
    !print*,'Working on it...'
    !call abort
    call wall_time(tb)
    print*,'Done'
    print*,'Time: ',tb-ta, ' s'
    print*,''
    write(*,'(A15,F18.12,A3)') ' E(CCSD(T))  = ', uncorr_energy + energy + t_corr, ' Ha'
    write(*,'(A15,F18.12,A3)') ' E(T)        = ', t_corr, ' Ha'
    write(*,'(A15,F18.12,A3)') ' Correlation = ', energy + t_corr, ' Ha'
    print*,''
  endif

  call save_energy(uncorr_energy + energy, t_corr)
  
  deallocate(f_oo,f_ov,f_vv,f_o,f_v)
  deallocate(v_ooov,v_vvoo,t1,t2)
  !deallocate(v_ovvv,v_vvvo,v_vovv)
  !deallocate(v_vvvv)
  
end

! Energy

subroutine ccsd_energy_spin(nO,nV,t1,t2,Fov,v_oovv,energy)

  implicit none

  BEGIN_DOC
  ! CCSD energy in spin orbitals
  END_DOC

  integer,intent(in)            :: nO,nV
  double precision,intent(in)   :: t1(nO,nV)
  double precision,intent(in)   :: t2(nO,nO,nV,nV)
  double precision,intent(in)   :: Fov(nO,nV)
  double precision,intent(in)   :: v_oovv(nO,nO,nV,nV)

  double precision,intent(out)  :: energy

  integer                       :: i,j,a,b


  energy = 0d0

  do i=1,nO
      do a=1,nV
      energy = energy + Fov(i,a) * t1(i,a)
    end do
  end do

  do i=1,nO
    do j=1,nO
      do a=1,nV
        do b=1,nV
          energy = energy                                & 
                 + 0.5d0 * v_oovv(i,j,a,b) * t1(i,a) * t1(j,b) &
                 + 0.25d0 * v_oovv(i,j,a,b) * t2(i,j,a,b)
        end do
      end do
    end do
  end do

end

! Tau

subroutine compute_tau_spin(nO,nV,t1,t2,tau)

  implicit none

  integer,intent(in)            :: nO,nV
  double precision,intent(in)   :: t1(nO,nV)
  double precision,intent(in)   :: t2(nO,nO,nV,nV)

  double precision,intent(out)  :: tau(nO,nO,nV,nV)
  
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d

  !$OMP PARALLEL &
  !$OMP SHARED(tau,t1,t2,nO,nV) &
  !$OMP PRIVATE(i,j,a,b) &
  !$OMP DEFAULT(NONE)
  !$OMP DO collapse(3)
  do i=1,nO
    do j=1,nO
      do a=1,nV
        do b=1,nV
          tau(i,j,a,b) = t2(i,j,a,b) + t1(i,a)*t1(j,b) - t1(i,b)*t1(j,a)
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

end

! Tau_t

subroutine compute_tau_t_spin(nO,nV,t1,t2,tau_t)

  implicit none

  integer,intent(in)            :: nO,nV
  double precision,intent(in)   :: t1(nO,nV)
  double precision,intent(in)   :: t2(nO,nO,nV,nV)

  double precision,intent(out)  :: tau_t(nO,nO,nV,nV)

  integer                       :: i,j,k,l
  integer                       :: a,b,c,d

  !$OMP PARALLEL &
  !$OMP SHARED(tau_t,t1,t2,nO,nV) &
  !$OMP PRIVATE(i,j,a,b) &
  !$OMP DEFAULT(NONE)
  !$OMP DO collapse(3)
  do i=1,nO
    do j=1,nO
      do a=1,nV
        do b=1,nV
          tau_t(i,j,a,b) = t2(i,j,a,b) + 0.5d0*(t1(i,a)*t1(j,b) - t1(i,b)*t1(j,a))
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

end

! R1

subroutine compute_r1_spin(nO,nV,t1,t2,f_o,f_v,Fov,cF_oo,cF_ov,cF_vv,v_oovo,v_ovov,r1)

  implicit none

  integer,intent(in)            :: nO,nV
  double precision,intent(in)   :: t1(nO,nV)
  double precision,intent(in)   :: t2(nO,nO,nV,nV)
  double precision,intent(in)   :: f_o(nO), f_v(nV)
  double precision,intent(in)   :: Fov(nO,nV)
  double precision,intent(in)   :: cF_oo(nO,nO)
  double precision,intent(in)   :: cF_ov(nO,nV)
  double precision,intent(in)   :: cF_vv(nV,nV)
  double precision,intent(in)   :: v_oovo(nO,nO,nV,nO)
  double precision,intent(in)   :: v_ovov(nO,nV,nO,nV)
  !double precision,intent(in)   :: v_ovvv(nO,nV,nV,nV)

  double precision,intent(out)  :: r1(nO,nV)

  integer                       :: i,j,m,n
  integer                       :: a,b,e,f
  !double precision, allocatable :: X_vovv(:,:,:,:)
  double precision, allocatable :: X_oovv(:,:,:,:)
  double precision              :: accu

  !$OMP PARALLEL &
  !$OMP SHARED(r1,t1,t2,Fov,cF_vv,cF_ov, &
  !$OMP v_ovov,nO,nV) &
  !$OMP PRIVATE(i,j,a,b,e,f,m,n) &
  !$OMP DEFAULT(NONE)
  
  !$OMP DO collapse(1)
  do a=1,nV
    do i=1,nO
      r1(i,a) = Fov(i,a)
      do e=1,nV
        do m=1,nO
          r1(i,a) = r1(i,a) + t2(i,m,a,e)*cF_ov(m,e)
        end do
      end do
      do f=1,nV
        do n=1,nO
          r1(i,a) = r1(i,a) - t1(n,f)*v_ovov(n,a,i,f)
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  !do a=1,nV
  !  do i=1,nO
  !    do e=1,nV
  !      r1(i,a) = r1(i,a) + t1(i,e)*cF_vv(a,e)
  !    end do
  !  end do
  !end do
  call dgemm('N','T', nO, nV, nV, &
             1d0, t1   , size(t1,1), &
                  cF_vv, size(cF_vv,1), &
             1d0, r1   , size(r1,1))
  
  !do a=1,nV
  !  do i=1,nO
  !    do m=1,nO
  !      r1(i,a) = r1(i,a) - t1(m,a)*cF_oo(m,i)
  !    end do
  !  end do
  !end do
  call dgemm('T','N', nO, nV, nO, &
             -1d0, cF_oo, size(cF_oo,1), &
                   t1   , size(t1,1), &
              1d0, r1   , size(r1,1))

  !do a=1,nV
  !  do i=1,nO
  !    do f=1,nV
  !      do e=1,nV
  !        do m=1,nO
  !          r1(i,a) = r1(i,a) - 0.5d0*t2(i,m,e,f)*v_ovvv(m,a,e,f)
  !        end do
  !      end do
  !    end do
  !  end do
  !end do

  !allocate(X_vovv(nV,nO,nV,nV))
  double precision, allocatable :: v_ovvf(:,:,:), X_vovf(:,:,:)
  allocate(v_ovvf(nO,nV,nV),X_vovf(nV,nO,nV))

  do f = 1, nV
    call gen_v_spin_3idx(cc_nO_m,cc_nV_m,cc_nV_m,cc_nV_m, f, cc_nO_S,cc_nV_S,cc_nV_S,cc_nV_S, &
           cc_list_occ_spin,cc_list_vir_spin,cc_list_vir_spin,cc_list_vir_spin, &
           nO,nV,nV, v_ovvf)
    !$OMP PARALLEL &
    !$OMP SHARED(r1,t1,t2,X_vovf,v_ovvf,nO,nV) &
    !$OMP PRIVATE(i,j,a,b,e,f,m,n) &
    !$OMP DEFAULT(NONE)
    

    !$OMP DO collapse(3)
    !do f = 1, nV
      do e = 1, nV
         do m = 1, nO
           do a = 1, nV
             !X_vovv(a,m,e,f) = v_ovvv(m,a,e,f)
             X_vovf(a,m,e) = v_ovvf(m,a,e)
          enddo
        enddo
      enddo
    !enddo
    !$OMP END DO
    !$OMP END PARALLEL
      
    call dgemm('N','T', nO, nV, nO*nV, &
             -0.5d0, t2(1,1,1,f), size(t2,1), &
                     X_vovf, size(X_vovf,1), &
              1d0  , r1    , size(r1,1))
  enddo
  
  !call dgemm('N','T', nO, nV, nO*nV*nV, &
  !           -0.5d0, t2    , size(t2,1), &
  !                   X_vovv, size(X_vovv,1), &
  !            1d0  , r1    , size(r1,1))
  
  deallocate(X_vovf)
  !deallocate(X_vovv)
  allocate(X_oovv(nO,nO,nV,nV))
  
  !$OMP PARALLEL &
  !$OMP SHARED(r1,t1,t2,X_oovv, &
  !$OMP f_o,f_v,v_oovo,nO,nV) &
  !$OMP PRIVATE(i,j,a,b,e,f,m,n) &
  !$OMP DEFAULT(NONE)
  
  !do a=1,nV
  !  do i=1,nO
  !    do e=1,nV
  !      do m=1,nO
  !        do n=1,nO
  !          r1(i,a) = r1(i,a) - 0.5d0*t2(m,n,a,e)*v_oovo(n,m,e,i)
  !        end do
  !      end do
  !    end do
  !  end do
  !end do
  
  !$OMP DO collapse(3)
  do a = 1, nV
    do e = 1, nV
      do m = 1, nO
        do n = 1, nO
          X_oovv(n,m,e,a) = t2(m,n,a,e)
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  
  call dgemm('T','N', nO, nV, nO*nO*nV, &
             -0.5d0, v_oovo, size(v_oovo,1) * size(v_oovo,2) * size(v_oovo,3), &
                     X_oovv, size(X_oovv,1) * size(X_oovv,2) * size(X_oovv,3), &
             1d0   , r1    , size(r1,1))
  
  !$OMP PARALLEL &
  !$OMP SHARED(r1,t1,X_oovv,f_o,f_v,nO,nV) &
  !$OMP PRIVATE(i,j,a,b,e,f,m,n) &
  !$OMP DEFAULT(NONE)
  !$OMP DO collapse(1)
  do a = 1, nV
    do i = 1, nO
      r1(i,a) = (f_o(i)-f_v(a)) * t1(i,a) - r1(i,a)
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  
  deallocate(X_oovv)

end

! R2

subroutine compute_r2_spin(nO,nV,t1,t2,tau,f_o,f_v,cF_oo,cF_ov,cF_vv,cW_oooo,cW_ovvo,v_ovoo,v_oovv,v_ovvo,r2)

  implicit none

  integer,intent(in)            :: nO,nV
  double precision,intent(in)   :: cF_oo(nO,nO)
  double precision,intent(in)   :: cF_ov(nO,nV)
  double precision,intent(in)   :: cF_vv(nV,nV)
  double precision,intent(in)   :: f_o(nO), f_v(nV)
  double precision,intent(in)   :: cW_oooo(nO,nO,nO,nO)
  !double precision,intent(in)   :: cW_vvvv(nV,nV,nV,nV)
  double precision,intent(in)   :: cW_ovvo(nO,nV,nV,nO)
  double precision,intent(in)   :: t1(nO,nV)
  double precision,intent(in)   :: t2(nO,nO,nV,nV)
  double precision,intent(in)   :: tau(nO,nO,nV,nV)
  double precision,intent(in)   :: v_ovoo(nO,nV,nO,nO)
  double precision,intent(in)   :: v_oovv(nO,nO,nV,nV)
  double precision,intent(in)   :: v_ovvo(nO,nV,nV,nO)
  !double precision,intent(in)   :: v_vvvo(nV,nV,nV,nO)!, v_vovv(nV,nO,nV,nV)

  double precision,intent(out)  :: r2(nO,nO,nV,nV)

  integer                       :: i,j,m,n
  integer                       :: a,b,e,f
  double precision, allocatable :: X_vvoo(:,:,:,:)
  !double precision, allocatable :: A_vvov(:,:,:,:)
  double precision, allocatable :: X_oovv(:,:,:,:), Y_oovv(:,:,:,:)
  double precision, allocatable :: A_vvoo(:,:,:,:), B_ovoo(:,:,:,:), C_ovov(:,:,:,:)
  double precision, allocatable :: A_ovov(:,:,:,:), B_ovvo(:,:,:,:), X_ovvo(:,:,:,:)
  double precision, allocatable :: A_vv(:,:)
  double precision, allocatable :: A_oo(:,:), B_oovv(:,:,:,:)
  double precision, allocatable :: A_vbov(:,:,:), X_vboo(:,:,:), v_vbvo(:,:,:)

  !do b=1,nV
  !  do a=1,nV
  !    do j=1,nO
  !      do i=1,nO
  !        r2(i,j,a,b) = v_oovv(i,j,a,b)
  !      end do
  !    end do
  !  end do
  !end do

  !do b=1,nV
  !  do a=1,nV
  !    do j=1,nO
  !      do i=1,nO

  !        do e=1,nV
  !          r2(i,j,a,b) = r2(i,j,a,b) + t2(i,j,a,e)*cF_vv(b,e)
  !          r2(i,j,a,b) = r2(i,j,a,b) - t2(i,j,b,e)*cF_vv(a,e)
  !        end do

  !      end do
  !    end do
  !  end do
  !end do
  allocate(X_oovv(nO,nO,nV,nV))
  call dgemm('N','T',nO*nO*nV, nV, nV, &
             1d0, t2    , size(t2,1) * size(t2,2) * size(t2,3), &
                  cF_VV , size(cF_vv,1), &
             0d0, X_oovv, size(X_oovv,1) * size(X_oovv,2) * size(X_oovv,3))

  !$OMP PARALLEL &
  !$OMP SHARED(r2,v_oovv,X_oovv,nO,nV) &
  !$OMP PRIVATE(i,j,a,b) &
  !$OMP DEFAULT(NONE)
  
  !$OMP DO collapse(3)
  do b=1,nV
    do a=1,nV
      do j=1,nO
        do i=1,nO
          r2(i,j,a,b) = v_oovv(i,j,a,b) + X_oovv(i,j,a,b) - X_oovv(i,j,b,a)
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  
  !deallocate(X_oovv)

  !do b=1,nV
  !  do a=1,nV
  !    do j=1,nO
  !      do i=1,nO

  !        do e=1,nV
  !          do m=1,nO
  !            r2(i,j,a,b) = r2(i,j,a,b) - 0.5d0*t2(i,j,a,e)*t1(m,b)*cF_ov(m,e)
  !            r2(i,j,a,b) = r2(i,j,a,b) + 0.5d0*t2(i,j,b,e)*t1(m,a)*cF_ov(m,e)
  !          end do
  !        end do

  !      end do
  !    end do
  !  end do
  !end do
  allocate(A_vv(nV,nV))!, X_oovv(nO,nO,nV,nV))
  call dgemm('T','N', nV, nV, nO, &
             1d0, t1   , size(t1,1), &
                  cF_ov, size(cF_ov,1), &
             0d0, A_vv , size(A_vv,1))

  call dgemm('N','T', nO*nO*nV, nV, nV, &
             0.5d0, t2    , size(t2,1) * size(t2,2) * size(t2,3), &
                    A_vv  , size(A_vv,1), &
             0d0  , X_oovv, size(X_oovv,1) * size(X_oovv,2) * size(X_oovv,3))
  
  !$OMP PARALLEL &
  !$OMP SHARED(r2,v_oovv,X_oovv,nO,nV) &
  !$OMP PRIVATE(i,j,a,b) &
  !$OMP DEFAULT(NONE)
  
  !$OMP DO collapse(3)
  do b=1,nV
    do a=1,nV
      do j=1,nO
        do i=1,nO
          r2(i,j,a,b) = r2(i,j,a,b) - X_oovv(i,j,a,b) + X_oovv(i,j,b,a) 
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
             
  deallocate(A_vv)!,X_oovv)

  !do b=1,nV
  !  do a=1,nV
  !    do j=1,nO
  !      do i=1,nO

  !        do m=1,nO
  !          r2(i,j,a,b) = r2(i,j,a,b) - t2(i,m,a,b)*cF_oo(m,j)
  !          r2(i,j,a,b) = r2(i,j,a,b) + t2(j,m,a,b)*cF_oo(m,i)
  !        end do

  !      end do
  !    end do
  !  end do
  !end do
  allocate(Y_oovv(nO,nO,nV,nV))!,X_oovv(nO,nO,nV,nV))
  !$OMP PARALLEL &
  !$OMP SHARED(t2,v_oovv,X_oovv,nO,nV) &
  !$OMP PRIVATE(i,m,a,b) &
  !$OMP DEFAULT(NONE)
  
  !$OMP DO collapse(3)
  do b=1,nV
    do a=1,nV
      do i=1,nO
        do m=1,nO
          X_oovv(m,i,a,b) = t2(i,m,a,b)
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  call dgemm('T','N', nO, nO*nV*nV, nO, &
             1d0, cF_oo , size(cF_oo,1), &
                  X_oovv, size(X_oovv,1), &
             0d0, Y_oovv, size(Y_oovv,1))

  !$OMP PARALLEL &
  !$OMP SHARED(r2,v_oovv,Y_oovv,nO,nV) &
  !$OMP PRIVATE(i,j,a,b) &
  !$OMP DEFAULT(NONE)
  
  !$OMP DO collapse(3)
  do b=1,nV
    do a=1,nV
      do j=1,nO
        do i=1,nO
          r2(i,j,a,b) = r2(i,j,a,b) - Y_oovv(j,i,a,b) + Y_oovv(i,j,a,b) 
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  deallocate(Y_oovv)!,X_oovv)

  !do b=1,nV
  !  do a=1,nV
  !    do j=1,nO
  !      do i=1,nO

  !        do e=1,nV
  !          do m=1,nO
  !            r2(i,j,a,b) = r2(i,j,a,b) - 0.5d0*t2(i,m,a,b)*t1(j,e)*cF_ov(m,e)
  !            r2(i,j,a,b) = r2(i,j,a,b) + 0.5d0*t2(j,m,a,b)*t1(i,e)*cF_ov(m,e)
  !          end do
  !        end do

  !      end do
  !    end do
  !  end do
  !end do
  allocate(A_oo(nO,nO),B_oovv(nO,nO,nV,nV))!,X_oovv(nO,nO,nV,nV))
  
  call dgemm('N','T', nO, nO, nV, &
        1d0, t1   , size(t1,1), &
             cF_ov, size(cF_ov,1), &
        0d0, A_oo , size(A_oo,1))
  
  !$OMP PARALLEL &
  !$OMP SHARED(t2,B_oovv,nO,nV) &
  !$OMP PRIVATE(i,m,a,b) &
  !$OMP DEFAULT(NONE)
  
  !$OMP DO collapse(3)
  do b = 1, nV
    do a = 1, nV
      do i = 1, nO
        do m = 1, nO
          B_oovv(m,i,a,b) = t2(i,m,a,b)
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  
  call dgemm('N','N', nO, nO*nV*nV, nO, &
             0.5d0, A_oo, size(A_oo,1), &
                    B_oovv, size(B_oovv,1), &
             0d0  , X_oovv, size(X_oovv,1))
  
  !$OMP PARALLEL &
  !$OMP SHARED(r2,X_oovv,nO,nV) &
  !$OMP PRIVATE(i,j,a,b) &
  !$OMP DEFAULT(NONE)
  
  !$OMP DO collapse(3)
  do b=1,nV
    do a=1,nV
      do j=1,nO
        do i=1,nO
          r2(i,j,a,b) = r2(i,j,a,b) - X_oovv(j,i,a,b) + X_oovv(i,j,a,b)
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  deallocate(A_oo,B_oovv,X_oovv)

  !do b=1,nV
  !  do a=1,nV
  !    do j=1,nO
  !      do i=1,nO

  !        do n=1,nO
  !          do m=1,nO
  !            r2(i,j,a,b) = r2(i,j,a,b) + 0.5d0*tau(m,n,a,b)*cW_oooo(m,n,i,j)
  !          end do
  !        end do

  !      end do
  !    end do
  !  end do
  !end do
  call dgemm('T','N', nO*nO, nV*nV, nO*nO, &
             0.5d0, cW_oooo, size(cW_oooo,1) * size(cW_oooo,2), &
                    tau    , size(tau,1) * size(tau,2), &
             1d0  , r2     , size(r2,1) * size(r2,2))
  
  !do b=1,nV
  !  do a=1,nV
  !    do j=1,nO
  !      do i=1,nO

  !        do f=1,nV
  !          do e=1,nV
  !            r2(i,j,a,b) = r2(i,j,a,b) + 0.5d0*tau(i,j,e,f)*cW_vvvv(a,b,e,f)
  !          end do
  !        end do

  !      end do
  !    end do
  !  end do
  !end do
  !call dgemm('N','T', nO*nO, nV*nV, nV*nV, &
  !           0.5d0, tau    , size(tau,1) * size(tau,2), &
  !                  cW_vvvv, size(cW_vvvv,1) * size(cW_vvvv,2), &
  !           1d0  , r2     , size(r2,1) * size(r2,2))
  double precision :: ti,tf
  call wall_time(ti)
  call use_cW_vvvf(nO,nV,t1,t2,tau,v_oovv,r2)
  call wall_time(tf)
  if (cc_dev) then
    print*,'cW_vvvv:',tf-ti,'s'
  endif
  
  !do b=1,nV
  !  do a=1,nV
  !    do j=1,nO
  !      do i=1,nO

  !        do e=1,nV
  !          do m=1,nO
  !            r2(i,j,a,b) = r2(i,j,a,b)                                                 & 
  !                        + t2(i,m,a,e)*cW_ovvo(m,b,e,j) &
  !                        - t2(j,m,a,e)*cW_ovvo(m,b,e,i) &
  !                        - t2(i,m,b,e)*cW_ovvo(m,a,e,j) &
  !                        + t2(j,m,b,e)*cW_ovvo(m,a,e,i) &
  !                        - t1(i,e)*t1(m,a)*v_ovvo(m,b,e,j) &
  !                        + t1(j,e)*t1(m,a)*v_ovvo(m,b,e,i) &
  !                        + t1(i,e)*t1(m,b)*v_ovvo(m,a,e,j) &
  !                        - t1(j,e)*t1(m,b)*v_ovvo(m,a,e,i)
  !          end do
  !        end do

  !      end do
  !    end do
  !  end do
  !end do
  allocate(A_ovov(nO,nV,nO,nV), B_ovvo(nO,nV,nV,nO), X_ovvo(nO,nV,nV,nO))
  !$OMP PARALLEL &
  !$OMP SHARED(t2,A_ovov,B_ovvo,cW_ovvo,nO,nV) &
  !$OMP PRIVATE(i,j,a,b,e,m) &
  !$OMP DEFAULT(NONE)
  
  !$OMP DO collapse(3)
  do a = 1, nV
    do i = 1, nO
      do e = 1, nV
        do m = 1, nO
          A_ovov(m,e,i,a) = t2(i,m,a,e)
        end do
      end do
    end do
  end do
  !$OMP END DO NOWAIT
  !$OMP DO collapse(3)
  do j = 1, nO
    do b = 1, nV
      do e = 1, nV
        do m = 1, nO
          B_ovvo(m,e,b,j) = cW_ovvo(m,b,e,j) 
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  
  call dgemm('T','N', nO*nV, nV*nO, nO*nV, &
             1d0, A_ovov, size(A_ovov,1) * size(A_ovov,2), &
                  B_ovvo, size(B_ovvo,1) * size(B_ovvo,2), &
             0d0, X_ovvo, size(X_ovvo,1) * size(X_ovvo,2))
  
  !$OMP PARALLEL &
  !$OMP SHARED(r2,X_ovvo,nO,nV) &
  !$OMP PRIVATE(i,j,a,b,e,m) &
  !$OMP DEFAULT(NONE)
  
  !$OMP DO collapse(3)
  do b = 1, nV
    do a = 1, nV
      do j = 1, nO
        do i = 1, nO
          r2(i,j,a,b) = r2(i,j,a,b) + X_ovvo(i,a,b,j) - X_ovvo(j,a,b,i) &
                                    - X_ovvo(i,b,a,j) + X_ovvo(j,b,a,i)
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  
  deallocate(A_ovov,B_ovvo,X_ovvo)
  allocate(A_vvoo(nV,nV,nO,nO), B_ovoo(nO,nV,nO,nO), C_ovov(nO,nV,nO,nV))
  
  !$OMP PARALLEL &
  !$OMP SHARED(A_vvoo,v_ovvo,nO,nV) &
  !$OMP PRIVATE(i,j,a,b,e,m) &
  !$OMP DEFAULT(NONE)
  
  !$OMP DO collapse(3)
  do m = 1, nO
    do j = 1, nO
      do b = 1, nV
        do e = 1, nV
          A_vvoo(e,b,j,m) = v_ovvo(m,b,e,j)
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  
  call dgemm('N','N', nO, nV*nO*nO, nV, &
             1d0, t1    , size(t1,1), &
                  A_vvoo, size(A_vvoo,1), &
             0d0, B_ovoo, size(B_ovoo,1))
  
  call dgemm('N','N', nO*nV*nO, nV, nO, &
             1d0, B_ovoo, size(B_ovoo,1) * size(B_ovoo,2) * size(B_ovoo,3), &
                  t1    , size(t1,1), &
             0d0, C_ovov, size(C_ovov,1) * size(C_ovov,2) * size(C_ovov,3))
  
  !$OMP PARALLEL &
  !$OMP SHARED(r2,C_ovov,nO,nV) &
  !$OMP PRIVATE(i,j,a,b,e,m) &
  !$OMP DEFAULT(NONE)
  
  !$OMP DO collapse(3)
  do b=1,nV
    do a=1,nV
      do j=1,nO
        do i=1,nO
          r2(i,j,a,b) = r2(i,j,a,b) - C_ovov(i,b,j,a) + C_ovov(j,b,i,a) &
                                    + C_ovov(i,a,j,b) - C_ovov(j,a,i,b)
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  
  deallocate(A_vvoo, B_ovoo, C_ovov)
                  
  !do b=1,nV
  !  do a=1,nV
  !    do j=1,nO
  !      do i=1,nO

  !        do e=1,nV
  !          r2(i,j,a,b) = r2(i,j,a,b) + t1(i,e)*v_vvvo(a,b,e,j) - t1(j,e)*v_vvvo(a,b,e,i)
  !        end do

  !      end do
  !    end do
  !  end do
  !end do
  !allocate(A_vvov(nV,nV,nO,nV), X_vvoo(nV,nV,nO,nO))
  allocate(A_vbov(nV,nO,nV), X_vboo(nV,nO,nO), v_vbvo(nV,nV,nO))
  do b = 1, nV

    call gen_v_spin_3idx_i_kl(cc_nV_m,cc_nV_m,cc_nV_m,cc_nO_m, b, cc_nV_S,cc_nV_S,cc_nV_S,cc_nO_S, &
         cc_list_vir_spin,cc_list_vir_spin,cc_list_vir_spin,cc_list_occ_spin, &
         nV,nV,nO, v_vbvo)
    
    !$OMP PARALLEL &
    !$OMP SHARED(b,A_vbov,v_vbvo,nO,nV) &
    !$OMP PRIVATE(i,j,a,e,m) &
    !$OMP DEFAULT(NONE)
    
    !$OMP DO collapse(2)
    do e = 1, nV
      do j = 1, nO
        !do b = 1, nV
          do a = 1, nV
            !A_vvov(a,b,j,e) = v_vvvo(a,b,e,j)
            A_vbov(a,j,e) = v_vbvo(a,e,j)
          enddo
        !enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    call dgemm('N','T', nV*nO, nO, nV, &
               1d0, A_vbov, size(A_vbov,1) * size(A_vbov,2), &
                    t1    , size(t1,1), &
               0d0, X_vboo, size(X_vboo,1) * size(X_vboo,2))
    !call dgemm('N','T', nV*nV*nO, nO, nV, &
    !           1d0, A_vvov, size(A_vvov,1) * size(A_vvov,2) * size(A_vvov,3), &
    !                t1    , size(t1,1), &
    !           0d0, X_vvoo, size(X_vvoo,1) * size(X_vvoo,2) * size(X_vvoo,3))
    
    !$OMP PARALLEL &
    !$OMP SHARED(b,r2,X_vboo,nO,nV) &
    !$OMP PRIVATE(i,j,a,e,m) &
    !$OMP DEFAULT(NONE)
    
    !$OMP DO collapse(2)
    !do b = 1, nV
      do a = 1, nV
        do j = 1, nO
          do i = 1, NO
             !r2(i,j,a,b ) = r2(i,j,a,b) + X_vvoo(a,b,j,i) - X_vvoo(a,b,i,j)
             r2(i,j,a,b) = r2(i,j,a,b) + X_vboo(a,j,i) - X_vboo(a,i,j)
          enddo
        enddo
      enddo
    !enddo
    !$OMP END DO
    !$OMP END PARALLEL
  enddo
  
  !deallocate(A_vvov)!,X_vvoo)
  deallocate(A_vbov, X_vboo, v_vbvo)
  allocate(X_vvoo(nV,nV,nO,nO))

  !do b=1,nV
  !  do a=1,nV
  !    do j=1,nO
  !      do i=1,nO

  !        do m=1,nO
  !          r2(i,j,a,b) = r2(i,j,a,b) - t1(m,a)*v_ovoo(m,b,i,j) + t1(m,b)*v_ovoo(m,a,i,j)
  !        end do

  !      end do
  !    end do
  !  end do
  !end do
  !allocate(X_vvoo(nV,nV,nO,nO))
  
  call dgemm('T','N', nV, nV*nO*nO, nO, &
             1d0, t1    , size(t1,1), &
                  v_ovoo, size(v_ovoo,1), &
             0d0, X_vvoo, size(X_vvoo,1))

  !$OMP PARALLEL &
  !$OMP SHARED(r2,X_vvoo,f_o,f_v,t2,nO,nV) &
  !$OMP PRIVATE(i,j,a,b,e,m) &
  !$OMP DEFAULT(NONE)
  
  !$OMP DO collapse(3)
  do b=1,nV
    do a=1,nV
      do j=1,nO
        do i=1,nO
          r2(i,j,a,b) = r2(i,j,a,b) - X_vvoo(a,b,i,j) + X_vvoo(b,a,i,j)
        end do
      end do
    end do
  end do
  !$OMP END DO
  
  !$OMP DO collapse(3)
  do b=1,nV
    do a=1,nV
      do j=1,nO
        do i=1,nO
          r2(i,j,a,b) = (f_o(i)+f_o(j)-f_v(a)-f_v(b)) * t2(i,j,a,b) - r2(i,j,a,b)
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  
  deallocate(X_vvoo)

end

! Use cF_oo

subroutine use_cF_oo(nO,nV,t1,t2,tau_t,F_oo,F_ov,v_ooov,v_oovv,r1,r2)

  implicit none

  integer,intent(in)              :: nO,nV
  double precision, intent(in)    :: t1(nO,nV), t2(nO,nO,nV,nV), tau_t(nO,nO,nV,nV)
  double precision, intent(in)    :: F_oo(nO,nV), F_ov(nO,nV)
  double precision, intent(in)    :: v_ooov(nO,nO,nO,nV), v_oovv(nO,nO,nV,nV)
  
  double precision, intent(inout) :: r1(nO,nV), r2(nO,nO,nV,nV)
  
  double precision, allocatable   :: cF_oo(:,:), X_oovv(:,:,:,:),Y_oovv(:,:,:,:)
  integer                         :: i,j,m,a,b

  allocate(cF_oo(nO,nO))
  
  call compute_cF_oo(nO,nV,t1,tau_t,F_oo,F_ov,v_ooov,v_oovv,cF_oo)
  
  !do a=1,nV
  !  do i=1,nO
  !    do m=1,nO
  !      r1(i,a) = r1(i,a) - t1(m,a)*cF_oo(m,i)
  !    end do
  !  end do
  !end do
  call dgemm('T','N', nO, nV, nO, &
             -1d0, cF_oo, size(cF_oo,1), &
                   t1   , size(t1,1), &
              1d0, r1   , size(r1,1))

  !do b=1,nV
  !  do a=1,nV
  !    do j=1,nO
  !      do i=1,nO

  !        do m=1,nO
  !          r2(i,j,a,b) = r2(i,j,a,b) - t2(i,m,a,b)*cF_oo(m,j)
  !          r2(i,j,a,b) = r2(i,j,a,b) + t2(j,m,a,b)*cF_oo(m,i)
  !        end do

  !      end do
  !    end do
  !  end do
  !end do
  
  allocate(Y_oovv(nO,nO,nV,nV),X_oovv(nO,nO,nV,nV))
  !$OMP PARALLEL &
  !$OMP SHARED(t2,v_oovv,X_oovv,nO,nV) &
  !$OMP PRIVATE(i,m,a,b) &
  !$OMP DEFAULT(NONE)
  
  !$OMP DO collapse(3)
  do b=1,nV
    do a=1,nV
      do i=1,nO
        do m=1,nO
          X_oovv(m,i,a,b) = t2(i,m,a,b)
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  call dgemm('T','N', nO, nO*nV*nV, nO, &
             1d0, cF_oo , size(cF_oo,1), &
                  X_oovv, size(X_oovv,1), &
             0d0, Y_oovv, size(Y_oovv,1))

  !$OMP PARALLEL &
  !$OMP SHARED(r2,v_oovv,Y_oovv,nO,nV) &
  !$OMP PRIVATE(i,j,a,b) &
  !$OMP DEFAULT(NONE)
  
  !$OMP DO collapse(3)
  do b=1,nV
    do a=1,nV
      do j=1,nO
        do i=1,nO
          r2(i,j,a,b) = r2(i,j,a,b) - Y_oovv(j,i,a,b) + Y_oovv(i,j,a,b) 
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  
  deallocate(cF_oo,X_oovv,Y_oovv)

end

! Use cF_ov

subroutine use_cF_ov(nO,nV,t1,t2,F_ov,v_oovv,r1,r2)

  implicit none

  integer, intent(in)             :: nO,nV
  double precision, intent(in)    :: t1(nO,nV), t2(nO,nO,nV,nV)
  double precision, intent(in)    :: F_ov(nO,nV), v_oovv(nO,nO,nV,nV)
  
  double precision, intent(inout) :: r1(nO,nV), r2(nO,nO,nV,nV)

  double precision, allocatable   :: cF_ov(:,:), A_oo(:,:), A_vv(:,:)
  double precision, allocatable   :: X_oovv(:,:,:,:), B_oovv(:,:,:,:)
  integer                         :: i,j,a,b,e,m

  allocate(cF_ov(nO,nV))
  
  call compute_cF_ov(nO,nV,t1,F_ov,v_oovv,cF_ov)

  !$OMP PARALLEL &
  !$OMP SHARED(r1,t2,cF_ov,nO,nV) &
  !$OMP PRIVATE(i,a,e,m) &
  !$OMP DEFAULT(NONE)
  
  !$OMP DO collapse(1)
  do a=1,nV
    do i=1,nO
      do e=1,nV
        do m=1,nO
          r1(i,a) = r1(i,a) + t2(i,m,a,e)*cF_ov(m,e)
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  !do b=1,nV
  !  do a=1,nV
  !    do j=1,nO
  !      do i=1,nO

  !        do e=1,nV
  !          do m=1,nO
  !            r2(i,j,a,b) = r2(i,j,a,b) - 0.5d0*t2(i,j,a,e)*t1(m,b)*cF_ov(m,e)
  !            r2(i,j,a,b) = r2(i,j,a,b) + 0.5d0*t2(i,j,b,e)*t1(m,a)*cF_ov(m,e)
  !          end do
  !        end do

  !      end do
  !    end do
  !  end do
  !end do
  allocate(A_vv(nV,nV), X_oovv(nO,nO,nV,nV))
  call dgemm('T','N', nV, nV, nO, &
             1d0, t1   , size(t1,1), &
                  cF_ov, size(cF_ov,1), &
             0d0, A_vv , size(A_vv,1))

  call dgemm('N','T', nO*nO*nV, nV, nV, &
             0.5d0, t2    , size(t2,1) * size(t2,2) * size(t2,3), &
                    A_vv  , size(A_vv,1), &
             0d0  , X_oovv, size(X_oovv,1) * size(X_oovv,2) * size(X_oovv,3))

  !$OMP PARALLEL &
  !$OMP SHARED(nO,nV,r2,X_oovv) &
  !$OMP PRIVATE(i,j,a,b) &
  !$OMP DEFAULT(NONE)
  
  !$OMP DO collapse(3)
  do b=1,nV
    do a=1,nV
      do j=1,nO
        do i=1,nO
          r2(i,j,a,b) = r2(i,j,a,b) - X_oovv(i,j,a,b) + X_oovv(i,j,b,a) 
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
             
  deallocate(A_vv)
  
  !do b=1,nV
  !  do a=1,nV
  !    do j=1,nO
  !      do i=1,nO

  !        do e=1,nV
  !          do m=1,nO
  !            r2(i,j,a,b) = r2(i,j,a,b) - 0.5d0*t2(i,m,a,b)*t1(j,e)*cF_ov(m,e)
  !            r2(i,j,a,b) = r2(i,j,a,b) + 0.5d0*t2(j,m,a,b)*t1(i,e)*cF_ov(m,e)
  !          end do
  !        end do

  !      end do
  !    end do
  !  end do
  !end do
  allocate(A_oo(nO,nO),B_oovv(nO,nO,nV,nV))!,X_oovv(nO,nO,nV,nV))
  
  call dgemm('N','T', nO, nO, nV, &
        1d0, t1   , size(t1,1), &
             cF_ov, size(cF_ov,1), &
        0d0, A_oo , size(A_oo,1))
  
  !$OMP PARALLEL &
  !$OMP SHARED(t2,B_oovv,nO,nV) &
  !$OMP PRIVATE(i,m,a,b) &
  !$OMP DEFAULT(NONE)
  
  !$OMP DO collapse(3)
  do b = 1, nV
    do a = 1, nV
      do i = 1, nO
        do m = 1, nO
          B_oovv(m,i,a,b) = t2(i,m,a,b)
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  
  call dgemm('N','N', nO, nO*nV*nV, nO, &
             0.5d0, A_oo, size(A_oo,1), &
                    B_oovv, size(B_oovv,1), &
             0d0  , X_oovv, size(X_oovv,1))

  !$OMP PARALLEL &
  !$OMP SHARED(r2,X_oovv,nO,nV) &
  !$OMP PRIVATE(i,j,a,b) &
  !$OMP DEFAULT(NONE)
  
  !$OMP DO collapse(3)
  do b=1,nV
    do a=1,nV
      do j=1,nO
        do i=1,nO
          r2(i,j,a,b) = r2(i,j,a,b) - X_oovv(j,i,a,b) + X_oovv(i,j,a,b)
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  
  deallocate(cF_ov,A_oo,B_oovv,X_oovv)
  
end

! Use cF_vv

subroutine use_cF_vv(nO,nV,t1,t2,r1,r2)

  implicit none

  integer, intent(in)             :: nO,nV
  double precision, intent(in)    :: t1(nO,nV), t2(nO,nO,nV,nV)
  
  double precision, intent(inout) :: r1(nO,nV), r2(nO,nO,nV,nV)

  double precision, allocatable   :: cF_vv(:,:)
  integer                         :: i,j,a,b,e,m

  allocate(cF_vv(nV,nV))
  
  !call compute_cF_vv(nO,nV,t1,tau_t,F_ov,F_vv,v_oovv,v_ovvv,cF_vv)

  deallocate(cF_vv)
  
end

! Use cW_vvvd

subroutine use_cW_vvvf(nO,nV,t1,t2,tau,v_oovv,r2)

  implicit none

  integer, intent(in)             :: nO,nV
  double precision, intent(in)    :: t1(nO,nV), t2(nO,nO,nV,nV), tau(nO,nO,nV,nV)
  double precision, intent(in)    :: v_oovv(nO,nO,nV,nV)
  !double precision, intent(in)    :: v_vovv(nV,nO,nV,nV)
  
  double precision, intent(inout) :: r2(nO,nO,nV,nV)

  double precision, allocatable   :: cW_vvvf(:,:,:), v_vvvf(:,:,:), tau_f(:,:,:), v_vovf(:,:,:)
  integer                         :: i,j,e,f
  double precision                :: ti,tf

  allocate(cW_vvvf(nV,nV,nV),v_vvvf(nV,nV,nV),tau_f(nO,nO,nV),v_vovf(nV,nO,nV))

  !PROVIDE cc_nVab
  
  !do b=1,nV
  !  do a=1,nV
  !    do j=1,nO
  !      do i=1,nO

  !        do f=1,nV
  !          do e=1,nV
  !            r2(i,j,a,b) = r2(i,j,a,b) + 0.5d0*tau(i,j,e,f)*cW_vvvv(a,b,e,f)
  !          end do
  !        end do

  !      end do
  !    end do
  !  end do
  !end do
  
  do f = 1, nV
    call wall_time(ti)
    !$OMP PARALLEL &
    !$OMP SHARED(tau,tau_f,f,nO,nV) &
    !$OMP PRIVATE(i,j,e) &
    !$OMP DEFAULT(NONE)
    
    !$OMP DO collapse(2)
    do e = 1, nV
      do j = 1, nO
        do i = 1, nO
          tau_f(i,j,e) = tau(i,j,e,f)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    call wall_time(tf)
    if (cc_dev .and. f == 1) then
      print*,'1st transpo', tf-ti
    endif

    call wall_time(ti)
    call gen_v_spin_3idx(cc_nV_m,cc_nV_m,cc_nV_m,cc_nV_m, f, cc_nV_S,cc_nV_S,cc_nV_S,cc_nV_S, &
         cc_list_vir_spin,cc_list_vir_spin,cc_list_vir_spin,cc_list_vir_spin, &
         nV,nV,nV, v_vvvf)
    call wall_time(tf)
    if (cc_dev .and. f == 1) then
      print*,'vvvf', tf-ti
    endif
    call wall_time(ti)
    call gen_v_spin_3idx(cc_nV_m,cc_nO_m,cc_nV_m,cc_nV_m, f, cc_nV_S,cc_nO_S,cc_nV_S,cc_nV_S, &
         cc_list_vir_spin,cc_list_occ_spin,cc_list_vir_spin,cc_list_vir_spin, &
         nV,nO,nV, v_vovf)
    call wall_time(tf)
    if (cc_dev .and. f == 1) then
      print*,'vovf', tf-ti
    endif
    
    call wall_time(ti)
    call compute_cW_vvvf(nO,nV,t1,t2,tau,f,v_vvvf,v_vovf,v_oovv,cW_vvvf)
    call wall_time(tf)
    if (cc_dev .and. f == 1) then
      print*,'cW_vvvf', tf-ti
    endif

    call wall_time(ti)
    call dgemm('N','T', nO*nO, nV*nV, nV, &
               0.5d0, tau_f    , size(tau_f,1) * size(tau_f,2), &
                      cW_vvvf, size(cW_vvvf,1) * size(cW_vvvf,2), &
               1d0  , r2     , size(r2,1) * size(r2,2))
    call wall_time(tf)
    if (cc_dev .and. f == 1) then
      print*,'last dgemm', tf-ti
    endif
  enddo

  deallocate(cW_vvvf,v_vvvf,v_vovf)
  
end

! cF_oo

subroutine compute_cF_oo(nO,nV,t1,tau_t,Foo,Fov,v_ooov,v_oovv,cF_oo)

  implicit none

  integer,intent(in)            :: nO,nV
  double precision,intent(in)   :: t1(nO,nV)
  double precision,intent(in)   :: tau_t(nO,nO,nV,nV)
  double precision,intent(in)   :: Foo(nO,nO)
  double precision,intent(in)   :: Fov(nO,nV)
  double precision,intent(in)   :: v_ooov(nO,nO,nO,nV)
  double precision,intent(in)   :: v_oovv(nO,nO,nV,nV)

  double precision,intent(out)  :: cF_oo(nO,nO)

  integer                       :: i,j,m,n
  integer                       :: a,b,e,f
  double precision,external     :: Kronecker_Delta

  !$OMP PARALLEL &
  !$OMP SHARED(cF_oo,Foo,t1,v_ooov,nO,nV) &
  !$OMP PRIVATE(i,m,n,e) &
  !$OMP DEFAULT(NONE)
  
  !do i=1,nO
  !  do m=1,nO
  !    cF_oo(m,i) = (1d0 - Kronecker_delta(m,i))*Foo(m,i)
  !  end do
  !end do
  !$OMP DO collapse(1)
  do i=1,nO
    do m=1,nO
        cF_oo(m,i) = Foo(m,i)
    end do
  end do
  !$OMP END DO
  !$OMP DO
  do i = 1, nO
    cF_oo(i,i) = 0d0
  end do
  !$OMP END DO
  
  do e=1,nV
    do n=1,nO
      !$OMP DO collapse(1)
      do i=1,nO
        do m=1,nO
          cF_oo(m,i) = cF_oo(m,i) + t1(n,e)*v_ooov(m,n,i,e)
        end do
      end do
      !$OMP END DO
    end do
  end do
  !$OMP END PARALLEL

  !do i=1,nO
  !  do m=1,nO
  !    do e=1,nV
  !      cF_oo(m,i) = cF_oo(m,i) + 0.5d0*t1(i,e)*Fov(m,e)
  !    end do
  !  end do
  !end do
  call dgemm('N','T', nO, nO, nV,&
             0.5d0, Fov  , size(Fov,1), &
                    t1   , size(t1,1), &
             1d0  , cF_oo, size(cF_oo,1))

  !do i=1,nO
  !  do m=1,nO
  !    do f=1,nV
  !      do e=1,nV
  !        do n=1,nO
  !          cF_oo(m,i) = cF_oo(m,i) + 0.5d0*tau_t(i,n,e,f)*v_oovv(m,n,e,f)
  !        end do
  !      end do
  !    end do
  !  end do
  !end do
  call dgemm('N','T', nO, nO, nO*nV*nV, &
             0.5d0, v_oovv, size(v_oovv,1), &
                    tau_t , size(tau_t,1), &
             1d0  , cF_oo , size(cF_oo,1)) 
  
end

! cF_ov

subroutine compute_cF_ov(nO,nV,t1,Fov,v_oovv,cF_ov)

  implicit none

  integer,intent(in)            :: nO,nV
  double precision,intent(in)   :: t1(nO,nV)
  double precision,intent(in)   :: Fov(nO,nV),v_oovv(nO,nO,nV,nV)

  double precision,intent(out)  :: cF_ov(nO,nV)

  integer                       :: i,j,m,n
  integer                       :: a,b,e,f

  !$OMP PARALLEL &
  !$OMP SHARED(cF_ov,Fov,t1,v_oovv,nO,nV) &
  !$OMP PRIVATE(i,a,m,n,e,f) &
  !$OMP DEFAULT(NONE)
  
  !cF_ov = Fov

  !$OMP DO collapse(1)
  do e=1,nV
    do m=1,nO
      cF_ov(m,e) = Fov(m,e)
      do f=1,nV
        do n=1,nO
          cF_ov(m,e) = cF_ov(m,e) + t1(n,f)*v_oovv(m,n,e,f)
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  
end

! cF_vv

subroutine compute_cF_vv(nO,nV,t1,tau_t,Fov,Fvv,v_oovv,cF_vv)

  implicit none

  integer,intent(in)            :: nO,nV
  double precision,intent(in)   :: t1(nO,nV)
  double precision,intent(in)   :: tau_t(nO,nO,nV,nV)
  double precision,intent(in)   :: Fov(nO,nV)
  double precision,intent(in)   :: Fvv(nV,nV)
  double precision,intent(in)   :: v_oovv(nO,nO,nV,nV)
  !double precision,intent(in)   :: v_ovvv(nO,nV,nV,nV)

  double precision,intent(out)  :: cF_vv(nV,nV)
  
  double precision, allocatable :: v_ovfv(:,:,:),X_ovfv(:,:,:)
  integer                       :: i,j,m,n
  integer                       :: a,b,e,f

  !$OMP PARALLEL &
  !$OMP SHARED(cF_vv,Fvv,nO,nV) &
  !$OMP PRIVATE(e,a) &
  !$OMP DEFAULT(NONE)
  !$OMP DO collapse(1)
  do e=1,nV
    do a=1,nV
      cF_vv(a,e) = Fvv(a,e)
    end do
  end do
  !$OMP END DO
  !$OMP DO
  do e = 1, nV
    cF_vv(e,e) = 0d0
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
 
  !do e=1,nV
  !  do a=1,nV
  !    do m=1,nO
  !      cF_vv(a,e) = cF_vv(a,e) - 0.5d0*t1(m,a)*Fov(m,e)
  !    end do
  !  end do
  !end do
  call dgemm('T','N', nV, nV, nO, &
             -0.5d0, t1   , size(t1,1), &
                     Fov  , size(Fov,1), &
              1d0  , cF_vv, size(cF_vv,1))
  
  !do e=1,nV
  !  do a=1,nV
  !    do m=1,nO
  !      do f=1,nV
  !        cF_vv(a,e) = cF_vv(a,e) + t1(m,f)*v_ovvv(m,a,f,e)
  !      end do
  !    end do
  !  end do
  !end do
  allocate(v_ovfv(nO,nV,nV),X_ovfv(nO,nV,nV))
  do f = 1, nV

     call gen_v_spin_3idx_ij_l(cc_nO_m,cc_nV_m,cc_nV_m,cc_nV_m, f, cc_nO_S,cc_nV_S,cc_nV_S,cc_nV_S, &
                              cc_list_occ_spin,cc_list_vir_spin,cc_list_vir_spin,cc_list_vir_spin, &
                              nO,nV,nV, v_ovfv)

    !$OMP PARALLEL &
    !$OMP SHARED(nO,nV,v_ovfv,X_ovfv,f) &
    !$OMP PRIVATE(m,a,e) &
    !$OMP DEFAULT(NONE)
    !$OMP DO collapse(2)
    do e = 1, nV
      do a = 1, nV
        do m = 1, nO
          !X_ovfv(m,a,e) = v_ovvv(m,a,f,e)
          X_ovfv(m,a,e) = v_ovfv(m,a,e)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    call dgemv('T', nO, nV*nV, &
               !1d0, v_ovvv(:,:,f,:), size(v_ovvv,1), &
               1d0, X_ovfv, size(X_ovfv,1), &
                    t1(1,f), 1, &
               1d0, cF_vv, 1)
  enddo
  deallocate(v_ovfv,X_ovfv)

  !do e=1,nV
  !  do a=1,nV
  !    do f=1,nV
  !      do n=1,nO
  !        do m=1,nO
  !          cF_vv(a,e) = cF_vv(a,e) - 0.5d0*tau_t(m,n,a,f)*v_oovv(m,n,e,f)
  !        end do
  !      end do
  !    end do
  !  end do
  !end do
  do f = 1, nV
     call dgemm('T','N', nV, nV, nO*nO,&
                -0.5d0, tau_t(1,1,1,f) , size(tau_t,1) * size(tau_t,2), &
                        v_oovv(1,1,1,f), size(v_oovv,1) * size(v_oovv,2), &
                1d0   , cF_vv, size(cF_vv,1))
  enddo

end

! cW_oooo

subroutine compute_cW_oooo(nO,nV,t1,t2,tau,v_oooo,v_ooov,v_oovv,cW_oooo)

  implicit none

  integer,intent(in)            :: nO,nV
  double precision,intent(in)   :: t1(nO,nV)
  double precision,intent(in)   :: t2(nO,nO,nV,nV)
  double precision,intent(in)   :: tau(nO,nO,nV,nV)
  double precision,intent(in)   :: v_oooo(nO,nO,nO,nO)
  double precision,intent(in)   :: v_ooov(nO,nO,nO,nV)
  double precision,intent(in)   :: v_oovv(nO,nO,nV,nV)

  double precision,intent(out)  :: cW_oooo(nO,nO,nO,nO)

  integer                       :: i,j,m,n
  integer                       :: a,b,e,f
  double precision, allocatable :: X_oooo(:,:,:,:)

  ! oooo block  

  !cW_oooo = v_oooo

  !do j=1,nO
  !  do i=1,nO
  !    do n=1,nO
  !      do m=1,nO

  !        do e=1,nV
  !          cW_oooo(m,n,i,j) = cW_oooo(m,n,i,j) + t1(j,e)*v_ooov(m,n,i,e) - t1(i,e)*v_ooov(m,n,j,e)
  !        end do

  !      end do
  !    end do
  !  end do
  !end do
  allocate(X_oooo(nO,nO,nO,nO))
  
  call dgemm('N','T', nO*nO*nO, nO, nV, &
             1d0, v_ooov, size(v_ooov,1) * size(v_ooov,2) * size(v_ooov,3), &
                  t1    , size(t1,1), &
             0d0, X_oooo, size(X_oooo,1) * size(X_oooo,1) * size(X_oooo,3))
  !$OMP PARALLEL &
  !$OMP SHARED(cW_oooo,v_oooo,X_oooo,nO,nV) &
  !$OMP PRIVATE(i,j,m,n) &
  !$OMP DEFAULT(NONE)
  !$OMP DO collapse(3)
  do j=1,nO
    do i=1,nO
      do n=1,nO
        do m=1,nO
          cW_oooo(m,n,i,j) = v_oooo(m,n,i,j) + X_oooo(m,n,i,j) - X_oooo(m,n,j,i)
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  
  deallocate(X_oooo)
  
  !do m=1,nO
  !  do n=1,nO
  !    do i=1,nO
  !      do j=1,nO
  !         
  !        do e=1,nV
  !          do f=1,nV
  !            cW_oooo(m,n,i,j) = cW_oooo(m,n,i,j) + 0.25d0*tau(i,j,e,f)*v_oovv(m,n,e,f)
  !          end do
  !        end do

  !      end do
  !    end do
  !  end do
  !end do

  call dgemm('N','T', nO*nO, nO*nO, nV*nV, &
             0.25d0, v_oovv , size(v_oovv,1) * size(v_oovv,2), &
                     tau    , size(tau,1) * size(tau,2), &
             1.d0  , cW_oooo, size(cW_oooo,1) * size(cW_oooo,2))
  
end

! cW_ovvo

subroutine compute_cW_ovvo(nO,nV,t1,t2,tau,v_ovvo,v_oovo,v_oovv,cW_ovvo)

  implicit none

  integer,intent(in)            :: nO,nV
  double precision,intent(in)   :: t1(nO,nV)
  double precision,intent(in)   :: t2(nO,nO,nV,nV)
  double precision,intent(in)   :: tau(nO,nO,nV,nV)
  double precision,intent(in)   :: v_oovo(nO,nO,nV,nO)
  double precision,intent(in)   :: v_oovv(nO,nO,nV,nV)
  double precision,intent(in)   :: v_ovvo(nO,nV,nV,nO)
  !double precision,intent(in)   :: v_ovvv(nO,nV,nV,nV)

  double precision,intent(out)  :: cW_ovvo(nO,nV,nV,nO)

  integer                       :: i,j,m,n
  integer                       :: a,b,e,f
  double precision, allocatable :: A_oovo(:,:,:,:), B_vovo(:,:,:,:)
  double precision, allocatable :: A_voov(:,:,:,:), B_voov(:,:,:,:), C_ovov(:,:,:,:)
  double precision, allocatable :: v_ovev(:,:,:), cW_oveo(:,:,:)

  !$OMP PARALLEL &
  !$OMP SHARED(cW_ovvo,v_ovvo,nO,nV) &
  !$OMP PRIVATE(i,j,a,b) &
  !$OMP DEFAULT(NONE)
  !$OMP DO collapse(3)
  do j = 1, nO
    do b = 1, nV
      do a = 1, nV
        do i = 1, nO
          cW_ovvo(i,a,b,j) = v_ovvo(i,a,b,j)
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  !do m=1,nO
  !  do b=1,nV
  !    do e=1,nV
  !      do j=1,nO
  !        do f=1,nV
  !          cW_ovvo(m,b,e,j) = cW_ovvo(m,b,e,j) + t1(j,f)*v_ovvv(m,b,e,f)
  !        end do
  !      end do
  !    end do
  !  end do
  !end do
  allocate(v_ovev(nO,nV,nV),cW_oveo(nO,nV,nO))
  do e = 1, nV

    call gen_v_spin_3idx_ij_l(cc_nO_m,cc_nV_m,cc_nV_m,cc_nV_m, e, cc_nO_S,cc_nV_S,cc_nV_S,cc_nV_S, &
                              cc_list_occ_spin,cc_list_vir_spin,cc_list_vir_spin,cc_list_vir_spin, &
                              nO,nV,nV, v_ovev)
     
    call dgemm('N','T', nO*nV, nO, nV, &
               1.d0, v_ovev , size(v_ovev,1) * size(v_ovev,2), &
                     t1     , size(t1,1), &
               0.d0, cW_oveo, size(cW_oveo,1) * size(cW_oveo,2))
    !$OMP PARALLEL &
    !$OMP SHARED(e,cW_ovvo,cW_oveo,nO,nV) &
    !$OMP PRIVATE(m,b,j) &
    !$OMP DEFAULT(NONE)
    !$OMP DO collapse(2)
    do j = 1, nO
      do b = 1, nV
        do m = 1, nO
          cW_ovvo(m,b,e,j) = cW_ovvo(m,b,e,j) + cW_oveo(m,b,j)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
  enddo
  deallocate(v_ovev,cW_oveo)
  !call dgemm('N','T', nO*nV*nV, nO, nV, &
  !           1.d0, v_ovvv , size(v_ovvv,1) * size(v_ovvv,2) * size(v_ovvv,3), &
  !                 t1     , size(t1,1), &
  !           1.d0, cW_ovvo, size(cW_ovvo,1) * size(cW_ovvo,2) * size(cW_ovvo,3))

  !do j=1,nO
  !  do e=1,nV
  !    do b=1,nV
  !      do m=1,nO
  !        do n=1,nO
  !          cW_ovvo(m,b,e,j) = cW_ovvo(m,b,e,j) - t1(n,b)*v_oovo(m,n,e,j)
  !        end do
  !      end do
  !    end do
  !  end do
  !end do
  
  allocate(A_oovo(nO,nO,nV,nO), B_vovo(nV,nO,nV,nO))
  
  !$OMP PARALLEL &
  !$OMP SHARED(A_oovo,v_oovo,nO,nV) &
  !$OMP PRIVATE(j,e,m,n) &
  !$OMP DEFAULT(NONE)
  
  !$OMP DO collapse(3)
  do j=1,nO
    do e=1,nV
      do m=1,nO
        do n=1,nO
          A_oovo(n,m,e,j) = v_oovo(m,n,e,j)
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  
  call dgemm('T','N', nV, nO*nV*nO, nO, &
             1d0, t1    , size(t1,1), &
                  A_oovo, size(A_oovo,1), &
             0d0, B_vovo, size(B_vovo,1))
  
  !$OMP PARALLEL &
  !$OMP SHARED(cW_ovvo,B_vovo,nO,nV) &
  !$OMP PRIVATE(j,e,m,b) &
  !$OMP DEFAULT(NONE)
  
  !$OMP DO collapse(3)
  do j=1,nO
    do e=1,nV
      do b=1,nV
        do m=1,nO
          cW_ovvo(m,b,e,j) = cW_ovvo(m,b,e,j) - B_vovo(b,m,e,j)
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  deallocate(A_oovo,B_vovo)

  !do j=1,nO
  !  do e=1,nV
  !    do b=1,nV
  !      do m=1,nO
  !        do f=1,nV
  !          do n=1,nO
  !            cW_ovvo(m,b,e,j) = cW_ovvo(m,b,e,j) &
  !                            - ( 0.5d0*t2(j,n,f,b) + t1(j,f)*t1(n,b) )*v_oovv(m,n,e,f)
  !          end do
  !        end do
  !      end do
  !    end do
  !  end do
  !end do
  allocate(A_voov(nV,nO,nO,nV), B_voov(nV,nO,nO,nV), C_ovov(nO,nV,nO,nV))

  !$OMP PARALLEL &
  !$OMP SHARED(nO,nV,A_voov,B_voov,v_oovv,t2,t1) &
  !$OMP PRIVATE(f,n,m,e,j,b) &
  !$OMP DEFAULT(NONE)
  
  !$OMP DO collapse(3)
  do b = 1, nV
    do j = 1, nO
      do n = 1, nO
        do f = 1, nV
          A_voov(f,n,j,b) = 0.5d0*t2(j,n,f,b) + t1(j,f)*t1(n,b)
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP DO collapse(3)
  do e = 1, nV
    do m = 1, nO
      do n = 1, nO
        do f = 1, nV
          B_voov(f,n,m,e) = v_oovv(m,n,e,f)
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  
  call dgemm('T','N', nO*nV, nV*nO, nV*nO, &
             1d0, A_voov, size(A_voov,1) * size(A_voov,2), &
                  B_voov, size(B_voov,1) * size(B_voov,2), &
             0d0, C_ovov, size(C_ovov,1) * size(C_ovov,2))
  
  deallocate(A_voov,B_voov)

  !$OMP PARALLEL &
  !$OMP SHARED(cW_ovvo,C_ovov,nO,nV) &
  !$OMP PRIVATE(j,e,m,b) &
  !$OMP DEFAULT(NONE)
  
  !$OMP DO collapse(3)
  do j = 1, nO
    do e = 1, nV
      do b = 1, nV
        do m = 1, nO
          cW_ovvo(m,b,e,j) = cW_ovvo(m,b,e,j) - C_ovov(j,b,m,e)
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  
  deallocate(C_ovov)

end

! cW_vvvv

subroutine compute_cW_vvvv(nO,nV,t1,t2,tau,v_vvvv,v_vovv,v_oovv,cW_vvvv)
 
  implicit none

  integer,intent(in)            :: nO,nV
  double precision,intent(in)   :: t1(nO,nV)
  double precision,intent(in)   :: t2(nO,nO,nV,nV)
  double precision,intent(in)   :: tau(nO,nO,nV,nV)
  double precision,intent(in)   :: v_oovv(nO,nO,nV,nV)
  double precision,intent(in)   :: v_vovv(nV,nO,nV,nV)
  double precision,intent(in)   :: v_vvvv(nV,nV,nV,nV)

  double precision,intent(out)  :: cW_vvvv(nV,nV,nV,nV)

  integer                       :: i,j,m,n
  integer                       :: a,b,c,d,e,f
  double precision, allocatable :: A_ovvv(:,:,:,:), B_vvvv(:,:,:,:)

  allocate(A_ovvv(nO,nV,nV,nV), B_vvvv(nV,nV,nV,nV))

  !$OMP PARALLEL &
  !$OMP SHARED(nO,nV,cW_vvvv,A_ovvv,v_vovv,v_vvvv) &
  !$OMP PRIVATE(a,b,c,d,e,f,m) &
  !$OMP DEFAULT(NONE)

  !$OMP DO collapse(3)
  do d = 1, nV
    do c = 1, nV
      do b = 1, nV
        do a = 1, nV
          cW_vvvv(a,b,c,d) = v_vvvv(a,b,c,d)
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO NOWAIT

  !do f=1,nV
  !  do e=1,nV
  !    do b=1,nV
  !      do a=1,nV
  !        do m=1,nO
  !          cW_vvvv(a,b,e,f) = cW_vvvv(a,b,e,f) - t1(m,b)*v_vovv(a,m,e,f) + t1(m,a)*v_vovv(b,m,e,f)
  !        end do
  !      end do
  !    end do
  !  end do
  !end do
  !$OMP DO collapse(3)
  do f=1,nV
    do e=1,nV
      do a=1,nV
        do m=1,nO
          A_ovvv(m,a,e,f) = v_vovv(a,m,e,f)
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  call dgemm('T','N', nV, nV*nV*nV, nO, &
             1d0, t1    , size(t1,1), &
                  A_ovvv, size(A_ovvv,1), &
             0d0, B_vvvv, size(B_vvvv,1))

  !$OMP PARALLEL &
  !$OMP SHARED(nO,nV,cW_vvvv,B_vvvv) &
  !$OMP PRIVATE(a,b,c,d,e,f,m) &
  !$OMP DEFAULT(NONE)

  !$OMP DO collapse(3)
  do f=1,nV
    do e=1,nV
      do b=1,nV
        do a=1,nV
          cW_vvvv(a,b,e,f) = cW_vvvv(a,b,e,f) - B_vvvv(b,a,e,f) + B_vvvv(a,b,e,f)
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  
  deallocate(A_ovvv,B_vvvv)

  !do a=1,nV
  !  do b=1,nV
  !    do e=1,nV
  !      do f=1,nV
  !         
  !        do m=1,nO
  !          do n=1,nO
  !            cW_vvvv(a,b,e,f) = cW_vvvv(a,b,e,f) + 0.25d0*tau(m,n,a,b)*v_oovv(m,n,e,f)
  !          end do
  !        end do

  !      end do
  !    end do
  !  end do
  !end do
  call dgemm('T','N', nV*nV, nV*nV, nO*nO, &
             0.25d0, tau    , size(tau,1) * size(tau,2), &
                     v_oovv , size(v_oovv,1) * size(v_oovv,2), &
             1.d0  , cW_vvvv, size(cW_vvvv,1) * size(cW_vvvv,2))

end

! cW_vvvf

subroutine compute_cW_vvvf(nO,nV,t1,t2,tau,f,v_vvvf,v_vovf,v_oovv,cW_vvvf)
 
  implicit none

  integer,intent(in)            :: nO,nV,f
  double precision,intent(in)   :: t1(nO,nV)
  double precision,intent(in)   :: t2(nO,nO,nV,nV)
  double precision,intent(in)   :: tau(nO,nO,nV,nV)
  double precision,intent(in)   :: v_oovv(nO,nO,nV,nV)
  double precision,intent(in)   :: v_vovf(nV,nO,nV)
  double precision,intent(in)   :: v_vvvf(nV,nV,nV)

  double precision,intent(out)  :: cW_vvvf(nV,nV,nV)

  integer                       :: i,j,m,n
  integer                       :: a,b,c,d,e
  double precision, allocatable :: A_ovvf(:,:,:), B_vvvf(:,:,:), v_oovf(:,:,:)
  double precision :: ti,tf

  allocate(A_ovvf(nO,nV,nV), B_vvvf(nV,nV,nV))
  allocate(v_oovf(nO,nO,nV))

  !$OMP PARALLEL &
  !$OMP SHARED(nO,nV,cW_vvvf,A_ovvf,v_vovf,v_vvvf,f) &
  !$OMP PRIVATE(a,b,c,d,e,m) &
  !$OMP DEFAULT(NONE)
  
  !$OMP DO collapse(2)
  do c = 1, nV
    do b = 1, nV
      do a = 1, nV
        !cW_vvvv(a,b,c,d) = v_vvvv(a,b,c,d)
        cW_vvvf(a,b,c) = v_vvvf(a,b,c)
      enddo
    enddo
  enddo
  !$OMP END DO NOWAIT

  !do f=1,nV
  !  do e=1,nV
  !    do b=1,nV
  !      do a=1,nV
  !        do m=1,nO
  !          cW_vvvv(a,b,e,f) = cW_vvvv(a,b,e,f) - t1(m,b)*v_vovv(a,m,e,f) + t1(m,a)*v_vovv(b,m,e,f)
  !        end do
  !      end do
  !    end do
  !  end do
  !end do

  !$OMP DO collapse(2)
  do e=1,nV
    do a=1,nV
      do m=1,nO
        !A_ovvv(m,a,e,f) = v_vovv(a,m,e,f)
        !A_ovvf(m,a,e) = v_vovv(a,m,e,f)
        A_ovvf(m,a,e) = v_vovf(a,m,e)
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  call dgemm('T','N', nV, nV*nV, nO, &
             1d0, t1    , size(t1,1), &
                  A_ovvf, size(A_ovvf,1), &
             0d0, B_vvvf, size(B_vvvf,1))
  
  !$OMP PARALLEL &
  !$OMP SHARED(nO,nV,cW_vvvf,B_vvvf,v_oovf,v_oovv,f) &
  !$OMP PRIVATE(a,b,c,d,e,m,n) &
  !$OMP DEFAULT(NONE)

  !$OMP DO collapse(3)
  do e=1,nV
    do b=1,nV
      do a=1,nV
        !cW_vvvv(a,b,e,f) = cW_vvvv(a,b,e,f) - B_vvvv(b,a,e,f) + B_vvvv(a,b,e,f)
        cW_vvvf(a,b,e) = cW_vvvf(a,b,e) - B_vvvf(b,a,e) + B_vvvf(a,b,e)
      end do
    end do
  end do
  !$OMP END DO NOWAIT
  
  !deallocate(A_ovvf,B_vvvf)

  !do a=1,nV
  !  do b=1,nV
  !    do e=1,nV
  !      do f=1,nV
  !         
  !        do m=1,nO
  !          do n=1,nO
  !            cW_vvvv(a,b,e,f) = cW_vvvv(a,b,e,f) + 0.25d0*tau(m,n,a,b)*v_oovv(m,n,e,f)
  !          end do
  !        end do

  !      end do
  !    end do
  !  end do
  !end do

  !$OMP DO collapse(2)
  do e = 1, nV
    do n = 1, nO
      do m = 1, nO
        v_oovf(m,n,e) = v_oovv(m,n,e,f)
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL  
  
  call dgemm('T','N', nV*nV, nV, nO*nO, &
             0.25d0, tau    , size(tau,1) * size(tau,2), &
                     v_oovf , size(v_oovf,1) * size(v_oovf,2), &
             1.d0  , cW_vvvf, size(cW_vvvf,1) * size(cW_vvvf,2))
  
  deallocate(v_oovf)
  deallocate(A_ovvf,B_vvvf)

end
