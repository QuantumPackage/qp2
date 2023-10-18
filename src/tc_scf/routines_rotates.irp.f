

! ---

subroutine LTxSxR(n, m, L, S, R, C)

  implicit none
  integer,          intent(in)  :: n, m
  double precision, intent(in)  :: L(n,m), S(n,n), R(n,m)
  double precision, intent(out) :: C(m,m)
  integer                       :: i, j
  double precision              :: accu_d, accu_nd
  double precision, allocatable :: tmp(:,:)

  ! L.T x S x R
  allocate(tmp(m,n))
  call dgemm( 'T', 'N', m, n, n, 1.d0      &
            , L, size(L, 1), S, size(S, 1) &
            , 0.d0, tmp, size(tmp, 1) )
  call dgemm( 'N', 'N', m, m, n, 1.d0          &
            , tmp, size(tmp, 1), R, size(R, 1) &
            , 0.d0, C, size(C, 1) )
  deallocate(tmp)

  accu_d  = 0.d0
  accu_nd = 0.d0
  do i = 1, m
    do j = 1, m
      if(j.eq.i) then
        accu_d += dabs(C(j,i))
      else
        accu_nd += C(j,i) * C(j,i)
      endif
    enddo
  enddo
  accu_nd = dsqrt(accu_nd)

  print*, ' accu_d  = ', accu_d
  print*, ' accu_nd = ', accu_nd

end subroutine LTxR

! ---


! ---

subroutine minimize_tc_orb_angles()

  BEGIN_DOC
  ! routine that minimizes the angle between left- and right-orbitals when degeneracies are found
  END_DOC

  implicit none
  logical          :: good_angles
  integer          :: i
  double precision :: thr_deg

  good_angles = .False.
  thr_deg = thr_degen_tc

  call print_energy_and_mos(good_angles)

  print *, ' Minimizing the angles between the TC orbitals'
  i = 1
  do while (.not. good_angles)
    print *, ' iteration = ', i
    call routine_save_rotated_mos(thr_deg, good_angles)
    thr_deg *= 10.d0
    i += 1
    if(i .gt. 100) then
      print *, ' minimize_tc_orb_angles does not seem to converge ..'
      print *, ' Something is weird in the tc orbitals ...'
      print *, ' STOPPING'
      stop
    endif
  enddo
  print *, ' Converged ANGLES MINIMIZATION !!'

  call print_angles_tc()
  call print_energy_and_mos(good_angles)

end

! ---

subroutine routine_save_rotated_mos(thr_deg, good_angles)

  implicit none

  double precision, intent(in)  :: thr_deg
  logical,          intent(out) :: good_angles

  integer                       :: i, j, k, n_degen_list, m, n, n_degen, ilast, ifirst
  double precision              :: max_angle, norm
  double precision              :: E_old, E_new, E_thr
  integer,          allocatable :: list_degen(:,:)
  double precision, allocatable :: new_angles(:)
  double precision, allocatable :: mo_r_coef_old(:,:), mo_l_coef_old(:,:)
  double precision, allocatable :: mo_r_coef_good(:,:), mo_l_coef_good(:,:)
  double precision, allocatable :: mo_r_coef_new(:,:)
  double precision, allocatable :: fock_diag(:), s_mat(:,:)
  double precision, allocatable :: stmp(:,:), T(:,:), Snew(:,:), smat2(:,:)
  double precision, allocatable :: mo_l_coef_tmp(:,:), mo_r_coef_tmp(:,:), mo_l_coef_new(:,:)

  E_thr = 1d-8
  E_old = TC_HF_energy
  allocate(mo_l_coef_old(ao_num,mo_num), mo_r_coef_old(ao_num,mo_num))
  mo_r_coef_old = mo_r_coef
  mo_l_coef_old = mo_l_coef

  good_angles = .False.

  allocate(mo_l_coef_good(ao_num, mo_num), mo_r_coef_good(ao_num,mo_num))

  print *, ' ***************************************'
  print *, ' ***************************************'
  print *, ' THRESHOLD FOR DEGENERACIES ::: ', thr_deg
  print *, ' ***************************************'
  print *, ' ***************************************'
  print *, ' Starting with the following TC energy gradient :', grad_non_hermit

  mo_r_coef_good = mo_r_coef
  mo_l_coef_good = mo_l_coef

  allocate(mo_r_coef_new(ao_num, mo_num))
  mo_r_coef_new = mo_r_coef
  do i = 1, mo_num
    norm = 1.d0/dsqrt(overlap_mo_r(i,i))
    do j = 1, ao_num
      mo_r_coef_new(j,i) *= norm
    enddo
  enddo

  allocate(list_degen(mo_num,0:mo_num), s_mat(mo_num,mo_num), fock_diag(mo_num))
  do i = 1, mo_num
    fock_diag(i) = Fock_matrix_tc_mo_tot(i,i)
  enddo

 ! compute the overlap between the left and rescaled right
  call build_s_matrix(ao_num, mo_num, mo_r_coef_new, mo_r_coef_new, ao_overlap, s_mat)
! call give_degen(fock_diag,mo_num,thr_deg,list_degen,n_degen_list)
  if(n_core_orb.ne.0)then
   call give_degen_full_listcore(fock_diag, mo_num, list_core, n_core_orb, thr_deg, list_degen, n_degen_list)
  else
   call give_degen_full_list(fock_diag, mo_num, thr_deg, list_degen, n_degen_list)
  endif
  print *, ' fock_matrix_mo'
  do i = 1, mo_num
    print *, i, fock_diag(i), angle_left_right(i)
  enddo
   
  do i = 1, n_degen_list
!  ifirst = list_degen(1,i)
!  ilast  = list_degen(2,i)
!  n_degen = ilast - ifirst +1

    n_degen = list_degen(i,0)
    if(n_degen .ge. 1000)n_degen = 1 ! convention for core orbitals
     
    if(n_degen .eq. 1) cycle

    allocate(stmp(n_degen,n_degen), smat2(n_degen,n_degen))
    allocate(mo_r_coef_tmp(ao_num,n_degen), mo_l_coef_tmp(ao_num,n_degen), mo_l_coef_new(ao_num,n_degen))
    allocate(T(n_degen,n_degen), Snew(n_degen,n_degen))

    do j = 1, n_degen
      mo_r_coef_tmp(1:ao_num,j) = mo_r_coef_new(1:ao_num,list_degen(i,j))
      mo_l_coef_tmp(1:ao_num,j) = mo_l_coef(1:ao_num,list_degen(i,j))
    enddo
    ! Orthogonalization of right functions
    print *, ' Orthogonalization of RIGHT functions'
    print *, ' ------------------------------------'
    call orthog_functions(ao_num, n_degen, mo_r_coef_tmp, ao_overlap)
  
    ! Orthogonalization of left functions
    print *, ' Orthogonalization of LEFT functions'
    print *, ' ------------------------------------'
    call orthog_functions(ao_num, n_degen, mo_l_coef_tmp, ao_overlap)

    print *, ' Overlap left-right '
    call build_s_matrix(ao_num, n_degen, mo_r_coef_tmp, mo_l_coef_tmp, ao_overlap, stmp)
    do j = 1, n_degen
     write(*,'(100(F8.4,X))') stmp(:,j)
    enddo
    call build_s_matrix(ao_num, n_degen, mo_l_coef_tmp, mo_l_coef_tmp, ao_overlap, stmp)

    !print*,'LEFT/LEFT OVERLAP '
    !do j = 1, n_degen
    ! write(*,'(100(F16.10,X))')stmp(:,j)
    !enddo
    call build_s_matrix(ao_num, n_degen, mo_r_coef_tmp, mo_r_coef_tmp, ao_overlap, stmp)
    !print*,'RIGHT/RIGHT OVERLAP '
    !do j = 1, n_degen
    ! write(*,'(100(F16.10,X))')stmp(:,j)
    !enddo

    if(maxovl_tc) then
      T    = 0.d0
      Snew = 0.d0
      call maxovl(n_degen, n_degen, stmp, T, Snew)
    !print*,'overlap after'
    !do j = 1, n_degen
    ! write(*,'(100(F16.10,X))')Snew(:,j)
    !enddo
      call dgemm( 'N', 'N', ao_num, n_degen, n_degen, 1.d0                  &
                , mo_l_coef_tmp, size(mo_l_coef_tmp, 1), T(1,1), size(T, 1) &
                , 0.d0, mo_l_coef_new, size(mo_l_coef_new, 1) )
     call build_s_matrix(ao_num, n_degen, mo_l_coef_new, mo_r_coef_tmp, ao_overlap, stmp)
    !print*,'Overlap test'
    !do j = 1, n_degen
    ! write(*,'(100(F16.10,X))')stmp(:,j)
    !enddo
    else 
      mo_l_coef_new = mo_l_coef_tmp
    endif

    call impose_weighted_biorthog_svd(ao_num, n_degen, ao_overlap, mo_l_coef_new, mo_r_coef_tmp)

    !call build_s_matrix(ao_num, n_degen, mo_l_coef_new, mo_r_coef_tmp, ao_overlap, stmp)
    !print*,'LAST OVERLAP '
    !do j = 1, n_degen
    ! write(*,'(100(F16.10,X))')stmp(:,j)
    !enddo
    !call build_s_matrix(ao_num, n_degen, mo_l_coef_new, mo_l_coef_new, ao_overlap, stmp)
    !print*,'LEFT OVERLAP '
    !do j = 1, n_degen
    ! write(*,'(100(F16.10,X))')stmp(:,j)
    !enddo
    !call build_s_matrix(ao_num, n_degen, mo_r_coef_tmp, mo_r_coef_tmp, ao_overlap, stmp)
    !print*,'RIGHT OVERLAP '
    !do j = 1, n_degen
    ! write(*,'(100(F16.10,X))')stmp(:,j)
    !enddo
    do j = 1, n_degen
!!!   mo_l_coef_good(1:ao_num,j+ifirst-1) = mo_l_coef_new(1:ao_num,j)
!!!   mo_r_coef_good(1:ao_num,j+ifirst-1) = mo_r_coef_tmp(1:ao_num,j)
      mo_l_coef_good(1:ao_num,list_degen(i,j)) = mo_l_coef_new(1:ao_num,j)
      mo_r_coef_good(1:ao_num,list_degen(i,j)) = mo_r_coef_tmp(1:ao_num,j)
    enddo

    deallocate(stmp, smat2)
    deallocate(mo_r_coef_tmp, mo_l_coef_tmp, mo_l_coef_new)
    deallocate(T, Snew)
  enddo

  !allocate(stmp(mo_num, mo_num))
  !call build_s_matrix(ao_num, mo_num, mo_l_coef_good, mo_r_coef_good, ao_overlap, stmp)
  !print*,'LEFT/RIGHT OVERLAP '
  !do j = 1, mo_num
  ! write(*,'(100(F16.10,X))')stmp(:,j)
  !enddo
  !call build_s_matrix(ao_num, mo_num, mo_l_coef_good, mo_l_coef_good, ao_overlap, stmp)
  !print*,'LEFT/LEFT OVERLAP '
  !do j = 1, mo_num
  ! write(*,'(100(F16.10,X))')stmp(:,j)
  !enddo
  !call build_s_matrix(ao_num, mo_num, mo_r_coef_good, mo_r_coef_good, ao_overlap, stmp)
  !print*,'RIGHT/RIGHT OVERLAP '
  !do j = 1, mo_num
  ! write(*,'(100(F16.10,X))')stmp(:,j)
  !enddo

  mo_r_coef = mo_r_coef_good
  mo_l_coef = mo_l_coef_good
  call ezfio_set_bi_ortho_mos_mo_l_coef(mo_l_coef)
  call ezfio_set_bi_ortho_mos_mo_r_coef(mo_r_coef)
  TOUCH mo_l_coef mo_r_coef

  ! check if TC energy has changed
  E_new = TC_HF_energy
  if(dabs(E_new - E_old) .gt. E_thr) then
    mo_r_coef = mo_r_coef_old
    mo_l_coef = mo_l_coef_old
    deallocate(mo_l_coef_old, mo_r_coef_old)
    call ezfio_set_bi_ortho_mos_mo_l_coef(mo_l_coef)
    call ezfio_set_bi_ortho_mos_mo_r_coef(mo_r_coef)
    TOUCH mo_l_coef mo_r_coef
    print*, ' TC energy bef rotation = ', E_old
    print*, ' TC energy aft rotation = ', E_new
    print*, ' the rotation is refused'
    stop
  endif

  allocate(new_angles(mo_num))
  new_angles(1:mo_num) = dabs(angle_left_right(1:mo_num))
  max_angle = maxval(new_angles)
  good_angles = max_angle.lt.thresh_lr_angle
  print *, ' max_angle = ', max_angle
  deallocate(new_angles)


  deallocate(mo_l_coef_old, mo_r_coef_old)
  deallocate(mo_l_coef_good, mo_r_coef_good)
  deallocate(mo_r_coef_new)
  
end

! ---

subroutine build_s_matrix(m, n, C1, C2, overlap, smat)

  implicit none
  integer,          intent(in)  :: m, n
  double precision, intent(in)  :: C1(m,n), C2(m,n), overlap(m,m)
  double precision, intent(out) :: smat(n,n)
  integer                       :: i, j, k, l
  double precision, allocatable :: S_tmp(:,:)

  smat = 0.d0

  !do i = 1, n
  !  do j = 1, n
  !    do k = 1, m
  !      do l = 1, m
  !        smat(i,j) += C1(k,i) * overlap(l,k) * C2(l,j) 
  !      enddo
  !    enddo
  !  enddo
  !enddo

  ! C1.T x overlap
  allocate(S_tmp(n,m))
  call dgemm( 'T', 'N', n, m, m, 1.d0                    &
            , C1, size(C1, 1), overlap, size(overlap, 1) &
            , 0.d0, S_tmp, size(S_tmp, 1) )
  ! C1.T x overlap x C2
  call dgemm( 'N', 'N', n, n, m, 1.d0                     &
            , S_tmp, size(S_tmp, 1), C2(1,1), size(C2, 1) &
            , 0.d0, smat, size(smat, 1) )
  deallocate(S_tmp)

end

! ---

subroutine orthog_functions(m, n, coef, overlap)

  implicit none

  integer,          intent(in)    :: m, n
  double precision, intent(in)    :: overlap(m,m)
  double precision, intent(inout) :: coef(m,n)
  double precision, allocatable   :: stmp(:,:)
  integer                         :: j, k

  allocate(stmp(n,n))
  call build_s_matrix(m, n, coef, coef, overlap, stmp)
! print*,'overlap before'
! do j = 1, n
!  write(*,'(100(F16.10,X))')stmp(:,j)
! enddo
  call impose_orthog_svd_overlap(m, n, coef, overlap)
  call build_s_matrix(m, n, coef, coef, overlap, stmp)
  do j = 1, n
    ! ---
    ! TODO: MANU check ici
    !coef(1,:m) *= 1.d0/dsqrt(stmp(j,j))
    do k = 1, m
      coef(k,j) *= 1.d0/dsqrt(stmp(j,j))
    enddo
    ! ---
  enddo
  call build_s_matrix(m, n, coef, coef, overlap, stmp)

 !print*,'overlap after'
 !do j = 1, n
 ! write(*,'(100(F16.10,X))')stmp(:,j)
 !enddo

 deallocate(stmp)

end

! ---

subroutine print_angles_tc()

  implicit none
  integer          :: i, j
  double precision :: left, right

  print *, ' product of norms, angle between vectors'                                                                  
  do i = 1, mo_num
    left  = overlap_mo_l(i,i)
    right = overlap_mo_r(i,i)
!  print*,Fock_matrix_tc_mo_tot(i,i),left*right,angle_left_right(i)
    print *, left*right, angle_left_right(i)
  enddo

end

! ---

subroutine print_energy_and_mos(good_angles)

  implicit none
  logical, intent(out) :: good_angles
  integer              :: i

  print *, ' '
  print *, ' TC energy              = ', TC_HF_energy
  print *, ' TC SCF energy gradient = ', grad_non_hermit
  print *, ' Max angle Left/right   = ', max_angle_left_right
  call print_angles_tc()

  if(max_angle_left_right .lt. thresh_lr_angle) then
    print *, ' Maximum angle BELOW 45 degrees, everthing is OK !'
    good_angles = .true.
  else if(max_angle_left_right .gt. thresh_lr_angle .and. max_angle_left_right .lt. 75.d0) then
    print *, ' Maximum angle between thresh_lr_angle and 75 degrees, this is not the best for TC-CI calculations ...'
    good_angles = .false.
  else if(max_angle_left_right .gt. 75.d0) then
    print *, ' Maximum angle between ABOVE 75 degrees, YOU WILL CERTAINLY FIND TROUBLES IN TC-CI calculations ...'
    good_angles = .false.
  endif

  print *, ' Diag Fock elem, product of left/right norm, angle left/right '
  do i = 1, mo_num
    write(*, '(I3,X,100(F16.10,X))') i, Fock_matrix_tc_mo_tot(i,i), overlap_mo_l(i,i)*overlap_mo_r(i,i), angle_left_right(i)
  enddo

end

! ---

subroutine sort_by_tc_fock
 implicit none 
 integer, allocatable :: iorder(:)
 double precision, allocatable :: mo_l_tmp(:,:), mo_r_tmp(:,:),fock(:)
 allocate(iorder(mo_num),fock(mo_num),mo_l_tmp(ao_num, mo_num),mo_r_tmp(ao_num,mo_num))
 integer :: i
 mo_l_tmp = mo_l_coef
 mo_r_tmp = mo_r_coef
 do i = 1, mo_num
  iorder(i) = i
  fock(i) = Fock_matrix_tc_mo_tot(i,i)
 enddo
 call dsort(fock,iorder,mo_num)
 do i = 1, mo_num
  mo_l_coef(1:ao_num,i) = mo_l_tmp(1:ao_num,iorder(i))
  mo_r_coef(1:ao_num,i) = mo_r_tmp(1:ao_num,iorder(i))
 enddo
 touch mo_l_coef mo_r_coef

end

