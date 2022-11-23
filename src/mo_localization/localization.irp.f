program localization
 call run_localization
end




! Variables:
! | pre_rot(mo_num, mo_num)   | double precision | Matrix for the pre rotation                       |
! | R(mo_num,mo_num)          | double precision | Rotation matrix                                   |
! | tmp_R(:,:)                | double precision | Rottation matrix in a subsapce                    |
! | prev_mos(ao_num, mo_num)  | double precision | Previous mo_coef                                  |
! | spatial_extent(mo_num)    | double precision | Spatial extent of the orbitals                    |
! | criterion                 | double precision | Localization criterion                            |
! | prev_criterion            | double precision | Previous criterion                                |
! | criterion_model           | double precision | Estimated next criterion                          |
! | rho                       | double precision | Ratio to measure the agreement between the model  |
! |                           |                  | and the reality                                   |
! | delta                     | double precision | Radisu of the trust region                        |
! | norm_grad                 | double precision | Norm of the gradient                              |
! | info                      | integer          | for dsyev from Lapack                             |
! | max_elem                  | double precision | maximal element in the gradient                   |
! | v_grad(:)                 | double precision | Gradient                                          |
! | H(:,:)                    | double precision | Hessian (diagonal)                                |
! | e_val(:)                  | double precision | Eigenvalues of the hessian                        |
! | W(:,:)                    | double precision | Eigenvectors of the hessian                       |
! | tmp_x(:)                  | double precision | Step in 1D (in a subaspace)                       |
! | tmp_m_x(:,:)              | double precision | Step in 2D (in a subaspace)                       |
! | tmp_list(:)               | double precision | List of MOs in a mo_class                         |
! | i,j,k                     | integer          | Indexes in the full MO space                      |
! | tmp_i, tmp_j, tmp_k       | integer          | Indexes in a subspace                             |
! | l                         | integer          | Index for the mo_class                            |
! | key(:)                    | integer          | Key to sort the eigenvalues of the hessian        |
! | nb_iter                   | integer          | Number of iterations                              |
! | must_exit                 | logical          | To exit the trust region loop                     |
! | cancel_step               | logical          | To cancel a step                                  |
! | not_*converged            | logical          | To localize the different mo classes              |
! | t*                        | double precision | To measure the time                               |
! | n                         | integer          | mo_num*(mo_num-1)/2, number of orbital parameters |
! | tmp_n                     | integer          | dim_subspace*(dim_subspace-1)/2                   |
! |                           |                  | Number of dimension in the subspace               |

! Variables in qp_edit for the localization:
! | localization_method      |
! | localization_max_nb_iter |
! | default_mo_class         |
! | thresh_loc_max_elem_grad |
! | kick_in_mos              |
! | angle_pre_rot            |

! + all the variables for the trust region

! Cf. qp_edit orbital optimization


subroutine run_localization

  include 'pi.h'

  BEGIN_DOC
  ! Orbital localization
  END_DOC

  implicit none

  ! Variables
  double precision, allocatable :: pre_rot(:,:), R(:,:)
  double precision, allocatable :: prev_mos(:,:), spatial_extent(:), tmp_R(:,:)
  double precision              :: criterion, norm_grad
  integer                       :: i,j,k,l,p, tmp_i, tmp_j, tmp_k
  integer                       :: info
  integer                       :: n, tmp_n, tmp_list_size
  double precision, allocatable :: v_grad(:), H(:,:), tmp_m_x(:,:), tmp_x(:),W(:,:),e_val(:)
  double precision              :: max_elem, t1, t2, t3, t4, t5, t6
  integer, allocatable          :: tmp_list(:), key(:)
  double precision              :: prev_criterion, rho, delta, criterion_model
  integer                       :: nb_iter, nb_sub_iter
  logical                       :: not_converged, not_core_converged
  logical                       :: not_act_converged, not_inact_converged, not_virt_converged
  logical                       :: use_trust_region, must_exit, cancel_step,enforce_step_cancellation

  n = mo_num*(mo_num-1)/2

  ! Allocation
  allocate(spatial_extent(mo_num))
  allocate(pre_rot(mo_num, mo_num), R(mo_num, mo_num))
  allocate(prev_mos(ao_num, mo_num))

  ! Locality before the localization
  call compute_spatial_extent(spatial_extent)

  ! Choice of the method (with qp_edit)
  print*,''
  print*,'Localization method:',localization_method
  if (localization_method == 'boys') then
    print*,'Foster-Boys localization'
  elseif (localization_method == 'pipek') then
    print*,'Pipek-Mezey localization'
  else
    print*,'Unknown localization_method, please select boys or pipek'
    call abort
  endif
  print*,''

  ! Localization criterion (FB, PM, ...) for each mo_class
  print*,'### Before the pre rotation'
  
  ! Debug
  if (debug_hf) then
    print*,'HF energy:', HF_energy
  endif
  
  do l = 1, 4
    if (l==1) then ! core
      tmp_list_size = dim_list_core_orb
    elseif (l==2) then ! act
      tmp_list_size = dim_list_act_orb
    elseif (l==3) then ! inact
      tmp_list_size = dim_list_inact_orb
    else ! virt
      tmp_list_size = dim_list_virt_orb
    endif

     ! Allocation tmp array  
    allocate(tmp_list(tmp_list_size))

    ! To give the list of MOs in a mo_class
    if (l==1) then ! core
      tmp_list = list_core
    elseif (l==2) then
      tmp_list = list_act
    elseif (l==3) then
      tmp_list = list_inact
    else
      tmp_list = list_virt
    endif

    if (tmp_list_size >= 2) then
      call criterion_localization(tmp_list_size, tmp_list,criterion)
      print*,'Criterion:', criterion, mo_class(tmp_list(1))
    endif

    deallocate(tmp_list)    

  enddo

  ! Debug
  !print*,'HF', HF_energy

  print*, 'Security mo_class:', security_mo_class

  ! The default mo_classes are setted only if the MOs to localize are not specified
  if (security_mo_class .and. (n_act_orb == mo_num .or. &
      n_core_orb + n_act_orb == mo_num)) then

    print*, 'WARNING'
    print*, 'You must set different mo_class with qp set_mo_class'
    print*, 'If you want to localize all the orbitals:'
    print*, 'qp set Orbital_optimization security_mo_class false'
    print*, ''
    print*, 'abort'

    call abort
  
  endif

! Loc

  ! Pre rotation, to give a little kick in the MOs
  call apply_pre_rotation()

  ! Criterion after the pre rotation
  ! Localization criterion (FB, PM, ...) for each mo_class
  print*,'### After the pre rotation'
  
  ! Debug
  if (debug_hf) then
    touch mo_coef
    print*,'HF energy:', HF_energy
  endif
  
  do l = 1, 4
    if (l==1) then ! core
      tmp_list_size = dim_list_core_orb
    elseif (l==2) then ! act
      tmp_list_size = dim_list_act_orb
    elseif (l==3) then ! inact
      tmp_list_size = dim_list_inact_orb
    else ! virt
      tmp_list_size = dim_list_virt_orb
    endif
    
    if (tmp_list_size >= 2) then
      ! Allocation tmp array  
      allocate(tmp_list(tmp_list_size))

      ! To give the list of MOs in a mo_class
      if (l==1) then ! core
        tmp_list = list_core
      elseif (l==2) then
        tmp_list = list_act
      elseif (l==3) then
        tmp_list = list_inact
      else
        tmp_list = list_virt
      endif

      call criterion_localization(tmp_list_size, tmp_list,criterion)
      print*,'Criterion:', criterion, trim(mo_class(tmp_list(1)))
      
      deallocate(tmp_list)    
    endif
    
  enddo

  ! Debug
  !print*,'HF', HF_energy

  print*,''
  print*,'========================'
  print*,'  Orbital localization'
  print*,'========================'
  print*,'' 

  !Initialization
  not_converged = .TRUE.

  ! To do the localization only if there is at least 2 MOs
  if (dim_list_core_orb >= 2) then
    not_core_converged = .TRUE.
  else
    not_core_converged = .FALSE. 
  endif
  
  if (dim_list_act_orb >= 2) then
    not_act_converged = .TRUE.
  else
    not_act_converged = .FALSE.
  endif

  if (dim_list_inact_orb >= 2) then
    not_inact_converged = .TRUE.
  else
    not_inact_converged = .FALSE.
  endif

  if (dim_list_virt_orb >= 2) then
    not_virt_converged = .TRUE.
  else
    not_virt_converged = .FALSE.
  endif
 
  ! Loop over the mo_classes
  do l = 1, 4

    if (l==1) then ! core
      not_converged = not_core_converged
      tmp_list_size = dim_list_core_orb
    elseif (l==2) then ! act
      not_converged = not_act_converged
      tmp_list_size = dim_list_act_orb
    elseif (l==3) then ! inact
      not_converged = not_inact_converged
      tmp_list_size = dim_list_inact_orb
    else ! virt
      not_converged = not_virt_converged
      tmp_list_size = dim_list_virt_orb
    endif

    ! Next iteration if converged = true 
    if (.not. not_converged) then
      cycle
    endif
 
    ! Allocation tmp array  
    allocate(tmp_list(tmp_list_size))

    ! To give the list of MOs in a mo_class
    if (l==1) then ! core
      tmp_list = list_core
    elseif (l==2) then
      tmp_list = list_act
    elseif (l==3) then
      tmp_list = list_inact
    else
      tmp_list = list_virt
    endif

    ! Display
    if (not_converged) then
      print*,''
      print*,'###', trim(mo_class(tmp_list(1))), 'MOs ###'
      print*,''
    endif

    ! Size for the 2D -> 1D transformation 
    tmp_n = tmp_list_size * (tmp_list_size - 1)/2

    ! Without hessian + trust region 
    if (.not. localization_use_hessian) then
       
      ! Allocation of temporary arrays
      allocate(v_grad(tmp_n), tmp_m_x(tmp_list_size, tmp_list_size))
      allocate(tmp_R(tmp_list_size, tmp_list_size), tmp_x(tmp_n))

      ! Criterion
      call criterion_localization(tmp_list_size, tmp_list, prev_criterion)

      ! Init
      nb_iter = 0
      delta = 1d0

      !Loop
      do while (not_converged)
         
        print*,''
        print*,'***********************'
        print*,'Iteration', nb_iter
        print*,'***********************'
        print*,''
        
        ! Angles of rotation
        call theta_localization(tmp_list, tmp_list_size, tmp_m_x, max_elem)
        tmp_m_x = - tmp_m_x * delta

        ! Rotation submatrix
        call rotation_matrix(tmp_m_x, tmp_list_size, tmp_R, tmp_list_size, tmp_list_size, &
             info, enforce_step_cancellation)

        ! To ensure that the rotation matrix is unitary
        if (enforce_step_cancellation) then
          print*, 'Step cancellation, too large error in the rotation matrix'
          delta = delta * 0.5d0
          cycle
        else
          delta = min(delta * 2d0, 1d0)
        endif

        ! Full rotation matrix and application of the rotation
        call sub_to_full_rotation_matrix(tmp_list_size, tmp_list, tmp_R, R)
        call apply_mo_rotation(R, prev_mos)

        ! Update the needed data
        call update_data_localization()

        ! New criterion
        call criterion_localization(tmp_list_size, tmp_list, criterion)
        print*,'Criterion:', trim(mo_class(tmp_list(1))), nb_iter, criterion
        print*,'Max elem :', max_elem
        print*,'Delta    :', delta
        
        nb_iter = nb_iter + 1

        ! Exit
        if (nb_iter >= localization_max_nb_iter .or. dabs(max_elem) < thresh_loc_max_elem_grad) then
           not_converged = .False.
        endif
      enddo

      ! Save the changes
      call update_data_localization()
      call save_mos()
      TOUCH mo_coef

      ! Deallocate
      deallocate(v_grad, tmp_m_x, tmp_list)
      deallocate(tmp_R, tmp_x)

    ! Trust region
    else
    
      ! Allocation of temporary arrays
      allocate(v_grad(tmp_n), H(tmp_n, tmp_n), tmp_m_x(tmp_list_size, tmp_list_size))
      allocate(tmp_R(tmp_list_size, tmp_list_size))
      allocate(tmp_x(tmp_n), W(tmp_n,tmp_n), e_val(tmp_n), key(tmp_n))

      ! ### Initialization ###
      delta = 0d0 ! can be deleted (normally)
      nb_iter = 0 ! Must start at 0 !!!
      rho = 0.5d0 ! Must be 0.5

      ! Compute the criterion before the loop
      call criterion_localization(tmp_list_size, tmp_list, prev_criterion)

      ! Loop until the convergence
      do while (not_converged)

        print*,''
        print*,'***********************'
        print*,'Iteration', nb_iter
        print*,'***********************'
        print*,''
  
        ! Gradient
        call gradient_localization(tmp_n, tmp_list_size, tmp_list, v_grad, max_elem, norm_grad)
        ! Diagonal hessian
        call hessian_localization(tmp_n, tmp_list_size, tmp_list, H)
        
        ! Diagonalization of the diagonal hessian by hands
        !call diagonalization_hessian(tmp_n,H,e_val,w)
        do i = 1, tmp_n
          e_val(i) = H(i,i)
        enddo

        ! Key list for dsort
        do i = 1, tmp_n 
          key(i) = i
        enddo

        ! Sort of the eigenvalues
        call dsort(e_val, key, tmp_n)

        ! Eigenvectors
        W = 0d0
        do i = 1, tmp_n
          j = key(i)
          W(j,i) = 1d0
        enddo

        ! To enter in the loop just after
        cancel_step = .True.
        nb_sub_iter = 0 

        ! Loop to reduce the trust radius until the criterion decreases and rho >= thresh_rho
        do while (cancel_step)
          print*,'-----------------------------'
          print*, mo_class(tmp_list(1))
          print*,'Iteration:', nb_iter
          print*,'Sub iteration:', nb_sub_iter
          print*,'-----------------------------'

          ! Hessian,gradient,Criterion -> x 
          call trust_region_step_w_expected_e(tmp_n, H, W, e_val, v_grad, prev_criterion, &
               rho, nb_iter, delta, criterion_model, tmp_x, must_exit)

          ! Internal loop exit condition
          if (must_exit) then
            print*,'trust_region_step_w_expected_e sent: Exit'
            exit 
          endif

          ! 1D tmp -> 2D tmp 
          call vec_to_mat_v2(tmp_n, tmp_list_size, tmp_x, tmp_m_x)

          ! Rotation submatrix (square matrix tmp_list_size by tmp_list_size)
          call rotation_matrix(tmp_m_x, tmp_list_size, tmp_R, tmp_list_size, tmp_list_size, &
               info, enforce_step_cancellation)

          if (enforce_step_cancellation) then
            print*, 'Step cancellation, too large error in the rotation matrix'
            rho = 0d0
            cycle
          endif

          ! tmp_R to R, subspace to full space
          call sub_to_full_rotation_matrix(tmp_list_size, tmp_list, tmp_R, R)
        
          ! Rotation of the MOs
          call apply_mo_rotation(R, prev_mos)   

          ! Update the things related to mo_coef
          call update_data_localization()

          ! Update the criterion
          call criterion_localization(tmp_list_size, tmp_list, criterion)
          print*,'Criterion:', trim(mo_class(tmp_list(1))), nb_iter, criterion

          ! Criterion -> step accepted or rejected 
          call trust_region_is_step_cancelled(nb_iter, prev_criterion, criterion, &
               criterion_model, rho, cancel_step)

          ! Cancellation of the step, previous MOs
          if (cancel_step) then
            mo_coef = prev_mos
          endif

          nb_sub_iter = nb_sub_iter + 1
        enddo
        !call save_mos() !### depend of the time for 1 iteration

        ! To exit the external loop if must_exti = .True.
        if (must_exit) then
          exit
        endif 

        ! Step accepted, nb iteration + 1
        nb_iter = nb_iter + 1

        ! External loop exit conditions
        if (DABS(max_elem) < thresh_loc_max_elem_grad) then
          not_converged = .False.
        endif
        if (nb_iter > localization_max_nb_iter) then
          not_converged = .False.
        endif
      enddo

      ! Deallocation of temporary arrays
      deallocate(v_grad, H, tmp_m_x, tmp_R, tmp_list, tmp_x, W, e_val, key)
      
      ! Save the MOs
      call save_mos()
      TOUCH mo_coef
 
      ! Debug
      if (debug_hf) then
        touch mo_coef
        print*,'HF energy:', HF_energy
      endif
      
    endif
  enddo


  TOUCH mo_coef 

  ! To sort the MOs using the diagonal elements of the Fock matrix
  if (sort_mos_by_e) then
    call run_sort_by_fock_energies()
  endif

  ! Debug
  if (debug_hf) then
    touch mo_coef
    print*,'HF energy:', HF_energy
  endif

  ! Locality after the localization
  call compute_spatial_extent(spatial_extent)

end
