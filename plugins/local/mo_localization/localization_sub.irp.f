! Gathering
! Gradient/hessian/criterion for the localization:
! They are chosen in function of the localization method

! Gradient:

! qp_edit : 
! | localization_method | method for the localization |

! Input: 
! | tmp_n                   | integer          | Number of parameters in the MO subspace           |
! | tmp_list_size           | integer          | Number of MOs in the mo_class we want to localize |
! | tmp_list(tmp_list_size) | integer          | MOs in the mo_class                               |

! Output:
! | v_grad(tmp_n)           | double precision | Gradient in the subspace                          |
! | max_elem                | double precision | Maximal element in the gradient                   |
! | norm_grad               | double precision | Norm of the gradient                              |



subroutine gradient_localization(tmp_n, tmp_list_size, tmp_list, v_grad, max_elem, norm_grad)
  
  include 'pi.h'

  implicit none

  BEGIN_DOC
  ! Compute the gradient of the chosen localization method
  END_DOC
  
  integer, intent(in)           :: tmp_n, tmp_list_size, tmp_list(tmp_list_size)
  double precision, intent(out) :: v_grad(tmp_n), max_elem, norm_grad

  if (localization_method == 'boys') then
    call gradient_FB_omp(tmp_n, tmp_list_size, tmp_list, v_grad, max_elem, norm_grad)
    !call gradient_FB(tmp_n, tmp_list_size, tmp_list, v_grad, max_elem, norm_grad)
  elseif (localization_method== 'pipek') then
    call gradient_PM(tmp_n, tmp_list_size, tmp_list, v_grad, max_elem, norm_grad)
  else
    print*,'Unkown method:'//localization_method
    call abort
  endif

end



! Hessian:

! Output:
! | H(tmp_n,tmp_n) | double precision | Gradient in the subspace        |
! | max_elem       | double precision | Maximal element in the gradient |
! | norm_grad      | double precision | Norm of the gradient            |


subroutine hessian_localization(tmp_n, tmp_list_size, tmp_list, H)

  include 'pi.h'

  implicit none

  BEGIN_DOC
  ! Compute the diagonal hessian of the chosen localization method
  END_DOC

  integer, intent(in)           :: tmp_n, tmp_list_size, tmp_list(tmp_list_size)
  double precision, intent(out) :: H(tmp_n)

  if (localization_method == 'boys') then
    call hessian_FB_omp(tmp_n, tmp_list_size, tmp_list, H)
    !call hessian_FB(tmp_n, tmp_list_size, tmp_list, H) ! non OMP for debugging
  elseif (localization_method == 'pipek') then
    call hessian_PM(tmp_n, tmp_list_size, tmp_list, H)
  else
    print*,'Unkown method: '//localization_method
    call abort
  endif

end



! Criterion:

! Output:
! | criterion | double precision | Criterion for the orbital localization |


subroutine criterion_localization(tmp_list_size, tmp_list,criterion)

  include 'pi.h'
  
  implicit none

  BEGIN_DOC
  ! Compute the localization criterion of the chosen localization method
  END_DOC

  integer, intent(in)           :: tmp_list_size, tmp_list(tmp_list_size)
  double precision, intent(out) :: criterion

  if (localization_method == 'boys') then
    call criterion_FB(tmp_list_size, tmp_list, criterion)
  elseif (localization_method == 'pipek') then
    !call criterion_PM(tmp_list_size, tmp_list,criterion)
    call criterion_PM_v3(tmp_list_size, tmp_list, criterion)
  else
    print*,'Unkown method: '//localization_method
    call abort
  endif

end



! Subroutine to update the datas needed for the localization

subroutine update_data_localization()

  include 'pi.h'

  implicit none

  if (localization_method == 'boys') then
    ! Update the dipoles
    call ao_to_mo_no_sym(ao_dipole_x, ao_num, mo_dipole_x, mo_num)
    call ao_to_mo_no_sym(ao_dipole_y, ao_num, mo_dipole_y, mo_num)
    call ao_to_mo_no_sym(ao_dipole_z, ao_num, mo_dipole_z, mo_num)
  elseif (localization_method == 'pipek') then
    ! Nothing required
  else
    print*,'Unkown method: '//localization_method
    call abort
  endif
end



! Angles:

! Output:
! | tmp_m_x(tmp_list_size, tmp_list_size) | double precision | Angles for the rotations in the subspace |
! | max_elem                              | double precision | Maximal angle                            |



subroutine theta_localization(tmp_list, tmp_list_size, tmp_m_x, max_elem)
  
  include 'pi.h'

  implicit none

  BEGIN_DOC
  ! Compute the rotation angles between the MOs for the chosen localization method
  END_DOC
  
  integer, intent(in)           :: tmp_list_size, tmp_list(tmp_list_size)
  double precision, intent(out) :: tmp_m_x(tmp_list_size,tmp_list_size), max_elem

  if (localization_method == 'boys') then
    call theta_FB(tmp_list, tmp_list_size, tmp_m_x, max_elem)
  elseif (localization_method== 'pipek') then
    call theta_PM(tmp_list, tmp_list_size, tmp_m_x, max_elem)
  else
    print*,'Unkown method: '//localization_method
    call abort
  endif

end

! Gradient
! Input: 
! | tmp_n                   | integer          | Number of parameters in the MO subspace           |
! | tmp_list_size           | integer          | Number of MOs in the mo_class we want to localize |
! | tmp_list(tmp_list_size) | integer          | MOs in the mo_class                               |

! Output:
! | v_grad(tmp_n)           | double precision | Gradient in the subspace                          |
! | max_elem                | double precision | Maximal element in the gradient                   |
! | norm_grad               | double precision | Norm of the gradient                              |

! Internal:
! | m_grad(tmp_n,tmp_n) | double precision | Gradient in the matrix form |
! | i,j,k               | integer          | indexes in the full space   |
! | tmp_i,tmp_j,tmp_k   | integer          | indexes in the subspace     |
! | t*                  | double precision | to compute the time         |


subroutine gradient_FB(tmp_n, tmp_list_size, tmp_list, v_grad, max_elem, norm_grad)
  
  implicit none
  
  BEGIN_DOC
  ! Compute the gradient for the Foster-Boys localization
  END_DOC
  
  integer, intent(in)           :: tmp_n, tmp_list_size, tmp_list(tmp_list_size)
  double precision, intent(out) :: v_grad(tmp_n), max_elem, norm_grad
  double precision, allocatable :: m_grad(:,:)
  integer                       :: i,j,k,tmp_i,tmp_j,tmp_k
  double precision              :: t1, t2, t3

  print*,''
  print*,'---gradient_FB---'

  call wall_time(t1)

  ! Allocation
  allocate(m_grad(tmp_list_size, tmp_list_size))

  ! Calculation
  do tmp_j = 1, tmp_list_size
    j = tmp_list(tmp_j)
    do tmp_i = 1, tmp_list_size
      i = tmp_list(tmp_i)
      m_grad(tmp_i,tmp_j) = 4d0 * mo_dipole_x(i,j) * (mo_dipole_x(i,i) - mo_dipole_x(j,j)) &
                           +4d0 * mo_dipole_y(i,j) * (mo_dipole_y(i,i) - mo_dipole_y(j,j)) &
                           +4d0 * mo_dipole_z(i,j) * (mo_dipole_z(i,i) - mo_dipole_z(j,j))
    enddo
  enddo
  
  ! 2D -> 1D
  do tmp_k = 1, tmp_n
    call vec_to_mat_index(tmp_k,tmp_i,tmp_j)
    v_grad(tmp_k) = m_grad(tmp_i,tmp_j) 
  enddo

  ! Maximum element in the gradient
  max_elem = 0d0
  do tmp_k = 1, tmp_n
    if (ABS(v_grad(tmp_k)) > max_elem) then
      max_elem = ABS(v_grad(tmp_k))
    endif
  enddo 
 
  ! Norm of the gradient
  norm_grad = 0d0
  do tmp_k = 1, tmp_n
    norm_grad = norm_grad + v_grad(tmp_k)**2
  enddo
  norm_grad = dsqrt(norm_grad)

  print*, 'Maximal element in the gradient:', max_elem
  print*, 'Norm of the gradient:', norm_grad  

  ! Deallocation
  deallocate(m_grad)

  call wall_time(t2)
  t3 = t2 - t1
  print*,'Time in gradient_FB:', t3

  print*,'---End gradient_FB---'

end subroutine

! Gradient (OMP)

subroutine gradient_FB_omp(tmp_n, tmp_list_size, tmp_list, v_grad, max_elem, norm_grad)
  
  use omp_lib

  implicit none

  BEGIN_DOC
  ! Compute the gradient for the Foster-Boys localization
  END_DOC
  
  integer, intent(in)           :: tmp_n, tmp_list_size, tmp_list(tmp_list_size)
  double precision, intent(out) :: v_grad(tmp_n), max_elem, norm_grad
  double precision, allocatable :: m_grad(:,:)
  integer                       :: i,j,k,tmp_i,tmp_j,tmp_k
  double precision              :: t1, t2, t3

  print*,''
  print*,'---gradient_FB_omp---'

  call wall_time(t1)

  ! Allocation
  allocate(m_grad(tmp_list_size, tmp_list_size))

  ! Initialization omp
  call omp_set_max_active_levels(1)

  !$OMP PARALLEL                                                             &
      !$OMP PRIVATE(i,j,tmp_i,tmp_j,tmp_k)                                         &
      !$OMP SHARED(tmp_n,tmp_list_size,m_grad,v_grad,mo_dipole_x,mo_dipole_y,mo_dipole_z,tmp_list) &
      !$OMP DEFAULT(NONE)

  ! Calculation
  !$OMP DO
  do tmp_j = 1, tmp_list_size
    j = tmp_list(tmp_j)
    do tmp_i = 1, tmp_list_size
      i = tmp_list(tmp_i)
      m_grad(tmp_i,tmp_j) = 4d0 * mo_dipole_x(i,j) * (mo_dipole_x(i,i) - mo_dipole_x(j,j)) &
                           +4d0 * mo_dipole_y(i,j) * (mo_dipole_y(i,i) - mo_dipole_y(j,j)) &
                           +4d0 * mo_dipole_z(i,j) * (mo_dipole_z(i,i) - mo_dipole_z(j,j))
    enddo
  enddo
  !$OMP END DO

  ! 2D -> 1D
  !$OMP DO
  do tmp_k = 1, tmp_n
    call vec_to_mat_index(tmp_k,tmp_i,tmp_j)
    v_grad(tmp_k) = m_grad(tmp_i,tmp_j) 
  enddo 
  !$OMP END DO

  !$OMP END PARALLEL

  call omp_set_max_active_levels(4)

  ! Maximum element in the gradient
  max_elem = 0d0
  do tmp_k = 1, tmp_n
    if (ABS(v_grad(tmp_k)) > max_elem) then
      max_elem = ABS(v_grad(tmp_k))
    endif
  enddo 

  ! Norm of the gradient
  norm_grad = 0d0
  do tmp_k = 1, tmp_n
    norm_grad = norm_grad + v_grad(tmp_k)**2
  enddo
  norm_grad = dsqrt(norm_grad)

  print*, 'Maximal element in the gradient:', max_elem
  print*, 'Norm of the gradient:', norm_grad  

  ! Deallocation
  deallocate(m_grad)

  call wall_time(t2)
  t3 = t2 - t1
  print*,'Time in gradient_FB_omp:', t3

  print*,'---End gradient_FB_omp---'

end subroutine

! Hessian 

! Output:
! | H(tmp_n,tmp_n) | double precision | Gradient in the subspace        |
! | max_elem       | double precision | Maximal element in the gradient |
! | norm_grad      | double precision | Norm of the gradient            |

! Internal:
! Internal:
! | beta(tmp_n,tmp_n) | double precision | beta in the documentation below to compute the hesian |
! | i,j,k             | integer          | indexes in the full space                             |
! | tmp_i,tmp_j,tmp_k | integer          | indexes in the subspace                               |
! | t*                | double precision | to compute the time                                   |


subroutine hessian_FB(tmp_n, tmp_list_size, tmp_list, H)

  implicit none
  
  BEGIN_DOC
  ! Compute the diagonal hessian for the Foster-Boys localization
  END_DOC

  integer, intent(in)           :: tmp_n, tmp_list_size, tmp_list(tmp_list_size)
  double precision, intent(out) :: H(tmp_n)
  double precision, allocatable :: beta(:,:)
  integer                       :: i,j,tmp_k,tmp_i, tmp_j
  double precision              :: max_elem, t1,t2,t3
   
  print*,''
  print*,'---hessian_FB---'

  call wall_time(t1)


  ! Allocation
  allocate(beta(tmp_list_size,tmp_list_size))
  
  ! Calculation
  do tmp_j = 1, tmp_list_size
    j = tmp_list(tmp_j)
    do tmp_i = 1, tmp_list_size
      i = tmp_list(tmp_i)
      beta(tmp_i,tmp_j) = (mo_dipole_x(i,i) - mo_dipole_x(j,j))**2 - 4d0 * mo_dipole_x(i,j)**2 &
                         +(mo_dipole_y(i,i) - mo_dipole_y(j,j))**2 - 4d0 * mo_dipole_y(i,j)**2 &
                         +(mo_dipole_z(i,i) - mo_dipole_z(j,j))**2 - 4d0 * mo_dipole_z(i,j)**2
    enddo
  enddo

  ! Diagonal of the hessian
  H = 0d0
  do tmp_k = 1, tmp_n
    call vec_to_mat_index(tmp_k,tmp_i,tmp_j)
    H(tmp_k) = 4d0 * beta(tmp_i, tmp_j)
  enddo
  
  ! Deallocation
  deallocate(beta)
 
  call wall_time(t2)
  t3 = t2 - t1
  print*,'Time in hessian_FB:', t3

  print*,'---End hessian_FB---'

end subroutine

! Hessian (OMP)

subroutine hessian_FB_omp(tmp_n, tmp_list_size, tmp_list, H)

  implicit none

  BEGIN_DOC
  ! Compute the diagonal hessian for the Foster-Boys localization
  END_DOC

  integer, intent(in)           :: tmp_n, tmp_list_size, tmp_list(tmp_list_size)
  double precision, intent(out) :: H(tmp_n)
  double precision, allocatable :: beta(:,:)
  integer                       :: i,j,tmp_k,tmp_i,tmp_j
  double precision              :: max_elem, t1,t2,t3
   
  print*,''
  print*,'---hessian_FB_omp---'

  call wall_time(t1)

  ! Allocation
  allocate(beta(tmp_list_size,tmp_list_size))

  ! Initialization omp
  call omp_set_max_active_levels(1)

  !$OMP PARALLEL                                                             &
      !$OMP PRIVATE(i,j,tmp_i,tmp_j,tmp_k)                                         &
      !$OMP SHARED(tmp_n,tmp_list_size,beta,H,mo_dipole_x,mo_dipole_y,mo_dipole_z,tmp_list) &
      !$OMP DEFAULT(NONE)

  
  ! Calculation
  !$OMP DO
  do tmp_j = 1, tmp_list_size
    j = tmp_list(tmp_j)
    do tmp_i = 1, tmp_list_size
      i = tmp_list(tmp_i)
      beta(tmp_i,tmp_j) = (mo_dipole_x(i,i) - mo_dipole_x(j,j))**2 - 4d0 * mo_dipole_x(i,j)**2 &
                         +(mo_dipole_y(i,i) - mo_dipole_y(j,j))**2 - 4d0 * mo_dipole_y(i,j)**2 &
                         +(mo_dipole_z(i,i) - mo_dipole_z(j,j))**2 - 4d0 * mo_dipole_z(i,j)**2
    enddo
  enddo
  !$OMP END DO

  ! Initialization
  !$OMP DO
  do i = 1, tmp_n
    H(i) = 0d0 
  enddo
  !$OMP END DO
  
  ! Diagonalm of the hessian
  !$OMP DO
  do tmp_k = 1, tmp_n
    call vec_to_mat_index(tmp_k,tmp_i,tmp_j)
    H(tmp_k) = 4d0 * beta(tmp_i, tmp_j)
  enddo
  !$OMP END DO
  
  !$OMP END PARALLEL

  call omp_set_max_active_levels(4)

  ! Deallocation
  deallocate(beta)
 
  call wall_time(t2)
  t3 = t2 - t1
  print*,'Time in hessian_FB_omp:', t3

  print*,'---End hessian_FB_omp---'

end subroutine

! Spatial extent

! The spatial extent of an orbital $i$ is computed as 
! \begin{align*}
! \sum_{\lambda=x,y,z}\sqrt{<i|\lambda^2|i> - <i|\lambda|i>^2}
! \end{align*}

! From that we can also compute the average and the standard deviation



subroutine compute_spatial_extent(spatial_extent)

  implicit none

  BEGIN_DOC
  ! Compute the spatial extent of the MOs
  END_DOC
 
  double precision, intent(out) :: spatial_extent(mo_num)
  double precision              :: average_core, average_act, average_inact, average_virt
  double precision              :: std_var_core, std_var_act, std_var_inact, std_var_virt
  integer                       :: i,j,k,l

  spatial_extent = 0d0
  
  do i = 1, mo_num
    spatial_extent(i) = mo_spread_x(i,i) - mo_dipole_x(i,i)**2
  enddo
  do i = 1, mo_num
    spatial_extent(i) = spatial_extent(i) + mo_spread_y(i,i) - mo_dipole_y(i,i)**2
  enddo
  do i = 1, mo_num
    spatial_extent(i) = spatial_extent(i) + mo_spread_z(i,i) - mo_dipole_z(i,i)**2
  enddo

  do i = 1, mo_num
    spatial_extent(i) = dsqrt(spatial_extent(i))
  enddo

  average_core = 0d0
  std_var_core = 0d0
  if (dim_list_core_orb >= 2) then
    call compute_average_sp_ext(spatial_extent, list_core, dim_list_core_orb, average_core)
    call compute_std_var_sp_ext(spatial_extent, list_core, dim_list_core_orb, average_core, std_var_core)
  endif

  average_act = 0d0
  std_var_act = 0d0
  if (dim_list_act_orb >= 2) then
    call compute_average_sp_ext(spatial_extent, list_act, dim_list_act_orb, average_act)
    call compute_std_var_sp_ext(spatial_extent, list_act, dim_list_act_orb, average_act, std_var_act)
  endif
  
  average_inact = 0d0
  std_var_inact = 0d0
  if (dim_list_inact_orb >= 2) then
    call compute_average_sp_ext(spatial_extent, list_inact, dim_list_inact_orb, average_inact)
    call compute_std_var_sp_ext(spatial_extent, list_inact, dim_list_inact_orb, average_inact, std_var_inact)
  endif

  average_virt = 0d0
  std_var_virt = 0d0
  if (dim_list_virt_orb >= 2) then
    call compute_average_sp_ext(spatial_extent, list_virt, dim_list_virt_orb, average_virt)
    call compute_std_var_sp_ext(spatial_extent, list_virt, dim_list_virt_orb, average_virt, std_var_virt)
  endif

  print*,''
  print*,'============================='
  print*,'  Spatial extent of the MOs'
  print*,'============================='
  print*,''

  print*, 'elec_num:', elec_num
  print*, 'elec_alpha_num:', elec_alpha_num
  print*, 'elec_beta_num:', elec_beta_num
  print*, 'core:', dim_list_core_orb
  print*, 'act:', dim_list_act_orb
  print*, 'inact:', dim_list_inact_orb
  print*, 'virt:', dim_list_virt_orb
  print*, 'mo_num:', mo_num
  print*,''
   
  print*,'-- Core MOs --'
  print*,'Average:', average_core
  print*,'Std var:', std_var_core
  print*,''
  
  print*,'-- Active MOs --'
  print*,'Average:', average_act
  print*,'Std var:', std_var_act
  print*,''

  print*,'-- Inactive MOs --'
  print*,'Average:', average_inact
  print*,'Std var:', std_var_inact
  print*,''

  print*,'-- Virtual MOs --'
  print*,'Average:', average_virt
  print*,'Std var:', std_var_virt
  print*,''

  print*,'Spatial extent:'
  do i = 1, mo_num
    print*, i, spatial_extent(i)
  enddo

end

subroutine compute_average_sp_ext(spatial_extent, list, list_size, average)

  implicit none
  
  BEGIN_DOC
  ! Compute the average spatial extent of the MOs
  END_DOC

  integer, intent(in) :: list_size, list(list_size)
  double precision, intent(in) :: spatial_extent(mo_num)
  double precision, intent(out) :: average
  integer :: i, tmp_i
  
  average = 0d0
  do tmp_i = 1, list_size
    i = list(tmp_i)
    average = average + spatial_extent(i)
  enddo

  average = average / DBLE(list_size)

end

subroutine compute_std_var_sp_ext(spatial_extent, list, list_size, average, std_var)

  implicit none

  BEGIN_DOC
  ! Compute the standard deviation of the spatial extent of the MOs
  END_DOC

  integer, intent(in)           :: list_size, list(list_size)
  double precision, intent(in)  :: spatial_extent(mo_num)
  double precision, intent(in)  :: average
  double precision, intent(out) :: std_var
  integer                       :: i, tmp_i

  std_var = 0d0

  do tmp_i = 1, list_size
    i = list(tmp_i)
    std_var = std_var + (spatial_extent(i) - average)**2
  enddo
  
  std_var = dsqrt(1d0/DBLE(list_size) * std_var)

end

! Utils


subroutine apply_pre_rotation()

  implicit none

  BEGIN_DOC
  ! Apply a rotation between the MOs
  END_DOC

  double precision, allocatable :: pre_rot(:,:), prev_mos(:,:), R(:,:)
  double precision              :: t1,t2,t3
  integer                       :: i,j,tmp_i,tmp_j
  integer                       :: info
  logical                       :: enforce_step_cancellation

  print*,'---apply_pre_rotation---'
  call wall_time(t1)

  allocate(pre_rot(mo_num,mo_num), prev_mos(ao_num,mo_num), R(mo_num,mo_num))

  ! Initialization of the matrix
  pre_rot = 0d0

  if (kick_in_mos) then
    ! Pre rotation for core MOs
    if (dim_list_core_orb >= 2) then
      do tmp_j = 1, dim_list_core_orb
        j = list_core(tmp_j)
        do tmp_i = 1, dim_list_core_orb
          i = list_core(tmp_i)
          if (i > j) then
            pre_rot(i,j) = angle_pre_rot
          elseif (i < j) then
            pre_rot(i,j) = - angle_pre_rot
          else
            pre_rot(i,j) = 0d0
          endif
        enddo
      enddo
    endif
    
    ! Pre rotation for active MOs
    if (dim_list_act_orb >= 2) then
      do tmp_j = 1, dim_list_act_orb
        j = list_act(tmp_j)
        do tmp_i = 1, dim_list_act_orb
          i = list_act(tmp_i)
          if (i > j) then
            pre_rot(i,j) = angle_pre_rot
          elseif (i < j) then
            pre_rot(i,j) = - angle_pre_rot
          else
            pre_rot(i,j) = 0d0
          endif
        enddo
      enddo
    endif
  
    ! Pre rotation for inactive MOs
    if (dim_list_inact_orb >= 2) then
      do tmp_j = 1, dim_list_inact_orb
        j = list_inact(tmp_j)
        do tmp_i = 1, dim_list_inact_orb
          i = list_inact(tmp_i)
          if (i > j) then
            pre_rot(i,j) = angle_pre_rot
          elseif (i < j) then
            pre_rot(i,j) = - angle_pre_rot
          else
            pre_rot(i,j) = 0d0
          endif
        enddo
      enddo
    endif
  
    ! Pre rotation for virtual MOs
    if (dim_list_virt_orb >= 2) then
      do tmp_j = 1, dim_list_virt_orb
        j = list_virt(tmp_j)
        do tmp_i = 1, dim_list_virt_orb
          i = list_virt(tmp_i)
          if (i > j) then
            pre_rot(i,j) = angle_pre_rot
          elseif (i < j) then
            pre_rot(i,j) = - angle_pre_rot
          else
            pre_rot(i,j) = 0d0
          endif
        enddo
      enddo
    endif
  
    ! Nothing for deleted ones
  
    ! Compute pre rotation matrix from pre_rot
    call rotation_matrix(pre_rot,mo_num,R,mo_num,mo_num,info,enforce_step_cancellation)
   
    if (enforce_step_cancellation) then
      print*, 'Cancellation of the pre rotation, too big error in the rotation matrix'
      print*, 'Reduce the angle for the pre rotation, abort'
      call abort
    endif
  
    ! New Mos (we don't car eabout the previous MOs prev_mos)
    call apply_mo_rotation(R,prev_mos)
  
    ! Update the things related to mo_coef
    TOUCH mo_coef
    call save_mos
  endif

  deallocate(pre_rot, prev_mos, R)

  call wall_time(t2)
  t3 = t2-t1
  print*,'Time in apply_pre_rotation:', t3
  print*,'---End apply_pre_rotation---'

end

subroutine x_tmp_orb_loc_v2(tmp_n, tmp_list_size, tmp_list, v_grad, H,tmp_x, tmp_m_x)

  implicit none

  integer, intent(in)           :: tmp_n, tmp_list_size, tmp_list(tmp_list_size)
  double precision, intent(in)  :: v_grad(tmp_n)
  double precision, intent(in)  :: H(tmp_n, tmp_n)
  double precision, intent(out) :: tmp_m_x(tmp_list_size, tmp_list_size), tmp_x(tmp_list_size)
  !double precision, allocatable :: x(:)
  double precision              :: lambda , accu, max_elem
  integer                       :: i,j,tmp_i,tmp_j,tmp_k

  ! Allocation
  !allocate(x(tmp_n))

  ! Level shifted hessian
  lambda = 0d0
  do tmp_k = 1, tmp_n
    if (H(tmp_k,tmp_k) < lambda) then
      lambda = H(tmp_k,tmp_k)
    endif
  enddo

  ! min element in the hessian
  if (lambda < 0d0) then
    lambda = -lambda + 1d-6
  endif  
  
  print*, 'lambda', lambda
 
  ! Good
  do tmp_k = 1, tmp_n
    if (ABS(H(tmp_k,tmp_k)) > 1d-6) then
       tmp_x(tmp_k) = - 1d0/(ABS(H(tmp_k,tmp_k))+lambda) * v_grad(tmp_k)!(-v_grad(tmp_k))
      !x(tmp_k) = - 1d0/(ABS(H(tmp_k,tmp_k))+lambda) * (-v_grad(tmp_k)) 
    endif
  enddo

  ! 1D tmp -> 2D tmp 
  tmp_m_x = 0d0
  do tmp_j = 1, tmp_list_size - 1
    do tmp_i = tmp_j + 1, tmp_list_size
      call mat_to_vec_index(tmp_i,tmp_j,tmp_k)
      tmp_m_x(tmp_i, tmp_j) = tmp_x(tmp_k)!x(tmp_k)
    enddo
  enddo

  ! Antisym
  do tmp_i = 1, tmp_list_size - 1
    do tmp_j = tmp_i + 1, tmp_list_size
      tmp_m_x(tmp_i,tmp_j) = - tmp_m_x(tmp_j,tmp_i) 
    enddo
  enddo

  ! Deallocation
  !deallocate(x)

end subroutine

subroutine ao_to_mo_no_sym(A_ao,LDA_ao,A_mo,LDA_mo)
  implicit none
  BEGIN_DOC
  ! Transform A from the |AO| basis to the |MO| basis
  !
  ! $C^\dagger.A_{ao}.C$
  END_DOC
  integer, intent(in)            :: LDA_ao,LDA_mo
  double precision, intent(in)   :: A_ao(LDA_ao,ao_num)
  double precision, intent(out)  :: A_mo(LDA_mo,mo_num)
  double precision, allocatable  :: T(:,:)

  allocate ( T(ao_num,mo_num) )
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: T

  call dgemm('N','N', ao_num, mo_num, ao_num,                    &
      1.d0, A_ao,LDA_ao,                                             &
      mo_coef, size(mo_coef,1),                                      &
      0.d0, T, size(T,1))

  call dgemm('T','N', mo_num, mo_num, ao_num,                &
      1.d0, mo_coef,size(mo_coef,1),                                 &
      T, ao_num,                                                     &
      0.d0, A_mo, size(A_mo,1))

  deallocate(T)
end

subroutine run_sort_by_fock_energies()
  
  implicit none
  
  BEGIN_DOC
  ! Saves the current MOs ordered by diagonal element of the Fock operator.
  END_DOC
  
  integer                        :: i,j,k,l,tmp_i,tmp_k,tmp_list_size
  integer, allocatable           :: iorder(:), tmp_list(:)
  double precision, allocatable  :: fock_energies_tmp(:), tmp_mo_coef(:,:)

  ! Test
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
      print*,'MO class: ',trim(mo_class(tmp_list(1)))

      allocate(iorder(tmp_list_size), fock_energies_tmp(tmp_list_size), tmp_mo_coef(ao_num,tmp_list_size))
      !print*,'MOs before sorting them by f_p^p energies:'
      do i = 1, tmp_list_size
        tmp_i = tmp_list(i)
        fock_energies_tmp(i) = Fock_matrix_diag_mo(tmp_i)
        iorder(i) = i
        !print*, tmp_i, fock_energies_tmp(i)
      enddo
  
      call dsort(fock_energies_tmp, iorder, tmp_list_size)
      
      print*,'MOs after sorting them by f_p^p energies:'
      do i = 1, tmp_list_size
        k = iorder(i)
        tmp_k = tmp_list(k)
        print*, tmp_k, fock_energies_tmp(k)
        do j = 1, ao_num
          tmp_mo_coef(j,k) = mo_coef(j,tmp_k)
        enddo
      enddo

      ! Update the MOs after sorting them by energies
      do i = 1, tmp_list_size
        tmp_i = tmp_list(i)
        do j = 1, ao_num
          mo_coef(j,tmp_i) = tmp_mo_coef(j,i)
        enddo
      enddo

      if (debug_hf) then
        touch mo_coef
        print*,'HF energy:', HF_energy
      endif
      print*,''
      
      deallocate(iorder, fock_energies_tmp, tmp_list, tmp_mo_coef)
    endif

  enddo

  touch mo_coef
  call save_mos
  
end

function is_core(i)

  implicit none

  BEGIN_DOC
  ! True if the orbital i is a core orbital
  END_DOC

  integer, intent(in) :: i
  logical             :: is_core

  integer             :: j

  ! Init
  is_core = .False.

  ! Search
  do j = 1, dim_list_core_orb
    if (list_core(j) == i) then
      is_core = .True.
      exit
    endif
  enddo

end

function is_del(i)

  implicit none

  BEGIN_DOC
  ! True if the orbital i is a deleted orbital
  END_DOC

  integer, intent(in) :: i
  logical             :: is_del

  integer             :: j

  ! Init
  is_del = .False.

  ! Search
  do j = 1, dim_list_core_orb
    if (list_core(j) == i) then
      is_del = .True.
      exit
    endif
  enddo

end

subroutine set_classes_loc()

  implicit none

  integer :: i
  logical :: ok1, ok2
  logical :: is_core, is_del
  integer(bit_kind) :: res(N_int,2)

  if (auto_mo_class) then
    do i = 1, mo_num
      if (is_core(i)) cycle
      if (is_del(i)) cycle
      call apply_hole(psi_det(1,1,1), 1, i, res, ok1, N_int)
      call apply_hole(psi_det(1,1,1), 2, i, res, ok2, N_int)
      if (ok1 .and. ok2) then
        mo_class(i) = 'Inactive'
      else if (.not. ok1 .and. .not. ok2) then
        mo_class(i) = 'Virtual'
      else
        mo_class(i) = 'Active'
      endif
    enddo
    touch mo_class
  endif
  
end

subroutine unset_classes_loc()

  implicit none

  integer :: i
  logical :: ok1, ok2
  logical :: is_core, is_del
  integer(bit_kind) :: res(N_int,2)

  if (auto_mo_class) then
    do i = 1, mo_num
      if (is_core(i)) cycle
      if (is_del(i)) cycle
      mo_class(i) = 'Active'
    enddo
    touch mo_class
  endif

end
