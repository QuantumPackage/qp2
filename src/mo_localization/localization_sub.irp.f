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
  double precision, intent(out) :: H(tmp_n, tmp_n)

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
  print*,''

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

  print*,''
  print*,'---End gradient_FB---'
  print*,''

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
  print*,''

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

  print*,''
  print*,'---End gradient_FB_omp---'
  print*,''

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
  double precision, intent(out) :: H(tmp_n, tmp_n)
  double precision, allocatable :: beta(:,:)
  integer                       :: i,j,tmp_k,tmp_i, tmp_j
  double precision              :: max_elem, t1,t2,t3

  print*,''
  print*,'---hessian_FB---'
  print*,''

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
    H(tmp_k,tmp_k) = 4d0 * beta(tmp_i, tmp_j)
  enddo

  ! Min elem
  max_elem = 0d0
  do tmp_k = 1, tmp_n
    if (H(tmp_k,tmp_k) < max_elem) then
      max_elem = H(tmp_k,tmp_k)
    endif
  enddo
  print*, 'Min elem H:', max_elem

  ! Max elem
  max_elem = 0d0
  do tmp_k = 1, tmp_n
    if (H(tmp_k,tmp_k) > max_elem) then
      max_elem = H(tmp_k,tmp_k)
    endif
  enddo
  print*, 'Max elem H:', max_elem

  ! Near 0
  max_elem = 1d10
  do tmp_k = 1, tmp_n
    if (ABS(H(tmp_k,tmp_k)) < ABS(max_elem)) then
      max_elem = H(tmp_k,tmp_k)
    endif
  enddo
  print*, 'Near 0 elem H:', max_elem

  ! Deallocation
  deallocate(beta)

  call wall_time(t2)
  t3 = t2 - t1
  print*,'Time in hessian_FB:', t3

  print*,''
  print*,'---End hessian_FB---'
  print*,''

end subroutine

! Hessian (OMP)

subroutine hessian_FB_omp(tmp_n, tmp_list_size, tmp_list, H)

  implicit none

  BEGIN_DOC
  ! Compute the diagonal hessian for the Foster-Boys localization
  END_DOC

  integer, intent(in)           :: tmp_n, tmp_list_size, tmp_list(tmp_list_size)
  double precision, intent(out) :: H(tmp_n, tmp_n)
  double precision, allocatable :: beta(:,:)
  integer                       :: i,j,tmp_k,tmp_i,tmp_j
  double precision              :: max_elem, t1,t2,t3

  print*,''
  print*,'---hessian_FB_omp---'
  print*,''

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
  do j = 1, tmp_n
    do i = 1, tmp_n
      H(i,j) = 0d0
    enddo
  enddo
  !$OMP END DO

  ! Diagonalm of the hessian
  !$OMP DO
  do tmp_k = 1, tmp_n
    call vec_to_mat_index(tmp_k,tmp_i,tmp_j)
    H(tmp_k,tmp_k) = 4d0 * beta(tmp_i, tmp_j)
  enddo
  !$OMP END DO

  !$OMP END PARALLEL

  call omp_set_max_active_levels(4)

  ! Min elem
  max_elem = 0d0
  do tmp_k = 1, tmp_n
    if (H(tmp_k,tmp_k) < max_elem) then
      max_elem = H(tmp_k,tmp_k)
    endif
  enddo
  print*, 'Min elem H:', max_elem

  ! Max elem
  max_elem = 0d0
  do tmp_k = 1, tmp_n
    if (H(tmp_k,tmp_k) > max_elem) then
      max_elem = H(tmp_k,tmp_k)
    endif
  enddo
  print*, 'Max elem H:', max_elem

  ! Near 0
  max_elem = 1d10
  do tmp_k = 1, tmp_n
    if (ABS(H(tmp_k,tmp_k)) < ABS(max_elem)) then
      max_elem = H(tmp_k,tmp_k)
    endif
  enddo
  print*, 'Near 0 elem H:', max_elem

  ! Deallocation
  deallocate(beta)

  call wall_time(t2)
  t3 = t2 - t1
  print*,'Time in hessian_FB_omp:', t3

  print*,''
  print*,'---End hessian_FB_omp---'
  print*,''

end subroutine

! Gradient v1

subroutine grad_pipek(tmp_n, tmp_list_size, tmp_list, v_grad, max_elem, norm_grad)

  implicit none

  BEGIN_DOC
  ! Compute gradient for the Pipek-Mezey localization
  END_DOC

  integer, intent(in)           :: tmp_n, tmp_list_size, tmp_list(tmp_list_size)
  double precision, intent(out) :: v_grad(tmp_n), max_elem, norm_grad
  double precision, allocatable :: m_grad(:,:), tmp_int(:,:)
  integer                       :: i,j,k,tmp_i,tmp_j,tmp_k, a, b, mu ,rho

  ! Allocation
  allocate(m_grad(tmp_list_size, tmp_list_size), tmp_int(tmp_list_size, tmp_list_size))

  ! Initialization
  m_grad = 0d0

  do a = 1, nucl_num ! loop over the nuclei
     tmp_int = 0d0 ! Initialization for each nuclei

     ! Loop over the MOs of the a given mo_class to compute <i|P_a|j>
     do tmp_j = 1, tmp_list_size
        j = tmp_list(tmp_j)
        do tmp_i = 1, tmp_list_size
           i = tmp_list(tmp_i)
           do rho = 1, ao_num ! loop over all the AOs
              do b = 1, nucl_n_aos(a) ! loop over the number of AOs which belongs to the nuclei a
                 mu = nucl_aos(a,b) ! AO centered on atom a

                 tmp_int(tmp_i,tmp_j) = tmp_int(tmp_i,tmp_j) + 0.5d0 * (mo_coef(rho,i) * ao_overlap(rho,mu) * mo_coef(mu,j) &
                      + mo_coef(mu,i) * ao_overlap(mu,rho) * mo_coef(rho,j))

              enddo
           enddo
        enddo
     enddo

     ! Gradient
     do tmp_j = 1, tmp_list_size
        do tmp_i = 1, tmp_list_size

           m_grad(tmp_i,tmp_j) = m_grad(tmp_i,tmp_j) +  4d0 * tmp_int(tmp_i,tmp_j) * (tmp_int(tmp_i,tmp_i) - tmp_int(tmp_j,tmp_j))

        enddo
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
  deallocate(m_grad,tmp_int)

end subroutine grad_pipek

! Gradient

! The gradient is

! \begin{align*}
! \left. \frac{\partial \mathcal{P} (\theta)}{\partial \theta} \right|_{\theta=0}= \gamma^{PM}
! \end{align*}
! with
! \begin{align*}
! \gamma_{st}^{PM} = \sum_{A=1}^N <s|P_A|t> \left[ <s| P_A |s> - <t|P_A|t> \right]
! \end{align*}

! \begin{align*}
! <s|P_A|t> = \frac{1}{2} \sum_{\rho} \sum_{\mu \in A} \left[ c_{\rho}^{s*} S_{\rho \nu} c_{\mu}^{t} +c_{\mu}^{s*} S_{\mu \rho} c_{\rho}^t \right]
! \end{align*}
! $\sum_{\rho}$ -> sum over all the AOs
! $\sum_{\mu \in A}$ -> sum over the AOs which belongs to atom A
! $c^t$ -> expansion coefficient of orbital |t>

! Input:
! | tmp_n                   | integer          | Number of parameters in the MO subspace           |
! | tmp_list_size           | integer          | Number of MOs in the mo_class we want to localize |
! | tmp_list(tmp_list_size) | integer          | MOs in the mo_class                               |

! Output:
! | v_grad(tmp_n)           | double precision | Gradient in the subspace                          |
! | max_elem                | double precision | Maximal element in the gradient                   |
! | norm_grad               | double precision | Norm of the gradient                              |

! Internal:
! | m_grad(tmp_list_size,tmp_list_size)       | double precision | Gradient in a 2D array                                   |
! | tmp_int(tmp_list_size,tmp_list_size)      |                  | Temporary array to store the integrals                   |
! | tmp_accu(tmp_list_size,tmp_list_size)     |                  | Temporary array to store a matrix                        |
! |                                           |                  | product and compute tmp_int                              |
! | CS(tmp_list_size,ao_num)                  |                  | Array to store the result of mo_coef * ao_overlap        |
! | tmp_mo_coef(ao_num,tmp_list_size)         |                  | Array to store just the useful MO coefficients           |
! |                                           |                  | depending of the mo_class                                |
! | tmp_mo_coef2(nucl_n_aos(a),tmp_list_size) |                  | Array to store just the useful MO coefficients           |
! |                                           |                  | depending of the nuclei                                  |
! | tmp_CS(tmp_list_size,nucl_n_aos(a))       |                  | Array to store just the useful mo_coef * ao_overlap      |
! |                                           |                  | values depending of the nuclei                           |
! | a                                         |                  | index to loop over the nuclei                            |
! | b                                         |                  | index to loop over the AOs which belongs to the nuclei a |
! | mu                                        |                  | index to refer to an AO which belongs to the nuclei a    |
! | rho                                       |                  | index to loop over all the AOs                           |


subroutine gradient_PM(tmp_n, tmp_list_size, tmp_list, v_grad, max_elem, norm_grad)

  implicit none

  BEGIN_DOC
  ! Compute gradient for the Pipek-Mezey localization
  END_DOC

  integer, intent(in)           :: tmp_n, tmp_list_size, tmp_list(tmp_list_size)
  double precision, intent(out) :: v_grad(tmp_n), max_elem, norm_grad
  double precision, allocatable :: m_grad(:,:), tmp_int(:,:), CS(:,:), tmp_mo_coef(:,:), tmp_mo_coef2(:,:),tmp_accu(:,:),tmp_CS(:,:)
  integer                       :: i,j,k,tmp_i,tmp_j,tmp_k, a, b, mu ,rho
  double precision              :: t1,t2,t3

  print*,''
  print*,'---gradient_PM---'
  print*,''

  call wall_time(t1)

  ! Allocation
  allocate(m_grad(tmp_list_size, tmp_list_size), tmp_int(tmp_list_size, tmp_list_size),tmp_accu(tmp_list_size, tmp_list_size))
  allocate(CS(tmp_list_size,ao_num),tmp_mo_coef(ao_num,tmp_list_size))


  ! submatrix of the mo_coef
  do tmp_i = 1, tmp_list_size
    i = tmp_list(tmp_i)
    do j = 1, ao_num

      tmp_mo_coef(j,tmp_i) = mo_coef(j,i)

    enddo
  enddo

  call dgemm('T','N',tmp_list_size,ao_num,ao_num,1d0,tmp_mo_coef,size(tmp_mo_coef,1),ao_overlap,size(ao_overlap,1),0d0,CS,size(CS,1))

  m_grad = 0d0

  do a = 1, nucl_num ! loop over the nuclei
    tmp_int = 0d0

    !do tmp_j = 1, tmp_list_size
    !  do tmp_i = 1, tmp_list_size
    !    do b = 1, nucl_n_aos(a) ! loop over the number of AOs which belongs to the nuclei a
    !      mu = nucl_aos(a,b)

    !      tmp_int(tmp_i,tmp_j) = tmp_int(tmp_i,tmp_j) + 0.5d0 * (CS(tmp_i,mu) * tmp_mo_coef(mu,tmp_j) + tmp_mo_coef(mu,tmp_i) * CS(tmp_j,mu))

    !                             !  (mo_coef(rho,i) * ao_overlap(rho,mu) * mo_coef(mu,j) &
    !                             !+ mo_coef(mu,i) * ao_overlap(mu,rho) * mo_coef(rho,j))

    !    enddo
    !  enddo
    !enddo

    allocate(tmp_mo_coef2(nucl_n_aos(a),tmp_list_size),tmp_CS(tmp_list_size,nucl_n_aos(a)))

    do tmp_i = 1, tmp_list_size
      do b = 1, nucl_n_aos(a)
        mu = nucl_aos(a,b)

        tmp_mo_coef2(b,tmp_i) = tmp_mo_coef(mu,tmp_i)

      enddo
    enddo

    do b = 1, nucl_n_aos(a)
      mu = nucl_aos(a,b)
      do tmp_i = 1, tmp_list_size

        tmp_CS(tmp_i,b) = CS(tmp_i,mu)

      enddo
    enddo

    call dgemm('N','N',tmp_list_size,tmp_list_size,nucl_n_aos(a),1d0,tmp_CS,size(tmp_CS,1),tmp_mo_coef2,size(tmp_mo_coef2,1),0d0,tmp_accu,size(tmp_accu,1))

    do tmp_j = 1, tmp_list_size
      do tmp_i = 1, tmp_list_size

        tmp_int(tmp_i,tmp_j) = 0.5d0 * (tmp_accu(tmp_i,tmp_j) + tmp_accu(tmp_j,tmp_i))

      enddo
    enddo

    deallocate(tmp_mo_coef2,tmp_CS)

    do tmp_j = 1, tmp_list_size
      do tmp_i = 1, tmp_list_size

        m_grad(tmp_i,tmp_j) = m_grad(tmp_i,tmp_j) +  4d0 * tmp_int(tmp_i,tmp_j) * (tmp_int(tmp_i,tmp_i) - tmp_int(tmp_j,tmp_j))

      enddo
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
  deallocate(m_grad,tmp_int,CS,tmp_mo_coef)

  call wall_time(t2)
  t3 = t2 - t1
  print*,'Time in gradient_PM:', t3

  print*,''
  print*,'---End gradient_PM---'
  print*,''

end

! Hessian v1

subroutine hess_pipek(tmp_n, tmp_list_size, tmp_list, H)

  implicit none

  BEGIN_DOC
  ! Compute diagonal hessian for the Pipek-Mezey localization
  END_DOC

  integer, intent(in)           :: tmp_n, tmp_list_size, tmp_list(tmp_list_size)
  double precision, intent(out) :: H(tmp_n, tmp_n)
  double precision, allocatable :: beta(:,:),tmp_int(:,:)
  integer                       :: i,j,tmp_k,tmp_i, tmp_j, a,b,rho,mu
  double precision              :: max_elem

  ! Allocation
  allocate(beta(tmp_list_size,tmp_list_size),tmp_int(tmp_list_size,tmp_list_size))

  beta = 0d0

  do a = 1, nucl_num
    tmp_int = 0d0

    do tmp_j = 1, tmp_list_size
      j = tmp_list(tmp_j)
      do tmp_i = 1, tmp_list_size
        i = tmp_list(tmp_i)
        do rho = 1, ao_num
          do b = 1, nucl_n_aos(a) ! loop over the number of AOs which belongs to the nuclei a
            mu = nucl_aos(a,b)

            tmp_int(tmp_i,tmp_j) = tmp_int(tmp_i,tmp_j) + 0.5d0 * (mo_coef(rho,i) * ao_overlap(rho,mu) * mo_coef(mu,j) &
                                   + mo_coef(mu,i) * ao_overlap(mu,rho) * mo_coef(rho,j))

          enddo
        enddo
      enddo
    enddo

    ! Calculation
    do tmp_j = 1, tmp_list_size
      do tmp_i = 1, tmp_list_size

        beta(tmp_i,tmp_j) = beta(tmp_i, tmp_j) +  (tmp_int(tmp_i,tmp_i) - tmp_int(tmp_j,tmp_j))**2 - 4d0 * tmp_int(tmp_i,tmp_j)**2

      enddo
    enddo

  enddo

  H = 0d0
  do tmp_k = 1, tmp_n
    call vec_to_mat_index(tmp_k,tmp_i,tmp_j)
    H(tmp_k,tmp_k) = 4d0 * beta(tmp_i, tmp_j)
  enddo

!  max_elem = 0d0
!  do tmp_k = 1, tmp_n
!    if (H(tmp_k,tmp_k) < max_elem) then
!      max_elem = H(tmp_k,tmp_k)
!    endif
!  enddo
!  print*, 'Min elem H:', max_elem
!
!  max_elem = 0d0
!  do tmp_k = 1, tmp_n
!    if (H(tmp_k,tmp_k) > max_elem) then
!      max_elem = H(tmp_k,tmp_k)
!    endif
!  enddo
!  print*, 'Max elem H:', max_elem
!
!  max_elem = 1d10
!  do tmp_k = 1, tmp_n
!    if (ABS(H(tmp_k,tmp_k)) < ABS(max_elem)) then
!      max_elem = H(tmp_k,tmp_k)
!    endif
!  enddo
!  print*, 'Near 0 elem H:', max_elem

  ! Deallocation
  deallocate(beta,tmp_int)

end

! Hessian

! The hessian is
! \begin{align*}
! \left. \frac{\partial^2 \mathcal{P} (\theta)}{\partial \theta^2}\right|_{\theta=0} = 4 \beta^{PM}
! \end{align*}
! \begin{align*}
! \beta_{st}^{PM} = \sum_{A=1}^N \left( <s|P_A|t>^2 - \frac{1}{4} \left[<s|P_A|s> - <t|P_A|t> \right]^2 \right)
! \end{align*}

! with
! \begin{align*}
! <s|P_A|t> = \frac{1}{2} \sum_{\rho} \sum_{\mu \in A} \left[ c_{\rho}^{s*} S_{\rho \nu} c_{\mu}^{t} +c_{\mu}^{s*} S_{\mu \rho} c_{\rho}^t \right]
! \end{align*}
! $\sum_{\rho}$ -> sum over all the AOs
! $\sum_{\mu \in A}$ -> sum over the AOs which belongs to atom A
! $c^t$ -> expansion coefficient of orbital |t>


subroutine hessian_PM(tmp_n, tmp_list_size, tmp_list, H)

  implicit none

  BEGIN_DOC
  ! Compute diagonal hessian for the Pipek-Mezey localization
  END_DOC

  integer, intent(in)           :: tmp_n, tmp_list_size, tmp_list(tmp_list_size)
  double precision, intent(out) :: H(tmp_n, tmp_n)
  double precision, allocatable :: beta(:,:),tmp_int(:,:),CS(:,:),tmp_mo_coef(:,:),tmp_mo_coef2(:,:),tmp_accu(:,:),tmp_CS(:,:)
  integer                       :: i,j,tmp_k,tmp_i, tmp_j, a,b,rho,mu
  double precision              :: max_elem, t1,t2,t3

  print*,''
  print*,'---hessian_PM---'
  print*,''

  call wall_time(t1)

  ! Allocation
  allocate(beta(tmp_list_size,tmp_list_size),tmp_int(tmp_list_size,tmp_list_size),tmp_accu(tmp_list_size,tmp_list_size))
  allocate(CS(tmp_list_size,ao_num),tmp_mo_coef(ao_num,tmp_list_size))

  beta = 0d0

  do tmp_i = 1, tmp_list_size
    i = tmp_list(tmp_i)
    do j = 1, ao_num

      tmp_mo_coef(j,tmp_i) = mo_coef(j,i)

    enddo
  enddo

  call dgemm('T','N',tmp_list_size,ao_num,ao_num,1d0,tmp_mo_coef,size(tmp_mo_coef,1),ao_overlap,size(ao_overlap,1),0d0,CS,size(CS,1))

  do a = 1, nucl_num ! loop over the nuclei
    tmp_int = 0d0

    !do tmp_j = 1, tmp_list_size
    !  do tmp_i = 1, tmp_list_size
    !    do b = 1, nucl_n_aos(a) ! loop over the number of AOs which belongs to the nuclei a
    !      mu = nucl_aos(a,b)

    !      tmp_int(tmp_i,tmp_j) = tmp_int(tmp_i,tmp_j) + 0.5d0 * (CS(tmp_i,mu) * tmp_mo_coef(mu,tmp_j) + tmp_mo_coef(mu,tmp_i) * CS(tmp_j,mu))

    !                             !  (mo_coef(rho,i) * ao_overlap(rho,mu) * mo_coef(mu,j) &
    !                             !+ mo_coef(mu,i) * ao_overlap(mu,rho) * mo_coef(rho,j))

    !    enddo
    !  enddo
    !enddo

    allocate(tmp_mo_coef2(nucl_n_aos(a),tmp_list_size),tmp_CS(tmp_list_size,nucl_n_aos(a)))

    do tmp_i = 1, tmp_list_size
      do b = 1, nucl_n_aos(a)
        mu = nucl_aos(a,b)

        tmp_mo_coef2(b,tmp_i) = tmp_mo_coef(mu,tmp_i)

      enddo
    enddo

    do b = 1, nucl_n_aos(a)
      mu = nucl_aos(a,b)
      do tmp_i = 1, tmp_list_size

        tmp_CS(tmp_i,b) = CS(tmp_i,mu)

      enddo
    enddo

    call dgemm('N','N',tmp_list_size,tmp_list_size,nucl_n_aos(a),1d0,tmp_CS,size(tmp_CS,1),tmp_mo_coef2,size(tmp_mo_coef2,1),0d0,tmp_accu,size(tmp_accu,1))

    do tmp_j = 1, tmp_list_size
      do tmp_i = 1, tmp_list_size

        tmp_int(tmp_i,tmp_j) = 0.5d0 * (tmp_accu(tmp_i,tmp_j) + tmp_accu(tmp_j,tmp_i))

      enddo
    enddo

    deallocate(tmp_mo_coef2,tmp_CS)

    ! Calculation
    do tmp_j = 1, tmp_list_size
      do tmp_i = 1, tmp_list_size

        beta(tmp_i,tmp_j) = beta(tmp_i, tmp_j) +  (tmp_int(tmp_i,tmp_i) - tmp_int(tmp_j,tmp_j))**2 - 4d0 * tmp_int(tmp_i,tmp_j)**2

      enddo
    enddo

  enddo

  H = 0d0
  do tmp_k = 1, tmp_n
    call vec_to_mat_index(tmp_k,tmp_i,tmp_j)
    H(tmp_k,tmp_k) = 4d0 * beta(tmp_i, tmp_j)
  enddo

  max_elem = 0d0
  do tmp_k = 1, tmp_n
    if (H(tmp_k,tmp_k) < max_elem) then
      max_elem = H(tmp_k,tmp_k)
    endif
  enddo
  print*, 'Min elem H:', max_elem

  max_elem = 0d0
  do tmp_k = 1, tmp_n
    if (H(tmp_k,tmp_k) > max_elem) then
      max_elem = H(tmp_k,tmp_k)
    endif
  enddo
  print*, 'Max elem H:', max_elem

  max_elem = 1d10
  do tmp_k = 1, tmp_n
    if (ABS(H(tmp_k,tmp_k)) < ABS(max_elem)) then
      max_elem = H(tmp_k,tmp_k)
    endif
  enddo
  print*, 'Near 0 elem H:', max_elem

  ! Deallocation
  deallocate(beta,tmp_int)

  call wall_time(t2)
  t3 = t2 - t1
  print*,'Time in hessian_PM:', t3

  print*,''
  print*,'---End hessian_PM---'
  print*,''

end

! Criterion PM (old)

subroutine compute_crit_pipek(criterion)

  implicit none

  BEGIN_DOC
  ! Compute the Pipek-Mezey localization criterion
  END_DOC

  double precision, intent(out) :: criterion
  double precision, allocatable :: tmp_int(:,:)
  integer                       :: i,j,k,tmp_i,tmp_j,tmp_k, a, b, mu ,rho

  ! Allocation
  allocate(tmp_int(mo_num, mo_num))

  criterion = 0d0

  do a = 1, nucl_num ! loop over the nuclei
    tmp_int = 0d0

    do i = 1, mo_num
      do rho = 1, ao_num ! loop over all the AOs
        do b = 1, nucl_n_aos(a) ! loop over the number of AOs which belongs to the nuclei a
          mu = nucl_aos(a,b)

          tmp_int(i,i) = tmp_int(i,i) + 0.5d0 * (mo_coef(rho,i) * ao_overlap(rho,mu) * mo_coef(mu,i) &
                                 + mo_coef(mu,i) * ao_overlap(mu,rho) * mo_coef(rho,i))

        enddo
      enddo
    enddo

    do i = 1, mo_num
      criterion = criterion + tmp_int(i,i)**2
    enddo

  enddo

  criterion = - criterion

  deallocate(tmp_int)

end

! Criterion PM

! The criterion is computed as
! \begin{align*}
! \mathcal{P} = \sum_{i=1}^n \sum_{A=1}^N \left[ <i|P_A|i> \right]^2
! \end{align*}
! with
! \begin{align*}
! <s|P_A|t> = \frac{1}{2} \sum_{\rho} \sum_{\mu \in A} \left[ c_{\rho}^{s*} S_{\rho \nu} c_{\mu}^{t} +c_{\mu}^{s*} S_{\mu \rho} c_{\rho}^t \right]
! \end{align*}


subroutine criterion_PM(tmp_list_size,tmp_list,criterion)

  implicit none

  BEGIN_DOC
  ! Compute the Pipek-Mezey localization criterion
  END_DOC

  integer, intent(in)           :: tmp_list_size, tmp_list(tmp_list_size)
  double precision, intent(out) :: criterion
  double precision, allocatable :: tmp_int(:,:),CS(:,:)
  integer                       :: i,j,k,tmp_i,tmp_j,tmp_k, a, b, mu ,rho

  print*,''
  print*,'---criterion_PM---'

  ! Allocation
  allocate(tmp_int(tmp_list_size, tmp_list_size),CS(mo_num,ao_num))

  ! Initialization
  criterion = 0d0

  call dgemm('T','N',mo_num,ao_num,ao_num,1d0,mo_coef,size(mo_coef,1),ao_overlap,size(ao_overlap,1),0d0,CS,size(CS,1))

  do a = 1, nucl_num ! loop over the nuclei
    tmp_int = 0d0

      do tmp_i = 1, tmp_list_size
        i = tmp_list(tmp_i)
        do b = 1, nucl_n_aos(a) ! loop over the number of AOs which belongs to the nuclei a
          mu = nucl_aos(a,b)

          tmp_int(tmp_i,tmp_i) = tmp_int(tmp_i,tmp_i) + 0.5d0 * (CS(i,mu) * mo_coef(mu,i) + mo_coef(mu,i) * CS(i,mu))

                                 !  (mo_coef(rho,i) * ao_overlap(rho,mu) * mo_coef(mu,j) &
                                 !+ mo_coef(mu,i) * ao_overlap(mu,rho) * mo_coef(rho,j))

      enddo
    enddo

    do tmp_i = 1, tmp_list_size
      criterion = criterion + tmp_int(tmp_i,tmp_i)**2
    enddo

  enddo

  criterion = - criterion

  deallocate(tmp_int,CS)

  print*,'---End criterion_PM---'
  print*,''

end

! Criterion PM v3

subroutine criterion_PM_v3(tmp_list_size,tmp_list,criterion)

  implicit none

  BEGIN_DOC
  ! Compute the Pipek-Mezey localization criterion
  END_DOC

  integer, intent(in)           :: tmp_list_size, tmp_list(tmp_list_size)
  double precision, intent(out) :: criterion
  double precision, allocatable :: tmp_int(:,:), CS(:,:), tmp_mo_coef(:,:), tmp_mo_coef2(:,:),tmp_accu(:,:),tmp_CS(:,:)
  integer                       :: i,j,k,tmp_i,tmp_j,tmp_k, a, b, mu ,rho,nu,c
  double precision              :: t1,t2,t3

  print*,''
  print*,'---criterion_PM_v3---'

  call wall_time(t1)

  ! Allocation
  allocate(tmp_int(tmp_list_size, tmp_list_size),tmp_accu(tmp_list_size, tmp_list_size))
  allocate(CS(tmp_list_size,ao_num),tmp_mo_coef(ao_num,tmp_list_size))

  criterion = 0d0

  ! submatrix of the mo_coef
  do tmp_i = 1, tmp_list_size
    i = tmp_list(tmp_i)
    do j = 1, ao_num

      tmp_mo_coef(j,tmp_i) = mo_coef(j,i)

    enddo
  enddo

  ! ao_overlap(ao_num,ao_num)
  ! mo_coef(ao_num,mo_num)
  call dgemm('T','N',tmp_list_size,ao_num,ao_num,1d0,tmp_mo_coef,size(tmp_mo_coef,1),ao_overlap,size(ao_overlap,1),0d0,CS,size(CS,1))

  do a = 1, nucl_num ! loop over the nuclei

    do j = 1, tmp_list_size
      do i = 1, tmp_list_size
        tmp_int(i,j) = 0d0
      enddo
    enddo

    !do tmp_j = 1, tmp_list_size
    !  do tmp_i = 1, tmp_list_size
    !    do b = 1, nucl_n_aos(a) ! loop over the number of AOs which belongs to the nuclei a
    !      mu = nucl_aos(a,b)

    !      tmp_int(tmp_i,tmp_j) = tmp_int(tmp_i,tmp_j) + 0.5d0 * (CS(tmp_i,mu) * tmp_mo_coef(mu,tmp_j) + tmp_mo_coef(mu,tmp_i) * CS(tmp_j,mu))

    !                             !  (mo_coef(rho,i) * ao_overlap(rho,mu) * mo_coef(mu,j) &
    !                             !+ mo_coef(mu,i) * ao_overlap(mu,rho) * mo_coef(rho,j))

    !    enddo
    !  enddo
    !enddo

    allocate(tmp_mo_coef2(nucl_n_aos(a),tmp_list_size),tmp_CS(tmp_list_size,nucl_n_aos(a)))

    do tmp_i = 1, tmp_list_size
      do b = 1, nucl_n_aos(a)
        mu = nucl_aos(a,b)

        tmp_mo_coef2(b,tmp_i) = tmp_mo_coef(mu,tmp_i)

      enddo
    enddo

    do b = 1, nucl_n_aos(a)
      mu = nucl_aos(a,b)
      do tmp_i = 1, tmp_list_size

         tmp_CS(tmp_i,b) = CS(tmp_i,mu)

      enddo
    enddo

    call dgemm('N','N',tmp_list_size,tmp_list_size,nucl_n_aos(a),1d0,tmp_CS,size(tmp_CS,1),tmp_mo_coef2,size(tmp_mo_coef2,1),0d0,tmp_accu,size(tmp_accu,1))

    ! Integrals
    do tmp_j = 1, tmp_list_size
      do tmp_i = 1, tmp_list_size

        tmp_int(tmp_i,tmp_j) = 0.5d0 * (tmp_accu(tmp_i,tmp_j) + tmp_accu(tmp_j,tmp_i))

      enddo
    enddo

    deallocate(tmp_mo_coef2,tmp_CS)

    ! Criterion
    do tmp_i = 1, tmp_list_size
      criterion = criterion + tmp_int(tmp_i,tmp_i)**2
    enddo

  enddo

  criterion = - criterion

  deallocate(tmp_int,CS,tmp_accu,tmp_mo_coef)

  call wall_time(t2)
  t3 = t2 - t1
  print*,'Time in criterion_PM_v3:', t3

  print*,'---End criterion_PM_v3---'
  print*,''

end

! Criterion FB (old)

! The criterion is just computed as

! \begin{align*}
! C = - \sum_i^{mo_{num}} (<i|x|i>^2 + <i|y|i>^2 + <i|z|i>^2)
! \end{align*}

! The minus sign is here in order to minimize this criterion

! Output:
! | criterion | double precision | criterion for the Foster-Boys localization |


subroutine criterion_FB_old(criterion)

  implicit none

  BEGIN_DOC
  ! Compute the Foster-Boys localization criterion
  END_DOC

  double precision, intent(out) :: criterion
  integer                       :: i

  ! Criterion (= \sum_i <i|r|i>^2 )
  criterion = 0d0
  do i = 1, mo_num
    criterion = criterion + mo_dipole_x(i,i)**2 + mo_dipole_y(i,i)**2 + mo_dipole_z(i,i)**2
  enddo
  criterion = - criterion

end subroutine

! Criterion FB

subroutine criterion_FB(tmp_list_size, tmp_list, criterion)

  implicit none

  BEGIN_DOC
  ! Compute the Foster-Boys localization criterion
  END_DOC

  integer, intent(in)           :: tmp_list_size, tmp_list(tmp_list_size)
  double precision, intent(out) :: criterion
  integer                       :: i, tmp_i

  ! Criterion (= - \sum_i <i|r|i>^2 )
  criterion = 0d0
  do tmp_i = 1, tmp_list_size
    i = tmp_list(tmp_i)
    criterion = criterion + mo_dipole_x(i,i)**2 + mo_dipole_y(i,i)**2 + mo_dipole_z(i,i)**2
  enddo
  criterion = - criterion

end subroutine

subroutine theta_FB(l, n, m_x, max_elem)

  include 'pi.h'

  BEGIN_DOC
  ! Compute the angles to minimize the Foster-Boys criterion by using pairwise rotations of the MOs
  ! Warning: you must give - the angles to build the rotation matrix...
  END_DOC

  implicit none

  integer, intent(in)           :: n, l(n)
  double precision, intent(out) :: m_x(n,n), max_elem

  integer                       :: i,j, tmp_i, tmp_j
  double precision, allocatable :: cos4theta(:,:), sin4theta(:,:)
  double precision, allocatable :: A(:,:), B(:,:), beta(:,:), gamma(:,:)
  integer                       :: idx_i,idx_j

  allocate(cos4theta(n, n), sin4theta(n, n))
  allocate(A(n,n), B(n,n), beta(n,n), gamma(n,n))

  do tmp_j = 1, n
    j = l(tmp_j)
    do tmp_i = 1, n
      i = l(tmp_i)
      A(tmp_i,tmp_j) = mo_dipole_x(i,j)**2 - 0.25d0 * (mo_dipole_x(i,i) - mo_dipole_x(j,j))**2 &
                     + mo_dipole_y(i,j)**2 - 0.25d0 * (mo_dipole_y(i,i) - mo_dipole_y(j,j))**2 &
                     + mo_dipole_z(i,j)**2 - 0.25d0 * (mo_dipole_z(i,i) - mo_dipole_z(j,j))**2
    enddo
    A(j,j) = 0d0
  enddo

  do tmp_j = 1, n
    j = l(tmp_j)
    do tmp_i = 1, n
      i = l(tmp_i)
      B(tmp_i,tmp_j) = mo_dipole_x(i,j) * (mo_dipole_x(i,i) - mo_dipole_x(j,j)) &
                     + mo_dipole_y(i,j) * (mo_dipole_y(i,i) - mo_dipole_y(j,j)) &
                     + mo_dipole_z(i,j) * (mo_dipole_z(i,i) - mo_dipole_z(j,j))
    enddo
  enddo

  !do tmp_j = 1, n
  !  j = l(tmp_j)
  !  do tmp_i = 1, n
  !    i = l(tmp_i)
  !    beta(tmp_i,tmp_j) =  (mo_dipole_x(i,i) - mo_dipole_x(j,j)) - 4d0 * mo_dipole_x(i,j)**2 &
  !               + (mo_dipole_y(i,i) - mo_dipole_y(j,j)) - 4d0 * mo_dipole_y(i,j)**2 &
  !               + (mo_dipole_z(i,i) - mo_dipole_z(j,j)) - 4d0 * mo_dipole_z(i,j)**2
  !  enddo
  !enddo

  !do tmp_j = 1, n
  !  j = l(tmp_j)
  !  do tmp_i = 1, n
  !    i = l(tmp_i)
  !    gamma(tmp_i,tmp_j) = 4d0 * ( mo_dipole_x(i,j) * (mo_dipole_x(i,i) - mo_dipole_x(j,j)) &
  !                       + mo_dipole_y(i,j) * (mo_dipole_y(i,i) - mo_dipole_y(j,j)) &
  !                       + mo_dipole_z(i,j) * (mo_dipole_z(i,i) - mo_dipole_z(j,j)))
  !  enddo
  !enddo

  !
  !do j = 1, n
  !  do i = 1, n
  !    cos4theta(i,j) = - A(i,j) / dsqrt(A(i,j)**2 + B(i,j)**2)
  !  enddo
  !enddo

  !do j = 1, n
  !  do i = 1, n
  !    sin4theta(i,j) = B(i,j) / dsqrt(A(i,j)**2 + B(i,j)**2)
  !  enddo
  !enddo

  ! Theta
  do j = 1, n
    do i = 1, n
      m_x(i,j) = 0.25d0 * atan2(B(i,j), -A(i,j))
      !m_x(i,j) = 0.25d0 * atan2(sin4theta(i,j), cos4theta(i,j))
    enddo
  enddo

  ! Enforce a perfect antisymmetry
  do j = 1, n-1
    do i = j+1, n
      m_x(j,i) = - m_x(i,j)
    enddo
  enddo
  do i = 1, n
    m_x(i,i) = 0d0
  enddo

  ! Max
  max_elem = 0d0
  do j = 1, n-1
    do i = j+1, n
      if (dabs(m_x(i,j)) > dabs(max_elem)) then
        max_elem = m_x(i,j)
        !idx_i = i
        !idx_j = j
      endif
    enddo
  enddo

  ! Debug
  !print*,''
  !print*,'sin/B'
  !do i = 1, n
  !  write(*,'(100F10.4)') sin4theta(i,:)
  !  !B(i,:)
  !enddo
  !print*,'cos/A'
  !do i = 1, n
  !  write(*,'(100F10.4)') cos4theta(i,:)
  !  !A(i,:)
  !enddo
  !print*,'X'
  !!m_x = 0d0
  !!m_x(idx_i,idx_j) = max_elem
  !!m_x(idx_j,idx_i) = -max_elem
  !do i = 1, n
  !  write(*,'(100F10.4)') m_x(i,:)
  !enddo
  !print*,idx_i,idx_j,max_elem

  max_elem = dabs(max_elem)

  deallocate(cos4theta, sin4theta)
  deallocate(A,B,beta,gamma)

end

subroutine theta_PM(l, n, m_x, max_elem)

  include 'pi.h'

  BEGIN_DOC
  ! Compute the angles to minimize the Foster-Boys criterion by using pairwise rotations of the MOs
  ! Warning: you must give - the angles to build the rotation matrix...
  END_DOC

  implicit none

  integer, intent(in)           :: n, l(n)
  double precision, intent(out) :: m_x(n,n), max_elem

  integer                       :: a,b,i,j,tmp_i,tmp_j,rho,mu,nu,idx_i,idx_j
  double precision, allocatable :: Aij(:,:), Bij(:,:), Pa(:,:)

  allocate(Aij(n,n), Bij(n,n), Pa(n,n))

  do a = 1, nucl_num ! loop over the nuclei
    Pa = 0d0 ! Initialization for each nuclei

    ! Loop over the MOs of the a given mo_class to compute <i|P_a|j>
    do tmp_j = 1, n
      j = l(tmp_j)
      do tmp_i = 1, n
         i = l(tmp_i)
        do rho = 1, ao_num ! loop over all the AOs
          do b = 1, nucl_n_aos(a) ! loop over the number of AOs which belongs to the nuclei a
            mu = nucl_aos(a,b) ! AO centered on atom a

            Pa(tmp_i,tmp_j) = Pa(tmp_i,tmp_j) + 0.5d0 * (mo_coef(rho,i) * ao_overlap(rho,mu) * mo_coef(mu,j) &
                                   + mo_coef(mu,i) * ao_overlap(mu,rho) * mo_coef(rho,j))

          enddo
        enddo
      enddo
    enddo

    ! A
    do j = 1, n
      do i = 1, n
        Aij(i,j) = Aij(i,j) + Pa(i,j)**2 - 0.25d0 * (Pa(i,i) - Pa(j,j))**2
      enddo
    enddo

    ! B
    do j = 1, n
      do i = 1, n
        Bij(i,j) = Bij(i,j) + Pa(i,j) * (Pa(i,i) - Pa(j,j))
      enddo
    enddo

  enddo

  ! Theta
  do j = 1, n
    do i = 1, n
      m_x(i,j) = 0.25d0 * atan2(Bij(i,j), -Aij(i,j))
    enddo
  enddo

  ! Enforce a perfect antisymmetry
  do j = 1, n-1
    do i = j+1, n
      m_x(j,i) = - m_x(i,j)
    enddo
  enddo
  do i = 1, n
    m_x(i,i) = 0d0
  enddo

  ! Max
  max_elem = 0d0
  do j = 1, n-1
    do i = j+1, n
      if (dabs(m_x(i,j)) > dabs(max_elem)) then
        max_elem = m_x(i,j)
        idx_i = i
        idx_j = j
      endif
    enddo
  enddo

  ! Debug
  !do i = 1, n
  !  write(*,'(100F10.4)') m_x(i,:)
  !enddo
  !print*,'Max',idx_i,idx_j,max_elem

  max_elem = dabs(max_elem)

  deallocate(Aij,Bij,Pa)

end

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
  integer, allocatable           :: iorder(:)
  double precision, allocatable  :: fock_energies_tmp(:), tmp_mo_coef(:,:)
  integer, allocatable           :: tmp_list(:)

!  allocate(iorder(mo_num), fock_energies_tmp(mo_num), new_mo_coef(ao_num, mo_num))
!
!  do i = 1, mo_num
!    fock_energies_tmp(i) = Fock_matrix_diag_mo(i)
!    print*,'fock_energies_tmp(i) = ',fock_energies_tmp(i)
!    iorder(i) = i
!  enddo
!
!  print*,''
!  print*,'Sorting by Fock energies'
!  print*,''
!
!  call dsort(fock_energies_tmp, iorder, mo_num)
!
!  do i = 1, mo_num
!    k = iorder(i)
!    print*,'fock_energies_new(i) = ',fock_energies_tmp(i)
!    do j = 1, ao_num
!      new_mo_coef(j,i) = mo_coef(j,k)
!    enddo
!  enddo

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
      print*,'MO class: ', trim(mo_class(tmp_list(1)))

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
