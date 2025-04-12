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

  print*,'---End gradient_PM---'

end

! Hessian v1

subroutine hess_pipek(tmp_n, tmp_list_size, tmp_list, H)

  implicit none

  BEGIN_DOC
  ! Compute diagonal hessian for the Pipek-Mezey localization
  END_DOC

  integer, intent(in)           :: tmp_n, tmp_list_size, tmp_list(tmp_list_size)
  double precision, intent(out) :: H(tmp_n)
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
    H(tmp_k) = 4d0 * beta(tmp_i, tmp_j)
  enddo
  
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
  double precision, intent(out) :: H(tmp_n)
  double precision, allocatable :: beta(:,:),tmp_int(:,:),CS(:,:),tmp_mo_coef(:,:),tmp_mo_coef2(:,:),tmp_accu(:,:),tmp_CS(:,:)
  integer                       :: i,j,tmp_k,tmp_i, tmp_j, a,b,rho,mu
  double precision              :: max_elem, t1,t2,t3
    
  print*,''
  print*,'---hessian_PM---'

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
    H(tmp_k) = 4d0 * beta(tmp_i, tmp_j)
  enddo
  
  ! Deallocation
  deallocate(beta,tmp_int)

  call wall_time(t2)
  t3 = t2 - t1
  print*,'Time in hessian_PM:', t3

  print*,'---End hessian_PM---'

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

end

subroutine theta_PM(l, n, m_x, max_elem)
  
  include 'pi.h'

  BEGIN_DOC
  ! Compute the angles to minimize the Pipek-Mezey criterion by using pairwise rotations of the MOs
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

