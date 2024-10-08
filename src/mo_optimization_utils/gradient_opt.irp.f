! Gradient

! The gradient of the CI energy with respects to the orbital rotation
! is:
! (C-c C-x C-l)
! $$
! G(p,q) = \mathcal{P}_{pq} \left[ \sum_r (h_p^r \gamma_r^q - h_r^q \gamma_p^r) +
! \sum_{rst}(v_{pt}^{rs} \Gamma_{rs}^{qt} - v_{rs}^{qt} \Gamma_{pt}^{rs})
! \right]
! $$


! $$
! \mathcal{P}_{pq}= 1 - (p \leftrightarrow q)
! $$

! $$
! G(p,q) = \left[
! \sum_r (h_p^r \gamma_r^q - h_r^q \gamma_p^r) +
! \sum_{rst}(v_{pt}^{rs} \Gamma_{rs}^{qt} - v_{rs}^{qt} \Gamma_{pt}^{rs})
! \right] - 
! \left[
! \sum_r (h_q^r \gamma_r^p - h_r^p \gamma_q^r) +
! \sum_{rst}(v_{qt}^{rs} \Gamma_{rs}^{pt} - v_{rs}^{pt}
! \Gamma_{qt}^{rs})
! \right]
! $$

! Where p,q,r,s,t are general spatial orbitals
! mo_num : the number of molecular orbitals
! $$h$$ : One electron integrals
! $$\gamma$$ : One body density matrix (state average in our case)
! $$v$$ : Two electron integrals
! $$\Gamma$$ : Two body density matrice (state average in our case)

! The gradient is a mo_num by mo_num matrix, p,q,r,s,t take all the
! values between 1 and mo_num (1 and mo_num include).

! To do that we compute $$G(p,q)$$ for all the pairs (p,q).

! Source :
! Seniority-based coupled cluster theory
! J. Chem. Phys. 141, 244104 (2014); https://doi.org/10.1063/1.4904384
! Thomas M. Henderson, Ireneusz W. Bulik, Tamar Stein, and Gustavo
! E. Scuseria

! *Compute the gradient of energy with respects to orbital rotations*

! Provided:
! | mo_num                                   | integer          | number of MOs             |
! | mo_one_e_integrals(mo_num,mo_num)        | double precision | mono_electronic integrals |
! | one_e_dm_mo(mo_num,mo_num)               | double precision | one e- density matrix     |
! | two_e_dm_mo(mo_num,mo_num,mo_num,mo_num) | double precision | two e- density matrix     |

! Input:
! | n | integer | mo_num*(mo_num-1)/2 |

! Output:
! | v_grad(n) | double precision | the gradient                    |
! | max_elem  | double precision | maximum element of the gradient |

! Internal:
! | grad(mo_num,mo_num)                | double precison   | gradient before the tranformation in a vector             |
! | A((mo_num,mo_num)                  | doubre precision  | gradient after the permutations                           |
! | norm                               | double precision  | norm of the gradient                                      |
! | p, q                               | integer           | indexes of the element in the matrix grad                 |
! | i                                  | integer           | index for the tranformation in a vector                   |
! | r, s, t                            | integer           | indexes dor the sums                                      |
! | t1, t2, t3                         | double precision  | t3 = t2 - t1, time to compute the gradient                |
! | t4, t5, t6                         | double precission | t6 = t5 - t4, time to compute each element                |
! | tmp_bi_int_3(mo_num,mo_num,mo_num) | double precision  | 3 indexes temporary array for the bi-electronic integrals |
! | tmp_2rdm_3(mo_num,mo_num,mo_num)   | double precision  | 3 indexes temporary array for the two e- density matrix   |
! | tmp_accu(mo_num,mo_num)            | double precision  | temporary array                                           |

! Function:
! | get_two_e_integral | double precision | bi-electronic integrals |
! | dnrm2              | double precision | (Lapack) norm           |


subroutine gradient_opt(n,v_grad,max_elem)
  use omp_lib
  include 'constants.h'

  implicit none
  
  ! Variables
  
  ! in
  integer, intent(in) :: n
  
  ! out
  double precision, intent(out) :: v_grad(n), max_elem

  ! internal
  double precision, allocatable :: grad(:,:),A(:,:)
  double precision              :: norm
  integer                       :: i,p,q,r,s,t
  double precision              :: t1,t2,t3,t4,t5,t6

  double precision, allocatable :: tmp_accu(:,:)
  double precision, allocatable :: tmp_bi_int_3(:,:,:), tmp_2rdm_3(:,:,:)

  ! Functions
  double precision :: get_two_e_integral, dnrm2

  
  print*,''
  print*,'---gradient---'

  ! Allocation of shared arrays
  allocate(grad(mo_num,mo_num),A(mo_num,mo_num)) 

  ! Initialization omp
  call omp_set_max_active_levels(1)

  !$OMP PARALLEL                                                 &
      !$OMP PRIVATE(                                             &
      !$OMP   p,q,r,s,t,                                         &
      !$OMP   tmp_accu, tmp_bi_int_3, tmp_2rdm_3)                &
      !$OMP SHARED(grad, one_e_dm_mo, mo_num,mo_one_e_integrals, &
      !$OMP mo_integrals_map,t4,t5,t6)                           &
      !$OMP DEFAULT(SHARED)
 
  ! Allocation of private arrays
  allocate(tmp_accu(mo_num,mo_num))
  allocate(tmp_bi_int_3(mo_num,mo_num,mo_num))
  allocate(tmp_2rdm_3(mo_num,mo_num,mo_num))

! Initialization

!$OMP DO
do q = 1, mo_num
  do p = 1,mo_num
    grad(p,q) = 0d0
  enddo
enddo
!$OMP END DO

! Term 1
  
! Without optimization the term 1 is :

! do p = 1, mo_num
!   do q = 1, mo_num
!      do r = 1, mo_num
!        grad(p,q) = grad(p,q) &
!                + mo_one_e_integrals(p,r) * one_e_dm_mo(r,q) &
!                - mo_one_e_integrals(r,q) * one_e_dm_mo(p,r)
!     enddo
!   enddo
! enddo
   
! Since the matrix multiplication A.B is defined like :
! \begin{equation}
! c_{ij} = \sum_k a_{ik}.b_{kj}
! \end{equation}
! The previous equation can be rewritten as a matrix multplication  
  

!****************
! Opt first term
!****************

!$OMP MASTER
CALL wall_TIME(t4)
!$OMP END MASTER

call dgemm('N','N',mo_num,mo_num,mo_num,1d0,mo_one_e_integrals,&
mo_num,one_e_dm_mo,mo_num,0d0,tmp_accu,mo_num)

!$OMP DO
do q = 1, mo_num
  do p = 1, mo_num

    grad(p,q) = grad(p,q) + (tmp_accu(p,q) - tmp_accu(q,p))

  enddo
enddo 
!$OMP END DO

!$OMP MASTER
CALL wall_TIME(t5)
t6 = t5-t4
print*,'Gradient, first term (s) :', t6 
!$OMP END MASTER

! Term 2
 
! Without optimization the second term is : 

! do p = 1, mo_num
!   do q = 1, mo_num 
!     do r = 1, mo_num
!       do s = 1, mo_num
!         do t= 1, mo_num

!         grad(p,q) = grad(p,q) &
!                 + get_two_e_integral(p,t,r,s,mo_integrals_map) * two_e_dm_mo(r,s,q,t) &
!                 - get_two_e_integral(r,s,q,t,mo_integrals_map) * two_e_dm_mo(p,t,r,s)
!        enddo
!       enddo
!     enddo
!   enddo
! enddo

! Using the bielectronic integral properties :
! get_two_e_integral(p,t,r,s,mo_integrals_map) = get_two_e_integral(r,s,p,t,mo_integrals_map)

! Using the two body matrix properties :
! two_e_dm_mo(p,t,r,s) = two_e_dm_mo(r,s,p,t)

! t is one the right, we can put it on the external loop and create 3
! indexes temporary array 
! r,s can be seen as one index

! By doing so, a matrix multiplication appears 


!*****************
! Opt second term  
!*****************

!$OMP MASTER
CALL wall_TIME(t4)
!$OMP END MASTER 

!$OMP DO
do t = 1, mo_num

  do p = 1, mo_num
    do s = 1, mo_num
      do r = 1, mo_num

        tmp_bi_int_3(r,s,p) = get_two_e_integral(r,s,p,t,mo_integrals_map)

      enddo
    enddo
  enddo

  do q = 1, mo_num
    do s = 1, mo_num
      do r = 1, mo_num

         tmp_2rdm_3(r,s,q) = two_e_dm_mo(r,s,q,t)

      enddo
    enddo
  enddo

  call dgemm('T','N',mo_num,mo_num,mo_num*mo_num,1d0,tmp_bi_int_3,&
    mo_num*mo_num,tmp_2rdm_3,mo_num*mo_num,0d0,tmp_accu,mo_num)

  !$OMP CRITICAL   
  do q = 1, mo_num
    do p = 1, mo_num

      grad(p,q) = grad(p,q) + tmp_accu(p,q) - tmp_accu(q,p)

    enddo
  enddo
  !$OMP END CRITICAL

enddo
!$OMP END DO

!$OMP MASTER
CALL wall_TIME(t5)
t6 = t5-t4
print*,'Gradient second term (s) : ', t6
!$OMP END MASTER

! Deallocation of private arrays

deallocate(tmp_bi_int_3,tmp_2rdm_3,tmp_accu)

!$OMP END PARALLEL

call omp_set_max_active_levels(4)

! Permutation, 2D matrix -> vector, transformation
! In addition there is a permutation in the gradient formula :
! \begin{equation}
! P_{pq} = 1 - (p <-> q) 
! \end{equation}

! We need a vector to use the gradient. Here the gradient is a
! antisymmetric matrix so we can transform it in a vector of length
! mo_num*(mo_num-1)/2.

! Here we do these two things at the same time.


do i=1,n
  call vec_to_mat_index(i,p,q)
  v_grad(i)=(grad(p,q) - grad(q,p))
enddo  

! Debug, diplay the vector containing the gradient elements 
if (debug) then  
  print*,'Vector containing the gradient :'
  write(*,'(100(F10.5))') v_grad(1:n)
endif

! Norm of the gradient
! The norm can be useful.

norm = dnrm2(n,v_grad,1)
print*, 'Gradient norm : ', norm

! Maximum element in the gradient
! The maximum element in the gradient is very important for the
! convergence criterion of the Newton method.


! Max element of the gradient
max_elem = 0d0
do i = 1, n
  if (ABS(v_grad(i)) > ABS(max_elem)) then
    max_elem = v_grad(i)
  endif
enddo

print*,'Max element in the gradient :', max_elem  

! Debug, display the matrix containting the gradient elements
if (debug) then
  ! Matrix gradient
  A = 0d0
  do q=1,mo_num
    do p=1,mo_num
      A(p,q) = grad(p,q) - grad(q,p)
    enddo
  enddo
  print*,'Matrix containing the gradient :'
  do i = 1, mo_num
    write(*,'(100(F10.5))') A(i,1:mo_num)
  enddo
endif

! Deallocation of shared arrays and end

deallocate(grad,A)

print*,'---End gradient---'

end subroutine
