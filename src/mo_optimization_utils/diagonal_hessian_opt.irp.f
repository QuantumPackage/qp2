! Diagonal hessian

! The hessian of the CI energy with respects to the orbital rotation is :
! (C-c C-x C-l)

! \begin{align*}
! H_{pq,rs} &= \dfrac{\partial^2 E(x)}{\partial x_{pq}^2} \\
!   &= \mathcal{P}_{pq} \mathcal{P}_{rs} [ \frac{1}{2} \sum_u [\delta_{qr}(h_p^u \gamma_u^s + h_u^s \gamma_p^u) 
!   + \delta_{ps}(h_r^u \gamma_u^q + h_u^q \gamma_r^u)]
!   -(h_p^s \gamma_r^q + h_r^q \gamma_p^s) \\
!   &+ \frac{1}{2} \sum_{tuv} [\delta_{qr}(v_{pt}^{uv} \Gamma_{uv}^{st} + v_{uv}^{st} \Gamma_{pt}^{uv})
!   + \delta_{ps}(v_{uv}^{qt} \Gamma_{rt}^{uv} + v_{rt}^{uv}\Gamma_{uv}^{qt})] \\
!   &+ \sum_{uv} (v_{pr}^{uv} \Gamma_{uv}^{qs} + v_{uv}^{qs}  \Gamma_{pr}^{uv}) 
!   - \sum_{tu} (v_{pu}^{st} \Gamma_{rt}^{qu}+v_{pu}^{tr} \Gamma_{tr}^{qu}+v_{rt}^{qu}\Gamma_{pu}^{st} + v_{tr}^{qu}\Gamma_{pu}^{ts}) 
! \end{align*}
! With pq a permutation operator :

! \begin{align*}
! \mathcal{P}_{pq}= 1 - (p \leftrightarrow q)
! \end{align*}
! \begin{align*}
! \mathcal{P}_{pq} \mathcal{P}_{rs} &= (1 - (p \leftrightarrow q))(1 - (r \leftrightarrow s)) \\
! &= 1 - (p \leftrightarrow q) - (r \leftrightarrow s) + (p \leftrightarrow q, r \leftrightarrow s)
! \end{align*}

! Where p,q,r,s,t,u,v are general spatial orbitals
! mo_num : the number of molecular orbitals
! $$h$$ : One electron integrals
! $$\gamma$$ : One body density matrix (state average in our case)
! $$v$$ : Two electron integrals
! $$\Gamma$$ : Two body density matrice (state average in our case)

! Source :
! Seniority-based coupled cluster theory
! J. Chem. Phys. 141, 244104 (2014); https://doi.org/10.1063/1.4904384
! Thomas M. Henderson, Ireneusz W. Bulik, Tamar Stein, and Gustavo E. Scuseria

! Here for the diagonal of the hessian it's a little more complicated
! than for the hessian. It's not just compute the diagonal terms of the
! hessian because of the permutations.

! The hessian is (p,q,r,s), so the diagonal terms are (p,q,p,q). But
! with the permutations : p <-> q, r <-> s, p <-> q and r <-> s, we have
! a diagonal term, if : 
! p = r and q = s, => (p,q,p,q)  
! or
! q = r and p = s, => (p,q,q,p) 

! For that reason, we will use 2D temporary arrays to store the
! elements. One for the terms (p,q,p,q) and an other for the terms of
! kind (p,q,q,p). We will also use a 1D temporary array to store the
! terms of the kind (p,p,p,p) due to the kronoecker delta.

! *Compute the diagonal hessian of energy with respects to orbital
! rotations*
! By diagonal hessian we mean, diagonal elements of the hessian

! Provided:
! | mo_num                            | integer          | number of MOs                         |
! | mo_one_e_integrals(mo_num,mo_num) | double precision | mono-electronic integrals             |
! | one_e_dm_mo(mo_num,mo_num)        | double precision | one e- density matrix (state average) |
! | two_e_dm_mo(mo_num,mo_num,mo_num) | double precision | two e- density matrix (state average) |

! Input:
! | n | integer | mo_num*(mo_num-1)/2 |

! Output:
! | H(n,n)                              | double precision | Hessian matrix                                   |
! | h_tmpr(mo_num,mo_num,mo_num,mo_num) | double precision | Complete hessian matrix before the tranformation |
! |                                     |                  | in n by n matrix                                 |

! Internal:
! | hessian(mo_num,mo_num,mo_num,mo_num)      | double precision | temporary array containing the hessian before                      |
! |                                           |                  | the permutations                                                   |
! | p, q, r, s                                | integer          | indexes of the hessian elements                                    |
! | t, u, v                                   | integer          | indexes for the sums                                               |
! | pq, rs                                    | integer          | indexes for the transformation of the hessian                      |
! |                                           |                  | (4D -> 2D)                                                         |
! | t1,t2,t3                                  | double precision | time to compute the hessian                                        |
! | t4,t5,t6                                  | double precision | time to compute the differ  each element                           |
! | tmp_bi_int_3(mo_num,mo_num,mo_num)        | double precision | 3 indexes temporary array for the bielectronic integrals (private) |
! | tmp_bi_int_3_shared(mo_num,mo_num,mo_num) | double precision | 3 indexes temporary array for the bielectronic integrals (shared)  |
! | tmp_2rdm_3(mo_num,mo_num,mo_num)          | double precision | 3 indexes temporary array for the 2 body density matrix (private)  |
! | tmp_2rdm_3_shared(mo_num,mo_num,mo_num)   | double precision | 3 indexes temporary array for the 2 body density matrix (shared)   |
! | tmp_accu(mo_num,mo_num)                   | double precision | temporary array (private)                                          |
! | tmp_accu_shared(mo_num,mo_num)            | double precision | temporary array (shared)                                           |
! | tmp_accu_1(mo_num)                        | double precision | temporary array (private)                                          |
! | tmp_accu_1_shared(mo_num)                 | double precision | temporary array (shared)                                           |
! | tmp_h_pppp(mo_num)                        | double precision | matrix containing the hessien elements hessian(p,p,p,p)            |
! | tmp_h_pqpq(mo_num,mo_num)                 | double precision | matrix containing the hessien elements hessian(p,q,p,q)            |
! | tmp_h_pqqp(mo_num,mo_num)                 | double precision | matrix containing the hessien elements hessian(p,q,q,p)            |

! Function:
! | get_two_e_integral | double precision | bi-electronic integrals |


subroutine diag_hessian_opt(n,H)!, h_tmpr)
 
  use omp_lib
  
  include 'constants.h' 

  implicit none

  ! Variables

  ! in
  integer, intent(in)           :: n 
 
  ! out
  double precision, intent(out) :: H(n)!,n),  h_tmpr(mo_num,mo_num,mo_num,mo_num)
  
  ! internal
  !double precision, allocatable :: hessian(:,:,:,:)!, h_tmpr(:,:,:,:)
  integer                       :: p,q,k
  integer                       :: r,s,t,u,v
  integer                       :: pq,rs
  integer                       :: istate
  double precision              :: t1,t2,t3,t4,t5,t6
  double precision, allocatable :: tmp_bi_int_3(:,:,:),tmp_bi_int_3_shared(:,:,:)
  double precision, allocatable :: tmp_2rdm_3(:,:,:),tmp_2rdm_3_shared(:,:,:)
  double precision, allocatable :: tmp_accu(:,:)
  double precision, allocatable :: tmp_accu_shared(:,:), tmp_accu_1_shared(:)
  double precision, allocatable :: tmp_h_pppp(:), tmp_h_pqpq(:,:), tmp_h_pqqp(:,:)
 
  ! Function
  double precision :: get_two_e_integral

  print*,''
  print*,'--- diagonal hessian---'
  print*,'Use the diagonal hessian'

  ! Allocation of shared arrays
  !allocate(hessian(mo_num,mo_num,mo_num,mo_num))!,h_tmpr(mo_num,mo_num,mo_num,mo_num))
  allocate(tmp_h_pppp(mo_num),tmp_h_pqpq(mo_num,mo_num),tmp_h_pqqp(mo_num,mo_num))
  allocate(tmp_2rdm_3_shared(mo_num,mo_num,mo_num))
  allocate(tmp_bi_int_3_shared(mo_num,mo_num,mo_num))
  allocate(tmp_accu_1_shared(mo_num),tmp_accu_shared(mo_num,mo_num))

  ! OMP
  call omp_set_max_active_levels(1)

  !$OMP PARALLEL                                                               &
      !$OMP PRIVATE(                                                           &
      !$OMP   p,q,r,s, tmp_accu,k,                                               &
      !$OMP   u,v,t, tmp_bi_int_3, tmp_2rdm_3)                                 &
      !$OMP SHARED(H, tmp_h_pppp, tmp_h_pqpq, tmp_h_pqqp,      &
      !$OMP  mo_num,n, mo_one_e_integrals, one_e_dm_mo,                        &
      !$OMP  tmp_bi_int_3_shared, tmp_2rdm_3_shared,tmp_accu_shared,           &
      !$OMP  tmp_accu_1_shared,two_e_dm_mo,mo_integrals_map,t1,t2,t3,t4,t5,t6) &
      !$OMP DEFAULT(NONE)

  ! Allocation of the private arrays
  allocate(tmp_2rdm_3(mo_num,mo_num,mo_num),tmp_bi_int_3(mo_num,mo_num,mo_num))
  allocate(tmp_accu(mo_num,mo_num))

! Initialization of the arrays

!!$OMP DO
!do s = 1,mo_num
!  do r = 1, mo_num
!    do q = 1, mo_num
!      do p = 1, mo_num
!        hessian(p,q,r,s) = 0d0
!      enddo
!    enddo
!  enddo
!enddo
!!$OMP END DO

!$OMP DO
do p = 1, mo_num
  tmp_h_pppp(p) = 0d0
enddo
!$OMP END DO

!$OMP DO
do q = 1, mo_num
  do p = 1, mo_num
    tmp_h_pqpq(p,q) = 0d0
  enddo
enddo
!$OMP END DO

!$OMP DO
do q = 1, mo_num
  do p = 1, mo_num
    tmp_h_pqqp(p,q) = 0d0
  enddo
enddo
!$OMP END DO

!$OMP MASTER
CALL wall_TIME(t1)
!$OMP END MASTER

! Line 1, term 1

! \begin{align*}
! \frac{1}{2} \sum_u \delta_{qr}(h_p^u \gamma_u^s + h_u^s \gamma_p^u)
! \end{align*}

! Without optimization :

! do p = 1, mo_num
!   do q = 1, mo_num
!     do r = 1, mo_num
!       do s = 1, mo_num

!         ! Permutations 
!         if (((p==r) .and. (q==s)) .or. ((q==r) .and. (p==s))) then
       
!           if (q==r) then
!             do u = 1, mo_num

!               hessian(p,q,r,s) = hessian(p,q,r,s) + 0.5d0 * (  &
!                 mo_one_e_integrals(u,p) * one_e_dm_mo(u,s) &
!               + mo_one_e_integrals(s,u) * one_e_dm_mo(p,u))

!             enddo
!           endif
!         endif
!       enddo
!     enddo
!   enddo
! enddo

! With optimization :

! *Part 1 : p=r and q=s and q=r*

!  hessian(p,q,r,s) -> hessian(p,p,p,p)

!   0.5d0 * (  &
!   mo_one_e_integrals(u,p) * one_e_dm_mo(u,s) &
! + mo_one_e_integrals(s,u) * one_e_dm_mo(p,u))
!   =  
!   0.5d0 * (  &
!   mo_one_e_integrals(u,p) * one_e_dm_mo(u,p) &
! + mo_one_e_integrals(p,u) * one_e_dm_mo(p,u))
!  =  
!   mo_one_e_integrals(u,p) * one_e_dm_mo(u,p) 


!$OMP MASTER
CALL wall_TIME(t4) 
!$OMP END MASTER

!$OMP DO
do p = 1, mo_num
  tmp_accu_1_shared(p) = 0d0
enddo
!$OMP END DO

!$OMP DO
do p = 1, mo_num
  do u = 1, mo_num

    tmp_accu_1_shared(p) = tmp_accu_1_shared(p) + mo_one_e_integrals(u,p) * one_e_dm_mo(u,p)

  enddo
enddo
!$OMP END DO

!$OMP DO
do p = 1, mo_num
  tmp_h_pppp(p) = tmp_h_pppp(p) + tmp_accu_1_shared(p)
enddo
!$OMP END DO


  
! *Part 2 : q=r and p=s and q=r*

!  hessian(p,q,r,s) -> hessian(p,q,q,p)
   
!   0.5d0 * (  &
!   mo_one_e_integrals(u,p) * one_e_dm_mo(u,s) &
! + mo_one_e_integrals(s,u) * one_e_dm_mo(p,u))
!   =  
!   0.5d0 * (  &
!   mo_one_e_integrals(u,p) * one_e_dm_mo(u,p) &
! + mo_one_e_integrals(p,u) * one_e_dm_mo(p,u))
!  =  
!   mo_one_e_integrals(u,p) * one_e_dm_mo(u,p)    


!$OMP DO
do p = 1, mo_num
  tmp_accu_1_shared(p) = 0d0
enddo
!$OMP END DO

!$OMP DO
do p = 1, mo_num
  do u = 1, mo_num

    tmp_accu_1_shared(p) = tmp_accu_1_shared(p) + mo_one_e_integrals(u,p) * one_e_dm_mo(u,p)

  enddo
enddo
!$OMP END DO

!$OMP DO
do q = 1, mo_num
  do p = 1, mo_num

    tmp_h_pqqp(p,q) = tmp_h_pqqp(p,q) + tmp_accu_1_shared(p)

  enddo
enddo
!$OMP END DO

!$OMP MASTER
CALL wall_TIME(t5)
t6= t5-t4
print*,'l1 1',t6
!$OMP END MASTER

! Line 1, term 2

! \begin{align*}
! \frac{1}{2} \sum_u \delta_{ps}(h_r^u \gamma_u^q + h_u^q \gamma_r^u)
! \end{align*}

! Without optimization :

! do p = 1, mo_num
!   do q = 1, mo_num
!     do r = 1, mo_num
!       do s = 1, mo_num

!         ! Permutations 
!         if (((p==r) .and. (q==s)) .or. ((q==r) .and. (p==s))) then
       
!           if (p==s) then
!             do u = 1, mo_num

!                 hessian(p,q,r,s) = hessian(p,q,r,s) + 0.5d0 * ( &
!                   mo_one_e_integrals(u,r) * one_e_dm_mo(u,q) &
!                 + mo_one_e_integrals(q,u) * one_e_dm_mo(r,u))
!             enddo
!           endif
!         endif
!       enddo
!     enddo
!   enddo
! enddo

! *Part 1 : p=r and q=s and p=s*

!  hessian(p,q,r,s) -> hessian(p,p,p,p)

!  0.5d0 * (&
!   mo_one_e_integrals(u,r) * one_e_dm_mo(u,q) &
! + mo_one_e_integrals(q,u) * one_e_dm_mo(r,u))
!  =
!  0.5d0 * ( &
!   mo_one_e_integrals(u,p) * one_e_dm_mo(u,p) &
! + mo_one_e_integrals(p,u) * one_e_dm_mo(p,u))
!  = 
!   mo_one_e_integrals(u,p) * one_e_dm_mo(u,p)


!$OMP MASTER
CALL wall_TIME(t4)
!$OMP END MASTER  

!$OMP DO
do p = 1, mo_num
  tmp_accu_1_shared(p) = 0d0 
enddo
!$OMP END DO

!$OMP DO
do p = 1, mo_num
  do u = 1, mo_num

    tmp_accu_1_shared(p) = tmp_accu_1_shared(p) +  mo_one_e_integrals(u,p) * one_e_dm_mo(u,p) 

  enddo
enddo
!$OMP END DO

!$OMP DO
do p = 1, mo_num

  tmp_h_pppp(p) = tmp_h_pppp(p) + tmp_accu_1_shared(p)

enddo
!$OMP END DO



! *Part 2 : q=r and p=s and p=s*

!  hessian(p,q,r,s) -> hessian(p,q,q,p)

!  0.5d0 * (&
!   mo_one_e_integrals(u,r) * one_e_dm_mo(u,q) &
! + mo_one_e_integrals(q,u) * one_e_dm_mo(r,u))
!  =
!  0.5d0 * ( &
!   mo_one_e_integrals(u,q) * one_e_dm_mo(u,q) &
! + mo_one_e_integrals(q,u) * one_e_dm_mo(q,u))
!  = 
!   mo_one_e_integrals(u,q) * one_e_dm_mo(u,q)


!$OMP DO
do p = 1, mo_num
  tmp_accu_1_shared(p) = 0d0
enddo
!$OMP END DO

!$OMP DO
do q = 1, mo_num
  do u = 1, mo_num

    tmp_accu_1_shared(q) = tmp_accu_1_shared(q) + mo_one_e_integrals(u,q) * one_e_dm_mo(u,q)

  enddo
enddo
!$OMP END DO

!$OMP DO
do q = 1, mo_num
  do p = 1, mo_num

    tmp_h_pqqp(p,q) = tmp_h_pqqp(p,q) + tmp_accu_1_shared(q)

  enddo
enddo
!$OMP END DO

!$OMP MASTER
CALL wall_TIME(t5)
t6= t5-t4
print*,'l1 2',t6
!$OMP END MASTER

! Line 1, term 3

! \begin{align*}
! -(h_p^s \gamma_r^q + h_r^q \gamma_p^s)
! \end{align*}

! Without optimization :

! do p = 1, mo_num
!   do q = 1, mo_num
!     do r = 1, mo_num
!       do s = 1, mo_num

!         ! Permutations 
!         if (((p==r) .and. (q==s)) .or. ((q==r) .and. (p==s))) then

!           hessian(p,q,r,s) = hessian(p,q,r,s) &
!            - mo_one_e_integrals(s,p) * one_e_rdm_mo(r,q) &
!            - mo_one_e_integrals(q,r) * one_e_rdm_mo(p,s)

!         endif
!       enddo
!     enddo
!   enddo
! enddo
 
! With optimization :

! *Part 1 : p=r and q=s*

! hessian(p,q,r,s) -> hessian(p,q,p,q)

!  - mo_one_e_integrals(s,p) * one_e_dm_mo(r,q) &
!  - mo_one_e_integrals(q,r) * one_e_dm_mo(p,s)
! =
!  - mo_one_e_integrals(q,p) * one_e_dm_mo(p,q) &
!  - mo_one_e_integrals(q,p) * one_e_dm_mo(p,q) 
! = 
!  - 2d0 mo_one_e_integrals(q,p) * one_e_dm_mo(p,q)


!$OMP MASTER
CALL wall_TIME(t4)
!$OMP END MASTER 

!$OMP DO
do q = 1, mo_num
  do p = 1, mo_num

    tmp_h_pqpq(p,q) = tmp_h_pqpq(p,q) &
      - 2d0 * mo_one_e_integrals(q,p) * one_e_dm_mo(p,q)

  enddo
enddo
!$OMP END DO



! *Part 2 : q=r and p=s*

! hessian(p,q,r,s) -> hessian(p,q,p,q)
 
!  - mo_one_e_integrals(s,p) * one_e_dm_mo(r,q) &
!  - mo_one_e_integrals(q,r) * one_e_dm_mo(p,s)
! =
!  - mo_one_e_integrals(q,p) * one_e_dm_mo(p,q) &
!  - mo_one_e_integrals(q,p) * one_e_dm_mo(p,q) 
! = 
!  - 2d0 mo_one_e_integrals(q,p) * one_e_dm_mo(p,q)


!$OMP DO
do q = 1, mo_num
  do p = 1, mo_num

    tmp_h_pqqp(p,q) = tmp_h_pqqp(p,q) &
      - 2d0 * mo_one_e_integrals(p,p) * one_e_dm_mo(q,q)

  enddo
enddo
!$OMP END DO

!$OMP MASTER
CALL wall_TIME(t5)
t6= t5-t4
print*,'l1 3',t6
!$OMP END MASTER

! Line 2, term 1

! \begin{align*}
! \frac{1}{2} \sum_{tuv} \delta_{qr}(v_{pt}^{uv} \Gamma_{uv}^{st} + v_{uv}^{st} \Gamma_{pt}^{uv})
! \end{align*}

! Without optimization :

! do p = 1, mo_num
!   do q = 1, mo_num
!     do r = 1, mo_num
!       do s = 1, mo_num

!         ! Permutations 
!         if (((p==r) .and. (q==s)) .or. ((q==r) .and. (p==s))) then

!           if (q==r) then
!             do t = 1, mo_num
!               do u = 1, mo_num
!                 do v = 1, mo_num

!                    hessian(p,q,r,s) = hessian(p,q,r,s) + 0.5d0 * (  &
!                      get_two_e_integral(u,v,p,t,mo_integrals_map) * two_e_dm_mo(u,v,s,t) &
!                    + get_two_e_integral(s,t,u,v,mo_integrals_map) * two_e_dm_mo(p,t,u,v))

!                 enddo
!               enddo
!             enddo
!           endif
!         endif
!       enddo
!     enddo
!   enddo
! enddo

! With optimization :

! *Part 1 : p=r and q=s and q=r*

!  hessian(p,q,r,s) -> hessian(p,p,p,p)
 
!  0.5d0 * (  &
!   get_two_e_integral(u,v,p,t,mo_integrals_map) * two_e_dm_mo(u,v,s,t) &
! + get_two_e_integral(s,t,u,v,mo_integrals_map) * two_e_dm_mo(p,t,u,v))
!  =
!  0.5d0 * (  &
!   get_two_e_integral(u,v,p,t,mo_integrals_map) * two_e_dm_mo(u,v,p,t) &
! + get_two_e_integral(p,t,u,v,mo_integrals_map) * two_e_dm_mo(p,t,u,v))
!  = 
!  get_two_e_integral(u,v,p,t,mo_integrals_map) * two_e_dm_mo(u,v,p,t)

! Just re-order the index and use 3D temporary arrays for optimal memory
! accesses.


!$OMP MASTER
CALL wall_TIME(t4)
!$OMP END MASTER

!$OMP DO
do p = 1, mo_num
  tmp_accu_1_shared(p) = 0d0
enddo
!$OMP END DO

!$OMP DO
do t = 1, mo_num 

  do p = 1, mo_num
    do v = 1, mo_num
      do u = 1, mo_num

         tmp_bi_int_3(u,v,p) = get_two_e_integral(u,v,p,t,mo_integrals_map)

      enddo
    enddo
  enddo 

  do p = 1, mo_num
    do v = 1, mo_num
      do u = 1, mo_num

         tmp_2rdm_3(u,v,p) = two_e_dm_mo(u,v,p,t)

      enddo
    enddo
  enddo

  !$OMP CRITICAL 
  do p = 1, mo_num
    do v = 1, mo_num
      do u = 1, mo_num

        tmp_accu_1_shared(p) = tmp_accu_1_shared(p) &
        + tmp_bi_int_3(u,v,p) * tmp_2rdm_3(u,v,p) 

      enddo
    enddo
  enddo
  !$OMP END CRITICAL 

enddo
!$OMP END DO

!$OMP DO
do p =1, mo_num

  tmp_h_pppp(p) = tmp_h_pppp(p) + tmp_accu_1_shared(p)

enddo
!$OMP END DO



! *Part 2 : q=r and p=s and q=r*

!  hessian(p,q,r,s) -> hessian(p,q,q,p)
 
!  0.5d0 * (  &
!   get_two_e_integral(u,v,p,t,mo_integrals_map) * two_e_dm_mo(u,v,s,t) &
! + get_two_e_integral(s,t,u,v,mo_integrals_map) * two_e_dm_mo(p,t,u,v))
!  =
!  0.5d0 * (  &
!   get_two_e_integral(u,v,p,t,mo_integrals_map) * two_e_dm_mo(u,v,p,t) &
! + get_two_e_integral(p,t,u,v,mo_integrals_map) * two_e_dm_mo(p,t,u,v))
!  = 
!  get_two_e_integral(u,v,p,t,mo_integrals_map) * two_e_dm_mo(u,v,p,t)

! Just re-order the index and use 3D temporary arrays for optimal memory
! accesses.


!$OMP DO
do p = 1, mo_num
  tmp_accu_1_shared(p) = 0d0
enddo
!$OMP END DO

!$OMP DO
do t = 1, mo_num

  do p = 1, mo_num
    do v = 1, mo_num
      do u = 1, mo_num

         tmp_bi_int_3(u,v,p) = get_two_e_integral(u,v,p,t,mo_integrals_map)

      enddo
    enddo
  enddo

  do p = 1, mo_num
    do v = 1, mo_num
      do u = 1, mo_num

         tmp_2rdm_3(u,v,p) = two_e_dm_mo(u,v,p,t)

      enddo
    enddo
  enddo

  !$OMP CRITICAL
  do p = 1, mo_num
    do v = 1, mo_num
      do u = 1, mo_num

        tmp_accu_1_shared(p) = tmp_accu_1_shared(p) + & 
          tmp_bi_int_3(u,v,p) * tmp_2rdm_3(u,v,p)

      enddo
    enddo
  enddo
  !$OMP END CRITICAL

enddo
!$OMP END DO

!$OMP DO
do q = 1, mo_num
  do p = 1, mo_num

    tmp_h_pqqp(p,q) = tmp_h_pqqp(p,q) + tmp_accu_1_shared(p) 

  enddo
enddo
!$OMP END DO

!$OMP MASTER
CALL wall_TIME(t5)
t6 = t5-t4
print*,'l2 1',t6
!$OMP END MASTER

! Line 2, term 2

! \begin{align*}
! \frac{1}{2} \sum_{tuv} \delta_{ps}(v_{uv}^{qt} \Gamma_{rt}^{uv} + v_{rt}^{uv}\Gamma_{uv}^{qt})
! \end{align*}

! Without optimization :

! do p = 1, mo_num
!   do q = 1, mo_num
!     do r = 1, mo_num
!       do s = 1, mo_num

!         ! Permutations 
!         if (((p==r) .and. (q==s)) .or. ((q==r) .and. (p==s))) then 

!           if (p==s) then
!             do t = 1, mo_num
!               do u = 1, mo_num
!                 do v = 1, mo_num

!                  hessian(p,q,r,s) = hessian(p,q,r,s) + 0.5d0 * ( &
!                     get_two_e_integral(q,t,u,v,mo_integrals_map) * two_e_dm_mo(r,t,u,v) &
!                   + get_two_e_integral(u,v,r,t,mo_integrals_map) * two_e_dm_mo(u,v,q,t))

!                 enddo
!               enddo
!             enddo
!           endif
!         endif
!       enddo
!     enddo
!   enddo
! enddo

! With optimization : 

! *Part 1 : p=r and q=s and p=s*

!  hessian(p,q,r,s) -> hessian(p,p,p,p)

!  0.5d0 * ( &
!   get_two_e_integral(q,t,u,v,mo_integrals_map) * two_e_dm_mo(r,t,u,v) &
! + get_two_e_integral(u,v,r,t,mo_integrals_map) * two_e_dm_mo(u,v,q,t)) 
!  = 
!  0.5d0 * ( &
!   get_two_e_integral(p,t,u,v,mo_integrals_map) * two_e_dm_mo(p,t,u,v) &
! + get_two_e_integral(u,v,p,t,mo_integrals_map) * two_e_dm_mo(u,v,p,t))
!  =
!  get_two_e_integral(u,v,p,t,mo_integrals_map) * two_e_dm_mo(u,v,p,t)

! Just re-order the index and use 3D temporary arrays for optimal memory
! accesses.


!$OMP MASTER
CALL wall_TIME(t4)
!$OMP END MASTER

!$OMP DO
do p = 1, mo_num
  tmp_accu_1_shared(p) = 0d0
enddo
!$OMP END DO   

!$OMP DO
do t = 1, mo_num

  do p = 1, mo_num
    do v = 1, mo_num
      do u = 1, mo_num

         tmp_bi_int_3(u,v,p) = get_two_e_integral(u,v,p,t,mo_integrals_map)

      enddo
    enddo
  enddo

  do p = 1, mo_num
    do v = 1, mo_num
      do u = 1, mo_num

         tmp_2rdm_3(u,v,p) = two_e_dm_mo(u,v,p,t)

      enddo
    enddo
  enddo

  !$OMP CRITICAL
  do p = 1, mo_num
    do v = 1, mo_num
      do u = 1, mo_num

        tmp_accu_1_shared(p) = tmp_accu_1_shared(p) +&
          tmp_bi_int_3(u,v,p) * tmp_2rdm_3(u,v,p)

      enddo
    enddo
  enddo
  !$OMP END CRITICAL

enddo
!$OMP END DO

!$OMP DO
do p = 1, mo_num

  tmp_h_pppp(p) = tmp_h_pppp(p) + tmp_accu_1_shared(p)

enddo
!$OMP END DO



! *Part 2 : q=r and p=s and p=s*

!  hessian(p,q,r,s) -> hessian(p,q,q,p)

!  0.5d0 * ( &
!   get_two_e_integral(q,t,u,v,mo_integrals_map) * two_e_dm_mo(r,t,u,v) &
! + get_two_e_integral(u,v,r,t,mo_integrals_map) * two_e_dm_mo(u,v,q,t)) 
!  = 
!  0.5d0 * ( &
!   get_two_e_integral(q,t,u,v,mo_integrals_map) * two_e_dm_mo(q,t,u,v) &
! + get_two_e_integral(u,v,q,t,mo_integrals_map) * two_e_dm_mo(u,v,q,t))
!  =
!  get_two_e_integral(u,v,q,t,mo_integrals_map) * two_e_dm_mo(u,v,q,t)

! Just re-order the index and use 3D temporary arrays for optimal memory
! accesses.


!$OMP DO
do p = 1,mo_num
  tmp_accu_1_shared(p) = 0d0
enddo
!$OMP END DO

!$OMP DO
do t = 1, mo_num

  do q = 1, mo_num
    do v = 1, mo_num
      do u = 1, mo_num

         tmp_bi_int_3(u,v,q) = get_two_e_integral(u,v,q,t,mo_integrals_map)

      enddo
    enddo
  enddo

  do q = 1, mo_num
    do v = 1, mo_num
      do u = 1, mo_num

         tmp_2rdm_3(u,v,q) = two_e_dm_mo(u,v,q,t)

      enddo
    enddo
  enddo

  !$OMP CRITICAL
  do q = 1, mo_num
    do v = 1, mo_num
      do u = 1, mo_num

        tmp_accu_1_shared(q) = tmp_accu_1_shared(q) +&
         tmp_bi_int_3(u,v,q) * tmp_2rdm_3(u,v,q)

      enddo
    enddo
  enddo
  !$OMP END CRITICAL

enddo
!$OMP END DO

!$OMP DO
do q = 1, mo_num
  do p = 1, mo_num

    tmp_h_pqqp(p,q) = tmp_h_pqqp(p,q) + tmp_accu_1_shared(p)

  enddo
enddo
!$OMP END DO   

!$OMP MASTER
CALL wall_TIME(t5)
t6 = t5-t4
print*,'l2 2',t6
!$OMP END MASTER

! Line 3, term 1

! \begin{align*}
! \sum_{uv} (v_{pr}^{uv} \Gamma_{uv}^{qs} + v_{uv}^{qs}  \Gamma_{pr}^{uv})
! \end{align*}

! Without optimization :

! do p = 1, mo_num
!   do q = 1, mo_num
!     do r = 1, mo_num
!       do s = 1, mo_num

!         ! Permutations 
!         if (((p==r) .and. (q==s)) .or. ((q==r) .and. (p==s)))) then

!           do u = 1, mo_num
!             do v = 1, mo_num

!               hessian(p,q,r,s) = hessian(p,q,r,s) &
!                + get_two_e_integral(u,v,p,r,mo_integrals_map) * two_e_dm_mo(u,v,q,s) &
!                + get_two_e_integral(q,s,u,v,mo_integrals_map) * two_e_dm_mo(p,r,u,v)

!             enddo
!           enddo
!         endif

!       enddo
!     enddo
!   enddo
! enddo

! With optimization
  
! *Part 1 : p=r and q=s*
  
!  hessian(p,q,r,s) -> hessian(p,q,p,q)
 
!   get_two_e_integral(u,v,p,r,mo_integrals_map) * two_e_dm_mo(u,v,q,s) &
! + get_two_e_integral(q,s,u,v,mo_integrals_map) * two_e_dm_mo(p,r,u,v) 
!  = 
!   get_two_e_integral(u,v,p,p,mo_integrals_map) * two_e_dm_mo(u,v,q,q) &
! + get_two_e_integral(q,q,u,v,mo_integrals_map) * two_e_dm_mo(p,p,u,v)
!  = 
!  2d0 * get_two_e_integral(u,v,p,p,mo_integrals_map) * two_e_dm_mo(u,v,q,q)

! Arrays of the kind (u,v,p,p) can be transform in 4D arrays (u,v,p).
! Using u,v as one variable a matrix multiplication appears.
! $$c_{p,q} = \sum_{uv} a_{p,uv} b_{uv,q}$$


!$OMP MASTER
CALL wall_TIME(t4)
!$OMP END MASTER

!$OMP DO
do q = 1, mo_num
  do v = 1, mo_num
    do u = 1, mo_num

      tmp_2rdm_3_shared(u,v,q) = two_e_dm_mo(u,v,q,q)

    enddo
  enddo
enddo
!$OMP END DO 

!$OMP DO
do p = 1, mo_num
  do v = 1, mo_num
    do u = 1, mo_num

      tmp_bi_int_3_shared(u,v,p) = get_two_e_integral(u,v,p,p,mo_integrals_map)

    enddo
  enddo
enddo
!$OMP END DO

call dgemm('T','N', mo_num, mo_num, mo_num*mo_num, 1d0, tmp_bi_int_3_shared,&
           mo_num*mo_num, tmp_2rdm_3_shared, mo_num*mo_num, 0d0, tmp_accu, mo_num)

!$OMP DO
do q = 1, mo_num
  do p = 1, mo_num

    tmp_h_pqpq(p,q) = tmp_h_pqpq(p,q) + tmp_accu(p,q) + tmp_accu(q,p)

  enddo
enddo
!$OMP END DO



! *Part 2 : q=r and p=s*
  
!  hessian(p,q,r,s) -> hessian(p,q,q,p)
 
!   get_two_e_integral(u,v,p,r,mo_integrals_map) * two_e_dm_mo(u,v,q,s) &
! + get_two_e_integral(q,s,u,v,mo_integrals_map) * two_e_dm_mo(p,r,u,v) 
!  = 
!   get_two_e_integral(u,v,p,q,mo_integrals_map) * two_e_dm_mo(u,v,q,p) &
! + get_two_e_integral(q,p,u,v,mo_integrals_map) * two_e_dm_mo(p,q,u,v)
!  = 
!  2d0 * get_two_e_integral(u,v,p,q,mo_integrals_map) * two_e_dm_mo(u,v,q,p)

! Just re-order the indexes and use 3D temporary arrays for optimal
! memory accesses.


!$OMP MASTER
call wall_time(t4)
!$OMP END MASTER

!$OMP DO
do q = 1, mo_num

  do p = 1, mo_num
    do v = 1, mo_num
      do u = 1, mo_num

        tmp_bi_int_3(u,v,p) = 2d0 * get_two_e_integral(u,v,q,p,mo_integrals_map)

      enddo
    enddo 
  enddo

  do p = 1, mo_num
    do v = 1, mo_num
      do u = 1, mo_num

        tmp_2rdm_3(u,v,p) = two_e_dm_mo(u,v,p,q)

      enddo
    enddo
  enddo

  do p = 1, mo_num
    do v = 1, mo_num
      do u = 1, mo_num

        tmp_h_pqqp(p,q) = tmp_h_pqqp(p,q) &
          + tmp_bi_int_3(u,v,p) * tmp_2rdm_3(u,v,p)

      enddo
    enddo
  enddo

enddo 
!$OMP END DO

!$OMP MASTER
CALL wall_TIME(t5)
t6= t5-t4
print*,'l3 1',t6
!$OMP END MASTER

! Line 3, term 2

! \begin{align*}
! - \sum_{tu} (v_{pu}^{st} \Gamma_{rt}^{qu}+v_{pu}^{tr} \Gamma_{tr}^{qu}+v_{rt}^{qu}\Gamma_{pu}^{st} + v_{tr}^{qu}\Gamma_{pu}^{ts})
! \end{align*}

! Without optimization :

! do p = 1, mo_num
!   do q = 1, mo_num
!     do r = 1, mo_num
!       do s = 1, mo_num

!         ! Permutations 
!         if (((p==r) .and. (q==s)) .or. ((q==r) .and. (p==s)) &
!               .or. ((p==s) .and. (q==r))) then

!           do t = 1, mo_num
!             do u = 1, mo_num

!               hessian(p,q,r,s) = hessian(p,q,r,s) &
!                - get_two_e_integral(s,t,p,u,mo_integrals_map) * two_e_dm_mo(r,t,q,u) &
!                - get_two_e_integral(t,s,p,u,mo_integrals_map) * two_e_dm_mo(t,r,q,u) &
!                - get_two_e_integral(q,u,r,t,mo_integrals_map) * two_e_dm_mo(p,u,s,t) &
!                - get_two_e_integral(q,u,t,r,mo_integrals_map) * two_e_dm_mo(p,u,t,s)

!             enddo
!           enddo

!         endif     

!       enddo
!     enddo
!   enddo
! enddo

! With optimization :

! *Part 1 : p=r and q=s*

!  hessian(p,q,r,s) -> hessian(p,q,p,q)

!  - get_two_e_integral(s,t,p,u,mo_integrals_map) * two_e_dm_mo(r,t,q,u) &
!  - get_two_e_integral(t,s,p,u,mo_integrals_map) * two_e_dm_mo(t,r,q,u) &
!  - get_two_e_integral(q,u,r,t,mo_integrals_map) * two_e_dm_mo(p,u,s,t) &
!  - get_two_e_integral(q,u,t,r,mo_integrals_map) * two_e_dm_mo(p,u,t,s)
!  =
!  - get_two_e_integral(q,t,p,u,mo_integrals_map) * two_e_dm_mo(p,t,q,u) &
!  - get_two_e_integral(t,q,p,u,mo_integrals_map) * two_e_dm_mo(t,p,q,u) &
!  - get_two_e_integral(q,u,p,t,mo_integrals_map) * two_e_dm_mo(p,u,q,t) &
!  - get_two_e_integral(q,u,t,p,mo_integrals_map) * two_e_dm_mo(p,u,t,q)
!  =
!  - 2d0 * get_two_e_integral(q,t,p,u,mo_integrals_map) * two_e_dm_mo(p,t,q,u) &
!  - 2d0 * get_two_e_integral(t,q,p,u,mo_integrals_map) * two_e_dm_mo(t,p,q,u) 
!  =
!  - 2d0 * get_two_e_integral(q,u,p,t,mo_integrals_map) * two_e_dm_mo(q,u,p,t) &
!  - 2d0 * get_two_e_integral(t,q,p,u,mo_integrals_map) * two_e_dm_mo(t,p,q,u)
 
! Just re-order the indexes and use 3D temporary arrays for optimal
! memory accesses.


!$OMP MASTER
CALL wall_TIME(t4)
!$OMP END MASTER

!----------
! Part 1.1
!----------
! - 2d0 * get_two_e_integral(q,u,p,t,mo_integrals_map) * two_e_dm_mo(q,u,p,t)

!$OMP DO
do q = 1, mo_num
  do p = 1, mo_num
    tmp_accu_shared(p,q) = 0d0
  enddo
enddo 
!$OMP END DO

!$OMP DO
do t = 1, mo_num

  do p = 1, mo_num
    do u = 1, mo_num
      do q = 1, mo_num

        tmp_bi_int_3(q,u,p) = 2d0 * get_two_e_integral(q,u,p,t,mo_integrals_map)

      enddo
    enddo
  enddo

  do p = 1, mo_num
    do u = 1, mo_num
      do q = 1, mo_num

         tmp_2rdm_3(q,u,p) = two_e_dm_mo(q,u,p,t)

      enddo
    enddo
  enddo

  !$OMP CRITICAL
  do p = 1, mo_num
    do u = 1, mo_num
      do q = 1, mo_num

         tmp_accu_shared(p,q) = tmp_accu_shared(p,q) &
         - tmp_bi_int_3(q,u,p) * tmp_2rdm_3(q,u,p) 

      enddo
    enddo
  enddo
  !$OMP END CRITICAL

enddo
!$OMP END DO

!$OMP DO
do q = 1, mo_num
  do p = 1, mo_num

    tmp_h_pqpq(p,q) = tmp_h_pqpq(p,q) + tmp_accu_shared(p,q)

  enddo
enddo
!$OMP END DO



! Just re-order the indexes and use 3D temporary arrays for optimal
! memory accesses.


!--------
! Part 1.2
!-------- 
! - 2d0 * get_two_e_integral(t,q,p,u,mo_integrals_map) * two_e_dm_mo(t,p,q,u)

!$OMP DO
do q = 1, mo_num
  do p = 1, mo_num
    tmp_accu_shared(p,q) = 0d0
  enddo
enddo
!$OMP END DO

!$OMP DO
do u = 1, mo_num

  do p = 1, mo_num
    do q = 1, mo_num
       do t = 1, mo_num

         tmp_bi_int_3(t,q,p) = 2d0*get_two_e_integral(t,q,p,u,mo_integrals_map)

       enddo
    enddo
  enddo

  do p= 1, mo_num
    do q = 1, mo_num
       do t = 1, mo_num

          tmp_2rdm_3(t,q,p) = two_e_dm_mo(t,p,q,u)

       enddo
    enddo
  enddo

  !$OMP CRITICAL
  do q = 1, mo_num
    do p = 1, mo_num
      do t = 1, mo_num

         tmp_accu_shared(p,q) = tmp_accu_shared(p,q) &
         - tmp_bi_int_3(t,q,p) * tmp_2rdm_3(t,q,p)

      enddo
    enddo
  enddo
  !$OMP END CRITICAL

enddo
!$OMP END DO

!$OMP DO
do q = 1, mo_num
  do p = 1, mo_num

    tmp_h_pqpq(p,q) = tmp_h_pqpq(p,q) + tmp_accu_shared(p,q)

  enddo
enddo
!$OMP END DO



! *Part 2 : q=r and p=s*

!  hessian(p,q,r,s) -> hessian(p,q,p,q)

!  - get_two_e_integral(s,t,p,u,mo_integrals_map) * two_e_dm_mo(r,t,q,u) &
!  - get_two_e_integral(t,s,p,u,mo_integrals_map) * two_e_dm_mo(t,r,q,u) &
!  - get_two_e_integral(q,u,r,t,mo_integrals_map) * two_e_dm_mo(p,u,s,t) &
!  - get_two_e_integral(q,u,t,r,mo_integrals_map) * two_e_dm_mo(p,u,t,s)
!  =
!  - get_two_e_integral(p,t,p,u,mo_integrals_map) * two_e_dm_mo(q,t,q,u) &
!  - get_two_e_integral(t,p,p,u,mo_integrals_map) * two_e_dm_mo(t,q,q,u) &
!  - get_two_e_integral(q,u,q,t,mo_integrals_map) * two_e_dm_mo(p,u,p,t) &
!  - get_two_e_integral(q,u,t,q,mo_integrals_map) * two_e_dm_mo(p,u,t,p)
!  =
!  - get_two_e_integral(p,t,p,u,mo_integrals_map) * two_e_dm_mo(q,t,q,u) &
!  - get_two_e_integral(q,t,q,u,mo_integrals_map) * two_e_dm_mo(p,t,p,u) &

!  - get_two_e_integral(t,u,p,p,mo_integrals_map) * two_e_dm_mo(t,q,q,u) &
!  - get_two_e_integral(t,u,q,q,mo_integrals_map) * two_e_dm_mo(t,p,p,u)
!  =
!  - get_two_e_integral(t,p,u,p,mo_integrals_map) * two_e_dm_mo(t,q,u,q) &
!  - get_two_e_integral(t,q,u,q,mo_integrals_map) * two_e_dm_mo(p,t,p,u) &

!  - get_two_e_integral(t,u,p,p,mo_integrals_map) * two_e_dm_mo(q,u,t,q) &
!  - get_two_e_integral(t,u,q,q,mo_integrals_map) * two_e_dm_mo(p,u,t,p)

! Arrays of the kind (t,p,u,p) can be transformed in 3D arrays. By doing
! so and using t,u as one variable, a matrix multiplication appears :
! $$c_{p,q} = \sum_{tu} a_{p,tu} b_{tu,q}$$


!----------
! Part 2.1
!----------
! - get_two_e_integral(t,p,u,p,mo_integrals_map) * two_e_dm_mo(t,q,u,q) &
! - get_two_e_integral(t,q,u,q,mo_integrals_map) * two_e_dm_mo(p,t,p,u)

!$OMP DO
do q = 1, mo_num
  do u = 1, mo_num
    do t = 1, mo_num

      tmp_2rdm_3_shared(t,u,q) = two_e_dm_mo(t,q,u,q)

    enddo
  enddo
enddo
!$OMP END DO

!$OMP DO
do p = 1, mo_num
  do u = 1, mo_num
    do t = 1, mo_num

      tmp_bi_int_3_shared(t,u,p) = get_two_e_integral(t,p,u,p,mo_integrals_map)

    enddo
  enddo
enddo
!$OMP END DO

call dgemm('T','N', mo_num, mo_num, mo_num*mo_num, 1d0, tmp_bi_int_3_shared,&
           mo_num*mo_num, tmp_2rdm_3_shared, mo_num*mo_num, 0d0, tmp_accu, mo_num)

!$OMP DO
do p = 1, mo_num
  do q = 1, mo_num

    tmp_h_pqqp(q,p) = tmp_h_pqqp(q,p) - tmp_accu(q,p) - tmp_accu(p,q)

  enddo
enddo
!$OMP END DO



! Arrays of the kind (t,u,p,p) can be transformed in 3D arrays. By doing
! so and using t,u as one variable, a matrix multiplication appears :
! $$c_{p,q} = \sum_{tu} a_{p,tu} b_{tu,q}$$


!--------
! Part 2.2 
!--------
! - get_two_e_integral(t,u,p,p,mo_integrals_map) * two_e_dm_mo(q,u,t,q) &
! - get_two_e_integral(t,u,q,q,mo_integrals_map) * two_e_dm_mo(p,u,t,p)

!$OMP DO
do p = 1, mo_num
  do u = 1, mo_num
    do t = 1, mo_num

      tmp_bi_int_3_shared(t,u,p) = get_two_e_integral(t,u,p,p,mo_integrals_map)

    enddo
  enddo
enddo
!$OMP END DO

!$OMP DO
do q = 1, mo_num
  do t = 1, mo_num
    do u = 1, mo_num

      tmp_2rdm_3_shared(u,t,q) = two_e_dm_mo(q,u,t,q)

    enddo
  enddo
enddo
!$OMP END DO

call dgemm('T','N', mo_num, mo_num, mo_num*mo_num, 1d0, tmp_2rdm_3_shared,&
           mo_num*mo_num, tmp_bi_int_3_shared, mo_num*mo_num, 0d0, tmp_accu, mo_num)

!$OMP DO
do q = 1, mo_num
  do p = 1, mo_num

    tmp_h_pqqp(p,q) = tmp_h_pqqp(p,q) - tmp_accu(p,q) - tmp_accu(q,p)

   enddo
enddo
!$OMP END DO

!$OMP MASTER
CALL wall_TIME(t5)
t6= t5-t4
print*,'l3 2',t6
!$OMP END MASTER  

!$OMP MASTER
CALL wall_TIME(t2)
t2 = t2 - t1
print*, 'Time to compute the hessian :', t2
!$OMP END MASTER

! Deallocation of private arrays
! In the OMP section !

deallocate(tmp_2rdm_3,tmp_bi_int_3)
deallocate(tmp_accu)

! Permutations
! As we mentioned before there are two permutation operator in the
! formula :
! Hessian(p,q,r,s) = P_pq P_rs [...]
! => Hessian(p,q,r,s) = (p,q,r,s) - (q,p,r,s) - (p,q,s,r) + (q,p,s,r)


!!$OMP DO
!do p = 1, mo_num
!  hessian(p,p,p,p) = hessian(p,p,p,p) + tmp_h_pppp(p)
!enddo
!!$OMP END DO

!!$OMP DO
!do q = 1, mo_num
!  do p = 1, mo_num
!    hessian(p,q,p,q) = hessian(p,q,p,q) + tmp_h_pqpq(p,q)
!  enddo
!enddo
!!$OMP END DO
!
!!$OMP DO
!do q = 1, mo_num
!  do p = 1, mo_num
!    hessian(p,q,q,p) = hessian(p,q,q,p) + tmp_h_pqqp(p,q)
!  enddo
!enddo
!!$OMP END DO

!!$OMP DO
!do s = 1, mo_num
!  do r = 1, mo_num
!    do q = 1, mo_num
!      do p = 1, mo_num

!        h_tmpr(p,q,r,s) = (hessian(p,q,r,s) - hessian(q,p,r,s) - hessian(p,q,s,r) + hessian(q,p,s,r))

!      enddo
!    enddo
!  enddo
!enddo
!!$OMP END DO

! 4D -> 2D matrix
! We need a 2D matrix for the Newton method's. Since the Hessian is
! "antisymmetric" : $$H_{pq,rs} = -H_{rs,pq}$$
! We can write it as a 2D matrix, N by N, with N = mo_num(mo_num-1)/2
! with p<q and r<s


!$OMP MASTER
CALL wall_TIME(t4)
!$OMP END MASTER

!$OMP DO
do pq = 1, n
  call vec_to_mat_index(pq,p,q)
  do rs = 1, n
    call vec_to_mat_index(rs,r,s)
    !H(pq,rs) = h_tmpr(p,q,r,s)
    if (pq == rs) then
      k = pq
      if (p == r) then
        H(k) = tmp_h_pqpq(p,q) + tmp_h_pqpq(q,p) - tmp_h_pqqp(p,q) - tmp_h_pqqp(q,p)
      elseif (p == s) then
        H(k) = - tmp_h_pqpq(p,q) - tmp_h_pqpq(q,p) + tmp_h_pqqp(p,q) + tmp_h_pqqp(q,p)
      endif
    endif
  enddo
enddo
!$OMP END DO

!!$OMP MASTER
!call wall_TIME(t5)
!t6 = t5-t4
!print*,'4D -> 2D :',t6
!!$OMP END MASTER

!$OMP END PARALLEL
call omp_set_max_active_levels(4)

! Display
!if (debug) then 
!  print*,'2D diag Hessian matrix'
!  do pq = 1, n
!    write(*,'(100(F10.5))') H(pq,:)
!  enddo 
!endif

! Deallocation of shared arrays, end


!deallocate(hessian)!,h_tmpr)
  deallocate(tmp_h_pppp,tmp_h_pqpq,tmp_h_pqqp)
  deallocate(tmp_accu_1_shared, tmp_accu_shared) 
 
  print*,'---diagonal_hessian'

end subroutine
