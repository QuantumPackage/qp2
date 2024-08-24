! Hessian

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

! The hessian is a 4D matrix of size mo_num, p,q,r,s,t,u,v take all the
! values between 1 and mo_num (1 and mo_num include).

! To do that we compute all the pairs (pq,rs)

! Source :
! Seniority-based coupled cluster theory
! J. Chem. Phys. 141, 244104 (2014); https://doi.org/10.1063/1.4904384
! Thomas M. Henderson, Ireneusz W. Bulik, Tamar Stein, and Gustavo E. Scuseria

! *Compute the hessian of energy with respects to orbital rotations*

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
! | hessian(mo_num,mo_num,mo_num,mo_num) | double precision | temporary array containing the hessian before            |
! |                                      |                  | the permutations                                         |
! | p, q, r, s                           | integer          | indexes of the hessian elements                          |
! | t, u, v                              | integer          | indexes for the sums                                     |
! | pq, rs                               | integer          | indexes for the transformation of the hessian            |
! |                                      |                  | (4D -> 2D)                                               |
! | t1,t2,t3                             | double precision | t3 = t2 - t1, time to compute the hessian                |
! | t4,t5,t6                             | double precision | t6 = t5 - t4, time to compute each element               |
! | tmp_bi_int_3(mo_num,mo_num,mo_num)   | double precision | 3 indexes temporary array for the bielectronic integrals |
! | tmp_2rdm_3(mo_num,mo_num,mo_num)     | double precision | 3 indexes temporary array for the 2 body density matrix  |
! | ind_3(mo_num,mo_num,mo_num)          | double precision | 3 indexes temporary array for matrix multiplication      |
! | tmp_accu(mo_num,mo_num)              | double precision | temporary array                                          |
! | tmp_accu_sym(mo_num,mo_num)          | double precision | temporary array                                          |

! Function:
! | get_two_e_integral | double precision | bielectronic integrals |


subroutine hessian_list_opt(n,m,list,H,h_tmpr)
  use omp_lib
  include 'constants.h' 

  implicit none
  
  ! Variables

  ! in
  integer, intent(in)           :: n,m,list(m) 
  
  ! out
  double precision, intent(out) :: H(n,n),h_tmpr(m,m,m,m)
 
  ! internal
  double precision, allocatable :: hessian(:,:,:,:)!, h_tmpr(:,:,:,:)
  double precision, allocatable :: H_test(:,:)
  integer                       :: p,q,tmp_p,tmp_q,tmp_r,tmp_s
  integer                       :: r,s,t,u,v,k
  integer                       :: pq,rs
  double precision              :: t1,t2,t3,t4,t5,t6
  ! H_test   : monum**2 by mo_num**2 double precision matrix to debug the H matrix

  double precision, allocatable :: tmp_bi_int_3(:,:,:), tmp_2rdm_3(:,:,:), ind_3(:,:,:),ind_3_3(:,:,:)
  double precision, allocatable :: tmp_bi_int_3_3(:,:,:), tmp_2rdm_3_3(:,:,:)
  double precision, allocatable :: tmp_accu(:,:), tmp_accu_sym(:,:),tmp_one_e_dm_mo(:,:),tmp_mo_one_e_integrals(:,:)

  ! Function 
  double precision              :: get_two_e_integral

  print*,''
  print*,'---hessian---'
  print*,'Use the full hessian'

  ! Allocation of shared arrays
  allocate(hessian(m,m,m,m),tmp_one_e_dm_mo(mo_num,m),tmp_mo_one_e_integrals(mo_num,m))

  ! Calculations

  ! OMP 
  call omp_set_max_active_levels(1)

  !$OMP PARALLEL                                                     &
      !$OMP PRIVATE(                                                 &
      !$OMP   tmp_p,tmp_q,tmp_r,tmp_s,p,q,r,s, tmp_accu, tmp_accu_sym,                       &
      !$OMP   u,v,t, tmp_bi_int_3, tmp_2rdm_3, ind_3, tmp_bi_int_3_3,tmp_2rdm_3_3, ind_3_3 )                       &
      !$OMP SHARED(m,list,hessian,h_tmpr,H, mo_num,n, & 
      !$OMP mo_one_e_integrals, one_e_dm_mo, &
      !$OMP two_e_dm_mo,mo_integrals_map, &
      !$OMP t1,t2,t3,t4,t5,t6,&
      !$OMP tmp_mo_one_e_integrals,tmp_one_e_dm_mo)&
      !$OMP DEFAULT(NONE)
 
  ! Allocation of private arrays 
  allocate(tmp_bi_int_3(mo_num,mo_num,mo_num))
  allocate(tmp_2rdm_3(mo_num,mo_num,mo_num), ind_3(mo_num,mo_num,mo_num))
  allocate(tmp_accu(m,m), tmp_accu_sym(mo_num,mo_num))

! Initialization of the arrays

!$OMP MASTER
do tmp_q = 1, m
  do tmp_p = 1, m
    tmp_accu(tmp_p,tmp_q) = 0d0
  enddo
enddo
!$OMP END MASTER

!$OMP MASTER
do tmp_q = 1, m
  do tmp_p = 1, m
    tmp_accu_sym(tmp_p,tmp_q) = 0d0
  enddo
enddo
!$OMP END MASTER

!$OMP DO
do tmp_s = 1, m
  do tmp_r = 1, m
    do tmp_q = 1, m
      do tmp_p = 1, m
        hessian(tmp_p,tmp_q,tmp_r,tmp_s) = 0d0
      enddo
    enddo
  enddo
enddo
!$OMP ENDDO

!$OMP MASTER
CALL wall_TIME(t1)
!$OMP END MASTER

! Line 1, term 1
   
! Without optimization the term 1 of the line 1 is :

! do p = 1, mo_num
!   do q = 1, mo_num
!     do r = 1, mo_num
!       do s = 1, mo_num

!         if (q==r) then
!           do u = 1, mo_num

!             hessian(p,q,r,s) = hessian(p,q,r,s) + 0.5d0 * (  &
!               mo_one_e_integrals(u,p) * one_e_dm_mo(u,s) &
!             + mo_one_e_integrals(s,u) * one_e_dm_mo(p,u))

!           enddo
!         endif

!       enddo
!     enddo
!   enddo
! enddo

! We can write the formula as matrix multiplication.
! $$c_{p,s} = \sum_u a_{p,u} b_{u,s}$$


!$OMP MASTER    
CALL wall_TIME(t4)
!$OMP END MASTER

!$OMP DO
do tmp_p = 1, m
  p = list(tmp_p)
  do u = 1, mo_num
    tmp_mo_one_e_integrals(u,tmp_p) = mo_one_e_integrals(u,p)
  enddo
enddo
!$OMP END DO

!$OMP DO
do tmp_s = 1, m
  s = list(tmp_s)
  do u = 1, mo_num
    tmp_one_e_dm_mo(u,tmp_s) = one_e_dm_mo(u,s)
  enddo
enddo
!$OMP END DO


call dgemm('T','N', m, m, mo_num, 1d0, tmp_mo_one_e_integrals,&
           size(tmp_mo_one_e_integrals,1), tmp_one_e_dm_mo, size(tmp_one_e_dm_mo,1),&
           0d0, tmp_accu, size(tmp_accu,1))

!$OMP DO
do tmp_s = 1, m
  do tmp_p = 1, m

    tmp_accu_sym(tmp_p,tmp_s) = 0.5d0 * (tmp_accu(tmp_p,tmp_s) + tmp_accu(tmp_s,tmp_p))

  enddo
enddo 
!$OMP END DO

!$OMP DO
do tmp_s = 1, m
  do tmp_p = 1, m
    do tmp_r = 1, m

      hessian(tmp_p,tmp_r,tmp_r,tmp_s) = hessian(tmp_p,tmp_r,tmp_r,tmp_s) + tmp_accu_sym(tmp_p,tmp_s)

    enddo
  enddo
enddo
!$OMP END DO

!$OMP MASTER
CALL wall_TIME(t5)
t6=t5-t4
print*,'l1 1',t6
!$OMP END MASTER

! Line 1, term 2
! do p = 1, mo_num
!   do q = 1, mo_num
!     do r = 1, mo_num
!       do s = 1, mo_num

!         if (p==s) then
!           do u = 1, mo_num

!                 hessian(p,q,r,s) = hessian(p,q,r,s) + 0.5d0 * ( &
!                   mo_one_e_integrals(u,r) * (one_e_dm_mo(u,q) &
!                 + mo_one_e_integrals(q,u) * (one_e_dm_mo(r,u))
!           enddo
!         endif

!       enddo
!     enddo
!   enddo
! enddo

! We can write the formula as matrix multiplication.
! $$c_{r,q} = \sum_u a_{r,u} b_{u,q}$$


!$OMP MASTER
CALL wall_TIME(t4)
!$OMP END MASTER

call dgemm('T','N', m, m, mo_num, 1d0, tmp_mo_one_e_integrals,&
           size(tmp_mo_one_e_integrals,1), tmp_one_e_dm_mo, size(tmp_one_e_dm_mo,1),&
           0d0, tmp_accu, size(tmp_accu,1))

!$OMP DO
do tmp_r = 1, m
  do tmp_q = 1, m

    tmp_accu_sym(tmp_q,tmp_r) = 0.5d0 * (tmp_accu(tmp_q,tmp_r) + tmp_accu(tmp_r,tmp_q))

  enddo
enddo
!OMP END DO

!$OMP DO
do tmp_r = 1, m
  do tmp_q = 1, m
    do tmp_s = 1, m

      hessian(tmp_s,tmp_q,tmp_r,tmp_s) = hessian(tmp_s,tmp_q,tmp_r,tmp_s) + tmp_accu_sym(tmp_q,tmp_r)

    enddo
  enddo
enddo
!OMP END DO

!$OMP MASTER
CALL wall_TIME(t5)
t6=t5-t4
print*,'l1 2',t6
!$OMP END MASTER

! Line 1, term 3 

! Without optimization the third term is :

! do p = 1, mo_num
!   do q = 1, mo_num
!     do r = 1, mo_num
!       do s = 1, mo_num

!         hessian(p,q,r,s) = hessian(p,q,r,s) &
!         - mo_one_e_integrals(s,p) * one_e_dm_mo(r,q) &
!         - mo_one_e_integrals(q,r) * one_e_dm_mo(p,s))

!       enddo
!     enddo
!   enddo
! enddo

! We can just re-order the indexes
 

!$OMP MASTER
CALL wall_TIME(t4)
!$OMP END MASTER

!$OMP DO
do tmp_s = 1, m
  s = list(tmp_s)
  do tmp_r = 1, m
    r = list(tmp_r)
    do tmp_q = 1, m
      q = list(tmp_q)
      do tmp_p = 1, m
        p = list(tmp_p)

        hessian(tmp_p,tmp_q,tmp_r,tmp_s) = hessian(tmp_p,tmp_q,tmp_r,tmp_s) &
          - mo_one_e_integrals(s,p) * one_e_dm_mo(r,q)&
          - mo_one_e_integrals(q,r) * one_e_dm_mo(p,s)

      enddo
    enddo
  enddo
enddo
!$OMP END DO

!$OMP MASTER
CALL wall_TIME(t5)
t6=t5-t4
print*,'l1 3',t6
!$OMP END MASTER

! Line 2, term 1

! Without optimization the fourth term is :

! do p = 1, mo_num
!   do q = 1, mo_num
!     do r = 1, mo_num
!       do s = 1, mo_num

!          if (q==r) then
!            do t = 1, mo_num
!              do u = 1, mo_num
!                do v = 1, mo_num

!                  hessian(p,q,r,s) = hessian(p,q,r,s) + 0.5d0 * (  &
!                    get_two_e_integral(u,v,p,t,mo_integrals_map) * two_e_dm_mo(u,v,s,t) &
!                  + get_two_e_integral(s,t,u,v,mo_integrals_map) * two_e_dm_mo(p,t,u,v))

!                enddo
!              enddo
!            enddo
!          endif

!       enddo
!     enddo
!   enddo
! enddo

! Using bielectronic integral properties :
! get_two_e_integral(s,t,u,v,mo_integrals_map) =
! get_two_e_integral(u,v,s,t,mo_integrals_map)

! Using the two electron density matrix properties :
! two_e_dm_mo(p,t,u,v) = two_e_dm_mo(u,v,p,t)

! With t on the external loop, using temporary arrays for each t and by
! taking u,v as one variable a matrix multplication appears. 
! $$c_{p,s} = \sum_{uv} a_{p,uv} b_{uv,s}$$

! There is a kroenecker delta $$\delta_{qr}$$, so we juste compute the
! terms like : hessian(p,r,r,s)


!$OMP MASTER 
call wall_TIME(t4)
!$OMP END MASTER

allocate(tmp_bi_int_3_3(mo_num,mo_num,m), tmp_2rdm_3_3(mo_num,mo_num,m))

!$OMP DO
do t = 1, mo_num

  do tmp_p = 1, m
    p = list(tmp_p)
    do v = 1, mo_num
      do u = 1, mo_num

        tmp_bi_int_3_3(u,v,tmp_p) = get_two_e_integral(u,v,p,t,mo_integrals_map)              

      enddo
    enddo
  enddo

  do tmp_p = 1, m ! error, the p might be replace by a s
  ! it's a temporary array, the result by replacing p and s will be the same
    p = list(tmp_p)
    do v = 1, mo_num
      do u = 1, mo_num

        tmp_2rdm_3_3(u,v,tmp_p) = two_e_dm_mo(u,v,p,t)  

      enddo
    enddo
  enddo

  call dgemm('T','N', m, m, mo_num*mo_num, 1.d0, &
             tmp_bi_int_3_3, mo_num*mo_num, tmp_2rdm_3_3, mo_num*mo_num, &
             0.d0, tmp_accu, size(tmp_accu,1))

  do tmp_p = 1, m
    do tmp_s = 1, m

      tmp_accu_sym(tmp_s,tmp_p) = 0.5d0 * (tmp_accu(tmp_p,tmp_s)+tmp_accu(tmp_s,tmp_p))

    enddo
  enddo

  !$OMP CRITICAL 
  do tmp_s = 1, m
    do tmp_r = 1, m
      do tmp_p = 1, m

        hessian(tmp_p,tmp_r,tmp_r,tmp_s) = hessian(tmp_p,tmp_r,tmp_r,tmp_s) + tmp_accu_sym(tmp_p,tmp_s) 

      enddo
    enddo
  enddo
  !$OMP END CRITICAL

enddo
!$OMP END DO

!$OMP MASTER
call wall_TIME(t5)
t6=t5-t4
print*,'l2 1', t6 
!$OMP END MASTER

! Line 2, term 2

! do p = 1, mo_num
!   do q = 1, mo_num
!     do r = 1, mo_num
!       do s = 1, mo_num

!         if (p==s) then
!           do t = 1, mo_num
!             do u = 1, mo_num
!               do v = 1, mo_num

!                 hessian(p,q,r,s) = hessian(p,q,r,s) + 0.5d0 * ( &
!                   get_two_e_integral(q,t,u,v,mo_integrals_map) * two_e_dm_mo(r,t,u,v) &
!                 + get_two_e_integral(u,v,r,t,mo_integrals_map) * two_e_dm_mo(u,v,q,t))

!               enddo
!             enddo
!           enddo
!         endif

!       enddo
!     enddo
!   enddo
! enddo

! Using the two electron density matrix properties :
! get_two_e_integral(q,t,u,v,mo_integrals_map) =
! get_two_e_integral(u,v,q,t,mo_integrals_map)

! Using the two electron density matrix properties :
! two_e_dm_mo(r,t,u,v) = two_e_dm_mo(u,v,r,t)

! With t on the external loop, using temporary arrays for each t and by
! taking u,v as one variable a matrix multplication appears. 
! $$c_{q,r} = \sum_uv a_{q,uv} b_{uv,r}$$ 

! There is a kroenecker delta $$\delta_{ps}$$, so we juste compute the
! terms like : hessian(s,q,r,s)


!******************************
! Opt Second line, second term
!******************************

!$OMP MASTER 
CALL wall_TIME(t4)
!$OMP END MASTER



!$OMP DO
do t = 1, mo_num

  do tmp_q = 1, m
    q = list(tmp_q)
    do v = 1, mo_num
      do u = 1, mo_num

        tmp_bi_int_3_3(u,v,tmp_q) = get_two_e_integral(u,v,q,t,mo_integrals_map)

      enddo
    enddo
  enddo

  do tmp_r = 1, m
    r = list(tmp_r)
    do v = 1, mo_num
      do u = 1, mo_num

         tmp_2rdm_3_3(u,v,tmp_r) = two_e_dm_mo(u,v,r,t)

      enddo
    enddo
  enddo

  call dgemm('T','N', m, m, mo_num*mo_num, 1.d0, &
             tmp_bi_int_3_3 , mo_num*mo_num, tmp_2rdm_3_3, mo_num*mo_num, &
             0.d0, tmp_accu, size(tmp_accu,1))

  do tmp_r = 1, m
    do tmp_q = 1, m

      tmp_accu_sym(tmp_q,tmp_r) = 0.5d0 * (tmp_accu(tmp_q,tmp_r) + tmp_accu(tmp_r,tmp_q))

    enddo
  enddo

  !$OMP CRITICAL
  do tmp_r = 1, m
    do tmp_q = 1, m
      do tmp_s = 1, m

        hessian(tmp_s,tmp_q,tmp_r,tmp_s) = hessian(tmp_s,tmp_q,tmp_r,tmp_s) + tmp_accu_sym(tmp_q,tmp_r)

      enddo
    enddo
  enddo
  !$OMP END CRITICAL

enddo
!$OMP END DO

deallocate(tmp_bi_int_3_3, tmp_2rdm_3_3)

!$OMP MASTER
CALL wall_TIME(t5)
t6=t5-t4
print*,'l2 2',t6
!$OMP END MASTER

! Line 3, term 1

! do p = 1, mo_num
!   do q = 1, mo_num
!     do r = 1, mo_num
!       do s = 1, mo_num

!         do u = 1, mo_num
!           do v = 1, mo_num

!             hessian(p,q,r,s) = hessian(p,q,r,s) &
!              + get_two_e_integral(u,v,p,r,mo_integrals_map) * two_e_dm_mo(u,v,q,s) &
!              + get_two_e_integral(q,s,u,v,mo_integrals_map) * two_e_dm_mo(p,r,u,v)

!           enddo
!         enddo

!       enddo
!     enddo
!   enddo
! enddo

! Using the two electron density matrix properties :
! get_two_e_integral(u,v,p,r,mo_integrals_map) =
! get_two_e_integral(p,r,u,v,mo_integrals_map)

! Using the two electron density matrix properties :
! two_e_dm_mo(u,v,q,s) =  two_e_dm_mo(q,s,u,v)

! With v on the external loop, using temporary arrays for each v and by
! taking p,r and q,s as one dimension a matrix multplication
! appears. $$c_{pr,qs} = \sum_u a_{pr,u} b_{u,qs}$$ 

! Part 1

!$OMP MASTER 
call wall_TIME(t4)
!$OMP END MASTER 

!--------
! part 1
! get_two_e_integral(u,v,p,r,mo_integrals_map) * two_e_dm_mo(u,v,q,s)
!--------

allocate(tmp_bi_int_3_3(m,m,mo_num), tmp_2rdm_3_3(mo_num,m,m),ind_3_3(m,m,m))

!$OMP DO
do v = 1, mo_num

  do u = 1, mo_num 
    do tmp_r = 1, m
      r = list(tmp_r)
      do tmp_p = 1, m
        p = list(tmp_p)

          tmp_bi_int_3_3(tmp_p,tmp_r,u) = get_two_e_integral(p,r,u,v,mo_integrals_map)

      enddo
    enddo
  enddo

  do tmp_s = 1, m
    s = list(tmp_s)
    do tmp_q = 1, m
      q = list(tmp_q)
      do u = 1, mo_num

        tmp_2rdm_3_3(u,tmp_q,tmp_s) = two_e_dm_mo(q,s,u,v)

      enddo
    enddo
  enddo

  do tmp_s = 1, m

    call dgemm('N','N',m*m, m, mo_num, 1d0, tmp_bi_int_3_3,&
               size(tmp_bi_int_3_3,1)*size(tmp_bi_int_3_3,2), tmp_2rdm_3_3(1,1,tmp_s),&
               mo_num, 0d0, ind_3_3, size(ind_3_3,1) * size(ind_3_3,2))

    !$OMP CRITICAL
    do tmp_r = 1, m
      do tmp_q = 1, m
        do tmp_p = 1, m
          hessian(tmp_p,tmp_q,tmp_r,tmp_s) = hessian(tmp_p,tmp_q,tmp_r,tmp_s) + ind_3_3(tmp_p,tmp_r,tmp_q)
        enddo
      enddo
    enddo
    !$OMP END CRITICAL

  enddo

enddo
!$OMP END DO

deallocate(tmp_bi_int_3_3, tmp_2rdm_3_3,ind_3_3)



! With v on the external loop, using temporary arrays for each v and by
! taking q,s and p,r as one dimension a matrix multplication
! appears. $$c_{qs,pr} = \sum_u a_{qs,u}*b_{u,pr}$$ 

! Part 2

!--------
! part 2
! get_two_e_integral(q,s,u,v,mo_integrals_map) * two_e_dm_mo(p,r,u,v)
!--------

allocate(tmp_bi_int_3_3(m,m,mo_num), tmp_2rdm_3_3(mo_num,m,m),ind_3_3(m,m,m))

!$OMP DO
do v = 1, mo_num

  do u = 1, mo_num
    do tmp_s = 1, m
    s = list(tmp_s)  
      do tmp_q = 1, m
        q = list(tmp_q)

          tmp_bi_int_3_3(tmp_q,tmp_s,u) = get_two_e_integral(q,s,u,v,mo_integrals_map)

      enddo
    enddo
  enddo

  do tmp_r = 1, m
    r = list(tmp_r)
    do tmp_p = 1, m
      p = list(tmp_p)
      do u = 1, mo_num

        tmp_2rdm_3_3(u,tmp_p,tmp_r) = two_e_dm_mo(p,r,u,v)

      enddo
    enddo
  enddo

  do tmp_r = 1, m

    call dgemm('N','N', m*m, m, mo_num, 1d0, tmp_bi_int_3_3,& 
               size(tmp_bi_int_3_3,1)*size(tmp_bi_int_3_3,2), tmp_2rdm_3_3(1,1,tmp_r),&
               mo_num, 0d0, ind_3_3, size(ind_3_3,1) * size(ind_3_3,2))

    !$OMP CRITICAL
    do tmp_s = 1, m
      s = list(tmp_s)
      do tmp_q = 1, m
        q = list(tmp_q)
        do tmp_p = 1, m
          p = list(tmp_p)
          hessian(tmp_p,tmp_q,tmp_r,tmp_s) = hessian(tmp_p,tmp_q,tmp_r,tmp_s) + ind_3_3(tmp_q,tmp_s,tmp_p)
        enddo
      enddo
    enddo
    !$OMP END CRITICAL

  enddo

enddo
!$OMP END DO

deallocate(tmp_bi_int_3_3, tmp_2rdm_3_3,ind_3_3)

!$OMP MASTER
call wall_TIME(t5)
t6 = t5 - t4
print*,'l3 1', t6
!$OMP END MASTER

! Line 3, term 2

! do p = 1, mo_num
!   do q = 1, mo_num
!     do r = 1, mo_num
!       do s = 1, mo_num

!         do t = 1, mo_num
!           do u = 1, mo_num

!             hessian(p,q,r,s) = hessian(p,q,r,s) &
!              - get_two_e_integral(s,t,p,u,mo_integrals_map) * two_e_dm_mo(r,t,q,u) &
!              - get_two_e_integral(t,s,p,u,mo_integrals_map) * two_e_dm_mo(t,r,q,u) &
!              - get_two_e_integral(q,u,r,t,mo_integrals_map) * two_e_dm_mo(p,u,s,t) &
!              - get_two_e_integral(q,u,t,r,mo_integrals_map) * two_e_dm_mo(p,u,t,s)

!           enddo
!         enddo

!       enddo
!     enddo
!   enddo
! enddo

! With q on the external loop, using temporary arrays for each p and q,
! and taking u,v as one variable, a matrix multiplication appears:
! $$c_{r,s} = \sum_{ut} a_{r,ut} b_{ut,s}$$

! Part 1

!--------
! Part 1 
! - get_two_e_integral(s,t,p,u,mo_integrals_map) * two_e_dm_mo(r,t,q,u)
!--------

!$OMP MASTER
CALL wall_TIME(t4) 
!$OMP END MASTER

allocate(tmp_bi_int_3_3(mo_num,mo_num,m), tmp_2rdm_3_3(mo_num,mo_num,m))

!$OMP DO
do tmp_q = 1, m
  q = list(tmp_q)

  do tmp_r = 1, m
    r = list(tmp_r)
    do t = 1, mo_num
      do u = 1, mo_num

        tmp_2rdm_3_3(u,t,tmp_r) = two_e_dm_mo(q,u,r,t)

      enddo
    enddo
  enddo

  do tmp_p = 1, m
    p = list(tmp_p)

    do tmp_s = 1, m
      s = list(tmp_s)
      do t = 1, mo_num
        do u = 1, mo_num

          tmp_bi_int_3_3(u,t,tmp_s) = - get_two_e_integral(u,s,t,p,mo_integrals_map)

        enddo
      enddo
    enddo

    call dgemm('T','N', m, m, mo_num*mo_num, 1d0, tmp_bi_int_3_3,&
               mo_num*mo_num, tmp_2rdm_3_3, mo_num*mo_num, 0d0, tmp_accu, m)

    !$OMP CRITICAL
    do tmp_s = 1, m
      do tmp_r = 1, m

         hessian(tmp_p,tmp_q,tmp_r,tmp_s) = hessian(tmp_p,tmp_q,tmp_r,tmp_s) + tmp_accu(tmp_s,tmp_r)

      enddo
    enddo
    !$OMP END CRITICAL

  enddo

enddo
!$OMP END DO

deallocate(tmp_bi_int_3_3, tmp_2rdm_3_3)



! With q on the external loop, using temporary arrays for each p and q,
! and taking u,v as one variable, a matrix multiplication appears:
! $$c_{r,s} = \sum_{ut} a_{r,ut} b_{ut,s}$$

! Part 2

!--------
! Part 2
!- get_two_e_integral(t,s,p,u,mo_integrals_map) * two_e_dm_mo(t,r,q,u)
!--------

allocate(tmp_bi_int_3_3(mo_num,mo_num,m), tmp_2rdm_3_3(mo_num,mo_num,m))

!$OMP DO 
do tmp_q = 1, m
  q = list(tmp_q)

  do tmp_r = 1, m
    r = list(tmp_r)
    do t = 1, mo_num
      do u = 1, mo_num

        tmp_2rdm_3_3(u,t,tmp_r) = two_e_dm_mo(q,u,t,r)

      enddo
    enddo
  enddo

  do tmp_p = 1, m
    p = list(tmp_p)

    do tmp_s = 1, m
      s = list(tmp_s)
      do t = 1, mo_num
        do u = 1, mo_num

          tmp_bi_int_3_3(u,t,tmp_s) = - get_two_e_integral(u,t,s,p,mo_integrals_map)

        enddo
      enddo
    enddo   

    call dgemm('T','N', m, m, mo_num*mo_num, 1d0, tmp_bi_int_3_3,&
               mo_num*mo_num, tmp_2rdm_3_3, mo_num*mo_num, 0d0, tmp_accu, m)

    !$OMP CRITICAL
    do tmp_s = 1, m
      do tmp_r = 1, m

        hessian(tmp_p,tmp_q,tmp_r,tmp_s) = hessian(tmp_p,tmp_q,tmp_r,tmp_s) + tmp_accu(tmp_s,tmp_r)

      enddo
    enddo
    !$OMP END CRITICAL

  enddo

enddo
!$OMP END DO

deallocate(tmp_bi_int_3_3, tmp_2rdm_3_3)



! With q on the external loop, using temporary arrays for each p and q,
! and taking u,v as one variable, a matrix multiplication appears:
! $$c_{r,s} = \sum_{ut} a_{r,ut} b_{ut,s}$$

! Part 3

!--------
! Part 3
!- get_two_e_integral(q,u,r,t,mo_integrals_map) * two_e_dm_mo(p,u,s,t) 
!--------

allocate(tmp_bi_int_3_3(mo_num,mo_num,m), tmp_2rdm_3_3(mo_num,mo_num,m))

!$OMP DO 
do tmp_q = 1, m
  q = list(tmp_q)

  do tmp_r = 1, m
    r = list(tmp_r)
    do t = 1, mo_num
      do u = 1, mo_num

        tmp_bi_int_3_3(u,t,tmp_r) = - get_two_e_integral(u,q,t,r,mo_integrals_map)

      enddo
    enddo
  enddo

  do tmp_p = 1, m
    p = list(tmp_p)

    do tmp_s = 1, m
      s = list(tmp_s)
      do t = 1, mo_num
        do u = 1, mo_num

          tmp_2rdm_3_3(u,t,tmp_s) = two_e_dm_mo(p,u,s,t)

        enddo
      enddo
    enddo

    call dgemm('T','N', m, m, mo_num*mo_num, 1d0, tmp_2rdm_3_3,&
               mo_num*mo_num, tmp_bi_int_3_3, mo_num*mo_num, 0d0, tmp_accu, m)

    !$OMP CRITICAL
    do tmp_s = 1, m
      do tmp_r = 1, m

        hessian(tmp_p,tmp_q,tmp_r,tmp_s) = hessian(tmp_p,tmp_q,tmp_r,tmp_s) + tmp_accu(tmp_s,tmp_r)

      enddo
    enddo
    !$OMP END CRITICAL

  enddo

enddo
!$OMP END DO

deallocate(tmp_bi_int_3_3, tmp_2rdm_3_3)



! With q on the external loop, using temporary arrays for each p and q,
! and taking u,v as one variable, a matrix multiplication appears:
! $$c_{r,s} = \sum_{ut} a_{r,ut} b_{ut,s}$$

! Part 4

!--------
! Part 4
! - get_two_e_integral(q,u,t,r,mo_integrals_map) * two_e_dm_mo(p,u,t,s)
!--------

allocate(tmp_bi_int_3_3(mo_num,mo_num,m), tmp_2rdm_3_3(mo_num,mo_num,m))

!$OMP DO
do tmp_q = 1, m
  q = list(tmp_q)

  do tmp_r = 1, m
    r = list(tmp_r)
    do t = 1, mo_num
      do u = 1, mo_num

        tmp_bi_int_3_3(u,t,tmp_r) = - get_two_e_integral(u,t,r,q,mo_integrals_map)

      enddo
    enddo
  enddo

  do tmp_p = 1, m
    p = list(tmp_p)

    do tmp_s = 1, m
      s = list(tmp_s)
      do t = 1, mo_num
        do u = 1, mo_num

          tmp_2rdm_3_3(u,t,tmp_s) = two_e_dm_mo(p,u,t,s)

        enddo
      enddo
    enddo

    call dgemm('T','N', m, m, mo_num*mo_num, 1d0, tmp_2rdm_3_3,&
               mo_num*mo_num, tmp_bi_int_3_3, mo_num*mo_num, 0d0, tmp_accu, m)

    !$OMP CRITICAL
    do tmp_s = 1, m
      do tmp_r = 1, m

        hessian(tmp_p,tmp_q,tmp_r,tmp_s) = hessian(tmp_p,tmp_q,tmp_r,tmp_s) + tmp_accu(tmp_s,tmp_r)

      enddo
    enddo
    !$OMP END CRITICAL

  enddo

enddo
!$OMP END DO  

deallocate(tmp_bi_int_3_3, tmp_2rdm_3_3)

!$OMP MASTER
call wall_TIME(t5)
t6 = t5-t4
print*,'l3 2',t6
!$OMP END MASTER

!$OMP MASTER
CALL wall_TIME(t2)
t3 = t2 -t1
print*,'Time to compute the hessian : ', t3
!$OMP END MASTER

! Deallocation of private arrays
! In the omp section !

deallocate(tmp_bi_int_3, tmp_2rdm_3, tmp_accu, tmp_accu_sym, ind_3)

! Permutations
! As we mentioned before there are two permutation operator in the
! formula :
! Hessian(p,q,r,s) = P_pq P_rs [...]
! => Hessian(p,q,r,s) = (p,q,r,s) - (q,p,r,s) - (p,q,s,r) + (q,p,s,r)


!$OMP MASTER
CALL wall_TIME(t4)
!$OMP END MASTER

!$OMP DO
do s = 1, m
  do r = 1, m
    do q = 1, m
      do p = 1, m

        h_tmpr(p,q,r,s) = (hessian(p,q,r,s) - hessian(q,p,r,s) - hessian(p,q,s,r) + hessian(q,p,s,r))

      enddo
    enddo
  enddo
enddo
!$OMP END DO

!$OMP MASTER
call wall_TIME(t5)
t6 = t5-t4
print*,'Time for permutations :',t6
!$OMP END MASTER

! 4D -> 2D matrix
! We need a 2D matrix for the Newton method's. Since the Hessian is
! "antisymmetric" : $$H_{pq,rs} = -H_{rs,pq}$$
! We can write it as a 2D matrix, N by N, with N = mo_num(mo_num-1)/2
! with p<q and r<s


!$OMP MASTER
CALL wall_TIME(t4)
!$OMP END MASTER

!$OMP DO
do rs = 1, n
  call vec_to_mat_index(rs,r,s)
  do pq = 1, n
    call vec_to_mat_index(pq,p,q)
    H(pq,rs) = h_tmpr(p,q,r,s)   
  enddo
enddo
!$OMP END DO 

!$OMP MASTER
call wall_TIME(t5)
t6 = t5-t4
print*,'4D -> 2D :',t6
!$OMP END MASTER

!$OMP END PARALLEL
call omp_set_max_active_levels(4)

! Display
if (debug) then 
  print*,'2D Hessian matrix'
  do pq = 1, n
    write(*,'(100(F10.5))') H(pq,:)
  enddo 
endif

! Deallocation of shared arrays, end

deallocate(hessian,tmp_one_e_dm_mo,tmp_mo_one_e_integrals)!,h_tmpr)
! h_tmpr is intent out in order to debug the subroutine
! It's why we don't deallocate it

  print*,'---End hessian---'

end subroutine
