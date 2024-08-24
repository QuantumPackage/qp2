subroutine first_hessian_list_opt(tmp_n,m,list,H,h_tmpr)

  include 'constants.h' 

  implicit none

  !==================================================================
  ! Compute the hessian of energy with respects to orbital rotations
  !==================================================================

  !===========
  ! Variables  
  !===========

  ! in
  integer, intent(in)           :: tmp_n, m, list(m)
  !tmp_n         : integer, tmp_n = m*(m-1)/2
  
  ! out
  double precision, intent(out) :: H(tmp_n,tmp_n),h_tmpr(m,m,m,m)
  ! H        : n by n double precision matrix containing the 2D hessian
 
  ! internal
  double precision, allocatable :: hessian(:,:,:,:)
  integer                       :: p,q, tmp_p,tmp_q
  integer                       :: r,s,t,u,v,tmp_r,tmp_s,tmp_t,tmp_u,tmp_v
  integer                       :: pq,rs,tmp_pq,tmp_rs
  double precision              :: t1,t2,t3,t4,t5,t6
  ! hessian  : mo_num 4D double precision matrix containing the hessian before the permutations
  ! h_tmpr   : mo_num 4D double precision matrix containing the hessian after the permutations
  ! p,q,r,s  : integer, indexes of the 4D hessian matrix
  ! t,u,v    : integer, indexes to compute hessian elements
  ! pq,rs    : integer, indexes for the conversion from 4D to 2D hessian matrix
  ! t1,t2,t3 : double precision, t3 = t2 - t1, time to compute the hessian 

  ! Funtion 
  double precision              :: get_two_e_integral
  ! get_two_e_integral :  double precision function, two e integrals 

  ! Provided :
  ! mo_one_e_integrals : mono e- integrals
  ! get_two_e_integral : two e- integrals
  ! one_e_dm_mo_alpha, one_e_dm_mo_beta : one body density matrix
  ! two_e_dm_mo : two body density matrix

  !============
  ! Allocation
  !============

  allocate(hessian(m,m,m,m))

  !=============
  ! Calculation
  !=============

  print*,'---first_hess_list---'

  ! From Anderson et. al. (2014) 
  ! The Journal of Chemical Physics 141, 244104 (2014); doi: 10.1063/1.4904384

  CALL wall_time(t1)

  ! Initialization
  hessian = 0d0

  !========================
  ! First line, first term
  !========================

  CALL wall_time(t4)

  do tmp_p = 1, m
    p = list(tmp_p)
    do tmp_q = 1, m
      q = list(tmp_q)
      do tmp_r = 1, m
        r = list(tmp_r)
        do tmp_s = 1, m
          s = list(tmp_s)

            if (q==r) then
              do u = 1, mo_num

                hessian(tmp_p,tmp_q,tmp_r,tmp_s) = hessian(tmp_p,tmp_q,tmp_r,tmp_s) + 0.5d0 * (  &
                  mo_one_e_integrals(u,p) * one_e_dm_mo(u,s) &
                + mo_one_e_integrals(s,u) * one_e_dm_mo(p,u))

              enddo
            endif

        enddo
      enddo
    enddo
  enddo

  CALL wall_time(t5)
  t6 = t5-t4
  print*,'l1 1 :', t6

  !=========================
  ! First line, second term
  !=========================

  CALL wall_time(t4)

  do tmp_p = 1, m
    p = list(tmp_p)
    do tmp_q = 1, m
      q = list(tmp_q)
      do tmp_r = 1, m
        r = list(tmp_r)
        do tmp_s = 1, m
          s = list(tmp_s)

          if (p==s) then
            do u = 1, mo_num

                  hessian(tmp_p,tmp_q,tmp_r,tmp_s) = hessian(tmp_p,tmp_q,tmp_r,tmp_s) + 0.5d0 * ( &
                    mo_one_e_integrals(u,r) * one_e_dm_mo(u,q) &
                  + mo_one_e_integrals(q,u) * one_e_dm_mo(r,u))
            enddo
          endif

        enddo
      enddo
    enddo
  enddo
  
  CALL wall_time(t5)
  t6 = t5-t4
  print*,'l1 2 :', t6

  !========================
  ! First line, third term
  !========================

  CALL wall_time(t4)

  do tmp_p = 1, m
    p = list(tmp_p)
    do tmp_q = 1, m
      q = list(tmp_q)
      do tmp_r = 1, m
        r = list(tmp_r)
        do tmp_s = 1, m
          s = list(tmp_s)

          hessian(tmp_p,tmp_q,tmp_r,tmp_s) = hessian(tmp_p,tmp_q,tmp_r,tmp_s) &
          - mo_one_e_integrals(s,p) * one_e_dm_mo(r,q)&
          - mo_one_e_integrals(q,r) * one_e_dm_mo(p,s)

        enddo
      enddo
    enddo
  enddo

  CALL wall_time(t5)
  t6 = t5-t4
  print*,'l1 3 :', t6


  !=========================
  ! Second line, first term
  !=========================

  CALL wall_time(t4)

  do tmp_p = 1, m
    p = list(tmp_p)
    do tmp_q = 1, m
      q = list(tmp_q)
      do tmp_r = 1, m
        r = list(tmp_r)
        do tmp_s = 1, m
          s = list(tmp_s)

           if (q==r) then
             do t = 1, mo_num
               do u = 1, mo_num
                 do v = 1, mo_num

                   hessian(tmp_p,tmp_q,tmp_r,tmp_s) = hessian(tmp_p,tmp_q,tmp_r,tmp_s) + 0.5d0 * (  &
                     get_two_e_integral(u,v,p,t,mo_integrals_map) * two_e_dm_mo(u,v,s,t) &
                   + get_two_e_integral(s,t,u,v,mo_integrals_map) * two_e_dm_mo(p,t,u,v))

                 enddo
               enddo
             enddo
           endif

        enddo
      enddo
    enddo
  enddo

  CALL wall_time(t5)
  t6 = t5-t4
  print*,'l2 1 :', t6

  !==========================
  ! Second line, second term
  !==========================

  CALL wall_time(t4)

  do tmp_p = 1, m
    p = list(tmp_p)
    do tmp_q = 1, m
      q = list(tmp_q)
      do tmp_r = 1, m
        r = list(tmp_r)
        do tmp_s = 1, m
          s = list(tmp_s)

          if (p==s) then
            do t = 1, mo_num
              do u = 1, mo_num
                do v = 1, mo_num

                  hessian(tmp_p,tmp_q,tmp_r,tmp_s) = hessian(tmp_p,tmp_q,tmp_r,tmp_s) + 0.5d0 * ( &
                    get_two_e_integral(q,t,u,v,mo_integrals_map) * two_e_dm_mo(r,t,u,v) &
                  + get_two_e_integral(u,v,r,t,mo_integrals_map) * two_e_dm_mo(u,v,q,t))

                enddo
              enddo
            enddo
          endif

        enddo
      enddo
    enddo
  enddo

  CALL wall_time(t5)
  t6 = t5-t4
  print*,'l2 2 :', t6

  !========================
  ! Third line, first term
  !========================

  CALL wall_time(t4)

  do tmp_p = 1, m
    p = list(tmp_p)
    do tmp_q = 1, m
      q = list(tmp_q)
      do tmp_r = 1, m
        r = list(tmp_r)
        do tmp_s = 1, m
          s = list(tmp_s)

          do u = 1, mo_num
            do v = 1, mo_num

              hessian(tmp_p,tmp_q,tmp_r,tmp_s) = hessian(tmp_p,tmp_q,tmp_r,tmp_s) &
               + get_two_e_integral(u,v,p,r,mo_integrals_map) * two_e_dm_mo(u,v,q,s) &
               + get_two_e_integral(q,s,u,v,mo_integrals_map) * two_e_dm_mo(p,r,u,v)

            enddo
          enddo

        enddo
      enddo
    enddo
  enddo
 
  CALL wall_time(t5)
  t6 = t5-t4
  print*,'l3 1 :', t6

  !=========================
  ! Third line, second term
  !=========================

  CALL wall_time(t4)

  do tmp_p = 1, m
    p = list(tmp_p)
    do tmp_q = 1, m
      q = list(tmp_q)
      do tmp_r = 1, m
        r = list(tmp_r)
        do tmp_s = 1, m
          s = list(tmp_s)

          do t = 1, mo_num
            do u = 1, mo_num

              hessian(tmp_p,tmp_q,tmp_r,tmp_s) = hessian(tmp_p,tmp_q,tmp_r,tmp_s) &
               - get_two_e_integral(s,t,p,u,mo_integrals_map) * two_e_dm_mo(r,t,q,u) &
               - get_two_e_integral(t,s,p,u,mo_integrals_map) * two_e_dm_mo(t,r,q,u) &
               - get_two_e_integral(q,u,r,t,mo_integrals_map) * two_e_dm_mo(p,u,s,t) &
               - get_two_e_integral(q,u,t,r,mo_integrals_map) * two_e_dm_mo(p,u,t,s)

            enddo
          enddo

        enddo
      enddo
    enddo
  enddo

  CALL wall_time(t5)
  t6 = t5-t4
  print*,'l3 2 :', t6

  CALL wall_time(t2)
  t3 = t2 -t1
  print*,'Time to compute the hessian : ', t3

  !==============
  ! Permutations 
  !==============

  ! Hessian(p,q,r,s) = P_pq P_rs [ ...]
  ! => Hessian(p,q,r,s) = (p,q,r,s) - (q,p,r,s) - (p,q,s,r) + (q,p,s,r) 

  do tmp_s = 1, m
    do tmp_r = 1, m
      do tmp_q = 1, m
        do tmp_p = 1, m

          h_tmpr(tmp_p,tmp_q,tmp_r,tmp_s) = (hessian(tmp_p,tmp_q,tmp_r,tmp_s) - hessian(tmp_q,tmp_p,tmp_r,tmp_s) &
                                           - hessian(tmp_p,tmp_q,tmp_s,tmp_r) + hessian(tmp_q,tmp_p,tmp_s,tmp_r))

        enddo
      enddo
    enddo
  enddo

  !========================
  ! 4D matrix to 2D matrix
  !========================

  ! Convert the hessian mo_num * mo_num * mo_num * mo_num matrix in a
  ! 2D n * n matrix (n = mo_num*(mo_num-1)/2)
  ! H(pq,rs) : p<q and r<s

  ! 4D mo_num matrix to 2D n matrix
  do tmp_pq = 1, tmp_n
    call vec_to_mat_index(tmp_pq,tmp_p,tmp_q)
    do tmp_rs = 1, tmp_n
      call vec_to_mat_index(tmp_rs,tmp_r,tmp_s)
      H(tmp_pq,tmp_rs) = h_tmpr(tmp_p,tmp_q,tmp_r,tmp_s)   
    enddo
  enddo

  ! Display
  if (debug) then 
    print*,'2D Hessian matrix'
    do tmp_pq = 1, tmp_n
      write(*,'(100(F10.5))') H(tmp_pq,:)
    enddo 
  endif

  !==============
  ! Deallocation
  !==============

  deallocate(hessian)

  print*,'---End first_hess_list---'

end subroutine
