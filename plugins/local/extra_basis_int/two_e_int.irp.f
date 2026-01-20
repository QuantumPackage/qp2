double precision function ao_two_e_integral_mixed_direct(i, j, k, l)

  BEGIN_DOC
  !  integral of the AO basis <ik|jl> or (ij|kl)
  !     i(r1) j(r1) 1/r12 k(r2) l(r2)
  !     A     A           B     B
  !
  ! where i,j belong to the REGULAR AO basis (system A) and k,l to the EXTRA basis (system B)
  END_DOC

  implicit none
  include 'utils/constants.include.F'

  integer, intent(in)            :: i, j, k, l

  integer                        :: p, q, r, s
  integer                        :: num_i,num_j,num_k,num_l,dim1,I_power(3),J_power(3),K_power(3),L_power(3)
  integer                        :: iorder_p(3), iorder_q(3)
  double precision               :: I_center(3), J_center(3), K_center(3), L_center(3)
  double precision               :: integral
  double precision               :: P_new(0:max_dim,3),P_center(3),fact_p,pp
  double precision               :: Q_new(0:max_dim,3),Q_center(3),fact_q,qq
  double precision               :: general_primitive_integral

  dim1 = n_pt_max_integrals

  num_i = ao_nucl(i)
  num_j = ao_nucl(j)
  num_k = ao_extra_nucl(k)
  num_l = ao_extra_nucl(l)
  ao_two_e_integral_mixed_direct = 0.d0

  do p = 1, 3
    I_power(p) = ao_power(i,p)
    J_power(p) = ao_power(j,p)
    K_power(p) = ao_extra_power(k,p)
    L_power(p) = ao_extra_power(l,p)
    I_center(p) = nucl_coord(num_i,p)
    J_center(p) = nucl_coord(num_j,p)
    K_center(p) = extra_nucl_coord(num_k,p)
    L_center(p) = extra_nucl_coord(num_l,p)
  enddo

  double precision               :: coef1, coef2, coef3, coef4
  double precision               :: p_inv,q_inv

  do p = 1, ao_prim_num(i)
    coef1 = ao_coef_normalized_ordered_transp(p,i)
    do q = 1, ao_prim_num(j)
      coef2 = coef1*ao_coef_normalized_ordered_transp(q,j)
      call give_explicit_poly_and_gaussian(P_new,P_center,pp,fact_p,iorder_p,&
          ao_expo_ordered_transp(p,i),ao_expo_ordered_transp(q,j),                 &
          I_power,J_power,I_center,J_center,dim1)
      p_inv = 1.d0/pp
      do r = 1, ao_extra_prim_num(k)
        coef3 = coef2*ao_extra_coef_normalized_ordered_transp(r,k)
        do s = 1, ao_extra_prim_num(l)
          coef4 = coef3*ao_extra_coef_normalized_ordered_transp(s,l)
          call give_explicit_poly_and_gaussian(Q_new,Q_center,qq,fact_q,iorder_q,&
              ao_extra_expo_ordered_transp(r,k),ao_extra_expo_ordered_transp(s,l),             &
              K_power,L_power,K_center,L_center,dim1)
          q_inv = 1.d0/qq
          integral = general_primitive_integral(dim1,              &
              P_new,P_center,fact_p,pp,p_inv,iorder_p,             &
              Q_new,Q_center,fact_q,qq,q_inv,iorder_q)
          ao_two_e_integral_mixed_direct = ao_two_e_integral_mixed_direct +  coef4 * integral
        enddo ! s
      enddo  ! r
    enddo   ! q
  enddo    ! p

end

double precision function ao_two_e_integral_mixed_exchange(i, j, k, l)

  BEGIN_DOC
  !  integral of the AO basis <ik|jl> or (ij|kl)
  !     i(r1) j(r1) 1/r12 k(r2) l(r2)
  !     A     B           A     B
  !
  ! where i,k belong to the REGULAR AO basis (system A) and j,l to the EXTRA basis (system B)
  END_DOC

  implicit none
  include 'utils/constants.include.F'

  integer, intent(in)            :: i, j, k, l

  integer                        :: p, q, r, s
  integer                        :: num_i,num_j,num_k,num_l,dim1,I_power(3),J_power(3),K_power(3),L_power(3)
  integer                        :: iorder_p(3), iorder_q(3)
  double precision               :: I_center(3), J_center(3), K_center(3), L_center(3)
  double precision               :: integral
  double precision               :: P_new(0:max_dim,3),P_center(3),fact_p,pp
  double precision               :: Q_new(0:max_dim,3),Q_center(3),fact_q,qq
  double precision               :: general_primitive_integral

  dim1 = n_pt_max_integrals

  num_i = ao_nucl(i)
  num_j = ao_extra_nucl(j)
  num_k = ao_nucl(k)
  num_l = ao_extra_nucl(l)
  ao_two_e_integral_mixed_exchange = 0.d0

  do p = 1, 3
    I_power(p) = ao_power(i,p)
    J_power(p) = ao_extra_power(j,p)
    K_power(p) = ao_power(k,p)
    L_power(p) = ao_extra_power(l,p)
    I_center(p) = nucl_coord(num_i,p)
    J_center(p) = extra_nucl_coord(num_j,p)
    K_center(p) = nucl_coord(num_k,p)
    L_center(p) = extra_nucl_coord(num_l,p)
  enddo

  double precision               :: coef1, coef2, coef3, coef4
  double precision               :: p_inv,q_inv

  do p = 1, ao_prim_num(i)
    coef1 = ao_coef_normalized_ordered_transp(p,i)
    do q = 1, ao_extra_prim_num(j)
      coef2 = coef1*ao_extra_coef_normalized_ordered_transp(q,j)
      call give_explicit_poly_and_gaussian(P_new,P_center,pp,fact_p,iorder_p,&
          ao_expo_ordered_transp(p,i),ao_extra_expo_ordered_transp(q,j),                 &
          I_power,J_power,I_center,J_center,dim1)
      p_inv = 1.d0/pp
      do r = 1, ao_prim_num(k)
        coef3 = coef2*ao_coef_normalized_ordered_transp(r,k)
        do s = 1, ao_extra_prim_num(l)
          coef4 = coef3*ao_extra_coef_normalized_ordered_transp(s,l)
          call give_explicit_poly_and_gaussian(Q_new,Q_center,qq,fact_q,iorder_q,&
              ao_expo_ordered_transp(r,k),ao_extra_expo_ordered_transp(s,l),             &
              K_power,L_power,K_center,L_center,dim1)
          q_inv = 1.d0/qq
          integral = general_primitive_integral(dim1,              &
              P_new,P_center,fact_p,pp,p_inv,iorder_p,             &
              Q_new,Q_center,fact_q,qq,q_inv,iorder_q)
          ao_two_e_integral_mixed_exchange = ao_two_e_integral_mixed_exchange +  coef4 * integral
        enddo ! s
      enddo  ! r
    enddo   ! q
  enddo    ! p

end

double precision function ao_two_e_integral_mixed_3_a_1_b(i, j, k, l)

  BEGIN_DOC
  !  integral of the AO basis <ik|jl> or (ij|kl)
  !     i(r1) j(r1) 1/r12 k(r2) l(r2)
  !     A     A           A     B
  !
  ! where i,j belong to the REGULAR AO basis (system A) and k,l to the EXTRA basis (system B)
  END_DOC

  implicit none
  include 'utils/constants.include.F'

  integer, intent(in)            :: i, j, k, l

  integer                        :: p, q, r, s
  integer                        :: num_i,num_j,num_k,num_l,dim1,I_power(3),J_power(3),K_power(3),L_power(3)
  integer                        :: iorder_p(3), iorder_q(3)
  double precision               :: I_center(3), J_center(3), K_center(3), L_center(3)
  double precision               :: integral
  double precision               :: P_new(0:max_dim,3),P_center(3),fact_p,pp
  double precision               :: Q_new(0:max_dim,3),Q_center(3),fact_q,qq
  double precision               :: general_primitive_integral

  dim1 = n_pt_max_integrals

  num_i = ao_nucl(i)
  num_j = ao_nucl(j)
  num_k = ao_nucl(k)
  num_l = ao_extra_nucl(l)
  ao_two_e_integral_mixed_3_a_1_b = 0.d0

  do p = 1, 3
    I_power(p) = ao_power(i,p)
    J_power(p) = ao_power(j,p)
    K_power(p) = ao_power(k,p)
    L_power(p) = ao_extra_power(l,p)
    I_center(p) = nucl_coord(num_i,p)
    J_center(p) = nucl_coord(num_j,p)
    K_center(p) = nucl_coord(num_k,p)
    L_center(p) = extra_nucl_coord(num_l,p)
  enddo

  double precision               :: coef1, coef2, coef3, coef4
  double precision               :: p_inv,q_inv

  do p = 1, ao_prim_num(i)
    coef1 = ao_coef_normalized_ordered_transp(p,i)
    do q = 1, ao_prim_num(j)
      coef2 = coef1*ao_coef_normalized_ordered_transp(q,j)
      call give_explicit_poly_and_gaussian(P_new,P_center,pp,fact_p,iorder_p,&
          ao_expo_ordered_transp(p,i),ao_expo_ordered_transp(q,j),                 &
          I_power,J_power,I_center,J_center,dim1)
      p_inv = 1.d0/pp
      do r = 1, ao_prim_num(k)
        coef3 = coef2*ao_coef_normalized_ordered_transp(r,k)
        do s = 1, ao_extra_prim_num(l)
          coef4 = coef3*ao_extra_coef_normalized_ordered_transp(s,l)
          call give_explicit_poly_and_gaussian(Q_new,Q_center,qq,fact_q,iorder_q,&
              ao_expo_ordered_transp(r,k),ao_extra_expo_ordered_transp(s,l),             &
              K_power,L_power,K_center,L_center,dim1)
          q_inv = 1.d0/qq
          integral = general_primitive_integral(dim1,              &
              P_new,P_center,fact_p,pp,p_inv,iorder_p,             &
              Q_new,Q_center,fact_q,qq,q_inv,iorder_q)
          ao_two_e_integral_mixed_3_a_1_b = ao_two_e_integral_mixed_3_a_1_b +  coef4 * integral
        enddo ! s
      enddo  ! r
    enddo   ! q
  enddo    ! p

end

double precision function ao_two_e_integral_mixed_3_b_1_a(i, j, k, l)

  BEGIN_DOC
  !  integral of the AO basis <ik|jl> or (ij|kl)
  !     i(r1) j(r1) 1/r12 k(r2) l(r2)
  !     A     B           B     B
  !
  ! where i,j belong to the REGULAR AO basis (system A) and k,l to the EXTRA basis (system B)
  END_DOC

  implicit none
  include 'utils/constants.include.F'

  integer, intent(in)            :: i, j, k, l

  integer                        :: p, q, r, s
  integer                        :: num_i,num_j,num_k,num_l,dim1,I_power(3),J_power(3),K_power(3),L_power(3)
  integer                        :: iorder_p(3), iorder_q(3)
  double precision               :: I_center(3), J_center(3), K_center(3), L_center(3)
  double precision               :: integral
  double precision               :: P_new(0:max_dim,3),P_center(3),fact_p,pp
  double precision               :: Q_new(0:max_dim,3),Q_center(3),fact_q,qq
  double precision               :: general_primitive_integral

  dim1 = n_pt_max_integrals

  num_i = ao_nucl(i)
  num_j = ao_extra_nucl(j)
  num_k = ao_extra_nucl(k)
  num_l = ao_extra_nucl(l)
  ao_two_e_integral_mixed_3_b_1_a = 0.d0

  do p = 1, 3
    I_power(p) = ao_power(i,p)
    J_power(p) = ao_extra_power(j,p)
    K_power(p) = ao_extra_power(k,p)
    L_power(p) = ao_extra_power(l,p)
    I_center(p) = nucl_coord(num_i,p)
    J_center(p) = extra_nucl_coord(num_j,p)
    K_center(p) = extra_nucl_coord(num_k,p)
    L_center(p) = extra_nucl_coord(num_l,p)
  enddo

  double precision               :: coef1, coef2, coef3, coef4
  double precision               :: p_inv,q_inv

  do p = 1, ao_prim_num(i)
    coef1 = ao_coef_normalized_ordered_transp(p,i)
    do q = 1, ao_extra_prim_num(j)
      coef2 = coef1*ao_extra_coef_normalized_ordered_transp(q,j)
      call give_explicit_poly_and_gaussian(P_new,P_center,pp,fact_p,iorder_p,&
          ao_expo_ordered_transp(p,i),ao_extra_expo_ordered_transp(q,j),                 &
          I_power,J_power,I_center,J_center,dim1)
      p_inv = 1.d0/pp
      do r = 1, ao_extra_prim_num(k)
        coef3 = coef2*ao_extra_coef_normalized_ordered_transp(r,k)
        do s = 1, ao_extra_prim_num(l)
          coef4 = coef3*ao_extra_coef_normalized_ordered_transp(s,l)
          call give_explicit_poly_and_gaussian(Q_new,Q_center,qq,fact_q,iorder_q,&
              ao_extra_expo_ordered_transp(r,k),ao_extra_expo_ordered_transp(s,l),             &
              K_power,L_power,K_center,L_center,dim1)
          q_inv = 1.d0/qq
          integral = general_primitive_integral(dim1,              &
              P_new,P_center,fact_p,pp,p_inv,iorder_p,             &
              Q_new,Q_center,fact_q,qq,q_inv,iorder_q)
          ao_two_e_integral_mixed_3_b_1_a = ao_two_e_integral_mixed_3_b_1_a +  coef4 * integral
        enddo ! s
      enddo  ! r
    enddo   ! q
  enddo    ! p

end
