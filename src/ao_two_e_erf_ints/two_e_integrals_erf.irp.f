double precision function ao_two_e_integral_erf(i,j,k,l)
  implicit none
  BEGIN_DOC
  !  integral of the AO basis <ik|jl> or (ij|kl)
  !     i(r1) j(r1) 1/r12 k(r2) l(r2)
  END_DOC

  integer,intent(in)             :: i,j,k,l
  integer                        :: p,q,r,s
  double precision               :: I_center(3),J_center(3),K_center(3),L_center(3)
  integer                        :: num_i,num_j,num_k,num_l,dim1,I_power(3),J_power(3),K_power(3),L_power(3)
  double precision               :: integral
  include 'utils/constants.include.F'
  double precision               :: P_new(0:max_dim,3),P_center(3),fact_p,pp
  double precision               :: Q_new(0:max_dim,3),Q_center(3),fact_q,qq
  integer                        :: iorder_p(3), iorder_q(3)
  double precision               :: ao_two_e_integral_schwartz_accel_erf
  
  provide mu_erf

   if (ao_prim_num(i) * ao_prim_num(j) * ao_prim_num(k) * ao_prim_num(l) > 1024 ) then
     ao_two_e_integral_erf = ao_two_e_integral_schwartz_accel_erf(i,j,k,l)
     return
   endif

  dim1 = n_pt_max_integrals

  num_i = ao_nucl(i)
  num_j = ao_nucl(j)
  num_k = ao_nucl(k)
  num_l = ao_nucl(l)
  ao_two_e_integral_erf = 0.d0

  if (num_i /= num_j .or. num_k /= num_l .or. num_j /= num_k)then
    do p = 1, 3
      I_power(p) = ao_power(i,p)
      J_power(p) = ao_power(j,p)
      K_power(p) = ao_power(k,p)
      L_power(p) = ao_power(l,p)
      I_center(p) = nucl_coord(num_i,p)
      J_center(p) = nucl_coord(num_j,p)
      K_center(p) = nucl_coord(num_k,p)
      L_center(p) = nucl_coord(num_l,p)
    enddo

    double precision               :: coef1, coef2, coef3, coef4
    double precision               :: p_inv,q_inv
    double precision               :: general_primitive_integral_erf

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
          do s = 1, ao_prim_num(l)
            coef4 = coef3*ao_coef_normalized_ordered_transp(s,l)
            call give_explicit_poly_and_gaussian(Q_new,Q_center,qq,fact_q,iorder_q,&
                ao_expo_ordered_transp(r,k),ao_expo_ordered_transp(s,l),             &
                K_power,L_power,K_center,L_center,dim1)
            q_inv = 1.d0/qq
            integral = general_primitive_integral_erf(dim1,              &
                P_new,P_center,fact_p,pp,p_inv,iorder_p,             &
                Q_new,Q_center,fact_q,qq,q_inv,iorder_q)
            ao_two_e_integral_erf = ao_two_e_integral_erf +  coef4 * integral
          enddo ! s
        enddo  ! r
      enddo   ! q
    enddo    ! p

  else

    do p = 1, 3
      I_power(p) = ao_power(i,p)
      J_power(p) = ao_power(j,p)
      K_power(p) = ao_power(k,p)
      L_power(p) = ao_power(l,p)
    enddo
    double  precision              :: ERI_erf

    do p = 1, ao_prim_num(i)
      coef1 = ao_coef_normalized_ordered_transp(p,i)
      do q = 1, ao_prim_num(j)
        coef2 = coef1*ao_coef_normalized_ordered_transp(q,j)
        do r = 1, ao_prim_num(k)
          coef3 = coef2*ao_coef_normalized_ordered_transp(r,k)
          do s = 1, ao_prim_num(l)
            coef4 = coef3*ao_coef_normalized_ordered_transp(s,l)
            integral = ERI_erf(                                          &
                ao_expo_ordered_transp(p,i),ao_expo_ordered_transp(q,j),ao_expo_ordered_transp(r,k),ao_expo_ordered_transp(s,l),&
                I_power(1),J_power(1),K_power(1),L_power(1),         &
                I_power(2),J_power(2),K_power(2),L_power(2),         &
                I_power(3),J_power(3),K_power(3),L_power(3))
            ao_two_e_integral_erf = ao_two_e_integral_erf + coef4 * integral
          enddo ! s
        enddo  ! r
      enddo   ! q
    enddo    ! p

  endif

end

double precision function ao_two_e_integral_schwartz_accel_erf(i,j,k,l)
  implicit none
  BEGIN_DOC
  !  integral of the AO basis <ik|jl> or (ij|kl)
  !     i(r1) j(r1) 1/r12 k(r2) l(r2)
  END_DOC
  integer,intent(in)             :: i,j,k,l
  integer                        :: p,q,r,s
  double precision               :: I_center(3),J_center(3),K_center(3),L_center(3)
  integer                        :: num_i,num_j,num_k,num_l,dim1,I_power(3),J_power(3),K_power(3),L_power(3)
  double precision               :: integral
  include 'utils/constants.include.F'
  double precision               :: P_new(0:max_dim,3),P_center(3),fact_p,pp
  double precision               :: Q_new(0:max_dim,3),Q_center(3),fact_q,qq
  integer                        :: iorder_p(3), iorder_q(3)
  double precision, allocatable  :: schwartz_kl(:,:)
  double precision               :: schwartz_ij

  dim1 = n_pt_max_integrals

  num_i = ao_nucl(i)
  num_j = ao_nucl(j)
  num_k = ao_nucl(k)
  num_l = ao_nucl(l)
  ao_two_e_integral_schwartz_accel_erf = 0.d0
  double precision               :: thr
  thr = ao_integrals_threshold*ao_integrals_threshold

  allocate(schwartz_kl(0:ao_prim_num(l),0:ao_prim_num(k)))

      double precision               :: coef3
      double precision               :: coef2
      double precision               :: p_inv,q_inv
      double precision               :: coef1
      double precision               :: coef4

  if (num_i /= num_j .or. num_k /= num_l .or. num_j /= num_k)then
    do p = 1, 3
      I_power(p) = ao_power(i,p)
      J_power(p) = ao_power(j,p)
      K_power(p) = ao_power(k,p)
      L_power(p) = ao_power(l,p)
      I_center(p) = nucl_coord(num_i,p)
      J_center(p) = nucl_coord(num_j,p)
      K_center(p) = nucl_coord(num_k,p)
      L_center(p) = nucl_coord(num_l,p)
    enddo

    schwartz_kl(0,0) = 0.d0
    do r = 1, ao_prim_num(k)
      coef1 = ao_coef_normalized_ordered_transp(r,k)*ao_coef_normalized_ordered_transp(r,k)
      schwartz_kl(0,r) = 0.d0
      do s = 1, ao_prim_num(l)
        coef2 = coef1 * ao_coef_normalized_ordered_transp(s,l) * ao_coef_normalized_ordered_transp(s,l)
        call give_explicit_poly_and_gaussian(Q_new,Q_center,qq,fact_q,iorder_q,&
            ao_expo_ordered_transp(r,k),ao_expo_ordered_transp(s,l),                 &
            K_power,L_power,K_center,L_center,dim1)
        q_inv = 1.d0/qq
        schwartz_kl(s,r) = general_primitive_integral_erf(dim1,          &
            Q_new,Q_center,fact_q,qq,q_inv,iorder_q,                 &
            Q_new,Q_center,fact_q,qq,q_inv,iorder_q)      &
            * coef2
        schwartz_kl(0,r) = max(schwartz_kl(0,r),schwartz_kl(s,r))
      enddo
      schwartz_kl(0,0) = max(schwartz_kl(0,r),schwartz_kl(0,0))
    enddo

    do p = 1, ao_prim_num(i)
      coef1 = ao_coef_normalized_ordered_transp(p,i)
      do q = 1, ao_prim_num(j)
        coef2 = coef1*ao_coef_normalized_ordered_transp(q,j)
        call give_explicit_poly_and_gaussian(P_new,P_center,pp,fact_p,iorder_p,&
            ao_expo_ordered_transp(p,i),ao_expo_ordered_transp(q,j),                 &
            I_power,J_power,I_center,J_center,dim1)
        p_inv = 1.d0/pp
        schwartz_ij = general_primitive_integral_erf(dim1,               &
            P_new,P_center,fact_p,pp,p_inv,iorder_p,                 &
            P_new,P_center,fact_p,pp,p_inv,iorder_p) *               &
            coef2*coef2
        if (schwartz_kl(0,0)*schwartz_ij < thr) then
           cycle
        endif
        do r = 1, ao_prim_num(k)
          if (schwartz_kl(0,r)*schwartz_ij < thr) then
             cycle
          endif
          coef3 = coef2*ao_coef_normalized_ordered_transp(r,k)
          do s = 1, ao_prim_num(l)
            if (schwartz_kl(s,r)*schwartz_ij < thr) then
               cycle
            endif
            coef4 = coef3*ao_coef_normalized_ordered_transp(s,l)
            double precision               :: general_primitive_integral_erf
            call give_explicit_poly_and_gaussian(Q_new,Q_center,qq,fact_q,iorder_q,&
                ao_expo_ordered_transp(r,k),ao_expo_ordered_transp(s,l),             &
                K_power,L_power,K_center,L_center,dim1)
            q_inv = 1.d0/qq
            integral = general_primitive_integral_erf(dim1,              &
                P_new,P_center,fact_p,pp,p_inv,iorder_p,             &
                Q_new,Q_center,fact_q,qq,q_inv,iorder_q)
            ao_two_e_integral_schwartz_accel_erf = ao_two_e_integral_schwartz_accel_erf + coef4 * integral
          enddo ! s
        enddo  ! r
      enddo   ! q
    enddo    ! p

  else

    do p = 1, 3
      I_power(p) = ao_power(i,p)
      J_power(p) = ao_power(j,p)
      K_power(p) = ao_power(k,p)
      L_power(p) = ao_power(l,p)
    enddo
    double  precision              :: ERI_erf

    schwartz_kl(0,0) = 0.d0
    do r = 1, ao_prim_num(k)
      coef1 = ao_coef_normalized_ordered_transp(r,k)*ao_coef_normalized_ordered_transp(r,k)
      schwartz_kl(0,r) = 0.d0
      do s = 1, ao_prim_num(l)
        coef2 = coef1*ao_coef_normalized_ordered_transp(s,l)*ao_coef_normalized_ordered_transp(s,l)
        schwartz_kl(s,r) = ERI_erf(                                      &
            ao_expo_ordered_transp(r,k),ao_expo_ordered_transp(s,l),ao_expo_ordered_transp(r,k),ao_expo_ordered_transp(s,l),&
            K_power(1),L_power(1),K_power(1),L_power(1),             &
            K_power(2),L_power(2),K_power(2),L_power(2),             &
            K_power(3),L_power(3),K_power(3),L_power(3)) * &
            coef2
        schwartz_kl(0,r) = max(schwartz_kl(0,r),schwartz_kl(s,r))
      enddo
      schwartz_kl(0,0) = max(schwartz_kl(0,r),schwartz_kl(0,0))
    enddo

    do p = 1, ao_prim_num(i)
      coef1 = ao_coef_normalized_ordered_transp(p,i)
      do q = 1, ao_prim_num(j)
        coef2 = coef1*ao_coef_normalized_ordered_transp(q,j)
        schwartz_ij = ERI_erf(                                          &
                ao_expo_ordered_transp(p,i),ao_expo_ordered_transp(q,j),ao_expo_ordered_transp(p,i),ao_expo_ordered_transp(q,j),&
                I_power(1),J_power(1),I_power(1),J_power(1),         &
                I_power(2),J_power(2),I_power(2),J_power(2),         &
                I_power(3),J_power(3),I_power(3),J_power(3))*coef2*coef2
        if (schwartz_kl(0,0)*schwartz_ij < thr) then
           cycle
        endif
        do r = 1, ao_prim_num(k)
          if (schwartz_kl(0,r)*schwartz_ij < thr) then
             cycle
          endif
          coef3 = coef2*ao_coef_normalized_ordered_transp(r,k)
          do s = 1, ao_prim_num(l)
            if (schwartz_kl(s,r)*schwartz_ij < thr) then
               cycle
            endif
            coef4 = coef3*ao_coef_normalized_ordered_transp(s,l)
            integral = ERI_erf(                                          &
                ao_expo_ordered_transp(p,i),ao_expo_ordered_transp(q,j),ao_expo_ordered_transp(r,k),ao_expo_ordered_transp(s,l),&
                I_power(1),J_power(1),K_power(1),L_power(1),         &
                I_power(2),J_power(2),K_power(2),L_power(2),         &
                I_power(3),J_power(3),K_power(3),L_power(3))
            ao_two_e_integral_schwartz_accel_erf = ao_two_e_integral_schwartz_accel_erf +  coef4 * integral
          enddo ! s
        enddo  ! r
      enddo   ! q
    enddo    ! p

  endif
  deallocate (schwartz_kl)

end


subroutine compute_ao_two_e_integrals_erf(j,k,l,sze,buffer_value)
  implicit none
  use map_module

  BEGIN_DOC
  ! Compute AO 1/r12 integrals for all i and fixed j,k,l
  END_DOC

  include 'utils/constants.include.F'
  integer, intent(in)            :: j,k,l,sze
  real(integral_kind), intent(out) :: buffer_value(sze)
  double precision               :: ao_two_e_integral_erf

  integer                        :: i
  logical, external              :: ao_one_e_integral_zero
  logical, external              :: ao_two_e_integral_zero

  if (ao_one_e_integral_zero(j,l)) then
    buffer_value = 0._integral_kind
    return
  endif
  if (ao_two_e_integral_erf_schwartz(j,l) < thresh ) then
    buffer_value = 0._integral_kind
    return
  endif

  do i = 1, ao_num
    if (ao_two_e_integral_zero(i,j,k,l)) then
      buffer_value(i) = 0._integral_kind
      cycle
    endif
    if (ao_two_e_integral_erf_schwartz(i,k)*ao_two_e_integral_erf_schwartz(j,l) < thresh ) then
      buffer_value(i) = 0._integral_kind
      cycle
    endif
    !DIR$ FORCEINLINE
    buffer_value(i) = ao_two_e_integral_erf(i,k,j,l)
  enddo

end

double precision function general_primitive_integral_erf(dim,            &
      P_new,P_center,fact_p,p,p_inv,iorder_p,                        &
      Q_new,Q_center,fact_q,q,q_inv,iorder_q)
  implicit none
  BEGIN_DOC
  ! Computes the integral <pq|rs> where p,q,r,s are Gaussian primitives
  END_DOC
  integer,intent(in)             :: dim
  include 'utils/constants.include.F'
  double precision, intent(in)   :: P_new(0:max_dim,3),P_center(3),fact_p,p,p_inv
  double precision, intent(in)   :: Q_new(0:max_dim,3),Q_center(3),fact_q,q,q_inv
  integer, intent(in)            :: iorder_p(3)
  integer, intent(in)            :: iorder_q(3)

  double precision               :: r_cut,gama_r_cut,rho,dist
  double precision               :: dx(0:max_dim),Ix_pol(0:max_dim),dy(0:max_dim),Iy_pol(0:max_dim),dz(0:max_dim),Iz_pol(0:max_dim)
  integer                        :: n_Ix,n_Iy,n_Iz,nx,ny,nz
  double precision               :: bla
  integer                        :: ix,iy,iz,jx,jy,jz,i
  double precision               :: a,b,c,d,e,f,accu,pq,const
  double precision               :: pq_inv, p10_1, p10_2, p01_1, p01_2,pq_inv_2
  integer                        :: n_pt_tmp,n_pt_out, iorder
  double precision               :: d1(0:max_dim),d_poly(0:max_dim),rint,d1_screened(0:max_dim)

  general_primitive_integral_erf = 0.d0

  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: dx,Ix_pol,dy,Iy_pol,dz,Iz_pol
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: d1, d_poly

  ! Gaussian Product
  ! ----------------
  double precision :: p_plus_q
  p_plus_q = (p+q) * ((p*q)/(p+q) + mu_erf*mu_erf)/(mu_erf*mu_erf)
  pq = p_inv*0.5d0*q_inv

  pq_inv = 0.5d0/p_plus_q
  p10_1 = q*pq  ! 1/(2p)
  p01_1 = p*pq  ! 1/(2q)
  pq_inv_2 = pq_inv+pq_inv
  p10_2 = pq_inv_2 * p10_1*q !0.5d0*q/(pq + p*p)
  p01_2 = pq_inv_2 * p01_1*p !0.5d0*p/(q*q + pq)


  accu = 0.d0
  iorder = iorder_p(1)+iorder_q(1)+iorder_p(1)+iorder_q(1)
  !DIR$ VECTOR ALIGNED
  do ix=0,iorder
    Ix_pol(ix) = 0.d0
  enddo
  n_Ix = 0
  do ix = 0, iorder_p(1)
    if (abs(P_new(ix,1)) < thresh) cycle
    a = P_new(ix,1)
    do jx = 0, iorder_q(1)
      d = a*Q_new(jx,1)
      if (abs(d) < thresh) cycle
      !DEC$ FORCEINLINE
      call give_polynom_mult_center_x(P_center(1),Q_center(1),ix,jx,p,q,iorder,pq_inv,pq_inv_2,p10_1,p01_1,p10_2,p01_2,dx,nx)
      !DEC$ FORCEINLINE
      call add_poly_multiply(dx,nx,d,Ix_pol,n_Ix)
    enddo
  enddo
  if (n_Ix == -1) then
    return
  endif
  iorder = iorder_p(2)+iorder_q(2)+iorder_p(2)+iorder_q(2)
  !DIR$ VECTOR ALIGNED
  do ix=0, iorder
    Iy_pol(ix) = 0.d0
  enddo
  n_Iy = 0
  do iy = 0, iorder_p(2)
    if (abs(P_new(iy,2)) > thresh) then
      b = P_new(iy,2)
      do jy = 0, iorder_q(2)
        e = b*Q_new(jy,2)
        if (abs(e) < thresh) cycle
        !DEC$ FORCEINLINE
        call   give_polynom_mult_center_x(P_center(2),Q_center(2),iy,jy,p,q,iorder,pq_inv,pq_inv_2,p10_1,p01_1,p10_2,p01_2,dy,ny)
        !DEC$ FORCEINLINE
        call add_poly_multiply(dy,ny,e,Iy_pol,n_Iy)
      enddo
    endif
  enddo
  if (n_Iy == -1) then
    return
  endif

  iorder = iorder_p(3)+iorder_q(3)+iorder_p(3)+iorder_q(3)
  do ix=0,iorder
    Iz_pol(ix) = 0.d0
  enddo
  n_Iz = 0
  do iz = 0, iorder_p(3)
    if (abs(P_new(iz,3)) > thresh) then
      c = P_new(iz,3)
      do jz = 0, iorder_q(3)
        f = c*Q_new(jz,3)
        if (abs(f) < thresh) cycle
        !DEC$ FORCEINLINE
        call   give_polynom_mult_center_x(P_center(3),Q_center(3),iz,jz,p,q,iorder,pq_inv,pq_inv_2,p10_1,p01_1,p10_2,p01_2,dz,nz)
        !DEC$ FORCEINLINE
        call add_poly_multiply(dz,nz,f,Iz_pol,n_Iz)
      enddo
    endif
  enddo
  if (n_Iz == -1) then
    return
  endif

  rho = p*q *pq_inv_2  ! le rho qui va bien
  dist =  (P_center(1) - Q_center(1))*(P_center(1) - Q_center(1)) +  &
      (P_center(2) - Q_center(2))*(P_center(2) - Q_center(2)) +      &
      (P_center(3) - Q_center(3))*(P_center(3) - Q_center(3))
  const = dist*rho

  n_pt_tmp = n_Ix+n_Iy
  do i=0,n_pt_tmp
    d_poly(i)=0.d0
  enddo

  !DEC$ FORCEINLINE
  call multiply_poly(Ix_pol,n_Ix,Iy_pol,n_Iy,d_poly,n_pt_tmp)
  if (n_pt_tmp == -1) then
    return
  endif
  n_pt_out = n_pt_tmp+n_Iz
  do i=0,n_pt_out
    d1(i)=0.d0
  enddo

  !DEC$ FORCEINLINE
  call multiply_poly(d_poly ,n_pt_tmp ,Iz_pol,n_Iz,d1,n_pt_out)
  double precision               :: rint_sum
  accu = accu + rint_sum(n_pt_out,const,d1)

  ! change p+q in dsqrt
  general_primitive_integral_erf = fact_p * fact_q * accu *pi_5_2*p_inv*q_inv/dsqrt(p_plus_q)
end


double precision function ERI_erf(alpha,beta,delta,gama,a_x,b_x,c_x,d_x,a_y,b_y,c_y,d_y,a_z,b_z,c_z,d_z)
  implicit none
  BEGIN_DOC
  ! Atomic primtive two-electron integral between the 4 primitives :
  !
  ! * primitive 1 : $x_1^{a_x} y_1^{a_y} z_1^{a_z} \exp(-\alpha * r1^2)$
  ! * primitive 2 : $x_1^{b_x} y_1^{b_y} z_1^{b_z} \exp(- \beta * r1^2)$
  ! * primitive 3 : $x_2^{c_x} y_2^{c_y} z_2^{c_z} \exp(-\delta * r2^2)$
  ! * primitive 4 : $x_2^{d_x} y_2^{d_y} z_2^{d_z} \exp(-\gamma * r2^2)$
  !
  END_DOC
  double precision, intent(in)   :: delta,gama,alpha,beta
  integer, intent(in)            :: a_x,b_x,c_x,d_x,a_y,b_y,c_y,d_y,a_z,b_z,c_z,d_z
  integer                        :: a_x_2,b_x_2,c_x_2,d_x_2,a_y_2,b_y_2,c_y_2,d_y_2,a_z_2,b_z_2,c_z_2,d_z_2
  integer                        :: i,j,k,l,n_pt
  integer                        :: n_pt_sup
  double precision               :: p,q,denom,coeff
  double precision               :: I_f
  integer                        :: nx,ny,nz
  include 'utils/constants.include.F'
  nx = a_x+b_x+c_x+d_x
  if(iand(nx,1) == 1) then
    ERI_erf = 0.d0
    return
  endif

  ny = a_y+b_y+c_y+d_y
  if(iand(ny,1) == 1) then
    ERI_erf = 0.d0
    return
  endif

  nz = a_z+b_z+c_z+d_z
  if(iand(nz,1) == 1) then
    ERI_erf = 0.d0
    return
  endif

  ASSERT (alpha >= 0.d0)
  ASSERT (beta >= 0.d0)
  ASSERT (delta >= 0.d0)
  ASSERT (gama >= 0.d0)
  p = alpha + beta
  q = delta + gama
  double precision :: p_plus_q
  p_plus_q = (p+q) * ((p*q)/(p+q) + mu_erf*mu_erf)/(mu_erf*mu_erf)
  ASSERT (p+q >= 0.d0)
  n_pt =  ishft( nx+ny+nz,1 )

  coeff = pi_5_2 / (p * q * dsqrt(p_plus_q))
  if (n_pt == 0) then
    ERI_erf = coeff
    return
  endif

  call integrale_new_erf(I_f,a_x,b_x,c_x,d_x,a_y,b_y,c_y,d_y,a_z,b_z,c_z,d_z,p,q,n_pt)

  ERI_erf = I_f * coeff
end



subroutine integrale_new_erf(I_f,a_x,b_x,c_x,d_x,a_y,b_y,c_y,d_y,a_z,b_z,c_z,d_z,p,q,n_pt)
  BEGIN_DOC
  ! Calculate the integral of the polynomial :
  !
  ! $I_x1(a_x+b_x, c_x+d_x,p,q) \, I_x1(a_y+b_y, c_y+d_y,p,q) \, I_x1(a_z+b_z, c_z+d_z,p,q)$
  !
  ! between $( 0 ; 1)$
  END_DOC


  implicit none
  include 'utils/constants.include.F'
  double precision               :: p,q
  integer                        :: a_x,b_x,c_x,d_x,a_y,b_y,c_y,d_y,a_z,b_z,c_z,d_z
  integer                        :: i, n_pt, j
  double precision               :: I_f, pq_inv, p10_1, p10_2, p01_1, p01_2,rho,pq_inv_2
  integer :: ix,iy,iz, jx,jy,jz, sx,sy,sz

  j = ishft(n_pt,-1)
  ASSERT (n_pt > 1)
  double precision :: p_plus_q
  p_plus_q = (p+q) * ((p*q)/(p+q) + mu_erf*mu_erf)/(mu_erf*mu_erf)

  pq_inv = 0.5d0/(p_plus_q)
  pq_inv_2 = pq_inv + pq_inv
  p10_1 = 0.5d0/p
  p01_1 = 0.5d0/q
  p10_2 = 0.5d0 *  q /(p * p_plus_q)
  p01_2 = 0.5d0 *  p /(q * p_plus_q)
  double precision               :: B00(n_pt_max_integrals)
  double precision               :: B10(n_pt_max_integrals), B01(n_pt_max_integrals)
  double precision               :: t1(n_pt_max_integrals), t2(n_pt_max_integrals)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: t1, t2, B10, B01, B00
  ix = a_x+b_x
  jx = c_x+d_x
  iy = a_y+b_y
  jy = c_y+d_y
  iz = a_z+b_z
  jz = c_z+d_z
  sx = ix+jx
  sy = iy+jy
  sz = iz+jz

  !DIR$ VECTOR ALIGNED
  do i = 1,n_pt
    B10(i)  = p10_1 -  gauleg_t2(i,j)* p10_2
    B01(i)  = p01_1 -  gauleg_t2(i,j)* p01_2
    B00(i)  = gauleg_t2(i,j)*pq_inv
  enddo
  if (sx > 0) then
    call I_x1_new(ix,jx,B10,B01,B00,t1,n_pt)
  else
    !DIR$ VECTOR ALIGNED
    do i = 1,n_pt
      t1(i) = 1.d0
    enddo
  endif
  if (sy > 0) then
    call I_x1_new(iy,jy,B10,B01,B00,t2,n_pt)
    !DIR$ VECTOR ALIGNED
    do i = 1,n_pt
      t1(i) = t1(i)*t2(i)
    enddo
  endif
  if (sz > 0) then
    call I_x1_new(iz,jz,B10,B01,B00,t2,n_pt)
    !DIR$ VECTOR ALIGNED
    do i = 1,n_pt
      t1(i) = t1(i)*t2(i)
    enddo
  endif
  I_f= 0.d0
  !DIR$ VECTOR ALIGNED
  do i = 1,n_pt
    I_f += gauleg_w(i,j)*t1(i)
  enddo



end


subroutine compute_ao_integrals_erf_jl(j,l,n_integrals,buffer_i,buffer_value)
  implicit none
  use map_module
  BEGIN_DOC
  !  Parallel client for AO integrals
  END_DOC

  integer, intent(in)            :: j,l
  integer,intent(out)            :: n_integrals
  integer(key_kind),intent(out)  :: buffer_i(ao_num*ao_num)
  real(integral_kind),intent(out) :: buffer_value(ao_num*ao_num)

  integer                        :: i,k
  double precision               :: ao_two_e_integral_erf,cpu_1,cpu_2, wall_1, wall_2
  double precision               :: integral, wall_0
  double precision               :: thr
  integer                        :: kk, m, j1, i1
  logical, external              :: ao_two_e_integral_zero

  thr = ao_integrals_threshold

  n_integrals = 0

  j1 = j+ishft(l*l-l,-1)
  do k = 1, ao_num           ! r1
    i1 = ishft(k*k-k,-1)
    if (i1 > j1) then
      exit
    endif
    do i = 1, k
      i1 += 1
      if (i1 > j1) then
        exit
      endif
      if (ao_two_e_integral_zero(i,j,k,l)) then
        cycle
      endif
      if (ao_two_e_integral_erf_schwartz(i,k)*ao_two_e_integral_erf_schwartz(j,l) < thr ) then
        cycle
      endif
      !DIR$ FORCEINLINE
      integral = ao_two_e_integral_erf(i,k,j,l)  ! i,k : r1    j,l : r2
      if (abs(integral) < thr) then
        cycle
      endif
      n_integrals += 1
      !DIR$ FORCEINLINE
      call two_e_integrals_index(i,j,k,l,buffer_i(n_integrals))
      buffer_value(n_integrals) = integral
    enddo
  enddo

end
