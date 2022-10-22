double precision function ao_tc_sym_two_e_pot(i,j,k,l)
  implicit none
  BEGIN_DOC
  !  integral of the AO basis <ik|jl> or (ij|kl)
  !     i(r1) j(r1) (tc_pot(r12,mu)) k(r2) l(r2)
  !
  ! where (tc_pot(r12,mu)) is the scalar part of the potential EXCLUDING the term erf(mu r12)/r12. 
  !
  ! See Eq. (32) of JCP 154, 084119 (2021). 
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
  double precision :: scw_gauss_int,general_primitive_integral_gauss

  dim1 = n_pt_max_integrals

  num_i = ao_nucl(i)
  num_j = ao_nucl(j)
  num_k = ao_nucl(k)
  num_l = ao_nucl(l)
  ao_tc_sym_two_e_pot = 0.d0
  double precision               :: thr
  thr = ao_integrals_threshold*ao_integrals_threshold

  allocate(schwartz_kl(0:ao_prim_num(l),0:ao_prim_num(k)))

      double precision               :: coef3
      double precision               :: coef2
      double precision               :: p_inv,q_inv
      double precision               :: coef1
      double precision               :: coef4

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
        scw_gauss_int = general_primitive_integral_gauss(dim1,              &
                Q_new,Q_center,fact_q,qq,q_inv,iorder_q,             &
                Q_new,Q_center,fact_q,qq,q_inv,iorder_q)

        schwartz_kl(s,r) = dabs(scw_gauss_int * coef2)
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
        scw_gauss_int = general_primitive_integral_gauss(dim1,              &
                P_new,P_center,fact_p,pp,p_inv,iorder_p,             &
                P_new,P_center,fact_p,pp,p_inv,iorder_p)
        schwartz_ij = dabs(scw_gauss_int * coef2*coef2)
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
            call give_explicit_poly_and_gaussian(Q_new,Q_center,qq,fact_q,iorder_q, &
                ao_expo_ordered_transp(r,k),ao_expo_ordered_transp(s,l),            &
                K_power,L_power,K_center,L_center,dim1)
            q_inv = 1.d0/qq
            integral = general_primitive_integral_gauss(dim1,              &
                P_new,P_center,fact_p,pp,p_inv,iorder_p,             &
                Q_new,Q_center,fact_q,qq,q_inv,iorder_q)
            ao_tc_sym_two_e_pot = ao_tc_sym_two_e_pot + coef4 * integral
          enddo ! s
        enddo  ! r
      enddo   ! q
    enddo    ! p

  deallocate (schwartz_kl)

end


double precision function general_primitive_integral_gauss(dim,      &
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
  double precision :: thr

  thr = ao_integrals_threshold

  general_primitive_integral_gauss = 0.d0

  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: dx,Ix_pol,dy,Iy_pol,dz,Iz_pol
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: d1, d_poly

  ! Gaussian Product
  ! ----------------

  pq = p_inv*0.5d0*q_inv
  pq_inv = 0.5d0/(p+q)
  p10_1 = q*pq  ! 1/(2p)
  p01_1 = p*pq  ! 1/(2q)
  pq_inv_2 = pq_inv+pq_inv
  p10_2 = pq_inv_2 * p10_1*q !0.5d0*q/(pq + p*p)
  p01_2 = pq_inv_2 * p01_1*p !0.5d0*p/(q*q + pq)


  accu = 0.d0
  iorder = iorder_p(1)+iorder_q(1)+iorder_p(1)+iorder_q(1)
  do ix=0,iorder
    Ix_pol(ix) = 0.d0
  enddo
  n_Ix = 0
  do ix = 0, iorder_p(1)
    if (abs(P_new(ix,1)) < thr) cycle
    a = P_new(ix,1)
    do jx = 0, iorder_q(1)
      d = a*Q_new(jx,1)
      if (abs(d) < thr) cycle
      !DIR$ FORCEINLINE
      call give_polynom_mult_center_x(P_center(1),Q_center(1),ix,jx,p,q,iorder,pq_inv,pq_inv_2,p10_1,p01_1,p10_2,p01_2,dx,nx)
      !DIR$ FORCEINLINE
      call add_poly_multiply(dx,nx,d,Ix_pol,n_Ix)
    enddo
  enddo
  if (n_Ix == -1) then
    return
  endif
  iorder = iorder_p(2)+iorder_q(2)+iorder_p(2)+iorder_q(2)
  do ix=0, iorder
    Iy_pol(ix) = 0.d0
  enddo
  n_Iy = 0
  do iy = 0, iorder_p(2)
    if (abs(P_new(iy,2)) > thr) then
      b = P_new(iy,2)
      do jy = 0, iorder_q(2)
        e = b*Q_new(jy,2)
        if (abs(e) < thr) cycle
        !DIR$ FORCEINLINE
        call   give_polynom_mult_center_x(P_center(2),Q_center(2),iy,jy,p,q,iorder,pq_inv,pq_inv_2,p10_1,p01_1,p10_2,p01_2,dy,ny)
        !DIR$ FORCEINLINE
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
    if (abs(P_new(iz,3)) > thr) then
      c = P_new(iz,3)
      do jz = 0, iorder_q(3)
        f = c*Q_new(jz,3)
        if (abs(f) < thr) cycle
        !DIR$ FORCEINLINE
        call   give_polynom_mult_center_x(P_center(3),Q_center(3),iz,jz,p,q,iorder,pq_inv,pq_inv_2,p10_1,p01_1,p10_2,p01_2,dz,nz)
        !DIR$ FORCEINLINE
        call add_poly_multiply(dz,nz,f,Iz_pol,n_Iz)
      enddo
    endif
  enddo
  if (n_Iz == -1) then
    return
  endif

  rho = p*q *pq_inv_2
  dist =  (P_center(1) - Q_center(1))*(P_center(1) - Q_center(1)) +  &
      (P_center(2) - Q_center(2))*(P_center(2) - Q_center(2)) +      &
      (P_center(3) - Q_center(3))*(P_center(3) - Q_center(3))
  const = dist*rho

  n_pt_tmp = n_Ix+n_Iy
  do i=0,n_pt_tmp
    d_poly(i)=0.d0
  enddo

  !DIR$ FORCEINLINE
  call multiply_poly(Ix_pol,n_Ix,Iy_pol,n_Iy,d_poly,n_pt_tmp)
  if (n_pt_tmp == -1) then
    return
  endif
  n_pt_out = n_pt_tmp+n_Iz
  do i=0,n_pt_out
    d1(i)=0.d0
  enddo

  !DIR$ FORCEINLINE
  call multiply_poly(d_poly ,n_pt_tmp ,Iz_pol,n_Iz,d1,n_pt_out)

  double precision :: aa,c_a,t_a,rho_old,w_a,pi_3,prefactor,inv_pq_3_2
  double precision :: gauss_int
  integer :: m
  gauss_int = 0.d0
  pi_3 = pi*pi*pi
  inv_pq_3_2 = (p_inv * q_inv)**(1.5d0)
  rho_old = (p*q)/(p+q)
  prefactor = pi_3 * inv_pq_3_2 * fact_p * fact_q 
  do i = 1, n_gauss_eff_pot ! browse the gaussians with different expo/coef
  !do i = 1, n_gauss_eff_pot-1
   aa = expo_gauss_eff_pot(i) 
   c_a = coef_gauss_eff_pot(i)
   t_a = dsqrt( aa /(rho_old + aa) ) 
   w_a = dexp(-t_a*t_a*rho_old*dist)
   accu = 0.d0
   ! evaluation of the polynom Ix(t_a) * Iy(t_a) * Iz(t_a)
   do m = 0, n_pt_out,2
    accu += d1(m) * (t_a)**(dble(m)) 
   enddo
   ! equation A8 of PRA-70-062505 (2004) of Toul. Col. Sav. 
   gauss_int = gauss_int + c_a * prefactor * (1.d0 - t_a*t_a)**(1.5d0) * w_a * accu
  enddo

  general_primitive_integral_gauss = gauss_int
end

subroutine compute_ao_integrals_gauss_jl(j,l,n_integrals,buffer_i,buffer_value)
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
  double precision               :: cpu_1,cpu_2, wall_1, wall_2
  double precision               :: integral, wall_0
  double precision               :: thr,ao_tc_sym_two_e_pot
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
!      if (ao_two_e_integral_zero(i,j,k,l)) then
      if (.False.) then
        cycle
      endif
      if (ao_two_e_integral_erf_schwartz(i,k)*ao_two_e_integral_erf_schwartz(j,l) < thr ) then
        cycle
      endif
      !DIR$ FORCEINLINE
      integral = ao_tc_sym_two_e_pot(i,k,j,l)  ! i,k : r1    j,l : r2
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
