
! ---

BEGIN_PROVIDER [ double precision, ao_integrals_n_e, (ao_num,ao_num)]

  BEGIN_DOC
  !  Nucleus-electron interaction, in the |AO| basis set.
  !
  !  :math:`\langle \chi_i | -\sum_A \frac{1}{|r-R_A|} | \chi_j \rangle`
  !
  !  These integrals also contain the pseudopotential integrals.
  END_DOC

  implicit none
  integer          :: num_A, num_B, power_A(3), power_B(3)
  integer          :: i, j, k, l, n_pt_in, m
  double precision :: alpha, beta
  double precision :: A_center(3),B_center(3),C_center(3)
  double precision :: overlap_x,overlap_y,overlap_z,overlap,dx,NAI_pol_mult

  ao_integrals_n_e = 0.d0

  if (read_ao_integrals_n_e) then

    call ezfio_get_ao_one_e_ints_ao_integrals_n_e(ao_integrals_n_e)
    print *,  'AO N-e integrals read from disk'

  else

    if(use_cosgtos) then
      !print *, " use_cosgtos for ao_integrals_n_e ?", use_cosgtos

      do j = 1, ao_num
        do i = 1, ao_num
          ao_integrals_n_e(i,j) = ao_integrals_n_e_cosgtos(i,j)
        enddo
      enddo

    else

      !$OMP PARALLEL                                                   &
          !$OMP DEFAULT (NONE)                                         &
          !$OMP PRIVATE (i,j,k,l,m,alpha,beta,A_center,B_center,C_center,power_A,power_B,&
          !$OMP          num_A,num_B,Z,c,c1,n_pt_in)                      &
          !$OMP SHARED (ao_num,ao_prim_num,ao_expo_ordered_transp,ao_power,ao_nucl,nucl_coord,ao_coef_normalized_ordered_transp,&
          !$OMP         n_pt_max_integrals,ao_integrals_n_e,nucl_num,nucl_charge)

      n_pt_in = n_pt_max_integrals

      !$OMP DO SCHEDULE (dynamic)

      do j = 1, ao_num
        num_A = ao_nucl(j)
        power_A(1:3)= ao_power(j,1:3)
        A_center(1:3) = nucl_coord(num_A,1:3)

        do i = 1, ao_num

          num_B = ao_nucl(i)
          power_B(1:3)= ao_power(i,1:3)
          B_center(1:3) = nucl_coord(num_B,1:3)

          do l=1,ao_prim_num(j)
            alpha = ao_expo_ordered_transp(l,j)

            do m=1,ao_prim_num(i)
              beta = ao_expo_ordered_transp(m,i)

              double precision               :: c, c1
              c = 0.d0

              do  k = 1, nucl_num
                double precision               :: Z
                Z = nucl_charge(k)

                C_center(1:3) = nucl_coord(k,1:3)

                !print *, ' '
                !print *, A_center, B_center, C_center, power_A, power_B
                !print *, alpha, beta

                c1 = NAI_pol_mult( A_center, B_center, power_A, power_B &
                                 , alpha, beta, C_center, n_pt_in )

                !print *, ' c1 = ', c1

                c = c - Z * c1

              enddo
              ao_integrals_n_e(i,j) = ao_integrals_n_e(i,j)  &
                  + ao_coef_normalized_ordered_transp(l,j)             &
                  * ao_coef_normalized_ordered_transp(m,i) * c
            enddo
          enddo
        enddo
      enddo

    !$OMP END DO
    !$OMP END PARALLEL

    endif


    IF(do_pseudo) THEN
       ao_integrals_n_e += ao_pseudo_integrals
    ENDIF
    IF(point_charges) THEN
       ao_integrals_n_e += ao_integrals_pt_chrg
    ENDIF 

  endif


  if (write_ao_integrals_n_e) then
    call ezfio_set_ao_one_e_ints_ao_integrals_n_e(ao_integrals_n_e)
    print *,  'AO N-e integrals written to disk'
  endif

END_PROVIDER

BEGIN_PROVIDER [ double precision, ao_integrals_n_e_imag, (ao_num,ao_num)]
  BEGIN_DOC
  !  Nucleus-electron interaction, in the |AO| basis set.
  !
  !  :math:`\langle \chi_i | -\sum_A \frac{1}{|r-R_A|} | \chi_j \rangle`
  END_DOC
  implicit none
  double precision               :: alpha, beta
  integer                        :: num_A,num_B
  double precision               :: A_center(3),B_center(3),C_center(3)
  integer                        :: power_A(3),power_B(3)
  integer                        :: i,j,k,l,n_pt_in,m
  double precision               :: overlap_x,overlap_y,overlap_z,overlap,dx,NAI_pol_mult

  if (read_ao_integrals_n_e) then
    call ezfio_get_ao_one_e_ints_ao_integrals_n_e_imag(ao_integrals_n_e_imag)
    print *,  'AO N-e integrals read from disk'
  else
   print *,  irp_here, ': Not yet implemented'
  endif
END_PROVIDER


BEGIN_PROVIDER [ double precision, ao_integrals_n_e_per_atom, (ao_num,ao_num,nucl_num)]
  BEGIN_DOC
! Nucleus-electron interaction in the |AO| basis set, per atom A.
!
! :math:`\langle \chi_i | -\frac{1}{|r-R_A|} | \chi_j \rangle`
  END_DOC
  implicit none
  double precision               :: alpha, beta
  integer                        :: i_c,num_A,num_B
  double precision               :: A_center(3),B_center(3),C_center(3)
  integer                        :: power_A(3),power_B(3)
  integer                        :: i,j,k,l,n_pt_in,m
  double precision               :: overlap_x,overlap_y,overlap_z,overlap,dx,NAI_pol_mult

  ao_integrals_n_e_per_atom = 0.d0

  !$OMP PARALLEL                                                    &
      !$OMP DEFAULT (NONE)                                          &
      !$OMP PRIVATE (i,j,k,l,m,alpha,beta,A_center,B_center,power_A,power_B,&
      !$OMP  num_A,num_B,c,n_pt_in,C_center)                        &
      !$OMP SHARED (ao_num,ao_prim_num,ao_expo_ordered_transp,ao_power,ao_nucl,nucl_coord,ao_coef_normalized_ordered_transp,&
      !$OMP  n_pt_max_integrals,ao_integrals_n_e_per_atom,nucl_num)
  n_pt_in = n_pt_max_integrals
  !$OMP DO SCHEDULE (dynamic)

  double precision               :: c
  do j = 1, ao_num
    power_A(1)= ao_power(j,1)
    power_A(2)= ao_power(j,2)
    power_A(3)= ao_power(j,3)
    num_A = ao_nucl(j)
    A_center(1) = nucl_coord(num_A,1)
    A_center(2) = nucl_coord(num_A,2)
    A_center(3) = nucl_coord(num_A,3)
    do  k = 1, nucl_num
      C_center(1) = nucl_coord(k,1)
      C_center(2) = nucl_coord(k,2)
      C_center(3) = nucl_coord(k,3)
      do i = 1, ao_num
        power_B(1)= ao_power(i,1)
        power_B(2)= ao_power(i,2)
        power_B(3)= ao_power(i,3)
        num_B = ao_nucl(i)
        B_center(1) = nucl_coord(num_B,1)
        B_center(2) = nucl_coord(num_B,2)
        B_center(3) = nucl_coord(num_B,3)
        c = 0.d0
        do l=1,ao_prim_num(j)
          alpha = ao_expo_ordered_transp(l,j)
          do m=1,ao_prim_num(i)
            beta = ao_expo_ordered_transp(m,i)
            c = c + NAI_pol_mult(A_center,B_center,power_A,power_B,  &
                alpha,beta,C_center,n_pt_in)                         &
                * ao_coef_normalized_ordered_transp(l,j)             &
                * ao_coef_normalized_ordered_transp(m,i)
          enddo
        enddo
        ao_integrals_n_e_per_atom(i,j,k) = -c
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

END_PROVIDER



double precision function NAI_pol_mult(A_center,B_center,power_A,power_B,alpha,beta,C_center,n_pt_in)
  BEGIN_DOC
! Computes the electron-nucleus attraction with two primitves.
!
! :math:`\langle g_i | \frac{1}{|r-R_c|} | g_j \rangle`
  END_DOC

  implicit none
  integer, intent(in)            :: n_pt_in
  double precision,intent(in)    :: C_center(3),A_center(3),B_center(3),alpha,beta
  integer                        :: power_A(3),power_B(3)
  integer                        :: i,j,k,l,n_pt
  double precision               :: P_center(3)
  double precision               :: d(0:n_pt_in),pouet,coeff,rho,dist,const,pouet_2,p,p_inv,factor
  double precision               :: I_n_special_exact,integrate_bourrin,I_n_bibi
  double precision               :: V_n_e,const_factor,dist_integral,tmp
  double precision               :: accu,epsilo,rint
  integer                        :: n_pt_out,lmax
  include 'utils/constants.include.F'
  if ( (A_center(1)/=B_center(1)).or.                                &
        (A_center(2)/=B_center(2)).or.                               &
        (A_center(3)/=B_center(3)).or.                               &
        (A_center(1)/=C_center(1)).or.                               &
        (A_center(2)/=C_center(2)).or.                               &
        (A_center(3)/=C_center(3))) then
    continue
  else
    NAI_pol_mult = V_n_e(power_A(1),power_A(2),power_A(3),           &
        power_B(1),power_B(2),power_B(3),alpha,beta)
    return
  endif
  p = alpha + beta
  p_inv = 1.d0/p
  rho = alpha * beta * p_inv
  dist = 0.d0
  dist_integral = 0.d0
  do i = 1, 3
    P_center(i) = (alpha * A_center(i) + beta * B_center(i)) * p_inv
    dist += (A_center(i) - B_center(i))*(A_center(i) - B_center(i))
    dist_integral += (P_center(i) - C_center(i))*(P_center(i) - C_center(i))
  enddo
  const_factor = dist*rho
  const = p * dist_integral
  if(const_factor > 80.d0)then
    NAI_pol_mult = 0.d0
    return
  endif
  factor = dexp(-const_factor)
  coeff = dtwo_pi * factor * p_inv
  lmax = 20

  !  print*, "b"
  do i = 0, n_pt_in
    d(i) = 0.d0
  enddo
  n_pt =  2 * ( (power_A(1) + power_B(1)) +(power_A(2) + power_B(2)) +(power_A(3) + power_B(3)) )
  if (n_pt == 0) then
    epsilo = 1.d0
    pouet = rint(0,const)
    NAI_pol_mult = coeff * pouet
    return
  endif

  call give_polynomial_mult_center_one_e(A_center,B_center,alpha,beta,power_A,power_B,C_center,n_pt_in,d,n_pt_out)


  if(n_pt_out<0)then
    NAI_pol_mult = 0.d0
    return
  endif
  accu = 0.d0

  ! 1/r1 standard attraction integral
  epsilo = 1.d0
  ! sum of integrals of type : int {t,[0,1]}  exp-(rho.(P-Q)^2 * t^2) * t^i
  do i =0 ,n_pt_out,2
    accu +=  d(i) * rint(i/2,const)
  enddo
  NAI_pol_mult = accu * coeff

end

! ---

subroutine give_polynomial_mult_center_one_e(A_center,B_center,alpha,beta,power_A,power_B,C_center,n_pt_in,d,n_pt_out)
  implicit none
  BEGIN_DOC
  ! Returns the explicit polynomial in terms of the "t" variable of the following
  !
  ! $I_{x1}(a_x, d_x,p,q) \times I_{x1}(a_y, d_y,p,q) \times I_{x1}(a_z, d_z,p,q)$.
  END_DOC
  integer, intent(in)            :: n_pt_in
  integer,intent(out)            :: n_pt_out
  double precision, intent(in)   :: A_center(3), B_center(3),C_center(3)
  double precision, intent(in)   :: alpha,beta
  integer, intent(in)            :: power_A(3), power_B(3)
  integer                        :: a_x,b_x,a_y,b_y,a_z,b_z
  double precision               :: d(0:n_pt_in)
  double precision               :: d1(0:n_pt_in)
  double precision               :: d2(0:n_pt_in)
  double precision               :: d3(0:n_pt_in)
  double precision               :: accu,  pq_inv, p10_1, p10_2, p01_1, p01_2
  double precision               :: p,P_center(3),rho,p_inv,p_inv_2

  accu = 0.d0

  ASSERT (n_pt_in > 1)
  p = alpha+beta
  p_inv = 1.d0/p
  p_inv_2 = 0.5d0/p
  do i =1, 3
    P_center(i) = (alpha * A_center(i) + beta * B_center(i)) * p_inv
  enddo

  double precision               :: R1x(0:2), B01(0:2), R1xp(0:2),R2x(0:2)
  R1x(0)  = (P_center(1) - A_center(1))
  R1x(1)  = 0.d0
  R1x(2)  = -(P_center(1) - C_center(1))

  R1xp(0)  = (P_center(1) - B_center(1))
  R1xp(1)  = 0.d0
  R1xp(2)  =-(P_center(1) - C_center(1))

  R2x(0)  =  p_inv_2
  R2x(1)  = 0.d0
  R2x(2)  = -p_inv_2

  do i = 0,n_pt_in
    d(i) = 0.d0
  enddo
  do i = 0,n_pt_in
    d1(i) = 0.d0
  enddo
  do i = 0,n_pt_in
    d2(i) = 0.d0
  enddo
  do i = 0,n_pt_in
    d3(i) = 0.d0
  enddo
  integer                        :: n_pt1,n_pt2,n_pt3,dim,i
  n_pt1 = n_pt_in
  n_pt2 = n_pt_in
  n_pt3 = n_pt_in
  a_x = power_A(1)
  b_x = power_B(1)
  call I_x1_pol_mult_one_e(a_x,b_x,R1x,R1xp,R2x,d1,n_pt1,n_pt_in)

  if(n_pt1<0)then
    n_pt_out = -1
    do i = 0,n_pt_in
      d(i) = 0.d0
    enddo
    return
  endif

  R1x(0)  = (P_center(2) - A_center(2))
  R1x(1)  = 0.d0
  R1x(2)  = -(P_center(2) - C_center(2))

  R1xp(0)  = (P_center(2) - B_center(2))
  R1xp(1)  = 0.d0
  R1xp(2)  =-(P_center(2) - C_center(2))

  a_y = power_A(2)
  b_y = power_B(2)
  call I_x1_pol_mult_one_e(a_y,b_y,R1x,R1xp,R2x,d2,n_pt2,n_pt_in)

  if(n_pt2<0)then
    n_pt_out = -1
    do i = 0,n_pt_in
      d(i) = 0.d0
    enddo
    return
  endif


  R1x(0)  = (P_center(3) - A_center(3))
  R1x(1)  = 0.d0
  R1x(2)  = -(P_center(3) - C_center(3))

  R1xp(0)  = (P_center(3) - B_center(3))
  R1xp(1)  = 0.d0
  R1xp(2)  =-(P_center(3) - C_center(3))

  a_z = power_A(3)
  b_z = power_B(3)

  call I_x1_pol_mult_one_e(a_z,b_z,R1x,R1xp,R2x,d3,n_pt3,n_pt_in)

  if(n_pt3<0)then
    n_pt_out = -1
    do i = 0,n_pt_in
      d(i) = 0.d0
    enddo
    return
  endif
  integer                        :: n_pt_tmp
  n_pt_tmp = 0
  call multiply_poly(d1,n_pt1,d2,n_pt2,d,n_pt_tmp)
  do i = 0,n_pt_tmp
    d1(i) = 0.d0
  enddo
  n_pt_out = 0
  call multiply_poly(d ,n_pt_tmp ,d3,n_pt3,d1,n_pt_out)
  do i = 0, n_pt_out
    d(i) = d1(i)
  enddo

end


recursive subroutine I_x1_pol_mult_one_e(a,c,R1x,R1xp,R2x,d,nd,n_pt_in)
  implicit none
  BEGIN_DOC
!  Recursive routine involved in the electron-nucleus potential
  END_DOC
  integer , intent(in)           :: n_pt_in
  double precision,intent(inout) :: d(0:n_pt_in)
  integer,intent(inout)          :: nd
  integer, intent(in)            :: a,c
  double precision, intent(in)   :: R1x(0:2),R1xp(0:2),R2x(0:2)
  include 'utils/constants.include.F'
  double precision               :: X(0:max_dim)
  double precision               :: Y(0:max_dim)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: X, Y
  integer                        :: nx, ix,dim,iy,ny
  dim = n_pt_in
  ! print*,'a,c = ',a,c
  ! print*,'nd_in = ',nd

  if( (a==0) .and. (c==0))then
    nd = 0
    d(0) = 1.d0
    return
  elseif( (c<0).or.(nd<0) )then
    nd = -1
    return
  else if ((a==0).and.(c.ne.0)) then
    call I_x2_pol_mult_one_e(c,R1x,R1xp,R2x,d,nd,n_pt_in)
  else if (a==1) then
    nx = nd
    do ix=0,n_pt_in
      X(ix) = 0.d0
      Y(ix) = 0.d0
    enddo
    call I_x2_pol_mult_one_e(c-1,R1x,R1xp,R2x,X,nx,n_pt_in)
    do ix=0,nx
      X(ix) *= dble(c)
    enddo
    call multiply_poly(X,nx,R2x,2,d,nd)
    ny=0
    call I_x2_pol_mult_one_e(c,R1x,R1xp,R2x,Y,ny,n_pt_in)
    call multiply_poly(Y,ny,R1x,2,d,nd)
  else
    do ix=0,n_pt_in
      X(ix) = 0.d0
      Y(ix) = 0.d0
    enddo
    nx = 0
    call I_x1_pol_mult_one_e(a-2,c,R1x,R1xp,R2x,X,nx,n_pt_in)
    do ix=0,nx
      X(ix) *= dble(a-1)
    enddo
    call multiply_poly(X,nx,R2x,2,d,nd)

    nx = nd
    do ix=0,n_pt_in
      X(ix) = 0.d0
    enddo
    call I_x1_pol_mult_one_e(a-1,c-1,R1x,R1xp,R2x,X,nx,n_pt_in)
    do ix=0,nx
      X(ix) *= dble(c)
    enddo
    call multiply_poly(X,nx,R2x,2,d,nd)
    ny=0
    call I_x1_pol_mult_one_e(a-1,c,R1x,R1xp,R2x,Y,ny,n_pt_in)
    call multiply_poly(Y,ny,R1x,2,d,nd)
  endif
end

recursive subroutine I_x2_pol_mult_one_e(c,R1x,R1xp,R2x,d,nd,dim)
  implicit none
  BEGIN_DOC
!  Recursive routine involved in the electron-nucleus potential
  END_DOC
  integer , intent(in)           :: dim
  include 'utils/constants.include.F'
  double precision               :: d(0:max_dim)
  integer,intent(inout)          :: nd
  integer, intent(in)            :: c
  double precision, intent(in)   :: R1x(0:2),R1xp(0:2),R2x(0:2)
  integer                        :: i

  if(c==0) then
    nd = 0
    d(0) = 1.d0
    return
  elseif ((nd<0).or.(c<0))then
    nd = -1
    return
  else
    integer                        :: nx, ix,ny
    double precision               :: X(0:max_dim),Y(0:max_dim)
    !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: X, Y
    do ix=0,dim
      X(ix) = 0.d0
      Y(ix) = 0.d0
    enddo
    nx = 0
    call I_x1_pol_mult_one_e(0,c-2,R1x,R1xp,R2x,X,nx,dim)
    do ix=0,nx
      X(ix) *= dble(c-1)
    enddo
    call multiply_poly(X,nx,R2x,2,d,nd)
    ny = 0
    do ix=0,dim
      Y(ix) = 0.d0
    enddo

    call I_x1_pol_mult_one_e(0,c-1,R1x,R1xp,R2x,Y,ny,dim)
    if(ny.ge.0)then
      call multiply_poly(Y,ny,R1xp,2,d,nd)
    endif
  endif
end

double precision function V_n_e(a_x,a_y,a_z,b_x,b_y,b_z,alpha,beta)
  implicit none
  BEGIN_DOC
! Primitve nuclear attraction between the two primitves centered on the same atom.
!
! $p_1 = x^{a_x} y^{a_y} z^{a_z} \exp(-\alpha r^2)$
!
! $p_2 = x^{b_x} y^{b_y} z^{b_z} \exp(-\beta  r^2)$
  END_DOC
  integer                        :: a_x,a_y,a_z,b_x,b_y,b_z
  double precision               :: alpha,beta
  double precision               :: V_r, V_phi, V_theta
  if(iand((a_x+b_x),1)==1.or.iand(a_y+b_y,1)==1.or.iand((a_z+b_z),1)==1)then
    V_n_e = 0.d0
  else
    V_n_e =   V_r(a_x+b_x+a_y+b_y+a_z+b_z+1,alpha+beta)              &
        * V_phi(a_x+b_x,a_y+b_y)                                     &
        * V_theta(a_z+b_z,a_x+b_x+a_y+b_y+1)
  endif

end


double precision function int_gaus_pol(alpha,n)
  implicit none
  BEGIN_DOC
! Computes the integral:
!
! $\int_{-\infty}^{\infty} x^n \exp(-\alpha x^2) dx$.
  END_DOC
  double precision               :: alpha
  integer                        :: n
  double precision               :: dble_fact
  include 'utils/constants.include.F'

  int_gaus_pol = 0.d0
  if(iand(n,1).eq.0)then
    int_gaus_pol = dsqrt(alpha/pi)
    double precision               :: two_alpha
    two_alpha = alpha+alpha
    integer                        :: i
    do i=1,n,2
      int_gaus_pol = int_gaus_pol * two_alpha
    enddo
    int_gaus_pol = dble_fact(n -1) / int_gaus_pol
  endif

end

double precision function V_r(n,alpha)
  implicit none
  BEGIN_DOC
  ! Computes the radial part of the nuclear attraction integral:
  !
  ! $\int_{0}^{\infty} r^n  \exp(-\alpha  r^2)  dr$
  !
  END_DOC
  double precision               :: alpha, fact
  integer                        :: n
  include 'utils/constants.include.F'
  if(iand(n,1).eq.1)then
    V_r = 0.5d0 * fact(shiftr(n,1)) / (alpha ** (shiftr(n,1) + 1))
  else
    V_r = sqpi * fact(n) / fact(shiftr(n,1)) * (0.5d0/sqrt(alpha)) ** (n+1)
  endif
end


