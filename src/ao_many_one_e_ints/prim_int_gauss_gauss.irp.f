! ---

double precision function overlap_gauss_r12(D_center, delta, A_center, B_center, power_A, power_B, alpha, beta)

  BEGIN_DOC
  !
  ! Computes the following integral :
  !
  ! .. math                      ::
  !
  !   \int dr exp(-delta (r - D)^2 ) (x-A_x)^a (x-B_x)^b \exp(-\alpha (x-A_x)^2 - \beta (x-B_x)^2 )
  !
  END_DOC

  include 'constants.include.F'

  implicit none
  double precision, intent(in) :: D_center(3), delta  ! pure gaussian "D"
  double precision, intent(in) :: A_center(3),B_center(3),alpha,beta ! gaussian/polynoms "A" and "B"
  integer, intent(in)          :: power_A(3),power_B(3)

  double precision             :: overlap_x,overlap_y,overlap_z,overlap
  ! First you multiply the usual gaussian "A" with the gaussian exp(-delta (r - D)^2 )
  double precision             :: A_new(0:max_dim,3)! new polynom
  double precision             :: A_center_new(3)   ! new center
  integer                      :: iorder_a_new(3)   ! i_order(i) = order of the new polynom ==> should be equal to power_A
  double precision             :: alpha_new         ! new exponent
  double precision             :: fact_a_new        ! constant factor
  double precision             :: accu, coefx, coefy, coefz, coefxy, coefxyz, thr
  integer                      :: d(3), i, lx, ly, lz, iorder_tmp(3), dim1

  dim1 = 100
  thr  = 1.d-10
  d(:) = 0 ! order of the polynom for the gaussian exp(-delta (r - D)^2 )  == 0
  overlap_gauss_r12 = 0.d0

  ! New gaussian/polynom defined by :: new pol new center new expo   cst fact new order
  call give_explicit_poly_and_gaussian(A_new , A_center_new , alpha_new, fact_a_new , iorder_a_new ,&
      delta,alpha,d,power_A,D_center,A_center,n_pt_max_integrals)
  if(fact_a_new.lt.thr)return
  ! The new gaussian exp(-delta (r - D)^2 ) (x-A_x)^a \exp(-\alpha (x-A_x)^2
  accu = 0.d0
  do lx = 0, iorder_a_new(1)
    coefx = A_new(lx,1)*fact_a_new
    if(dabs(coefx).lt.thr)cycle
    iorder_tmp(1) = lx

    do ly = 0, iorder_a_new(2)
      coefy  = A_new(ly,2)
      coefxy = coefx * coefy
      if(dabs(coefxy) .lt. thr) cycle
      iorder_tmp(2) = ly

      do lz = 0, iorder_a_new(3)
        coefz   = A_new(lz,3)
        coefxyz = coefxy * coefz
        if(dabs(coefxyz) .lt. thr) cycle
        iorder_tmp(3) = lz

        call overlap_gaussian_xyz( A_center_new, B_center, alpha_new, beta, iorder_tmp, power_B &
                                 , overlap_x, overlap_y, overlap_z, overlap, dim1)

        accu += coefxyz * overlap
      enddo
    enddo
  enddo
  overlap_gauss_r12 = accu
end

!---
double precision function overlap_abs_gauss_r12(D_center,delta,A_center,B_center,power_A,power_B,alpha,beta)
  BEGIN_DOC
  ! Computes the following integral :
  !
  ! .. math                      ::
  !
  !   \int dr exp(-delta (r - D)^2 ) |(x-A_x)^a (x-B_x)^b \exp(-\alpha (x-A_x)^2 - \beta (x-B_x)^2 )|
  !
  END_DOC

  implicit none
  include 'constants.include.F'
  double precision, intent(in)   :: D_center(3), delta  ! pure gaussian "D"
  double precision, intent(in)   :: A_center(3),B_center(3),alpha,beta ! gaussian/polynoms "A" and "B"
  integer, intent(in)            :: power_A(3),power_B(3)

  double precision               :: overlap_x,overlap_y,overlap_z,overlap
  ! First you multiply the usual gaussian "A" with the gaussian exp(-delta (r - D)^2 )
  double precision               :: A_new(0:max_dim,3)! new polynom
  double precision               :: A_center_new(3)   ! new center
  integer                        :: iorder_a_new(3)   ! i_order(i) = order of the new polynom ==> should be equal to power_A
  double precision               :: alpha_new         ! new exponent
  double precision               :: fact_a_new        ! constant factor
  double precision               :: accu,coefx,coefy,coefz,coefxy,coefxyz,thr,dx,lower_exp_val
  integer                        :: d(3),i,lx,ly,lz,iorder_tmp(3),dim1
  dim1=50
  lower_exp_val = 40.d0
  thr = 1.d-12
  d(:) = 0 ! order of the polynom for the gaussian exp(-delta (r - D)^2 )  == 0
  overlap_abs_gauss_r12 = 0.d0

  ! New gaussian/polynom defined by :: new pol new center new expo   cst fact new order
  call give_explicit_poly_and_gaussian(A_new , A_center_new , alpha_new, fact_a_new , iorder_a_new ,&
      delta,alpha,d,power_A,D_center,A_center,n_pt_max_integrals)
  if(fact_a_new.lt.thr)return
  ! The new gaussian exp(-delta (r - D)^2 ) (x-A_x)^a \exp(-\alpha (x-A_x)^2
  accu = 0.d0
  do lx = 0, iorder_a_new(1)
    coefx = A_new(lx,1)*fact_a_new
!    if(dabs(coefx).lt.thr)cycle
    iorder_tmp(1) = lx
    do ly = 0, iorder_a_new(2)
      coefy = A_new(ly,2)
      coefxy = coefx * coefy
      if(dabs(coefxy).lt.thr)cycle
      iorder_tmp(2) = ly
      do lz = 0, iorder_a_new(3)
        coefz = A_new(lz,3)
        coefxyz = coefxy * coefz
        if(dabs(coefxyz).lt.thr)cycle
        iorder_tmp(3) = lz
        call overlap_x_abs(A_center_new(1),B_center(1),alpha_new,beta,iorder_tmp(1),power_B(1),overlap_x,lower_exp_val,dx,dim1)
        call overlap_x_abs(A_center_new(2),B_center(2),alpha_new,beta,iorder_tmp(2),power_B(2),overlap_y,lower_exp_val,dx,dim1)
        call overlap_x_abs(A_center_new(3),B_center(3),alpha_new,beta,iorder_tmp(3),power_B(3),overlap_z,lower_exp_val,dx,dim1)
        accu += dabs(coefxyz * overlap_x * overlap_y * overlap_z)
      enddo
    enddo
  enddo
  overlap_abs_gauss_r12= accu
end

!---

! TODO apply Gaussian product three times first
subroutine overlap_gauss_r12_v(D_center, LD_D, delta, A_center, B_center, power_A, power_B, alpha, beta, rvec, LD_rvec, n_points)

  BEGIN_DOC
  !
  ! Computes the following integral :
  !
  !   \int dr exp(-delta (r - D)^2) (x-A_x)^a (x-B_x)^b \exp(-\alpha (x-A_x)^2 - \beta (x-B_x)^2)
  !   using an array of D_centers
  !
  ! n_points: nb of integrals
  !
  END_DOC

  implicit none

  include 'constants.include.F'

  integer,          intent(in)  :: LD_D, LD_rvec, n_points
  integer,          intent(in)  :: power_A(3), power_B(3)
  double precision, intent(in)  :: D_center(LD_D,3), delta
  double precision, intent(in)  :: A_center(3), B_center(3), alpha, beta 
  double precision, intent(out) :: rvec(LD_rvec)

  integer                       :: maxab
  integer                       :: d(3), i, lx, ly, lz, iorder_tmp(3), ipoint
  double precision              :: overlap_x, overlap_y, overlap_z
  double precision              :: alpha_new
  double precision              :: accu, thr, coefxy
  integer,          allocatable :: iorder_a_new(:)
  double precision, allocatable :: overlap(:)
  double precision, allocatable :: A_new(:,:,:), A_center_new(:,:)
  double precision, allocatable :: fact_a_new(:)

  thr  = 1.d-10
  d(:) = 0

  maxab = maxval(power_A(1:3))

  allocate(A_new(n_points,0:maxab,3), A_center_new(n_points,3), fact_a_new(n_points), iorder_a_new(3), overlap(n_points))

  call give_explicit_poly_and_gaussian_v(A_new, maxab, A_center_new, alpha_new, fact_a_new, iorder_a_new, delta, alpha, d, power_A, D_center, LD_D, A_center, n_points)

  rvec(:) = 0.d0

  do lx = 0, iorder_a_new(1)
    iorder_tmp(1) = lx

    do ly = 0, iorder_a_new(2)
      iorder_tmp(2) = ly

      do lz = 0, iorder_a_new(3)
        iorder_tmp(3) = lz

        call overlap_gaussian_xyz_v(A_center_new, B_center, alpha_new, beta, iorder_tmp, power_B, overlap, n_points)

        do ipoint = 1, n_points
          rvec(ipoint) = rvec(ipoint) + A_new(ipoint,lx,1) * A_new(ipoint,ly,2) * A_new(ipoint,lz,3) * overlap(ipoint)
        enddo
      enddo
    enddo
  enddo

  do ipoint = 1, n_points
    rvec(ipoint) = rvec(ipoint) * fact_a_new(ipoint)
  enddo

  deallocate(A_new, A_center_new, fact_a_new, iorder_a_new, overlap)

end subroutine overlap_gauss_r12_v

!---

subroutine overlap_gauss_xyz_r12(D_center, delta, A_center, B_center, power_A, power_B, alpha, beta, gauss_ints)

  BEGIN_DOC
  ! Computes the following integral :
  !
  ! .. math::
  !
  !   gauss_ints(m) = \int dr exp(-delta (r - D)^2 ) * x/y/z (x-A_x)^a (x-B_x)^b \exp(-\alpha (x-A_x)^2 - \beta (x-B_x)^2 )
  !
  ! with m == 1 ==> x, m == 2 ==> y, m == 3 ==> z
  END_DOC

 implicit none
  include 'constants.include.F'
 double precision, intent(in)    :: D_center(3), delta  ! pure gaussian "D"
 double precision, intent(in)    :: A_center(3),B_center(3),alpha,beta ! gaussian/polynoms "A" and "B"
 integer, intent(in)             :: power_A(3),power_B(3)
 double precision, intent(out)   :: gauss_ints(3)

 double precision  :: overlap_x,overlap_y,overlap_z,overlap
 ! First you multiply the usual gaussian "A" with the gaussian exp(-delta (r - D)^2 )
 double precision  :: A_new(0:max_dim,3)! new polynom
 double precision  :: A_center_new(3)   ! new center
 integer           :: iorder_a_new(3)   ! i_order(i) = order of the new polynom ==> should be equal to power_A
 integer           :: power_B_new(3)
 double precision  :: alpha_new         ! new exponent
 double precision  :: fact_a_new        ! constant factor
 double precision  :: coefx,coefy,coefz,coefxy,coefxyz,thr
 integer           :: d(3),i,lx,ly,lz,iorder_tmp(3),dim1,m
 dim1=100
 thr = 1.d-10
 d = 0 ! order of the polynom for the gaussian exp(-delta (r - D)^2 )  == 0

 ! New gaussian/polynom defined by :: new pol new center new expo   cst fact new order
 call give_explicit_poly_and_gaussian(A_new , A_center_new , alpha_new, fact_a_new , iorder_a_new , &
                                      delta,alpha,d,power_A,D_center,A_center,n_pt_max_integrals)
 ! The new gaussian exp(-delta (r - D)^2 ) (x-A_x)^a \exp(-\alpha (x-A_x)^2
 gauss_ints = 0.d0
 do lx = 0, iorder_a_new(1)
  coefx = A_new(lx,1)
  if(dabs(coefx).lt.thr)cycle
  iorder_tmp(1) = lx
  do ly = 0, iorder_a_new(2)
   coefy = A_new(ly,2)
   coefxy = coefx * coefy
   if(dabs(coefxy).lt.thr)cycle
   iorder_tmp(2) = ly
   do lz = 0, iorder_a_new(3)
    coefz = A_new(lz,3)
    coefxyz = coefxy * coefz
    if(dabs(coefxyz).lt.thr)cycle
    iorder_tmp(3) = lz
    do m = 1, 3
     ! change (x-Bx)^bx --> (x-Bx)^(bx+1) + Bx(x-Bx)^bx
     power_B_new = power_B
     power_B_new(m) += 1 ! (x-Bx)^(bx+1)
     call overlap_gaussian_xyz(A_center_new,B_center,alpha_new,beta,iorder_tmp,power_B_new,overlap_x,overlap_y,overlap_z,overlap,dim1)
     gauss_ints(m) += coefxyz * overlap

     power_B_new = power_B
     call overlap_gaussian_xyz(A_center_new,B_center,alpha_new,beta,iorder_tmp,power_B_new,overlap_x,overlap_y,overlap_z,overlap,dim1)
     gauss_ints(m) += coefxyz * overlap * B_center(m) ! Bx (x-Bx)^(bx)
    enddo
   enddo
  enddo
 enddo
 gauss_ints *= fact_a_new
end

double precision function overlap_gauss_xyz_r12_specific(D_center,delta,A_center,B_center,power_A,power_B,alpha,beta,mx)
  BEGIN_DOC
  ! Computes the following integral :
  !
  ! .. math::
  !
  !    \int dr exp(-delta (r - D)^2 ) * x/y/z (x-A_x)^a (x-B_x)^b \exp(-\alpha (x-A_x)^2 - \beta (x-B_x)^2 )
  !
  ! with mx == 1 ==> x, mx == 2 ==> y, mx == 3 ==> z
  END_DOC

 implicit none
  include 'constants.include.F'
 double precision, intent(in)    :: D_center(3), delta  ! pure gaussian "D"
 double precision, intent(in)    :: A_center(3),B_center(3),alpha,beta ! gaussian/polynoms "A" and "B"
 integer, intent(in)             :: power_A(3),power_B(3),mx

 double precision  :: overlap_x,overlap_y,overlap_z,overlap
 ! First you multiply the usual gaussian "A" with the gaussian exp(-delta (r - D)^2 )
 double precision  :: A_new(0:max_dim,3)! new polynom
 double precision  :: A_center_new(3)   ! new center
 integer           :: iorder_a_new(3)   ! i_order(i) = order of the new polynom ==> should be equal to power_A
 integer           :: power_B_new(3)
 double precision  :: alpha_new         ! new exponent
 double precision  :: fact_a_new        ! constant factor
 double precision  :: coefx,coefy,coefz,coefxy,coefxyz,thr
 integer           :: d(3),i,lx,ly,lz,iorder_tmp(3),dim1,m
 dim1=100
 thr = 1.d-10
 d = 0 ! order of the polynom for the gaussian exp(-delta (r - D)^2 )  == 0

 ! New gaussian/polynom defined by :: new pol new center new expo   cst fact new order
 call give_explicit_poly_and_gaussian(A_new , A_center_new , alpha_new, fact_a_new , iorder_a_new , &
                                      delta,alpha,d,power_A,D_center,A_center,n_pt_max_integrals)
 ! The new gaussian exp(-delta (r - D)^2 ) (x-A_x)^a \exp(-\alpha (x-A_x)^2
 overlap_gauss_xyz_r12_specific = 0.d0
 do lx = 0, iorder_a_new(1)
  coefx = A_new(lx,1)
  if(dabs(coefx).lt.thr)cycle
  iorder_tmp(1) = lx
  do ly = 0, iorder_a_new(2)
   coefy = A_new(ly,2)
   coefxy = coefx * coefy
   if(dabs(coefxy).lt.thr)cycle
   iorder_tmp(2) = ly
   do lz = 0, iorder_a_new(3)
    coefz = A_new(lz,3)
    coefxyz = coefxy * coefz
    if(dabs(coefxyz).lt.thr)cycle
    iorder_tmp(3) = lz
    m = mx
    ! change (x-Bx)^bx --> (x-Bx)^(bx+1) + Bx(x-Bx)^bx
    power_B_new = power_B
    power_B_new(m) += 1 ! (x-Bx)^(bx+1)
    call overlap_gaussian_xyz(A_center_new,B_center,alpha_new,beta,iorder_tmp,power_B_new,overlap_x,overlap_y,overlap_z,overlap,dim1)
    overlap_gauss_xyz_r12_specific += coefxyz * overlap

    power_B_new = power_B
    call overlap_gaussian_xyz(A_center_new,B_center,alpha_new,beta,iorder_tmp,power_B_new,overlap_x,overlap_y,overlap_z,overlap,dim1)
    overlap_gauss_xyz_r12_specific += coefxyz * overlap * B_center(m) ! Bx (x-Bx)^(bx)
   enddo
  enddo
 enddo
 overlap_gauss_xyz_r12_specific *= fact_a_new
end
