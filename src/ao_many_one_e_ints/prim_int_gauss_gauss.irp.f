
double precision function overlap_gauss_r12(D_center,delta,A_center,B_center,power_A,power_B,alpha,beta)
  BEGIN_DOC
  ! Computes the following integral :
  !
  ! .. math::
  ! 
  !   \int dr exp(-delta (r - D)^2 ) (x-A_x)^a (x-B_x)^b \exp(-\alpha (x-A_x)^2 - \beta (x-B_x)^2 )
  !
  END_DOC

 implicit none
  include 'constants.include.F'
 double precision, intent(in)    :: D_center(3), delta  ! pure gaussian "D" 
 double precision, intent(in)    :: A_center(3),B_center(3),alpha,beta ! gaussian/polynoms "A" and "B"
 integer, intent(in)             :: power_A(3),power_B(3)

 double precision  :: overlap_x,overlap_y,overlap_z,overlap
 ! First you multiply the usual gaussian "A" with the gaussian exp(-delta (r - D)^2 )
 double precision  :: A_new(0:max_dim,3)! new polynom 
 double precision  :: A_center_new(3)   ! new center
 integer           :: iorder_a_new(3)   ! i_order(i) = order of the new polynom ==> should be equal to power_A
 double precision  :: alpha_new         ! new exponent
 double precision  :: fact_a_new        ! constant factor
 double precision  :: accu,coefx,coefy,coefz,coefxy,coefxyz,thr
 integer           :: d(3),i,lx,ly,lz,iorder_tmp(3),dim1
 dim1=100
 thr = 1.d-10
 d = 0 ! order of the polynom for the gaussian exp(-delta (r - D)^2 )  == 0

 ! New gaussian/polynom defined by :: new pol new center new expo   cst fact new order                                
 call give_explicit_poly_and_gaussian(A_new , A_center_new , alpha_new, fact_a_new , iorder_a_new , & 
                                      delta,alpha,d,power_A,D_center,A_center,n_pt_max_integrals)
 ! The new gaussian exp(-delta (r - D)^2 ) (x-A_x)^a \exp(-\alpha (x-A_x)^2
 accu = 0.d0
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
    call overlap_gaussian_xyz(A_center_new,B_center,alpha_new,beta,iorder_tmp,power_B,overlap_x,overlap_y,overlap_z,overlap,dim1)
    accu += coefxyz * overlap
   enddo
  enddo
 enddo
 overlap_gauss_r12 = fact_a_new * accu 
end


subroutine overlap_gauss_xyz_r12(D_center,delta,A_center,B_center,power_A,power_B,alpha,beta,gauss_ints)
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
