double precision function NAI_pol_mult_erf_gauss_r12(D_center,delta,A_center,B_center,power_A,power_B,alpha,beta,C_center,mu)
  BEGIN_DOC
  ! Computes the following integral R^3 :
  !
  ! .. math::
  ! 
  !   \int dr  (x-A_x)^a (x-B_x)^b \exp(-\alpha (x-A_x)^2 - \beta (x-B_x)^2 )
  !   \frac{\erf(\mu | r - R_C | )}{ | r - R_C | }$ exp(-delta (r - D)^2 ).
  !
  END_DOC

 implicit none
  include 'constants.include.F'
 double precision, intent(in)    :: D_center(3), delta  ! pure gaussian "D" 
 double precision, intent(in)    :: C_center(3),mu      ! coulomb center "C" and "mu" in the erf(mu*x)/x function
 double precision, intent(in)    :: A_center(3),B_center(3),alpha,beta ! gaussian/polynoms "A" and "B"
 integer, intent(in)             :: power_A(3),power_B(3)

 double precision  :: NAI_pol_mult_erf
 ! First you multiply the usual gaussian "A" with the gaussian exp(-delta (r - D)^2 )
 double precision  :: A_new(0:max_dim,3)! new polynom 
 double precision  :: A_center_new(3)   ! new center
 integer           :: iorder_a_new(3)   ! i_order(i) = order of the new polynom ==> should be equal to power_A
 double precision  :: alpha_new         ! new exponent
 double precision  :: fact_a_new        ! constant factor
 double precision  :: accu,coefx,coefy,coefz,coefxy,coefxyz,thr
 integer           :: d(3),i,lx,ly,lz,iorder_tmp(3)
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
    accu += coefxyz * NAI_pol_mult_erf(A_center_new,B_center,iorder_tmp,power_B,alpha_new,beta,C_center,n_pt_max_integrals,mu)
   enddo
  enddo
 enddo
 NAI_pol_mult_erf_gauss_r12 = fact_a_new * accu 
end

subroutine erfc_mu_gauss_xyz(D_center,delta,mu,A_center,B_center,power_A,power_B,alpha,beta,n_pt_in,xyz_ints)
  BEGIN_DOC
  ! Computes the following integral :
  !
  ! .. math::
  ! 
  !   \int dr exp(-delta (r - D)^2 ) x/y/z * (1 - erf(mu |r-r'|))/ |r-r'| * (x-A_x)^a (x-B_x)^b \exp(-\alpha (x-A_x)^2 - \beta (x-B_x)^2 )
  !
  !   xyz_ints(1) = x , xyz_ints(2) = y, xyz_ints(3) = z, xyz_ints(4) = x^0 
  END_DOC

 implicit none
  include 'constants.include.F'
 double precision, intent(in)    :: D_center(3), delta,mu  ! pure gaussian "D" and mu parameter
 double precision, intent(in)    :: A_center(3),B_center(3),alpha,beta ! gaussian/polynoms "A" and "B"
 integer, intent(in)             :: power_A(3),power_B(3),n_pt_in
 double precision, intent(out)   :: xyz_ints(4)

 double precision  :: NAI_pol_mult_erf
 ! First you multiply the usual gaussian "A" with the gaussian exp(-delta (r - D)^2 )
 double precision  :: A_new(0:max_dim,3)! new polynom 
 double precision  :: A_center_new(3)   ! new center
 integer           :: iorder_a_new(3)   ! i_order(i) = order of the new polynom ==> should be equal to power_A
 double precision  :: alpha_new         ! new exponent
 double precision  :: fact_a_new        ! constant factor
 double precision  :: accu,coefx,coefy,coefz,coefxy,coefxyz,thr,contrib,contrib_inf,mu_inf
 integer           :: d(3),i,lx,ly,lz,iorder_tmp(3),dim1,mm
 integer           :: power_B_tmp(3)
 dim1=100
 mu_inf = 1.d+10
 thr = 1.d-10
 d = 0 ! order of the polynom for the gaussian exp(-delta (r - D)^2 )  == 0

 ! New gaussian/polynom defined by :: new pol new center new expo   cst fact new order                                
 call give_explicit_poly_and_gaussian(A_new , A_center_new , alpha_new, fact_a_new , iorder_a_new , & 
                                      delta,alpha,d,power_A,D_center,A_center,n_pt_max_integrals)
 ! The new gaussian exp(-delta (r - D)^2 ) (x-A_x)^a \exp(-\alpha (x-A_x)^2
 xyz_ints = 0.d0
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
     power_B_tmp = power_B
     contrib = NAI_pol_mult_erf(A_center_new,B_center,iorder_tmp,power_B_tmp,alpha_new,beta,D_center,n_pt_in,mu)  
     contrib_inf = NAI_pol_mult_erf(A_center_new,B_center,iorder_tmp,power_B_tmp,alpha_new,beta,D_center,n_pt_in,mu_inf)  
     xyz_ints(4) += (contrib_inf - contrib) * coefxyz ! usual term with no x/y/z 
                                      
     do mm = 1, 3 
      ! (x phi_i ) * phi_j 
      ! x * (x - B_x)^b_x = B_x (x - B_x)^b_x + 1 * (x - B_x)^{b_x+1}
      
      !
      ! first contribution :: B_x (x - B_x)^b_x :: usual integral multiplied by B_x
      power_B_tmp = power_B
      contrib_inf = NAI_pol_mult_erf(A_center_new,B_center,iorder_tmp,power_B_tmp,alpha_new,beta,D_center,n_pt_in,mu_inf)  
      contrib = NAI_pol_mult_erf(A_center_new,B_center,iorder_tmp,power_B_tmp,alpha_new,beta,D_center,n_pt_in,mu)  
      xyz_ints(mm) += (contrib_inf - contrib) * B_center(mm) * coefxyz 
                                                       
      !
      ! second contribution :: (x - B_x)^(b_x+1) :: integral with b_x=>b_x+1 
      power_B_tmp(mm) += 1
      contrib = NAI_pol_mult_erf(A_center_new,B_center,iorder_tmp,power_B_tmp,alpha_new,beta,D_center,n_pt_in,mu)  
      contrib_inf = NAI_pol_mult_erf(A_center_new,B_center,iorder_tmp,power_B_tmp,alpha_new,beta,D_center,n_pt_in,mu_inf)  
      xyz_ints(mm) += (contrib_inf -  contrib) * coefxyz     
     enddo
   enddo
  enddo
 enddo
 xyz_ints *= fact_a_new 
end


double precision function erf_mu_gauss(D_center,delta,mu,A_center,B_center,power_A,power_B,alpha,beta,n_pt_in)
  BEGIN_DOC
  ! Computes the following integral :
  !
  ! .. math::
  ! 
  !   \int dr exp(-delta (r - D)^2 ) erf(mu*|r-r'|)/ |r-r'| * (x-A_x)^a (x-B_x)^b \exp(-\alpha (x-A_x)^2 - \beta (x-B_x)^2 )
  !
  END_DOC

 implicit none
  include 'constants.include.F'
 double precision, intent(in)    :: D_center(3), delta,mu  ! pure gaussian "D" and mu parameter
 double precision, intent(in)    :: A_center(3),B_center(3),alpha,beta ! gaussian/polynoms "A" and "B"
 integer, intent(in)             :: power_A(3),power_B(3),n_pt_in

 double precision  :: NAI_pol_mult_erf
 ! First you multiply the usual gaussian "A" with the gaussian exp(-delta (r - D)^2 )
 double precision  :: A_new(0:max_dim,3)! new polynom 
 double precision  :: A_center_new(3)   ! new center
 integer           :: iorder_a_new(3)   ! i_order(i) = order of the new polynom ==> should be equal to power_A
 double precision  :: alpha_new         ! new exponent
 double precision  :: fact_a_new        ! constant factor
 double precision  :: accu,coefx,coefy,coefz,coefxy,coefxyz,thr,contrib,contrib_inf,mu_inf
 integer           :: d(3),i,lx,ly,lz,iorder_tmp(3),dim1,mm
 dim1=100
 mu_inf = 1.d+10
 thr = 1.d-10
 d = 0 ! order of the polynom for the gaussian exp(-delta (r - D)^2 )  == 0

 ! New gaussian/polynom defined by :: new pol new center new expo   cst fact new order                                
 call give_explicit_poly_and_gaussian(A_new , A_center_new , alpha_new, fact_a_new , iorder_a_new , & 
                                      delta,alpha,d,power_A,D_center,A_center,n_pt_max_integrals)
 ! The new gaussian exp(-delta (r - D)^2 ) (x-A_x)^a \exp(-\alpha (x-A_x)^2
 erf_mu_gauss = 0.d0
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
    contrib = NAI_pol_mult_erf(A_center_new,B_center,iorder_tmp,power_B,alpha_new,beta,D_center,n_pt_in,mu)
    erf_mu_gauss += contrib * coefxyz     
   enddo
  enddo
 enddo
 erf_mu_gauss *= fact_a_new 
end

