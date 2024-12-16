subroutine rotate_nuclei(phi,theta,psi,nucl_centers,n_nucl,nucl_centers_after)
 implicit none
 BEGIN_DOC
 ! routine that rotates a set of nuclei according to three axis corresponding to angles phi, theta, psi.
 END_DOC
 double precision, intent(in) :: phi,theta,psi
 double precision, intent(in) :: nucl_centers(3,n_nucl)
 integer, intent(in) :: n_nucl
 double precision, intent(out):: nucl_centers_after(3,n_nucl)
 double precision :: r_mat(3,3)
 call r_phi_theta_psi_matrix(phi,theta,psi,r_mat)
 call get_AB_prod(r_mat,3,3,nucl_centers,n_nucl,nucl_centers_after)

end


subroutine r_phi_theta_psi_matrix(phi,theta,psi,r_mat)
 implicit none
 BEGIN_DOC
 ! routine that creates the rotation matrix corresponding to phi,theta,psi 
 !
 ! according to conventions in MDFT code
 END_DOC
 double precision, intent(in) :: phi,theta,psi
 double precision, intent(out):: r_mat(3,3)
 double precision :: ctheta, stheta 
 double precision :: cphi  , sphi 
 double precision :: cpsi  , spsi 
 ctheta = dcos(theta)
 cphi = dcos(phi)
 cpsi = dcos(psi)

 stheta = dsin(theta)
 sphi = dsin(phi)
 spsi = dsin(psi)

 r_mat(1,1) = ctheta*cphi*cpsi-sphi*spsi
 r_mat(1,2) = -ctheta*cphi*spsi-sphi*cpsi
 r_mat(1,3) = stheta*cphi
 
 r_mat(2,1) = ctheta*sphi*cpsi+cphi*spsi
 r_mat(2,2) = -ctheta*sphi*spsi+cphi*cpsi
 r_mat(2,3) = stheta*sphi

 r_mat(3,1) = -stheta*cpsi
 r_mat(3,2) = stheta*spsi
 r_mat(3,3) = ctheta

end
