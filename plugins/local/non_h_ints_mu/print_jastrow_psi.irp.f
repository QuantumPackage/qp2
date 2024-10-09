program print_j_psi
 implicit none
 integer :: i,j,a,b
 do i = 1, elec_beta_num ! r2
  do j = 1, elec_alpha_num ! r1
   do a = elec_beta_num+1, mo_num ! r2
    do b = elec_alpha_num+1, mo_num ! r1 
     print*,b,a,j,i
     print*,c_ij_ab_jastrow(b,a,j,i)
    enddo
   enddo
  enddo
 enddo

end
