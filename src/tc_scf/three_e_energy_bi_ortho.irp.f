
subroutine contrib_3e_diag_sss(i,j,k,integral)
 implicit none
 integer, intent(in) :: i,j,k
 BEGIN_DOC
 ! returns the pure same spin contribution to diagonal matrix element of 3e term
 END_DOC
 double precision, intent(out) :: integral
 double precision :: direct_int, exch_13_int, exch_23_int, exch_12_int, c_3_int, c_minus_3_int
 call  give_integrals_3_body_bi_ort(i, k, j, i, k, j, direct_int )!!! < i k j | i k j >
 call  give_integrals_3_body_bi_ort(i, k, j, j, i, k, c_3_int)      ! < i k j | j i k >
 call  give_integrals_3_body_bi_ort(i, k, j, k, j, i, c_minus_3_int)! < i k j | k j i >
 integral = direct_int + c_3_int + c_minus_3_int 
 ! negative terms :: exchange contrib
 call  give_integrals_3_body_bi_ort(i, k, j, j, k, i, exch_13_int)!!! < i k j | j k i > : E_13 
 call  give_integrals_3_body_bi_ort(i, k, j, i, j, k, exch_23_int)!!! < i k j | i j k > : E_23
 call  give_integrals_3_body_bi_ort(i, k, j, k, i, j, exch_12_int)!!! < i k j | k i j > : E_12
 integral += - exch_13_int - exch_23_int  - exch_12_int 
 integral = -integral
end

subroutine contrib_3e_diag_soo(i,j,k,integral)
 implicit none
 integer, intent(in) :: i,j,k
 BEGIN_DOC
 ! returns the pure same spin contribution to diagonal matrix element of 3e term
 END_DOC
 double precision, intent(out) :: integral
 double precision :: direct_int, exch_23_int
 call  give_integrals_3_body_bi_ort(i, k, j, i, k, j, direct_int) ! < i k j | i k j >
 call  give_integrals_3_body_bi_ort(i, k, j, i, j, k, exch_23_int)! < i k j | i j k > : E_23
 integral = direct_int - exch_23_int 
 integral = -integral
end


subroutine give_aaa_contrib_bis(integral_aaa)
 implicit none
 double precision, intent(out) :: integral_aaa
 double precision :: integral
 integer :: i,j,k
 integral_aaa = 0.d0
 do i = 1, elec_alpha_num
  do j = i+1, elec_alpha_num
   do k = j+1, elec_alpha_num
    call contrib_3e_diag_sss(i,j,k,integral)
    integral_aaa += integral
   enddo
  enddo
 enddo

end

subroutine give_aaa_contrib(integral_aaa)
 implicit none
 double precision, intent(out) :: integral_aaa
 double precision :: integral
 integer :: i,j,k
 integral_aaa = 0.d0
 do i = 1, elec_alpha_num
  do j = 1, elec_alpha_num
   do k = 1, elec_alpha_num
    call contrib_3e_diag_sss(i,j,k,integral)
    integral_aaa += integral
   enddo
  enddo
 enddo
 integral_aaa *= 1.d0/6.d0 
end


subroutine give_aab_contrib(integral_aab)
 implicit none
 double precision, intent(out) :: integral_aab
 double precision :: integral
 integer :: i,j,k
 integral_aab = 0.d0
 do i = 1, elec_beta_num
  do j = 1, elec_alpha_num
   do k = 1, elec_alpha_num
    call contrib_3e_diag_soo(i,j,k,integral)
    integral_aab += integral
   enddo
  enddo
 enddo
 integral_aab *= 0.5d0
end


subroutine give_aab_contrib_bis(integral_aab)
 implicit none
 double precision, intent(out) :: integral_aab
 double precision :: integral
 integer :: i,j,k
 integral_aab = 0.d0
 do i = 1, elec_beta_num
  do j = 1, elec_alpha_num
   do k = j+1, elec_alpha_num
    call contrib_3e_diag_soo(i,j,k,integral)
    integral_aab += integral
   enddo
  enddo
 enddo
end


subroutine give_abb_contrib(integral_abb)
 implicit none
 double precision, intent(out) :: integral_abb
 double precision :: integral
 integer :: i,j,k
 integral_abb = 0.d0
 do i = 1, elec_alpha_num
  do j = 1, elec_beta_num
   do k = 1, elec_beta_num
    call contrib_3e_diag_soo(i,j,k,integral)
    integral_abb += integral
   enddo
  enddo
 enddo
 integral_abb *= 0.5d0
end

subroutine give_abb_contrib_bis(integral_abb)
 implicit none
 double precision, intent(out) :: integral_abb
 double precision :: integral
 integer :: i,j,k
 integral_abb = 0.d0
 do i = 1, elec_alpha_num
  do j = 1, elec_beta_num
   do k = j+1, elec_beta_num
    call contrib_3e_diag_soo(i,j,k,integral)
    integral_abb += integral
   enddo
  enddo
 enddo
end

subroutine give_bbb_contrib_bis(integral_bbb)
 implicit none
 double precision, intent(out) :: integral_bbb
 double precision :: integral
 integer :: i,j,k
 integral_bbb = 0.d0
 do i = 1, elec_beta_num
  do j = i+1, elec_beta_num
   do k = j+1, elec_beta_num
    call contrib_3e_diag_sss(i,j,k,integral)
    integral_bbb += integral
   enddo
  enddo
 enddo

end

subroutine give_bbb_contrib(integral_bbb)
 implicit none
 double precision, intent(out) :: integral_bbb
 double precision :: integral
 integer :: i,j,k
 integral_bbb = 0.d0
 do i = 1, elec_beta_num
  do j = 1, elec_beta_num
   do k = 1, elec_beta_num
    call contrib_3e_diag_sss(i,j,k,integral)
    integral_bbb += integral
   enddo
  enddo
 enddo
 integral_bbb *= 1.d0/6.d0 
end


