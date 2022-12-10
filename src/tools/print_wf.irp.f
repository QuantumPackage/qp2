program print_wf
 implicit none
 BEGIN_DOC
 ! Print the ground state wave function stored in the |EZFIO| directory
 ! in the intermediate normalization.
 !
 ! It also prints a lot of information regarding the excitation
 ! operators from the reference determinant ! and a first-order
 ! perturbative analysis of the wave function.
 !
 ! If the wave function strongly deviates from the first-order analysis,
 ! something funny is going on :)
 END_DOC


 ! this has to be done in order to be sure that N_det, psi_det and
 ! psi_coef_sorted are the wave function stored in the |EZFIO| directory.
 read_wf = .True.
 touch read_wf
 call write_wf
 call routine
end

subroutine routine
 implicit none
 integer :: i
 integer :: degree
 double precision :: hij,hii,coef_1,h00
 integer          :: exc(0:2,2,2)
 double precision :: phase
 integer :: h1,p1,h2,p2,s1,s2
 double precision :: get_two_e_integral
 double precision :: norm_mono_a,norm_mono_b
 double precision :: norm_mono_a_2,norm_mono_b_2
 double precision :: norm_mono_a_pert_2,norm_mono_b_pert_2
 double precision :: norm_mono_a_pert,norm_mono_b_pert,norm_double_1
 double precision :: delta_e,coef_2_2

 norm_mono_a = 0.d0
 norm_mono_b = 0.d0
 norm_mono_a_2 = 0.d0
 norm_mono_b_2 = 0.d0
 norm_mono_a_pert = 0.d0
 norm_mono_b_pert = 0.d0
 norm_mono_a_pert_2 = 0.d0
 norm_mono_b_pert_2 = 0.d0
 norm_double_1 = 0.d0
 do i = 1, min(N_det_print_wf,N_det)
  print*,''
  print*,'i = ',i
  call debug_det(psi_det_sorted(1,1,i),N_int)
  call get_excitation_degree(psi_det_sorted(1,1,i),psi_det_sorted(1,1,1),degree,N_int)
  print*,'degree = ',degree
  if(degree == 0)then
   print*,'Reference determinant '
   call i_H_j(psi_det_sorted(1,1,i),psi_det_sorted(1,1,i),N_int,h00)
  else if(degree .le. 2)then
   call i_H_j(psi_det_sorted(1,1,i),psi_det_sorted(1,1,i),N_int,hii)
   call i_H_j(psi_det_sorted(1,1,1),psi_det_sorted(1,1,i),N_int,hij)
   delta_e = hii - h00
   coef_1 = hij/(h00-hii)
   if(hij.ne.0.d0)then
    if (delta_e > 0.d0) then
      coef_2_2 = 0.5d0 * (delta_e - dsqrt(delta_e * delta_e + 4.d0 * hij * hij ))/ hij
    else
      coef_2_2 = 0.5d0 * (delta_e + dsqrt(delta_e * delta_e + 4.d0 * hij * hij )) /hij
    endif
   else
    coef_2_2 = 0.d0
   endif
   call get_excitation(psi_det_sorted(1,1,1),psi_det_sorted(1,1,i),exc,degree,phase,N_int)
   call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
   print*,'phase = ',phase
   if(degree == 1)then
    print*,'s1',s1
    print*,'h1,p1 = ',h1,p1
    if(s1 == 1)then
     norm_mono_a += dabs(psi_coef_sorted(i,1)/psi_coef_sorted(1,1))
     norm_mono_a_2 += dabs(psi_coef_sorted(i,1)/psi_coef_sorted(1,1))**2
     norm_mono_a_pert += dabs(coef_1)
     norm_mono_a_pert_2 += dabs(coef_1)**2
    else
     norm_mono_b += dabs(psi_coef_sorted(i,1)/psi_coef_sorted(1,1))
     norm_mono_b_2 += dabs(psi_coef_sorted(i,1)/psi_coef_sorted(1,1))**2
     norm_mono_b_pert += dabs(coef_1)
     norm_mono_b_pert_2 += dabs(coef_1)**2
    endif
    double precision :: hmono,hdouble
    call  i_H_j_verbose(psi_det_sorted(1,1,1),psi_det_sorted(1,1,i),N_int,hij,hmono,hdouble,phase)
    print*,'hmono         = ',hmono
    print*,'hdouble       = ',hdouble
    print*,'hmono+hdouble = ',hmono+hdouble
    print*,'hij           = ',hij
   else if(degree ==2)then
    print*,'s1',s1
    print*,'h1,p1 = ',h1,p1
    print*,'s2',s2
    print*,'h2,p2 = ',h2,p2
    norm_double_1 += dabs(psi_coef_sorted(i,1)/psi_coef_sorted(1,1))
   endif

   print*,'<Ref| H |D_I> = ',hij
   print*,'Delta E       = ',h00-hii
   print*,'coef pert (1) = ',coef_1
   print*,'coef 2x2      = ',coef_2_2
   print*,'Delta E_corr  = ',psi_coef_sorted(i,1)/psi_coef_sorted(1,1) * hij
  endif
   print*,'amplitude     = ',psi_coef_sorted(i,1)/psi_coef_sorted(1,1)

 enddo


 print*,''
 print*,'L1 norm of mono alpha = ',norm_mono_a
 print*,'L1 norm of mono beta  = ',norm_mono_b
 print*,'L1 norm of double exc = ',norm_double_1
 print*, '---'
 print*,'L2 norm of mono alpha = ',norm_mono_a_2
 print*,'L2 norm of mono beta  = ',norm_mono_b_2
 print*, '-- perturbative mono'
 print*,''
 print*,'L1 norm of pert alpha = ',norm_mono_a_pert
 print*,'L1 norm of pert beta  = ',norm_mono_b_pert
 print*,'L2 norm of pert alpha = ',norm_mono_a_pert_2
 print*,'L2 norm of pert beta  = ',norm_mono_b_pert_2

end

subroutine write_wf
 implicit none
 character*(128) :: output
 integer :: i_unit_output,getUnitAndOpen                                                                                 
 output=trim(ezfio_filename)//'.wf'
 i_unit_output = getUnitAndOpen(output,'w')
 integer :: i
 print*,'Writing the sorted wf'
 do i = 1, N_det
  write(i_unit_output,*)i,psi_coef_sorted(i,1)/psi_coef_sorted(1,1)
 enddo


end

