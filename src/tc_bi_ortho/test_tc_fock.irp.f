program test_tc_fock
 implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  print *, 'Hello world'
  my_grid_becke = .True.
  my_n_pt_r_grid = 30
  my_n_pt_a_grid = 50
  read_wf = .True.
  touch read_wf
  touch  my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  !call routine_1
  !call routine_2
!  call routine_3()

! call test_3e
 call routine_tot
end

! ---

subroutine test_3e
 implicit none
 double precision :: integral_aaa,integral_aab,integral_abb,integral_bbb,accu
 double precision ::  hmono, htwoe, hthree, htot
 call htilde_mu_mat_bi_ortho(ref_bitmask, ref_bitmask, N_int, hmono, htwoe, hthree, htot)
! call diag_htilde_three_body_ints_bi_ort(N_int, ref_bitmask, hthree)
 print*,'hmono = ',hmono
 print*,'htwoe = ',htwoe
 print*,'hthree= ',hthree
 print*,'htot  = ',htot
 print*,''
 print*,''
 print*,'TC_one= ',tc_hf_one_e_energy
 print*,'TC_two= ',TC_HF_two_e_energy
 print*,'TC_3e = ',diag_three_elem_hf
 print*,'TC_tot= ',TC_HF_energy
 print*,''
 print*,''
 call give_aaa_contrib(integral_aaa)
 print*,'integral_aaa = ',integral_aaa
 call give_aab_contrib(integral_aab)
 print*,'integral_aab = ',integral_aab
 call give_abb_contrib(integral_abb)
 print*,'integral_abb = ',integral_abb
 call give_bbb_contrib(integral_bbb)
 print*,'integral_bbb = ',integral_bbb
 accu = integral_aaa + integral_aab + integral_abb + integral_bbb
 print*,'accu = ',accu
 print*,'delta = ',hthree - accu

end

subroutine routine_3()

  use bitmasks ! you need to include the bitmasks_module.f90 features

  implicit none
  integer                        :: i, a, i_ok, s1
  double precision               :: hmono, htwoe, hthree, htilde_ij
  double precision               :: err_ai, err_tot, ref, new
  integer(bit_kind), allocatable :: det_i(:,:)

  allocate(det_i(N_int,2))

  err_tot = 0.d0
 
  do s1 = 1, 2

   det_i = ref_bitmask
   call debug_det(det_i, N_int)
   print*, ' HF det'
   call debug_det(det_i, N_int)
 
   do i = 1, elec_num_tab(s1)
     do a = elec_num_tab(s1)+1, mo_num ! virtual 
 
 
       det_i = ref_bitmask
       call do_single_excitation(det_i, i, a, s1, i_ok)
       if(i_ok == -1) then
        print*, 'PB !!'
        print*, i, a
        stop
       endif
       print*, ' excited det'
       call debug_det(det_i, N_int)
 
       call htilde_mu_mat_bi_ortho(det_i, ref_bitmask, N_int, hmono, htwoe, hthree, htilde_ij)
       if(dabs(hthree).lt.1.d-10)cycle
       ref = hthree 
       if(s1 == 1)then
        new = fock_a_tot_3e_bi_orth(a,i)
       else if(s1 == 2)then
        new = fock_b_tot_3e_bi_orth(a,i)
       endif
       err_ai = dabs(dabs(ref) - dabs(new))
       if(err_ai .gt. 1d-7) then
         print*,'s1 = ',s1
         print*, ' warning on', i, a
         print*, ref,new,err_ai
       endif
       print*, ref,new,err_ai
       err_tot += err_ai
 
       write(22, *) htilde_ij
     enddo
   enddo
  enddo

  print *, ' err_tot = ', err_tot

  deallocate(det_i)

end subroutine routine_3

! ---
subroutine routine_tot()

  use bitmasks ! you need to include the bitmasks_module.f90 features

  implicit none
  integer                        :: i, a, i_ok, s1,other_spin(2)
  double precision               :: hmono, htwoe, hthree, htilde_ij
  double precision               :: err_ai, err_tot, ref, new
  integer(bit_kind), allocatable :: det_i(:,:)

  allocate(det_i(N_int,2))
  other_spin(1) = 2
  other_spin(2) = 1

  err_tot = 0.d0
 
!  do s1 = 1, 2
   s1 = 2
   det_i = ref_bitmask
   call debug_det(det_i, N_int)
   print*, ' HF det'
   call debug_det(det_i, N_int)
 
!   do i = 1, elec_num_tab(s1)
!     do a = elec_num_tab(s1)+1, mo_num ! virtual 
   do i = 1, elec_beta_num
     do a = elec_beta_num+1, elec_alpha_num! virtual 
!   do i = elec_beta_num+1, elec_alpha_num
!     do a = elec_alpha_num+1, mo_num! virtual 
       print*,i,a 
 
       det_i = ref_bitmask
       call do_single_excitation(det_i, i, a, s1, i_ok)
       if(i_ok == -1) then
        print*, 'PB !!'
        print*, i, a
        stop
       endif
 
       call htilde_mu_mat_bi_ortho(det_i, ref_bitmask, N_int, hmono, htwoe, hthree, htilde_ij)
       print*,htilde_ij
       if(dabs(htilde_ij).lt.1.d-10)cycle
       print*, ' excited det'
       call debug_det(det_i, N_int)

       if(s1 == 1)then
        new = Fock_matrix_tc_mo_alpha(a,i)
       else
        new = Fock_matrix_tc_mo_beta(a,i)
       endif
       ref = htilde_ij
!       if(s1 == 1)then
!        new = fock_a_tot_3e_bi_orth(a,i)
!       else if(s1 == 2)then
!        new = fock_b_tot_3e_bi_orth(a,i)
!       endif
       err_ai = dabs(dabs(ref) - dabs(new))
       if(err_ai .gt. 1d-7) then
         print*,'s1 = ',s1
         print*, ' warning on', i, a
         print*, ref,new,err_ai
       endif
       print*, ref,new,err_ai
       err_tot += err_ai
 
       write(22, *) htilde_ij
     enddo
   enddo
!  enddo

  print *, ' err_tot = ', err_tot

  deallocate(det_i)

end subroutine routine_3
