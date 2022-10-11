program test_normal_order
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
  call test
end

subroutine test
 implicit none
  use bitmasks ! you need to include the bitmasks_module.f90 features
 integer :: h1,h2,p1,p2,s1,s2,i_ok,degree,Ne(2)
 integer                        :: exc(0:2,2,2)
 integer(bit_kind), allocatable :: det_i(:,:)
 double precision :: hmono,htwoe,hthree,htilde_ij,accu,phase,normal
 integer,           allocatable :: occ(:,:)
 allocate( occ(N_int*bit_kind_size,2) )
 call bitstring_to_list_ab(ref_bitmask, occ, Ne, N_int)
 allocate(det_i(N_int,2))
 s1 = 1 
 s2 = 2
 accu = 0.d0
 do h1 = 1, elec_beta_num
  do p1 = elec_beta_num+1, mo_num
   do h2 = 1, elec_beta_num
    do p2 = elec_beta_num+1, mo_num
     det_i = ref_bitmask
     call do_single_excitation(det_i,h1,p1,s1,i_ok)
     call do_single_excitation(det_i,h2,p2,s2,i_ok)
     call htilde_mu_mat_bi_ortho(det_i,HF_bitmask,N_int,hmono,htwoe,hthree,htilde_ij)
     call get_excitation_degree(ref_bitmask,det_i,degree,N_int)
     call get_excitation(ref_bitmask,det_i,exc,degree,phase,N_int)
     hthree *= phase
     normal = normal_two_body_bi_orth_ab(p2,h2,p1,h1)
     accu += dabs(hthree-normal)
    enddo
   enddo
  enddo
 enddo
 print*,'accu opposite spin = ',accu

 s1 = 2 
 s2 = 2
 accu = 0.d0
 do h1 = 1, elec_beta_num
  do p1 = elec_beta_num+1, mo_num
   do h2 = h1+1, elec_beta_num
    do p2 = elec_beta_num+1, mo_num
     det_i = ref_bitmask
     call do_single_excitation(det_i,h1,p1,s1,i_ok)
     call do_single_excitation(det_i,h2,p2,s2,i_ok)
     if(i_ok.ne.1)cycle
     call htilde_mu_mat_bi_ortho(det_i,ref_bitmask,N_int,hmono,htwoe,hthree,htilde_ij)
     call get_excitation_degree(ref_bitmask,det_i,degree,N_int)
     call get_excitation(ref_bitmask,det_i,exc,degree,phase,N_int)
     hthree *= phase
     normal = normal_two_body_bi_orth_aa_bb(p2,h2,p1,h1)
     accu += dabs(hthree-normal)
    enddo
   enddo
  enddo
 enddo
 print*,'accu same spin = ',accu
end


