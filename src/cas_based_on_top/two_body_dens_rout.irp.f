
subroutine give_n2_ii_val_ab(r1,r2,two_bod_dens)
 implicit none
 BEGIN_DOC
! contribution from purely inactive orbitals to n2_{\Psi^B}(r_1,r_2) for a CAS wave function 
 END_DOC
 double precision, intent(in) :: r1(3),r2(3)
 double precision, intent(out):: two_bod_dens
 integer :: i,j,m,n,i_m,i_n
 integer :: i_i,i_j
 double precision, allocatable :: mos_array_inact_r1(:),mos_array_inact_r2(:)
 double precision, allocatable :: mos_array_r1(:) , mos_array_r2(:) 
 ! You get all orbitals in r_1 and r_2
 allocate(mos_array_r1(mo_num) , mos_array_r2(mo_num) )
 call give_all_mos_at_r(r1,mos_array_r1) 
 call give_all_mos_at_r(r2,mos_array_r2) 
 ! You extract the inactive orbitals 
 allocate(mos_array_inact_r1(n_inact_orb) , mos_array_inact_r2(n_inact_orb) )
 do i_m = 1, n_inact_orb
  mos_array_inact_r1(i_m) = mos_array_r1(list_inact(i_m))
 enddo
 do i_m = 1, n_inact_orb
  mos_array_inact_r2(i_m) = mos_array_r2(list_inact(i_m))
 enddo

 two_bod_dens = 0.d0
 ! You browse all OCCUPIED ALPHA electrons in the \mathcal{A} space
 do m = 1, n_inact_orb ! electron 1 
 ! You browse all OCCUPIED BETA  electrons in the \mathcal{A} space
  do n = 1, n_inact_orb ! electron 2 
   ! two_bod_dens(r_1,r_2) = n_alpha(r_1) * n_beta(r_2)
   two_bod_dens += mos_array_inact_r1(m) * mos_array_inact_r1(m) * mos_array_inact_r2(n) * mos_array_inact_r2(n) 
  enddo
 enddo
end


subroutine give_n2_ia_val_ab(r1,r2,two_bod_dens,istate)
 BEGIN_DOC
! contribution from inactive and active orbitals to n2_{\Psi^B}(r_1,r_2) for the "istate" state of a CAS wave function 
 END_DOC
 implicit none
 integer, intent(in) :: istate
 double precision, intent(in) :: r1(3),r2(3)
 double precision, intent(out):: two_bod_dens
 integer :: i,orb_i,a,orb_a,n,m,b
 double precision :: rho
 double precision, allocatable :: mos_array_r1(:) , mos_array_r2(:) 
 double precision, allocatable :: mos_array_inact_r1(:),mos_array_inact_r2(:)
 double precision, allocatable :: mos_array_act_r1(:),mos_array_act_r2(:)

 two_bod_dens = 0.d0
 ! You get all orbitals in r_1 and r_2
 allocate(mos_array_r1(mo_num) , mos_array_r2(mo_num) )
 call give_all_mos_at_r(r1,mos_array_r1) 
 call give_all_mos_at_r(r2,mos_array_r2) 

 ! You extract the inactive orbitals 
 allocate( mos_array_inact_r1(n_inact_orb) , mos_array_inact_r2(n_inact_orb) )
 do i = 1, n_inact_orb
  mos_array_inact_r1(i) = mos_array_r1(list_inact(i))
 enddo
 do i= 1, n_inact_orb
  mos_array_inact_r2(i) = mos_array_r2(list_inact(i))
 enddo

 ! You extract the active orbitals 
 allocate( mos_array_act_r1(n_act_orb) , mos_array_act_r2(n_act_orb) )
 do i= 1, n_act_orb
  mos_array_act_r1(i) = mos_array_r1(list_act(i))
 enddo
 do i= 1, n_act_orb
  mos_array_act_r2(i) = mos_array_r2(list_act(i))
 enddo

 ! Contracted density : intermediate quantity
 two_bod_dens = 0.d0
 do a = 1, n_act_orb 
  do i = 1, n_inact_orb
   do b = 1, n_act_orb 
    rho = one_e_act_dm_beta_mo_for_dft(b,a,istate) + one_e_act_dm_alpha_mo_for_dft(b,a,istate) 
    two_bod_dens += mos_array_inact_r1(i) * mos_array_inact_r1(i) * mos_array_act_r2(a) * mos_array_act_r2(b) * rho
   enddo
  enddo
 enddo
end


subroutine give_n2_aa_val_ab(r1,r2,two_bod_dens,istate)
 BEGIN_DOC
! contribution from purely active orbitals to n2_{\Psi^B}(r_1,r_2) for the "istate" state of a CAS wave function 
 END_DOC
 implicit none
 integer, intent(in) :: istate
 double precision, intent(in) :: r1(3),r2(3)
 double precision, intent(out):: two_bod_dens
 integer :: i,orb_i,a,orb_a,n,m,b,c,d
 double precision :: rho
 double precision, allocatable :: mos_array_r1(:) , mos_array_r2(:) 
 double precision, allocatable :: mos_array_act_r1(:),mos_array_act_r2(:)

 two_bod_dens = 0.d0
 ! You get all orbitals in r_1 and r_2
 allocate(mos_array_r1(mo_num) , mos_array_r2(mo_num) )
 call give_all_mos_at_r(r1,mos_array_r1) 
 call give_all_mos_at_r(r2,mos_array_r2) 

 ! You extract the active orbitals 
 allocate( mos_array_act_r1(n_act_orb) , mos_array_act_r2(n_act_orb) )
 do i= 1, n_act_orb
  mos_array_act_r1(i) = mos_array_r1(list_act(i))
 enddo
 do i= 1, n_act_orb
  mos_array_act_r2(i) = mos_array_r2(list_act(i))
 enddo

 ! Contracted density : intermediate quantity
 two_bod_dens = 0.d0
 do a = 1, n_act_orb  ! 1
  do b = 1, n_act_orb  ! 2
   do c = 1, n_act_orb  ! 1
    do d = 1, n_act_orb  ! 2
     rho = mos_array_act_r1(c) * mos_array_act_r2(d) * act_2_rdm_ab_mo(d,c,b,a,istate)
     two_bod_dens += rho * mos_array_act_r1(a) * mos_array_act_r2(b) 
    enddo
   enddo
  enddo
 enddo

end

subroutine give_n2_cas(r1,r2,istate,n2_psi)
 implicit none
 BEGIN_DOC
! returns n2_psi for a general cas wave function
 END_DOC
 integer, intent(in) :: istate
 double precision, intent(in)  :: r1(3),r2(3)
 double precision, intent(out) :: n2_psi
 double precision :: two_bod_dens_ii
 double precision :: two_bod_dens_ia
 double precision :: two_bod_dens_aa
 ! inactive-inactive part of n2_psi(r1,r2)
 call give_n2_ii_val_ab(r1,r2,two_bod_dens_ii)
 ! inactive-active part of n2_psi(r1,r2)
 call give_n2_ia_val_ab(r1,r2,two_bod_dens_ia,istate)
 ! active-active part of n2_psi(r1,r2)
 call give_n2_aa_val_ab(r1,r2,two_bod_dens_aa,istate)

 n2_psi = two_bod_dens_ii + two_bod_dens_ia + two_bod_dens_aa
 
end
