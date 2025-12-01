BEGIN_PROVIDER [double precision, Fock_matrix_tc_mo_core_eri, (mo_num, mo_num)]

  BEGIN_DOC
  ! Two-electron part with BARE COULOMB interaction of the Fock matrix 
  ! 
  ! with contributions coming only from the CORE ELECTRONS
  END_DOC

  implicit none
  double precision              :: t0, t1, tt0, tt1
  double precision, allocatable :: tmp(:,:)

  PROVIDE mo_l_coef mo_r_coef

  call ao_to_mo_bi_ortho( ao_two_e_core_integral, size(ao_two_e_core_integral, 1) &
                        ,  Fock_matrix_tc_mo_core_eri, size( Fock_matrix_tc_mo_core_eri, 1) )


END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, Fock_matrix_alpha_tc_mo_valence_eri, (mo_num, mo_num)]

  BEGIN_DOC
  ! Two-electron part with BARE COULOMB interaction of the Fock matrix 
  ! 
  ! with contributions coming only from the VALENCE ELECTRONS
  END_DOC

  implicit none
  double precision              :: t0, t1, tt0, tt1
  double precision, allocatable :: tmp(:,:)

  PROVIDE mo_l_coef mo_r_coef

  call ao_to_mo_bi_ortho( ao_two_e_valence_integral_alpha, size(ao_two_e_valence_integral_alpha, 1) &
                        ,  Fock_matrix_alpha_tc_mo_valence_eri, size(Fock_matrix_alpha_tc_mo_valence_eri, 1) )


END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, Fock_matrix_beta_tc_mo_valence_eri, (mo_num, mo_num)]

  BEGIN_DOC
  ! Two-electron part with BARE COULOMB interaction of the Fock matrix 
  ! 
  ! with contributions coming only from the VALENCE ELECTRONS
  END_DOC

  implicit none
  double precision              :: t0, t1, tt0, tt1
  double precision, allocatable :: tmp(:,:)

  PROVIDE mo_l_coef mo_r_coef

  call ao_to_mo_bi_ortho( ao_two_e_valence_integral_beta, size(ao_two_e_valence_integral_beta, 1) &
                        ,  Fock_matrix_beta_tc_mo_valence_eri, size(Fock_matrix_beta_tc_mo_valence_eri, 1) )


END_PROVIDER

! ---


BEGIN_PROVIDER [ double precision, Fock_matrix_tc_eri_mo_valence, (mo_num, mo_num) ]
 implicit none
 BEGIN_DOC
 ! Fock matrix containing only the 1/r12 interaction and only the VALENCE density matrix
 END_DOC
 integer :: i,j
  if(elec_alpha_num == elec_beta_num) then

    PROVIDE Fock_matrix_alpha_tc_mo_valence_eri

    Fock_matrix_tc_eri_mo_valence = Fock_matrix_alpha_tc_mo_valence_eri

  else

    PROVIDE Fock_matrix_beta_tc_mo_valence_eri Fock_matrix_alpha_tc_mo_valence_eri

    do j = 1, elec_beta_num
      ! F-K
      do i = 1, elec_beta_num !CC
        Fock_matrix_tc_eri_mo_valence(i,j) = 0.5d0*(Fock_matrix_alpha_tc_mo_valence_eri(i,j)+Fock_matrix_beta_tc_mo_valence_eri(i,j))&
             - (Fock_matrix_beta_tc_mo_valence_eri(i,j) - Fock_matrix_alpha_tc_mo_valence_eri(i,j))
      enddo
      ! F+K/2
      do i = elec_beta_num+1, elec_alpha_num  !CA
        Fock_matrix_tc_eri_mo_valence(i,j) = 0.5d0*(Fock_matrix_alpha_tc_mo_valence_eri(i,j)+Fock_matrix_beta_tc_mo_valence_eri(i,j))&
             + 0.5d0*(Fock_matrix_beta_tc_mo_valence_eri(i,j) - Fock_matrix_alpha_tc_mo_valence_eri(i,j))
      enddo
      ! F
      do i = elec_alpha_num+1, mo_num !CV
        Fock_matrix_tc_eri_mo_valence(i,j) = 0.5d0*(Fock_matrix_alpha_tc_mo_valence_eri(i,j)+Fock_matrix_beta_tc_mo_valence_eri(i,j))
      enddo
    enddo

    do j = elec_beta_num+1, elec_alpha_num
      ! F+K/2
      do i = 1, elec_beta_num !AC
        Fock_matrix_tc_eri_mo_valence(i,j) = 0.5d0*(Fock_matrix_alpha_tc_mo_valence_eri(i,j)+Fock_matrix_beta_tc_mo_valence_eri(i,j))&
             + 0.5d0*(Fock_matrix_beta_tc_mo_valence_eri(i,j) - Fock_matrix_alpha_tc_mo_valence_eri(i,j))
      enddo
      ! F
      do i = elec_beta_num+1, elec_alpha_num !AA
        Fock_matrix_tc_eri_mo_valence(i,j) = 0.5d0*(Fock_matrix_alpha_tc_mo_valence_eri(i,j)+Fock_matrix_beta_tc_mo_valence_eri(i,j))
      enddo
      ! F-K/2
      do i = elec_alpha_num+1, mo_num !AV
        Fock_matrix_tc_eri_mo_valence(i,j) = 0.5d0*(Fock_matrix_alpha_tc_mo_valence_eri(i,j)+Fock_matrix_beta_tc_mo_valence_eri(i,j))&
             - 0.5d0*(Fock_matrix_beta_tc_mo_valence_eri(i,j) - Fock_matrix_alpha_tc_mo_valence_eri(i,j))
      enddo
    enddo

    do j = elec_alpha_num+1, mo_num
      ! F
      do i = 1, elec_beta_num !VC
        Fock_matrix_tc_eri_mo_valence(i,j) = 0.5d0*(Fock_matrix_alpha_tc_mo_valence_eri(i,j)+Fock_matrix_beta_tc_mo_valence_eri(i,j))
      enddo
      ! F-K/2
      do i = elec_beta_num+1, elec_alpha_num !VA
        Fock_matrix_tc_eri_mo_valence(i,j) = 0.5d0*(Fock_matrix_alpha_tc_mo_valence_eri(i,j)+Fock_matrix_beta_tc_mo_valence_eri(i,j))&
             - 0.5d0*(Fock_matrix_beta_tc_mo_valence_eri(i,j) - Fock_matrix_alpha_tc_mo_valence_eri(i,j))
      enddo
      ! F+K
      do i = elec_alpha_num+1, mo_num !VV
        Fock_matrix_tc_eri_mo_valence(i,j) = 0.5d0*(Fock_matrix_alpha_tc_mo_valence_eri(i,j)+Fock_matrix_beta_tc_mo_valence_eri(i,j)) &
             + (Fock_matrix_beta_tc_mo_valence_eri(i,j) - Fock_matrix_alpha_tc_mo_valence_eri(i,j))
      enddo
    enddo
  endif
END_PROVIDER 


 BEGIN_PROVIDER [ double precision, ao_two_e_core_integral, (ao_num, ao_num) ]
 use map_module
 implicit none
 BEGIN_DOC
 ! Fock matrices in AO basis set for only the core orbitals 
 ! 
 ! WARNING : use the TCSCF_bi_ort_core_dm_ao matrix  
 END_DOC

 integer :: m,n,l,s,j
 double precision :: integral
 double precision, allocatable :: X(:), X2(:,:,:,:), X3(:,:,:,:)

 allocate (X(cholesky_ao_num))

 ! X(j) = \sum_{mn} TCSCF_bi_ort_core_dm_ao(m,n) * cholesky_ao(m,n,j)
 call dgemm('T','N',cholesky_ao_num,1,ao_num*ao_num,1.d0,            &
     cholesky_ao, ao_num*ao_num,                                     &
     TCSCF_bi_ort_core_dm_ao, ao_num*ao_num,0.d0,                      &
     X, cholesky_ao_num)
!

 ! ao_two_e_core_integral(m,n) = \sum_{j} cholesky_ao(m,n,j) * X(j)
 call dgemm('N','N',ao_num*ao_num,1,cholesky_ao_num, 1.d0,           &
     cholesky_ao, ao_num*ao_num,                                     &
     X, cholesky_ao_num, 0.d0,                                       &
     ao_two_e_core_integral, ao_num*ao_num)

 deallocate(X)

 double precision :: rss, mem0, mem
 double precision :: memory_of_double

 integer :: iblock
 integer :: block_size

 call resident_memory(mem0)

 block_size = 1024

 rss = memory_of_double(2.d0*ao_num*ao_num)
 do
   mem = mem0 + block_size*rss
   if ( (block_size < 2).or.(mem < qp_max_mem) ) exit
   block_size = block_size/2
 enddo

 call check_mem(block_size*rss, irp_here)

 allocate(X2(ao_num,ao_num,block_size,2))
 allocate(X3(ao_num,block_size,ao_num,2))

! ao_two_e_core_integral (l,s) -= cholesky_ao(l,m,j) * TCSCF_bi_ort_core_dm_ao(m,n) * cholesky_ao(n,s,j)

 do iblock=1,cholesky_ao_num,block_size

   call dgemm('N','N',ao_num,ao_num*min(cholesky_ao_num-iblock+1,block_size),ao_num, 1.d0,      &
       TCSCF_bi_ort_core_dm_ao, ao_num,    &
       cholesky_ao(1,1,iblock), ao_num, 0.d0,  &
       X2(1,1,1,1), ao_num)

     do s=1,ao_num
      do j=1,min(cholesky_ao_num-iblock+1,block_size)
       do m=1,ao_num
        X3(m,j,s,1) = X2(m,s,j,1)
       enddo
      enddo
     enddo

   call dgemm('N','N',ao_num,ao_num,ao_num*min(cholesky_ao_num-iblock+1,block_size), -1.d0,     &
       cholesky_ao(1,1,iblock), ao_num,       &
       X3(1,1,1,1), ao_num*block_size, 1.d0,  &
       ao_two_e_core_integral, ao_num)

 enddo

 deallocate(X2,X3)
 ao_two_e_core_integral = 2.D0 * ao_two_e_core_integral ! count for alpha + beta

END_PROVIDER

 BEGIN_PROVIDER [ double precision, ao_two_e_valence_integral_alpha, (ao_num, ao_num) ]
&BEGIN_PROVIDER [ double precision, ao_two_e_valence_integral_beta ,  (ao_num, ao_num) ]
 use map_module
 implicit none
 BEGIN_DOC
 ! Alpha and Beta Fock matrices in AO basis set
 END_DOC

 integer :: m,n,l,s,j
 double precision :: integral
 double precision, allocatable :: X(:), X2(:,:,:,:), X3(:,:,:,:)

 allocate (X(cholesky_ao_num))

 ! X(j) = \sum_{mn} SCF_density_matrix_ao(m,n) * cholesky_ao(m,n,j)
 call dgemm('T','N',cholesky_ao_num,1,ao_num*ao_num,1.d0,            &
     cholesky_ao, ao_num*ao_num,                                     &
     TCSCF_bi_ort_valence_dm_ao_full, ao_num*ao_num,0.d0,                      &
     X, cholesky_ao_num)
!

 ! ao_two_e_integral_alpha(m,n) = \sum_{j} cholesky_ao(m,n,j) * X(j)
 call dgemm('N','N',ao_num*ao_num,1,cholesky_ao_num, 1.d0,           &
     cholesky_ao, ao_num*ao_num,                                     &
     X, cholesky_ao_num, 0.d0,                                       &
     ao_two_e_valence_integral_alpha, ao_num*ao_num)

 deallocate(X)

 if (elec_alpha_num > elec_beta_num) then
   ao_two_e_valence_integral_beta = ao_two_e_valence_integral_alpha
 endif


 double precision :: rss, mem0, mem
 double precision :: memory_of_double

 integer :: iblock
 integer :: block_size

 call resident_memory(mem0)

 block_size = 1024

 rss = memory_of_double(2.d0*ao_num*ao_num)
 do
   mem = mem0 + block_size*rss
   if ( (block_size < 2).or.(mem < qp_max_mem) ) exit
   block_size = block_size/2
 enddo

 call check_mem(block_size*rss, irp_here)

 allocate(X2(ao_num,ao_num,block_size,2))
 allocate(X3(ao_num,block_size,ao_num,2))

! ao_two_e_valence_integral_alpha (l,s) -= cholesky_ao(l,m,j) * SCF_density_matrix_ao_beta (m,n) * cholesky_ao(n,s,j)

 do iblock=1,cholesky_ao_num,block_size

   call dgemm('N','N',ao_num,ao_num*min(cholesky_ao_num-iblock+1,block_size),ao_num, 1.d0,      &
       TCSCF_bi_ort_valence_dm_ao_alpha, ao_num,    &
       cholesky_ao(1,1,iblock), ao_num, 0.d0,  &
       X2(1,1,1,1), ao_num)

   if (elec_alpha_num > elec_beta_num) then
     call dgemm('N','N',ao_num,ao_num*min(cholesky_ao_num-iblock+1,block_size),ao_num, 1.d0,      &
       TCSCF_bi_ort_valence_dm_ao_beta, ao_num,     &
       cholesky_ao(1,1,iblock), ao_num, 0.d0,  &
       X2(1,1,1,2), ao_num)

     do s=1,ao_num
      do j=1,min(cholesky_ao_num-iblock+1,block_size)
       do m=1,ao_num
        X3(m,j,s,1) = X2(m,s,j,1)
        X3(m,j,s,2) = X2(m,s,j,2)
       enddo
      enddo
     enddo

   else

     do s=1,ao_num
      do j=1,min(cholesky_ao_num-iblock+1,block_size)
       do m=1,ao_num
        X3(m,j,s,1) = X2(m,s,j,1)
       enddo
      enddo
     enddo
   endif

   call dgemm('N','N',ao_num,ao_num,ao_num*min(cholesky_ao_num-iblock+1,block_size), -1.d0,     &
       cholesky_ao(1,1,iblock), ao_num,       &
       X3(1,1,1,1), ao_num*block_size, 1.d0,  &
       ao_two_e_valence_integral_alpha, ao_num)

   if (elec_alpha_num > elec_beta_num) then
     call dgemm('N','N',ao_num,ao_num,ao_num*min(cholesky_ao_num-iblock+1,block_size), -1.d0,     &
       cholesky_ao(1,1,iblock), ao_num,        &
       X3(1,1,1,2), ao_num*block_size, 1.d0,   &
       ao_two_e_valence_integral_beta, ao_num)
   endif

 enddo

 if (elec_alpha_num == elec_beta_num) then
   ao_two_e_valence_integral_beta = ao_two_e_valence_integral_alpha
 endif
 deallocate(X2,X3)

END_PROVIDER
