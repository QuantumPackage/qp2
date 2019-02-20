BEGIN_PROVIDER [ double precision, eigenvectors_Fock_matrix_mo, (ao_num,mo_num) ]
   implicit none
   BEGIN_DOC
   ! Eigenvectors of the Fock matrix in the MO basis obtained with level shift.
   END_DOC

   integer                        :: i,j
   integer                        :: liwork, lwork, n, info
   integer, allocatable           :: iwork(:)
   double precision, allocatable  :: work(:), F(:,:), S(:,:)
   double precision, allocatable  :: diag(:)


   allocate( F(mo_num,mo_num) )
   allocate (diag(mo_num) )

   do j=1,mo_num
     do i=1,mo_num
       F(i,j) = Fock_matrix_mo(i,j)
     enddo
   enddo

   if(frozen_orb_scf)then
     integer                        :: iorb,jorb
     do i = 1, n_core_orb
      iorb = list_core(i)
      do j = 1, n_act_orb
       jorb = list_act(j)
       F(iorb,jorb) = 0.d0
       F(jorb,iorb) = 0.d0
      enddo
     enddo
   endif

   ! Insert level shift here
   do i = elec_beta_num+1, elec_alpha_num
     F(i,i) += 0.5d0*level_shift
   enddo

   do i = elec_alpha_num+1, mo_num
     F(i,i) += level_shift
   enddo

   n = mo_num
   lwork = 1+6*n + 2*n*n
   liwork = 3 + 5*n

   allocate(work(lwork))
   allocate(iwork(liwork) )

   lwork = -1
   liwork = -1

   call dsyevd( 'V', 'U', mo_num, F,                             &
       size(F,1), diag, work, lwork, iwork, liwork, info)

   if (info /= 0) then
     print *,  irp_here//' DSYEVD failed : ', info
     stop 1
   endif
   lwork = int(work(1))
   liwork = iwork(1)
   deallocate(iwork)
   deallocate(work)

   allocate(work(lwork))
   allocate(iwork(liwork) )
   call dsyevd( 'V', 'U', mo_num, F,                             &
       size(F,1), diag, work, lwork, iwork, liwork, info)
   deallocate(iwork)


   if (info /= 0) then
     call dsyev( 'V', 'L', mo_num, F,                            &
         size(F,1), diag, work, lwork, info)

     if (info /= 0) then
       print *,  irp_here//' DSYEV failed : ', info
       stop 1
     endif
   endif

   call dgemm('N','N',ao_num,mo_num,mo_num, 1.d0,            &
       mo_coef, size(mo_coef,1), F, size(F,1),                       &
       0.d0, eigenvectors_Fock_matrix_mo, size(eigenvectors_Fock_matrix_mo,1))
   deallocate(work, F, diag)


END_PROVIDER

