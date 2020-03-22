subroutine orb_range_all_states_2_rdm(big_array,dim1,norb,list_orb,list_orb_reverse,ispin,u_0,N_st,sze)
   use bitmasks
   implicit none
   BEGIN_DOC
   ! if ispin == 1 :: alpha/alpha 2rdm 
   !          == 2 :: beta /beta  2rdm 
   !          == 3 :: alpha/beta  2rdm 
   !          == 4 :: spin traced 2rdm :: aa + bb + 0.5 (ab + ba))
   !
   ! Assumes that the determinants are in psi_det
   !
   ! istart, iend, ishift, istep are used in ZMQ parallelization.
   END_DOC
   integer, intent(in)             :: N_st,sze
   integer, intent(in)             :: dim1,norb,list_orb(norb),ispin
   integer, intent(in)             :: list_orb_reverse(mo_num)
   double precision, intent(inout) :: big_array(dim1,dim1,dim1,dim1,N_st)
   double precision, intent(in)    :: u_0(sze,N_st)

   integer                        :: k
   double precision, allocatable  :: u_t(:,:)
   !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: u_t
   allocate(u_t(N_st,N_det))
   do k=1,N_st
     call dset_order(u_0(1,k),psi_bilinear_matrix_order,N_det)
   enddo
   call dtranspose(                                                  &
       u_0,                                                          &
       size(u_0, 1),                                                 &
       u_t,                                                          &
       size(u_t, 1),                                                 &
       N_det, N_st)
   
   call orb_range_all_states_2_rdm_work(big_array,dim1,norb,list_orb,list_orb_reverse,ispin,u_t,N_st,sze,1,N_det,0,1)
   deallocate(u_t)
   
   do k=1,N_st
     call dset_order(u_0(1,k),psi_bilinear_matrix_order_reverse,N_det)
   enddo
   
end

subroutine orb_range_all_states_2_rdm_work(big_array,dim1,norb,list_orb,list_orb_reverse,ispin,u_t,N_st,sze,istart,iend,ishift,istep)
   use bitmasks
   implicit none
   BEGIN_DOC
   ! Computes two-rdm
   !
   ! Default should be 1,N_det,0,1
   END_DOC
   integer, intent(in)            :: N_st,sze,istart,iend,ishift,istep
   integer, intent(in)             :: dim1,norb,list_orb(norb),ispin
   integer, intent(in)             :: list_orb_reverse(mo_num)
   double precision, intent(inout) :: big_array(dim1,dim1,dim1,dim1,N_st)
   double precision, intent(in)   :: u_t(N_st,N_det)
   
   integer :: k
   
   PROVIDE N_int
   
   select case (N_int)
     case (1)
       call orb_range_all_states_2_rdm_work_1(big_array,dim1,norb,list_orb,list_orb_reverse,ispin,u_t,N_st,sze,istart,iend,ishift,istep)
     case (2)
       call orb_range_all_states_2_rdm_work_2(big_array,dim1,norb,list_orb,list_orb_reverse,ispin,u_t,N_st,sze,istart,iend,ishift,istep)
     case (3)
       call orb_range_all_states_2_rdm_work_3(big_array,dim1,norb,list_orb,list_orb_reverse,ispin,u_t,N_st,sze,istart,iend,ishift,istep)
     case (4)
       call orb_range_all_states_2_rdm_work_4(big_array,dim1,norb,list_orb,list_orb_reverse,ispin,u_t,N_st,sze,istart,iend,ishift,istep)
       case default
       call orb_range_all_states_2_rdm_work_N_int(big_array,dim1,norb,list_orb,list_orb_reverse,ispin,u_t,N_st,sze,istart,iend,ishift,istep)
   end select
end
 
 
 

 BEGIN_TEMPLATE
subroutine orb_range_all_states_2_rdm_work_$N_int(big_array,dim1,norb,list_orb,list_orb_reverse,ispin,u_t,N_st,sze,istart,iend,ishift,istep)
   use bitmasks
   implicit none
   BEGIN_DOC
   ! Computes the two rdm for the N_st vectors |u_t>
   ! if ispin == 1 :: alpha/alpha 2rdm 
   !          == 2 :: beta /beta  2rdm 
   !          == 3 :: alpha/beta  2rdm 
   !          == 4 :: spin traced 2rdm :: aa + bb + 0.5 (ab + ba))
   ! The 2rdm will be computed only on the list of orbitals list_orb, which contains norb
   ! Default should be 1,N_det,0,1 for istart,iend,ishift,istep
   END_DOC
   integer, intent(in)            :: N_st,sze,istart,iend,ishift,istep
   double precision, intent(in)   :: u_t(N_st,N_det)
   integer, intent(in)            :: dim1,norb,list_orb(norb),ispin
   integer, intent(in)             :: list_orb_reverse(mo_num)
   double precision, intent(inout) :: big_array(dim1,dim1,dim1,dim1,N_st)
   
   integer                        :: i,j,k,l
   integer                        :: k_a, k_b, l_a, l_b, m_a, m_b
   integer                        :: istate
   integer                        :: krow, kcol, krow_b, kcol_b
   integer                        :: lrow, lcol
   integer                        :: mrow, mcol
   integer(bit_kind)              :: spindet($N_int)
   integer(bit_kind)              :: tmp_det($N_int,2)
   integer(bit_kind)              :: tmp_det2($N_int,2)
   integer(bit_kind)              :: tmp_det3($N_int,2)
   integer(bit_kind), allocatable :: buffer(:,:)
   integer                        :: n_doubles
   integer, allocatable           :: doubles(:)
   integer, allocatable           :: singles_a(:)
   integer, allocatable           :: singles_b(:)
   integer, allocatable           :: idx(:), idx0(:)
   integer                        :: maxab, n_singles_a, n_singles_b, kcol_prev
   integer*8                      :: k8
   double precision,allocatable   :: c_contrib(:)

   logical                        :: alpha_alpha,beta_beta,alpha_beta,spin_trace
   integer(bit_kind)              :: orb_bitmask($N_int)
   alpha_alpha = .False.   
   beta_beta   = .False.
   alpha_beta  = .False.
   spin_trace  = .False.
   if(     ispin == 1)then
    alpha_alpha = .True.
   else if(ispin == 2)then
    beta_beta   = .True.
   else if(ispin == 3)then
    alpha_beta  = .True.
   else if(ispin == 4)then
    spin_trace  = .True.
   else
    print*,'Wrong parameter for ispin in general_2_rdm_dm_nstates_work'
    print*,'ispin = ',ispin
    stop
   endif
   
   PROVIDE N_int

   call list_to_bitstring( orb_bitmask, list_orb, norb, N_int)
 
   maxab = max(N_det_alpha_unique, N_det_beta_unique)+1
   allocate(idx0(maxab))
   
   do i=1,maxab
     idx0(i) = i
   enddo
   
   ! Prepare the array of all alpha single excitations
   ! -------------------------------------------------
   
   PROVIDE N_int nthreads_davidson
   !!$OMP PARALLEL DEFAULT(NONE) NUM_THREADS(nthreads_davidson)      &
       !    !$OMP   SHARED(psi_bilinear_matrix_rows, N_det,          &
       !    !$OMP          psi_bilinear_matrix_columns,              &
       !    !$OMP          psi_det_alpha_unique, psi_det_beta_unique,&
       !    !$OMP          n_det_alpha_unique, n_det_beta_unique, N_int,&
       !    !$OMP          psi_bilinear_matrix_transp_rows,          &
       !    !$OMP          psi_bilinear_matrix_transp_columns,       &
       !    !$OMP          psi_bilinear_matrix_transp_order, N_st,   &
       !    !$OMP          psi_bilinear_matrix_order_transp_reverse, &
       !    !$OMP          psi_bilinear_matrix_columns_loc,          &
       !    !$OMP          psi_bilinear_matrix_transp_rows_loc,      &
       !    !$OMP          istart, iend, istep, irp_here, v_t, s_t,  &
       !    !$OMP          ishift, idx0, u_t, maxab)                 &
       !    !$OMP   PRIVATE(krow, kcol, tmp_det, spindet, k_a, k_b, i,&
       !    !$OMP          lcol, lrow, l_a, l_b,                     &
       !    !$OMP          buffer, doubles, n_doubles,               &
       !    !$OMP          tmp_det2, idx, l, kcol_prev,              &
       !    !$OMP          singles_a, n_singles_a, singles_b,        &
       !    !$OMP          n_singles_b, k8)
   
   ! Alpha/Beta double excitations
   ! =============================
   
   allocate( buffer($N_int,maxab),                                   &
       singles_a(maxab),                                             &
       singles_b(maxab),                                             &
       doubles(maxab),                                               &
       idx(maxab),c_contrib(N_st))
   
   kcol_prev=-1
   
   ASSERT (iend <= N_det)
   ASSERT (istart > 0)
   ASSERT (istep  > 0)
   
   !!$OMP DO SCHEDULE(dynamic,64)
   do k_a=istart+ishift,iend,istep
     
     krow = psi_bilinear_matrix_rows(k_a)
     ASSERT (krow <= N_det_alpha_unique)
     
     kcol = psi_bilinear_matrix_columns(k_a)
     ASSERT (kcol <= N_det_beta_unique)
     
     tmp_det(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, krow)
     tmp_det(1:$N_int,2) = psi_det_beta_unique (1:$N_int, kcol)
     
     if (kcol /= kcol_prev) then
       call get_all_spin_singles_$N_int(                             &
           psi_det_beta_unique, idx0,                                &
           tmp_det(1,2), N_det_beta_unique,                          &
           singles_b, n_singles_b)
     endif
     kcol_prev = kcol
     
     ! Loop over singly excited beta columns
     ! -------------------------------------
     
     do i=1,n_singles_b
       lcol = singles_b(i)
       
       tmp_det2(1:$N_int,2) = psi_det_beta_unique(1:$N_int, lcol)
       
       l_a = psi_bilinear_matrix_columns_loc(lcol)
       ASSERT (l_a <= N_det)
       
       do j=1,psi_bilinear_matrix_columns_loc(lcol+1) - l_a
         lrow = psi_bilinear_matrix_rows(l_a)
         ASSERT (lrow <= N_det_alpha_unique)
         
         buffer(1:$N_int,j) = psi_det_alpha_unique(1:$N_int, lrow)
         
         ASSERT (l_a <= N_det)
         idx(j) = l_a
         l_a = l_a+1
       enddo
       j = j-1
       
       call get_all_spin_singles_$N_int(                             &
           buffer, idx, tmp_det(1,1), j,                             &
           singles_a, n_singles_a )
       
       ! Loop over alpha singles
       ! -----------------------
       
      if(alpha_beta.or.spin_trace)then
       do k = 1,n_singles_a
         l_a = singles_a(k)
         ASSERT (l_a <= N_det)
         
         lrow = psi_bilinear_matrix_rows(l_a)
         ASSERT (lrow <= N_det_alpha_unique)
         
         tmp_det2(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, lrow)
         c_contrib = 0.d0
         do l= 1, N_st
           c_1(l) = u_t(l,l_a)
           c_2(l) = u_t(l,k_a)
           c_contrib(l)  = c_1(l) * c_2(l) 
         enddo
         call orb_range_off_diagonal_double_to_2_rdm_ab_dm_all_states(tmp_det,tmp_det2,c_contrib,N_st,big_array,dim1,orb_bitmask,list_orb_reverse,ispin)
       enddo
      endif
       
     enddo
     
   enddo
   ! !$OMP END DO
   
   ! !$OMP DO SCHEDULE(dynamic,64)
   do k_a=istart+ishift,iend,istep
     
     
     ! Single and double alpha exitations
     ! ===================================
     
     
     ! Initial determinant is at k_a in alpha-major representation
     ! -----------------------------------------------------------------------
     
     krow = psi_bilinear_matrix_rows(k_a)
     ASSERT (krow <= N_det_alpha_unique)
     
     kcol = psi_bilinear_matrix_columns(k_a)
     ASSERT (kcol <= N_det_beta_unique)
     
     tmp_det(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, krow)
     tmp_det(1:$N_int,2) = psi_det_beta_unique (1:$N_int, kcol)
     
     ! Initial determinant is at k_b in beta-major representation
     ! ----------------------------------------------------------------------
     
     k_b = psi_bilinear_matrix_order_transp_reverse(k_a)
     ASSERT (k_b <= N_det)
     
     spindet(1:$N_int) = tmp_det(1:$N_int,1)
     
     ! Loop inside the beta column to gather all the connected alphas
     lcol = psi_bilinear_matrix_columns(k_a)
     l_a = psi_bilinear_matrix_columns_loc(lcol)
     do i=1,N_det_alpha_unique
       if (l_a > N_det) exit
       lcol = psi_bilinear_matrix_columns(l_a)
       if (lcol /= kcol) exit
       lrow = psi_bilinear_matrix_rows(l_a)
       ASSERT (lrow <= N_det_alpha_unique)
       
       buffer(1:$N_int,i) = psi_det_alpha_unique(1:$N_int, lrow)
       idx(i) = l_a
       l_a = l_a+1
     enddo
     i = i-1
     
     call get_all_spin_singles_and_doubles_$N_int(                   &
         buffer, idx, spindet, i,                                    &
         singles_a, doubles, n_singles_a, n_doubles )
     
     ! Compute Hij for all alpha singles
     ! ----------------------------------
     
     tmp_det2(1:$N_int,2) = psi_det_beta_unique (1:$N_int, kcol)
     do i=1,n_singles_a
       l_a = singles_a(i)
       ASSERT (l_a <= N_det)
       
       lrow = psi_bilinear_matrix_rows(l_a)
       ASSERT (lrow <= N_det_alpha_unique)
       
       tmp_det2(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, lrow)
       c_contrib = 0.d0
       do l= 1, N_st
         c_1(l) = u_t(l,l_a)
         c_2(l) = u_t(l,k_a)
         c_contrib(l) = c_1(l) * c_2(l) 
       enddo
       if(alpha_beta.or.spin_trace.or.alpha_alpha)then
       ! increment the alpha/beta part for single excitations
       call orb_range_off_diagonal_single_to_2_rdm_ab_dm_all_states(tmp_det, tmp_det2,c_contrib,N_st,big_array,dim1,orb_bitmask,list_orb_reverse,ispin)
       ! increment the alpha/alpha part for single excitations
       call orb_range_off_diagonal_single_to_2_rdm_aa_dm_all_states(tmp_det,tmp_det2,c_contrib,N_st,big_array,dim1,orb_bitmask,list_orb_reverse,ispin)
       endif
       
     enddo
     
     
     ! Compute Hij for all alpha doubles
     ! ----------------------------------
     
     if(alpha_alpha.or.spin_trace)then
      do i=1,n_doubles
        l_a = doubles(i)
        ASSERT (l_a <= N_det)
        
        lrow = psi_bilinear_matrix_rows(l_a)
        ASSERT (lrow <= N_det_alpha_unique)
        
        c_contrib = 0.d0
        do l= 1, N_st
          c_1(l) = u_t(l,l_a)
          c_2(l) = u_t(l,k_a)
          c_contrib(l) += c_1(l) * c_2(l) 
        enddo
        call orb_range_off_diagonal_double_to_2_rdm_aa_dm_all_states(tmp_det(1,1),psi_det_alpha_unique(1, lrow),c_contrib,N_st,big_array,dim1,orb_bitmask,list_orb_reverse,ispin)
      enddo
     endif
     
     
     ! Single and double beta excitations
     ! ==================================
     
     
     ! Initial determinant is at k_a in alpha-major representation
     ! -----------------------------------------------------------------------
     
     krow = psi_bilinear_matrix_rows(k_a)
     kcol = psi_bilinear_matrix_columns(k_a)
     
     tmp_det(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, krow)
     tmp_det(1:$N_int,2) = psi_det_beta_unique (1:$N_int, kcol)
     
     spindet(1:$N_int) = tmp_det(1:$N_int,2)
     
     ! Initial determinant is at k_b in beta-major representation
     ! -----------------------------------------------------------------------
     
     k_b = psi_bilinear_matrix_order_transp_reverse(k_a)
     ASSERT (k_b <= N_det)
     
     ! Loop inside the alpha row to gather all the connected betas
     lrow = psi_bilinear_matrix_transp_rows(k_b)
     l_b = psi_bilinear_matrix_transp_rows_loc(lrow)
     do i=1,N_det_beta_unique
       if (l_b > N_det) exit
       lrow = psi_bilinear_matrix_transp_rows(l_b)
       if (lrow /= krow) exit
       lcol = psi_bilinear_matrix_transp_columns(l_b)
       ASSERT (lcol <= N_det_beta_unique)
       
       buffer(1:$N_int,i) = psi_det_beta_unique(1:$N_int, lcol)
       idx(i) = l_b
       l_b = l_b+1
     enddo
     i = i-1
     
     call get_all_spin_singles_and_doubles_$N_int(                   &
         buffer, idx, spindet, i,                                    &
         singles_b, doubles, n_singles_b, n_doubles )
     
     ! Compute Hij for all beta singles
     ! ----------------------------------
     
     tmp_det2(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, krow)
     do i=1,n_singles_b
       l_b = singles_b(i)
       ASSERT (l_b <= N_det)
       
       lcol = psi_bilinear_matrix_transp_columns(l_b)
       ASSERT (lcol <= N_det_beta_unique)
       
       tmp_det2(1:$N_int,2) = psi_det_beta_unique (1:$N_int, lcol)
       l_a = psi_bilinear_matrix_transp_order(l_b)
       c_contrib = 0.d0
       do l= 1, N_st
         c_1(l) = u_t(l,l_a)
         c_2(l) = u_t(l,k_a)
         c_contrib(l) = c_1(l) * c_2(l) 
       enddo
       if(alpha_beta.or.spin_trace.or.beta_beta)then
        ! increment the alpha/beta  part for single excitations
        call orb_range_off_diagonal_single_to_2_rdm_ab_dm_all_states(tmp_det, tmp_det2,c_contrib,N_st,big_array,dim1,orb_bitmask,list_orb_reverse,ispin)
        ! increment the beta /beta  part for single excitations
        call orb_range_off_diagonal_single_to_2_rdm_bb_dm_all_states(tmp_det, tmp_det2,c_contrib,N_st,big_array,dim1,orb_bitmask,list_orb_reverse,ispin)
       endif
     enddo
     
     ! Compute Hij for all beta doubles
     ! ----------------------------------
     
     if(beta_beta.or.spin_trace)then
      do i=1,n_doubles
        l_b = doubles(i)
        ASSERT (l_b <= N_det)
        
        lcol = psi_bilinear_matrix_transp_columns(l_b)
        ASSERT (lcol <= N_det_beta_unique)
        
        l_a = psi_bilinear_matrix_transp_order(l_b)
        c_contrib = 0.d0
        do l= 1, N_st
          c_1(l) = u_t(l,l_a)
          c_2(l) = u_t(l,k_a)
          c_contrib(l) = c_1(l) * c_2(l) 
        enddo
        call orb_range_off_diagonal_double_to_2_rdm_bb_dm_all_states(tmp_det(1,2),psi_det_beta_unique(1, lcol),c_contrib,N_st,big_array,dim1,orb_bitmask,list_orb_reverse,ispin)
        ASSERT (l_a <= N_det)
        
      enddo
     endif
     
     
     ! Diagonal contribution
     ! =====================
     
     
     ! Initial determinant is at k_a in alpha-major representation
     ! -----------------------------------------------------------------------
     
     krow = psi_bilinear_matrix_rows(k_a)
     ASSERT (krow <= N_det_alpha_unique)
     
     kcol = psi_bilinear_matrix_columns(k_a)
     ASSERT (kcol <= N_det_beta_unique)
     
     tmp_det(1:$N_int,1) = psi_det_alpha_unique(1:$N_int, krow)
     tmp_det(1:$N_int,2) = psi_det_beta_unique (1:$N_int, kcol)
     
     double precision, external     :: diag_wee_mat_elem, diag_S_mat_elem
     
     double precision               :: c_1(N_states),c_2(N_states)
     c_contrib = 0.d0
     do l = 1, N_st
       c_1(l) = u_t(l,k_a)
       c_contrib(l) = c_1(l) * c_1(l) 
     enddo
     

     call orb_range_diagonal_contrib_to_all_2_rdm_dm_all_states(tmp_det,c_contrib,N_st,big_array,dim1,orb_bitmask,list_orb_reverse,ispin)
     
   end do
   !!$OMP END DO
   deallocate(buffer, singles_a, singles_b, doubles, idx)
   !!$OMP END PARALLEL
   
end
 
 SUBST [ N_int ]
 
 1;;
 2;;
 3;;
 4;;
 N_int;;
 
 END_TEMPLATE
 
