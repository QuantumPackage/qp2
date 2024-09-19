 subroutine reorder_mo_max_overlap
  implicit none
   BEGIN_DOC
   ! routines that compute the projection of each MO of the current `mo_coef` on the space spanned by the occupied orbitals of `mo_coef_begin_iteration`
   END_DOC
   integer :: i,j,k,l
   double precision, allocatable :: overlap(:,:)
   double precision, allocatable :: proj(:)
   integer, allocatable :: iorder(:)
   double precision, allocatable :: mo_coef_tmp(:,:)
   double precision, allocatable :: tmp(:,:)
   allocate(overlap(mo_num,mo_num),proj(mo_num),iorder(mo_num),mo_coef_tmp(ao_num,mo_num),tmp(mo_num,ao_num))
   
   overlap(:,:) = 0d0
   mo_coef_tmp(:,:) = 0d0
   proj(:) = 0d0
   iorder(:) = 0d0
   tmp(:,:) = 0d0
   
   ! These matrix products compute the overlap bewteen the initial and the current MOs
   call dgemm('T','N', mo_num, ao_num, ao_num, 1.d0, &
        mo_coef_begin_iteration, size(mo_coef_begin_iteration,1), &
        ao_overlap, size(ao_overlap,1), 0.d0, &
        tmp, size(tmp,1))

   call dgemm('N','N', mo_num, mo_num, ao_num, 1.d0, &
        tmp, size(tmp,1), &
        mo_coef, size(mo_coef, 1), 0.d0, &
        overlap, size(overlap,1) )
   
   
   ! for each orbital compute the best overlap
   do i = 1, mo_num
      iorder(i) = i ! initialize the iorder list as we need it to sort later
      do j = 1, elec_alpha_num
         proj(i) += overlap(j,i)*overlap(j,i) ! compute the projection of current orbital i on the occupied space of the initial orbitals
      enddo
      proj(i) = dsqrt(proj(i))
   enddo
   ! sort the list of projection to find the mos with the largest overlap
   call dsort(proj(:),iorder(:),mo_num)
   ! reorder orbitals according to projection
   do i=1,mo_num
      mo_coef_tmp(:,i) = mo_coef(:,iorder(mo_num+1-i)) 
   enddo

   ! update the orbitals
   mo_coef(:,:) =  mo_coef_tmp(:,:)
   
   ! if the determinant is open-shell we need to make sure that the singly occupied orbital correspond to the initial ones
   if (elec_alpha_num > elec_beta_num) then
      double precision, allocatable :: overlap_alpha(:,:)
      double precision, allocatable :: proj_alpha(:)
      integer, allocatable :: iorder_alpha(:)
      allocate(overlap_alpha(mo_num,elec_alpha_num),proj_alpha(elec_alpha_num),iorder_alpha(elec_alpha_num))
      overlap_alpha(:,:) = 0d0
      mo_coef_tmp(:,:) = 0d0
      proj_alpha(:) = 0d0
      iorder_alpha(:) = 0d0
      tmp(:,:) = 0d0 
      ! These matrix products compute the overlap bewteen the initial and the current MOs
      call dgemm('T','N', mo_num, ao_num, ao_num, 1.d0, &
           mo_coef_begin_iteration, size(mo_coef_begin_iteration,1), &
           ao_overlap, size(ao_overlap,1), 0.d0, &
           tmp, size(tmp,1))

      call dgemm('N','N', mo_num, elec_alpha_num, ao_num, 1.d0, &
           tmp, size(tmp,1), &
           mo_coef, size(mo_coef, 1), 0.d0, &
           overlap_alpha, size(overlap_alpha,1) )

      do i = 1, elec_alpha_num
         iorder_alpha(i) = i ! initialize the iorder list as we need it to sort later
         do j = 1, elec_beta_num
            proj_alpha(i) += overlap_alpha(j,i)*overlap_alpha(j,i) ! compute the projection of current orbital i on the beta occupied space of the initial orbitals
         enddo
         proj_alpha(i) = dsqrt(proj_alpha(i))
      enddo
      
      ! sort the list of projection to find the mos with the largest overlap
      call dsort(proj_alpha(:),iorder_alpha(:),elec_alpha_num)
      ! reorder orbitals according to projection
      do i=1,elec_alpha_num
         mo_coef_tmp(:,i) = mo_coef(:,iorder_alpha(elec_alpha_num+1-i)) 
      enddo
      do i=1,elec_alpha_num
         mo_coef(:,i) = mo_coef_tmp(:,i) 
      enddo

      deallocate(overlap_alpha, proj_alpha, iorder_alpha)
   endif
   
   deallocate(overlap, proj, iorder, mo_coef_tmp, tmp)

 end

