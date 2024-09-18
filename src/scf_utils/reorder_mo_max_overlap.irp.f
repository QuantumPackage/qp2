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
   allocate(overlap(mo_num,mo_num),proj(mo_num),iorder(mo_num),mo_coef_tmp(ao_num,mo_num))
   
   overlap(:,:) = 0d0
   mo_coef_tmp(:,:) = 0d0
   proj(:) = 0d0
   iorder(:) = 0d0
   
   ! this loop compute the overlap between the initial and the current MOS
   do i = 1, mo_num ! old mo
      do j = 1, mo_num ! curent mo
         do k = 1, ao_num
            do l = 1, ao_num
               overlap(i,j) += mo_coef_begin_iteration(k,i)* ao_overlap(k,l) * mo_coef(l,j) 
            enddo
         enddo
      enddo
   enddo

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
      proj_alpha(:) = 0d0
      iorder_alpha(:) = 0d0
      do i = 1, mo_num ! old mo
         do j = 1, elec_alpha_num ! curent mo
            do k = 1, ao_num
               do l = 1, ao_num
                  overlap_alpha(i,j) += mo_coef_begin_iteration(k,i) * ao_overlap(k,l) * mo_coef(l,j) 
               enddo
            enddo
         enddo
      enddo
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
   endif
   
   deallocate(overlap, proj, iorder, mo_coef_tmp)

  
 end

