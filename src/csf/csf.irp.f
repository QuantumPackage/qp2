BEGIN_PROVIDER [ integer, csf_unique_vector_value_size ]
 implicit none
 BEGIN_DOC
 ! Is reset in the provider of csf_unique_vector_value
 END_DOC
 csf_unique_vector_value_size = 10
END_PROVIDER


 BEGIN_PROVIDER [ integer, N_csf ]
&BEGIN_PROVIDER [ integer, psi_configuration_to_psi_csf, (2,N_configuration) ]
&BEGIN_PROVIDER [ integer, psi_configuration_n_csf, (N_configuration) ]
&BEGIN_PROVIDER [ integer, csf_unique_vector_index, (N_det) ]
&BEGIN_PROVIDER [ double precision, csf_unique_vector_value, (csf_unique_vector_value_size) ]
 implicit none
 BEGIN_DOC
 ! To get the det-csf transformation matrix, call get_det_csf_transformation
 END_DOC

 integer :: i_det, j_det
 integer :: i_cfg, i_csf, i, j, k, lwork
 integer :: sze, istart, iend
 double precision, allocatable :: s2mat(:,:)
 double precision :: sij

 double precision, allocatable :: eigenvalues(:), work(:)
 integer :: info
 integer :: index, index_max, n_unique
 integer*2 :: compressed
 double precision, parameter :: undefined = 2.d0

 double precision, allocatable :: resize_array(:)
 logical :: found
 integer, allocatable :: unique_index(:)

 PROVIDE expected_s2 N_int psi_det 
 lwork = n_det_per_config_max**3

 allocate(eigenvalues(n_det_per_config_max), work(lwork))
 allocate(s2mat(n_det_per_config_max,n_det_per_config_max))
 allocate(unique_index(N_det))

 N_csf = 0

 n_unique = 1
 unique_index(1) = 1
 csf_unique_vector_value(1) = 1.d0
 index_max = 2

 do i_cfg=1,N_configuration

   istart = psi_configuration_to_psi_det(1,i_cfg)
   iend   = psi_configuration_to_psi_det(2,i_cfg)
   sze    = iend-istart+1

   if (sze > n_det_per_config_max) then
    call qp_bug(irp_here, -1, 'Configuration has more determinants than allowed by n_det_per_config_max')
   endif

   if (sze == 1) then

     N_csf = N_csf + 1
     csf_unique_vector_index(N_csf) = 1
     psi_configuration_to_psi_csf(1,i_cfg) = N_csf
     psi_configuration_to_psi_csf(2,i_cfg) = N_csf
     psi_configuration_n_csf(i_cfg) = 1

   else

     do j=1,sze
       j_det = psi_configuration_to_psi_det_data(j+istart-1)
       do i=1,j
         i_det = psi_configuration_to_psi_det_data(i+istart-1)
         call get_s2(psi_det(1,1,i_det),psi_det(1,1,j_det),N_int,sij)
         s2mat(i,j) = sij
         s2mat(j,i) = sij
       enddo
     enddo

     call dsyev('V','L', sze, s2mat, size(s2mat,1), eigenvalues, work, lwork, info)
     if (info /= 0) then
       call qp_bug(irp_here, info, 'Failed in dsyev')
     endif

     psi_configuration_to_psi_csf(1,i_cfg) = N_csf+1
     i_csf = 0

     do j=1,sze
       if (dabs(eigenvalues(j) - expected_s2) < 0.1d0) then
         i_csf = i_csf + 1
         N_csf = N_csf + 1

         found = .False.
         do k=2,n_unique
           index = unique_index(k)
           if (dabs(s2mat(1,j)-csf_unique_vector_value(index)) < 1.d-9) then
             ! We have found that the 1st component matches, we now check the dot product
             ! is 1.
             double precision :: dot
             dot = 0.d0
             do i=1,sze
               dot = dot + s2mat(i,j)*csf_unique_vector_value(index+i-1)
             enddo
             if (dabs(1.d0-dot) < 1.d-9) then
               csf_unique_vector_index(N_csf) = index
               found = .True.
               exit
             end if
           end if
         end do

         if (.not.found) then

           ! Reallocate csf_unique_vector_value 2x larger if needed
           if (index_max+sze > csf_unique_vector_value_size) then
             allocate(resize_array(csf_unique_vector_value_size))
             resize_array(:) = csf_unique_vector_value(:)
             deallocate(csf_unique_vector_value)
             allocate(csf_unique_vector_value(2*csf_unique_vector_value_size))
             csf_unique_vector_value(:csf_unique_vector_value_size) = resize_array(:)
             deallocate(resize_array)
             csf_unique_vector_value_size *= 2
!             SOFT_TOUCH csf_unique_vector_value_size
           endif

           n_unique = n_unique+1
           csf_unique_vector_index(N_csf) = index_max
           unique_index(n_unique) = index_max
           csf_unique_vector_value(index_max:index_max+sze-1) = s2mat(1:sze,j)
           index_max = index_max + sze
         end if

       end if
     enddo
     if (i_csf == 0) cycle
     psi_configuration_n_csf(i_cfg) = i_csf
     psi_configuration_to_psi_csf(2,i_cfg) = N_csf

   endif

 enddo

 deallocate(s2mat, work, eigenvalues)

END_PROVIDER





subroutine get_det_csf_transformation(matrix, dim1, i_cfg)
  implicit none
  double precision, intent(out) :: matrix(dim1,dim1)
  integer, intent(in) :: dim1, i_cfg

  integer :: ncsf, ndet
  ndet = psi_configuration_n_det(i_cfg)

  if (ndet == 1) then
    matrix(1,1) = 1.d0
    return
  endif

  ncsf = psi_configuration_n_csf(i_cfg)

  if (dim1 < ndet) call qp_bug(irp_here, 1, 'dimensions too small')
  if (dim1 < ncsf) call qp_bug(irp_here, 2, 'dimensions too small')

  integer :: i, j, i_csf
  integer :: index
  i_csf = psi_configuration_to_psi_csf(1,i_cfg)
  do j=1,ncsf
    i = csf_unique_vector_index(i_csf+j-1)
    matrix(1:ndet,j) = csf_unique_vector_value(i:i+ndet-1)
  enddo

end




subroutine convertWFfromDETtoCSF(N_st,psi_coef_det_in, psi_coef_csf_out)
 implicit none
 BEGIN_DOC
 ! Converts a wave functions from determinant space to CSF space
 END_DOC
 integer, intent(in)            :: N_st
 double precision, intent(in)   :: psi_coef_det_in(N_det,N_st)
 double precision, intent(out)  :: psi_coef_csf_out(N_csf,N_st)

 integer :: i_state, i, i_csf, i_cfg, j
 integer :: startdet, enddet, i_det, j_det

 double precision, allocatable :: coef(:,:), matrix(:,:)

 allocate (coef(n_det_per_config_max,N_st))
 allocate (matrix(n_det_per_config_max,n_det_per_config_max))

 do i_cfg=1,N_configuration

   ! Copy CI coefficients of the configuration in a temporary contiguous buffer

   startdet = psi_configuration_to_psi_det(1,i_cfg)
   do i_det=1,psi_configuration_n_det(i_cfg)
     j_det = psi_configuration_to_psi_det_data(startdet + i_det - 1)
     do i_state=1,N_st
       coef(i_det,i_state) = psi_coef_det_in(j_det, i_state)
     enddo
   enddo

   ! Matrix-multiplication to convert the vector into CSF space

   if (psi_configuration_n_det(i_cfg) == 1) then

     do i_state=1,N_st
       psi_coef_csf_out(psi_configuration_to_psi_csf(1,i_cfg), i_state) = coef(1,i_state)
     enddo

   else

      integer :: ncsf, ndet
      ncsf = psi_configuration_n_csf(i_cfg)
      ndet = psi_configuration_n_det(i_cfg)

      call get_det_csf_transformation(matrix, size(matrix,1), i_cfg)

      call dgemm('T','N', ncsf, N_st, ndet, 1.d0, &
                matrix, size(matrix,1), &
                coef, size(coef,1), &
                0.d0, psi_coef_csf_out( psi_configuration_to_psi_csf(1,i_cfg), 1), &
                size(psi_coef_csf_out,1))
   endif

 end do

 deallocate(coef)
end


subroutine convertWFfromCSFtoDET(N_st,psi_coef_csf_in, psi_coef_det_out)
 implicit none
 BEGIN_DOC
 ! Converts a wave functions from determinant space to CSF space
 END_DOC
 integer, intent(in)            :: N_st
 double precision, intent(in)   :: psi_coef_csf_in(N_csf,N_st)
 double precision, intent(out)  :: psi_coef_det_out(N_det,N_st)

 integer :: i_state, i, i_csf, i_cfg, j
 integer :: startdet, enddet, i_det, j_det

 double precision, allocatable :: coef(:,:), matrix(:,:)

 allocate (coef(n_det_per_config_max,N_st))
 allocate (matrix(n_det_per_config_max,n_det_per_config_max))

 do i_cfg=1,N_configuration

   ! Matrix-multiplication to convert the vector into CSF space

   if (psi_configuration_n_det(i_cfg) == 1) then

      do i_state = 1,N_st
        coef(1,i_state) = psi_coef_csf_in( psi_configuration_to_psi_csf(1,i_cfg), i_state)
      enddo

   else

      integer :: ncsf, ndet
      ncsf = psi_configuration_n_csf(i_cfg)
      ndet = psi_configuration_n_det(i_cfg)

      call get_det_csf_transformation(matrix, size(matrix,1), i_cfg)

      call dgemm('N','N', ndet, N_st, ncsf, 1.d0, &
                matrix, size(matrix,1), &
                psi_coef_csf_in( psi_configuration_to_psi_csf(1,i_cfg), 1), &
                size(psi_coef_csf_in,1), &
                0.d0, coef, size(coef,1))

    end if


   ! Copy CI coefficients of temporary contiguous buffer to the CI vector

   startdet = psi_configuration_to_psi_det(1,i_cfg)
   do i_det=1,psi_configuration_n_det(i_cfg)
     j_det = psi_configuration_to_psi_det_data(startdet + i_det - 1)
     do i_state = 1,N_st
       psi_coef_det_out(j_det, i_state) = coef(i_det,i_state)
     enddo
   enddo

 end do

 deallocate(coef)
end


BEGIN_PROVIDER [ double precision, psi_csf_coef, (N_csf, N_states) ]
 implicit none
 BEGIN_DOC
 ! Wafe function in CSF basis
 END_DOC

 call convertWFfromDETtoCSF(N_states, psi_coef, psi_csf_coef)
END_PROVIDER

