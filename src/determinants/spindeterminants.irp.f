!==============================================================================!
!                                                                              !
!                          Independent alpha/beta parts                        !
!                                                                              !
!==============================================================================!

use bitmasks

integer*8 function spin_det_search_key(det,Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
! Returns an integer(8) corresponding to a determinant index for searching
  END_DOC
  integer, intent(in) :: Nint
  integer(bit_kind), intent(in) :: det(Nint)
  integer(bit_kind), parameter :: unsigned_shift = -huge(1_bit_kind) ! 100...00
  integer :: i
  spin_det_search_key = det(1)
  do i=2,Nint
    spin_det_search_key = ieor(spin_det_search_key,det(i))
  enddo
  spin_det_search_key = spin_det_search_key+unsigned_shift
end


BEGIN_PROVIDER [ integer(bit_kind), psi_det_alpha, (N_int,psi_det_size) ]
  implicit none
  BEGIN_DOC
  ! List of $\alpha$ determinants of psi_det
  END_DOC
  integer                        :: i,k

  do i=1,N_det
    do k=1,N_int
      psi_det_alpha(k,i) = psi_det(k,1,i)
    enddo
  enddo

END_PROVIDER

BEGIN_PROVIDER [ integer(bit_kind), psi_det_beta, (N_int,psi_det_size) ]
  implicit none
  BEGIN_DOC
  ! List of $\beta$ determinants of psi_det
  END_DOC
  integer                        :: i,k

  do i=1,N_det
    do k=1,N_int
      psi_det_beta(k,i) = psi_det(k,2,i)
    enddo
  enddo
END_PROVIDER


BEGIN_TEMPLATE

 BEGIN_PROVIDER [ integer(bit_kind), psi_det_$alpha_unique, (N_int,psi_det_size) ]
&BEGIN_PROVIDER [ integer, N_det_$alpha_unique ]
  implicit none
  BEGIN_DOC
  ! Unique $\\$alpha$ determinants
  END_DOC

  integer                        :: i,j,k
  integer, allocatable           :: iorder(:)
  integer*8, allocatable         :: bit_tmp(:)
  integer*8                      :: last_key
  integer*8, external            :: spin_det_search_key
  logical,allocatable            :: duplicate(:)

  allocate ( iorder(N_det), bit_tmp(N_det), duplicate(N_det) )

  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
  do i=1,N_det
    iorder(i) = i
    bit_tmp(i) = spin_det_search_key(psi_det_$alpha(1,i),N_int)
  enddo
  !$OMP END PARALLEL DO


  call i8sort(bit_tmp,iorder,N_det)

  N_det_$alpha_unique = 0
  last_key = 0_8
  do i=1,N_det
    last_key = bit_tmp(i)
    N_det_$alpha_unique += 1
    do k=1,N_int
      psi_det_$alpha_unique(k,N_det_$alpha_unique) = psi_det_$alpha(k,iorder(i))
    enddo
    duplicate(i) = .False.
  enddo

  j=1
  do i=1,N_det_$alpha_unique-1
    if (duplicate(i)) then
      cycle
    endif
    j = i+1
    do while (bit_tmp(j)==bit_tmp(i))
      if (duplicate(j)) then
        j += 1
        cycle
      endif
      duplicate(j) = .True.
      do k=1,N_int
        if (psi_det_$alpha_unique(k,i) /= psi_det_$alpha_unique(k,j)) then
          duplicate(j) = .False.
          exit
        endif
      enddo
      j+=1
      if (j > N_det_$alpha_unique) then
        exit
      endif
    enddo
  enddo

  j=1
  do i=2,N_det_$alpha_unique
    if (duplicate(i)) then
      cycle
    else
      j += 1
      psi_det_$alpha_unique(:,j) = psi_det_$alpha_unique(:,i)
    endif
  enddo
  N_det_$alpha_unique = j

  call write_int(6,N_det_$alpha_unique,'Number of unique $alpha determinants')

  deallocate (iorder, bit_tmp, duplicate)

END_PROVIDER


SUBST [ alpha ]

alpha ;;
beta ;;

END_TEMPLATE




integer function get_index_in_psi_det_alpha_unique(key,Nint)
  use bitmasks
  BEGIN_DOC
  ! Returns the index of the determinant in the :c:data:`psi_det_alpha_unique` array
  END_DOC
  implicit none

  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key(Nint)

  integer                        :: i, ibegin, iend, istep, l
  integer(8)                     :: det_ref, det_search
  integer(8), external           :: spin_det_search_key
  logical                        :: in_wavefunction

  in_wavefunction = .False.
  get_index_in_psi_det_alpha_unique = 0
  ibegin = 1
  iend   = N_det_alpha_unique + 1

  !DIR$ FORCEINLINE
  det_ref = spin_det_search_key(key,Nint)

  !DIR$ FORCEINLINE
  det_search = spin_det_search_key(psi_det_alpha_unique(1,1),Nint)

  istep = shiftr(iend-ibegin,1)
  i=ibegin+istep
  do while (istep > 0)
    !DIR$ FORCEINLINE
    det_search = spin_det_search_key(psi_det_alpha_unique(1,i),Nint)
    if ( det_search > det_ref ) then
      iend = i
    else if ( det_search == det_ref ) then
      exit
    else
      ibegin = i
    endif
    istep = shiftr(iend-ibegin,1)
    i = ibegin + istep
  end do

  !DIR$ FORCEINLINE
  do while (spin_det_search_key(psi_det_alpha_unique(1,i),Nint) == det_ref)
    i = i-1
    if (i == 0) then
      exit
    endif
  enddo
  i += 1

  ASSERT (i <= N_det_alpha_unique)

  !DIR$ FORCEINLINE
  do while (spin_det_search_key(psi_det_alpha_unique(1,i),Nint) == det_ref)
    if (key(1) /= psi_det_alpha_unique(1,i)) then
      continue
    else
      in_wavefunction = .True.
      !DIR$ IVDEP
      !DIR$ LOOP COUNT MIN(3)
      do l=2,Nint
        if (key(l) /= psi_det_alpha_unique(l,i)) then
          in_wavefunction = .False.
        endif
      enddo
      if (in_wavefunction) then
        get_index_in_psi_det_alpha_unique = i
        return
      endif
    endif
    i += 1
    if (i > N_det_alpha_unique) then
      ASSERT (get_index_in_psi_det_alpha_unique > 0)
      return
    endif

  enddo

end

integer function get_index_in_psi_det_beta_unique(key,Nint)
  use bitmasks
  BEGIN_DOC
  ! Returns the index of the determinant in the :c:data:`psi_det_beta_unique` array
  END_DOC
  implicit none

  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key(Nint)

  integer                        :: i, ibegin, iend, istep, l
  integer(8)                     :: det_ref, det_search
  integer(8), external           :: spin_det_search_key
  logical                        :: in_wavefunction

  in_wavefunction = .False.
  get_index_in_psi_det_beta_unique = 0
  ibegin = 1
  iend   = N_det_beta_unique + 1

  !DIR$ FORCEINLINE
  det_ref = spin_det_search_key(key,Nint)
  !DIR$ FORCEINLINE
  det_search = spin_det_search_key(psi_det_beta_unique(1,1),Nint)

  istep = shiftr(iend-ibegin,1)
  i=ibegin+istep
  do while (istep > 0)
    !DIR$ FORCEINLINE
    det_search = spin_det_search_key(psi_det_beta_unique(1,i),Nint)
    if ( det_search > det_ref ) then
      iend = i
    else if ( det_search == det_ref ) then
      exit
    else
      ibegin = i
    endif
    istep = shiftr(iend-ibegin,1)
    i = ibegin + istep
  end do

  !DIR$ FORCEINLINE
  do while (spin_det_search_key(psi_det_beta_unique(1,i),Nint) == det_ref)
    i = i-1
    if (i == 0) then
      exit
    endif
  enddo
  i += 1

  ASSERT (i <= N_det_beta_unique)

  !DIR$ FORCEINLINE
  do while (spin_det_search_key(psi_det_beta_unique(1,i),Nint) == det_ref)
    if (key(1) /= psi_det_beta_unique(1,i)) then
      continue
    else
      in_wavefunction = .True.
      !DIR$ IVDEP
      !DIR$ LOOP COUNT MIN(3)
      do l=2,Nint
        if (key(l) /= psi_det_beta_unique(l,i)) then
          in_wavefunction = .False.
        endif
      enddo
      if (in_wavefunction) then
        get_index_in_psi_det_beta_unique = i
        return
      endif
    endif
    i += 1
    if (i > N_det_beta_unique) then
      ASSERT (get_index_in_psi_det_beta_unique > 0)
      return
    endif

  enddo

end


subroutine write_spindeterminants
  use bitmasks
  implicit none
  integer(8), allocatable        :: tmpdet(:,:)
  integer                        :: N_int2
  integer                        :: i,j,k
  integer(8)                     :: det_8(100)
  integer(bit_kind)              :: det_bk((100*8)/bit_kind)
  equivalence (det_8, det_bk)

  N_int2 = (N_int*bit_kind)/8
  call ezfio_set_spindeterminants_n_det_alpha(N_det_alpha_unique)
  call ezfio_set_spindeterminants_n_det_beta(N_det_beta_unique)
  call ezfio_set_spindeterminants_n_det(N_det)
  call ezfio_set_spindeterminants_n_int(N_int)
  call ezfio_set_spindeterminants_bit_kind(bit_kind)
  call ezfio_set_spindeterminants_n_states(N_states)

  allocate(tmpdet(N_int2,N_det_alpha_unique))
  do i=1,N_det_alpha_unique
    do k=1,N_int
      det_bk(k) = psi_det_alpha_unique(k,i)
    enddo
    do k=1,N_int2
      tmpdet(k,i) = det_8(k)
    enddo
  enddo
  call ezfio_set_spindeterminants_psi_det_alpha(psi_det_alpha_unique)
  deallocate(tmpdet)

  allocate(tmpdet(N_int2,N_det_beta_unique))
  do i=1,N_det_beta_unique
    do k=1,N_int
      det_bk(k) = psi_det_beta_unique(k,i)
    enddo
    do k=1,N_int2
      tmpdet(k,i) = det_8(k)
    enddo
  enddo
  call ezfio_set_spindeterminants_psi_det_beta(psi_det_beta_unique)
  deallocate(tmpdet)

  call ezfio_set_spindeterminants_psi_coef_matrix_values(psi_bilinear_matrix_values)
  call ezfio_set_spindeterminants_psi_coef_matrix_rows(psi_bilinear_matrix_rows)
  call ezfio_set_spindeterminants_psi_coef_matrix_columns(psi_bilinear_matrix_columns)

end

subroutine read_spindeterminants
  use bitmasks
  implicit none
  integer                        :: k

  call ezfio_get_spindeterminants_n_det(N_det)
  call ezfio_get_spindeterminants_n_states(N_states)
  TOUCH N_det N_states

  call ezfio_get_spindeterminants_n_det_alpha(N_det_alpha_unique)
  call ezfio_get_spindeterminants_n_det_beta(N_det_beta_unique)
  call ezfio_get_spindeterminants_psi_coef_matrix_values(psi_bilinear_matrix_values)
  call ezfio_get_spindeterminants_psi_coef_matrix_rows(psi_bilinear_matrix_rows)
  call ezfio_get_spindeterminants_psi_coef_matrix_columns(psi_bilinear_matrix_columns)
  call ezfio_get_spindeterminants_psi_det_alpha(psi_det_alpha_unique)
  call ezfio_get_spindeterminants_psi_det_beta(psi_det_beta_unique)
  do k=1,N_det
    psi_bilinear_matrix_order(k) = k
  enddo
  TOUCH psi_bilinear_matrix_values psi_bilinear_matrix_rows psi_bilinear_matrix_columns N_det_alpha_unique N_det_beta_unique psi_det_alpha_unique psi_det_beta_unique psi_bilinear_matrix_order

  call wf_of_psi_bilinear_matrix(.True.)

end

 BEGIN_PROVIDER [ double precision, det_alpha_norm, (N_det_alpha_unique) ]
&BEGIN_PROVIDER [ double precision, det_beta_norm, (N_det_beta_unique) ]
  implicit none
  BEGIN_DOC
  ! Norm of the $\alpha$ and $\beta$ spin determinants in the wave function:
  !
  ! $||D_\alpha||_i = \sum_j C_{ij}^2$
  END_DOC

  integer                        :: i,j,k,l
  double precision               :: f

  det_alpha_norm = 0.d0
  det_beta_norm  = 0.d0
  do k=1,N_det
    i = psi_bilinear_matrix_rows(k)
    j = psi_bilinear_matrix_columns(k)
    f = 0.d0
    do l=1,N_states
      f += psi_bilinear_matrix_values(k,l)*psi_bilinear_matrix_values(k,l) * state_average_weight(l)
    enddo
    det_alpha_norm(i) += f
    det_beta_norm(j)  += f
  enddo
  det_alpha_norm = det_alpha_norm
  det_beta_norm = det_beta_norm

END_PROVIDER


!==============================================================================!
!                                                                              !
!                               Alpha x Beta Matrix                            !
!                                                                              !
!==============================================================================!

 BEGIN_PROVIDER  [ double precision, psi_bilinear_matrix_values, (N_det,N_states) ]
&BEGIN_PROVIDER [ integer, psi_bilinear_matrix_rows   , (N_det) ]
&BEGIN_PROVIDER [ integer, psi_bilinear_matrix_columns, (N_det) ]
&BEGIN_PROVIDER [ integer, psi_bilinear_matrix_order  , (N_det) ]
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Sparse coefficient matrix if the wave function is expressed in a bilinear form :
  !  $D_\alpha^\dagger.C.D_\beta$
  !
  ! Rows are $\alpha$ determinants and columns are $\beta$.
  !
  ! Order refers to psi_det
  END_DOC
  integer                        :: i,j,k, l
  integer(bit_kind)              :: tmp_det(N_int,2)
  integer, external              :: get_index_in_psi_det_sorted_bit


  PROVIDE psi_coef_sorted_bit

  integer*8, allocatable         :: to_sort(:)
  integer, external              :: get_index_in_psi_det_alpha_unique
  integer, external              :: get_index_in_psi_det_beta_unique
  allocate(to_sort(N_det))
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,l)
  do k=1,N_det
    i = get_index_in_psi_det_alpha_unique(psi_det(1,1,k),N_int)
    ASSERT (i>0)
    ASSERT (i<=N_det_alpha_unique)

    j = get_index_in_psi_det_beta_unique (psi_det(1,2,k),N_int)
    ASSERT (j>0)
    ASSERT (j<=N_det_beta_unique)

    do l=1,N_states
      psi_bilinear_matrix_values(k,l) = psi_coef(k,l)
    enddo
    psi_bilinear_matrix_rows(k) = i
    psi_bilinear_matrix_columns(k) = j
    to_sort(k) = int(N_det_alpha_unique,8) * int(j-1,8) + int(i,8)
    ASSERT (to_sort(k) > 0_8)
    psi_bilinear_matrix_order(k) = k
  enddo
  !$OMP END PARALLEL DO
  call i8sort(to_sort, psi_bilinear_matrix_order, N_det)
  !$OMP PARALLEL
  !$OMP SINGLE
  call iset_order(psi_bilinear_matrix_rows,psi_bilinear_matrix_order,N_det)
  !$OMP END SINGLE
  !$OMP SINGLE
  call iset_order(psi_bilinear_matrix_columns,psi_bilinear_matrix_order,N_det)
  !$OMP END SINGLE
  !$OMP DO
  do l=1,N_states
    call dset_order(psi_bilinear_matrix_values(1,l),psi_bilinear_matrix_order,N_det)
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  deallocate(to_sort)
  ASSERT (minval(psi_bilinear_matrix_rows) == 1)
  ASSERT (minval(psi_bilinear_matrix_columns) == 1)
  ASSERT (minval(psi_bilinear_matrix_order) == 1)
  ASSERT (maxval(psi_bilinear_matrix_rows) == N_det_alpha_unique)
  ASSERT (maxval(psi_bilinear_matrix_columns) == N_det_beta_unique)
  ASSERT (maxval(psi_bilinear_matrix_order) == N_det)

END_PROVIDER


BEGIN_PROVIDER [ integer, psi_bilinear_matrix_order_reverse , (N_det) ]
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Order which allows to go from :c:data:`psi_bilinear_matrix` to :c:data:`psi_det`
  END_DOC
  integer                        :: k
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k)
  do k=1,N_det
    psi_bilinear_matrix_order_reverse(psi_bilinear_matrix_order(k)) = k
  enddo
  !$OMP END PARALLEL DO
  ASSERT (minval(psi_bilinear_matrix_order) == 1)
  ASSERT (maxval(psi_bilinear_matrix_order) == N_det)
END_PROVIDER


BEGIN_PROVIDER [ integer, psi_bilinear_matrix_columns_loc, (N_det_beta_unique+1) ]
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Sparse coefficient matrix if the wave function is expressed in a bilinear form :
  !
  ! $D_\alpha^\dagger.C.D_\beta$
  !
  ! Rows are $\alpha$ determinants and columns are $\beta$.
  !
  ! Order refers to :c:data:`psi_det`
  END_DOC
  integer                        :: i,j,k, l

  l = psi_bilinear_matrix_columns(1)
  psi_bilinear_matrix_columns_loc(l) = 1
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k,l)
  do k=2,N_det
    if (psi_bilinear_matrix_columns(k) == psi_bilinear_matrix_columns(k-1)) then
      cycle
    else
      l = psi_bilinear_matrix_columns(k)
      psi_bilinear_matrix_columns_loc(l) = k
    endif
    if (psi_bilinear_matrix_columns(k) < 1) then
      stop '(psi_bilinear_matrix_columns(k) < 1)'
    endif
  enddo
  !$OMP END PARALLEL DO
  psi_bilinear_matrix_columns_loc(N_det_beta_unique+1) = N_det+1
  ASSERT (minval(psi_bilinear_matrix_columns_loc) == 1)
  ASSERT (maxval(psi_bilinear_matrix_columns_loc) == N_det+1)

END_PROVIDER

 BEGIN_PROVIDER  [ double precision, psi_bilinear_matrix_transp_values, (N_det,N_states) ]
&BEGIN_PROVIDER [ integer, psi_bilinear_matrix_transp_rows   , (N_det) ]
&BEGIN_PROVIDER [ integer, psi_bilinear_matrix_transp_columns, (N_det) ]
&BEGIN_PROVIDER [ integer, psi_bilinear_matrix_transp_order  , (N_det) ]
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Transpose of :c:data:`psi_bilinear_matrix`
  !
  ! $D_\beta^\dagger.C^\dagger.D_\alpha$
  !
  ! Rows are $\alpha$ determinants and columns are $\beta$, but the matrix is stored in row major
  ! format.
  END_DOC
  integer                        :: i,j,k,l

  PROVIDE psi_coef_sorted_bit

  integer*8, allocatable         :: to_sort(:)
  allocate(to_sort(N_det))
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,l)
  do l=1,N_states
    !$OMP DO
    do k=1,N_det
      psi_bilinear_matrix_transp_values (k,l) = psi_bilinear_matrix_values (k,l)
    enddo
    !$OMP ENDDO NOWAIT
  enddo
  !$OMP DO
  do k=1,N_det
    psi_bilinear_matrix_transp_columns(k) = psi_bilinear_matrix_columns(k)
    ASSERT (psi_bilinear_matrix_transp_columns(k) > 0)
    ASSERT (psi_bilinear_matrix_transp_columns(k) <= N_det)

    psi_bilinear_matrix_transp_rows   (k) = psi_bilinear_matrix_rows   (k)
    ASSERT (psi_bilinear_matrix_transp_rows(k) > 0)
    ASSERT (psi_bilinear_matrix_transp_rows(k) <= N_det)

    i = psi_bilinear_matrix_transp_columns(k)
    j = psi_bilinear_matrix_transp_rows   (k)
    to_sort(k) = int(N_det_beta_unique,8) * int(j-1,8) + int(i,8)
    ASSERT (to_sort(k) > 0)
    psi_bilinear_matrix_transp_order(k) = k
  enddo
  !$OMP ENDDO
  !$OMP END PARALLEL
  call i8sort(to_sort, psi_bilinear_matrix_transp_order, N_det)
  call iset_order(psi_bilinear_matrix_transp_rows,psi_bilinear_matrix_transp_order,N_det)
  call iset_order(psi_bilinear_matrix_transp_columns,psi_bilinear_matrix_transp_order,N_det)
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(l)
  do l=1,N_states
    call dset_order(psi_bilinear_matrix_transp_values(1,l),psi_bilinear_matrix_transp_order,N_det)
  enddo
  !$OMP END PARALLEL DO
  deallocate(to_sort)
  ASSERT (minval(psi_bilinear_matrix_transp_columns) == 1)
  ASSERT (minval(psi_bilinear_matrix_transp_rows) == 1)
  ASSERT (minval(psi_bilinear_matrix_transp_order) == 1)
  ASSERT (maxval(psi_bilinear_matrix_transp_columns) == N_det_beta_unique)
  ASSERT (maxval(psi_bilinear_matrix_transp_rows) == N_det_alpha_unique)
  ASSERT (maxval(psi_bilinear_matrix_transp_order) == N_det)

END_PROVIDER

BEGIN_PROVIDER [ integer, psi_bilinear_matrix_transp_rows_loc, (N_det_alpha_unique+1) ]
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Location of the columns in the :c:data:`psi_bilinear_matrix`
  END_DOC
  integer                        :: i,j,k, l

  l = psi_bilinear_matrix_transp_rows(1)
  psi_bilinear_matrix_transp_rows_loc(l) = 1
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k,l)
  do k=2,N_det
    if (psi_bilinear_matrix_transp_rows(k) == psi_bilinear_matrix_transp_rows(k-1)) then
      cycle
    else
      l = psi_bilinear_matrix_transp_rows(k)
      psi_bilinear_matrix_transp_rows_loc(l) = k
    endif
  enddo
  !$OMP END PARALLEL DO
  psi_bilinear_matrix_transp_rows_loc(N_det_alpha_unique+1) = N_det+1
  ASSERT (minval(psi_bilinear_matrix_transp_rows_loc) == 1)
  ASSERT (maxval(psi_bilinear_matrix_transp_rows_loc) == N_det+1)

END_PROVIDER

BEGIN_PROVIDER [ integer, psi_bilinear_matrix_order_transp_reverse , (N_det) ]
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Order which allows to go from :c:data:`psi_bilinear_matrix_order_transp` to
  ! :c:data:`psi_bilinear_matrix`
  END_DOC
  integer                        :: k
  psi_bilinear_matrix_order_transp_reverse = -1
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(k)
  do k=1,N_det
    psi_bilinear_matrix_order_transp_reverse(psi_bilinear_matrix_transp_order(k)) = k
  enddo
  !$OMP END PARALLEL DO
  ASSERT (minval(psi_bilinear_matrix_order_transp_reverse) == 1)
  ASSERT (maxval(psi_bilinear_matrix_order_transp_reverse) == N_det)
END_PROVIDER


BEGIN_PROVIDER [ double precision, psi_bilinear_matrix, (N_det_alpha_unique,N_det_beta_unique,N_states) ]
  implicit none
  BEGIN_DOC
  ! Coefficient matrix if the wave function is expressed in a bilinear form :
  !
  ! $D_\alpha^\dagger.C.D_\beta$
  END_DOC
  integer                        :: i,j,k,istate
  psi_bilinear_matrix = 0.d0
  do k=1,N_det
    i = psi_bilinear_matrix_rows(k)
    j = psi_bilinear_matrix_columns(k)
    do istate=1,N_states
      psi_bilinear_matrix(i,j,istate) = psi_bilinear_matrix_values(k,istate)
    enddo
  enddo
END_PROVIDER

subroutine create_wf_of_psi_bilinear_matrix(truncate)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Generates a wave function containing all possible products
  ! of $\alpha$ and $\beta$ determinants
  END_DOC
  logical, intent(in)            :: truncate

  call generate_all_alpha_beta_det_products
  call update_wf_of_psi_bilinear_matrix(truncate)

end

subroutine update_wf_of_psi_bilinear_matrix(truncate)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Updates a wave function when psi_bilinear_matrix was modified
  END_DOC
  logical, intent(in)            :: truncate
  integer                        :: i,j,k
  integer(bit_kind)              :: tmp_det(N_int,2)
  integer                        :: idx
  integer, external              :: get_index_in_psi_det_sorted_bit
  double precision               :: norm(N_states)
  PROVIDE psi_bilinear_matrix

  norm = 0.d0
  !$OMP PARALLEL DO DEFAULT(NONE)                                    &
      !$OMP PRIVATE(i,j,k,idx,tmp_det)                               &
      !$OMP SHARED(N_det_alpha_unique, N_det_beta_unique, N_det,     &
      !$OMP N_int, N_states, norm, psi_det_beta_unique,              &
      !$OMP psi_det_alpha_unique, psi_bilinear_matrix,               &
      !$OMP psi_coef_sorted_bit)
  do j=1,N_det_beta_unique
    do k=1,N_int
      tmp_det(k,2) = psi_det_beta_unique(k,j)
    enddo
    do i=1,N_det_alpha_unique
      do k=1,N_int
        tmp_det(k,1) = psi_det_alpha_unique(k,i)
      enddo
      idx = get_index_in_psi_det_sorted_bit(tmp_det,N_int)
      if (idx > 0) then
        do k=1,N_states
          psi_coef_sorted_bit(idx,k) = psi_bilinear_matrix(i,j,k)
          !$OMP ATOMIC
          norm(k) += psi_bilinear_matrix(i,j,k)*psi_bilinear_matrix(i,j,k)
        enddo
      endif
    enddo
  enddo
  !$OMP END PARALLEL DO

  do k=1,N_states
    norm(k) = 1.d0/dsqrt(norm(k))
    do i=1,N_det
      psi_coef_sorted_bit(i,k) = psi_coef_sorted_bit(i,k)*norm(k)
    enddo
  enddo
  psi_det  = psi_det_sorted_bit
  psi_coef = psi_coef_sorted_bit
  TOUCH psi_det psi_coef N_det_beta_unique N_det_alpha_unique psi_det_beta_unique psi_det_alpha_unique
  psi_det  = psi_det_sorted
  psi_coef = psi_coef_sorted
  norm(1) = 0.d0
  do i=1,N_det
    norm(1) += psi_average_norm_contrib_sorted(i)
    if (truncate) then
      if (norm(1) >= 1.d0) then
        exit
      endif
      if (psi_average_norm_contrib_sorted(i) == 0.d0) then
        exit
      endif
    endif
  enddo
  N_det = min(i,N_det)
  SOFT_TOUCH psi_det psi_coef N_det N_det_beta_unique N_det_alpha_unique psi_det_beta_unique psi_det_alpha_unique

end

subroutine generate_all_alpha_beta_det_products
  implicit none
  BEGIN_DOC
  ! Creates a wave function from all possible $\alpha \times \beta$ determinants
  END_DOC
  integer                        :: i,j,k,l
  integer                        :: iproc
  integer, external              :: get_index_in_psi_det_sorted_bit
  integer(bit_kind), allocatable :: tmp_det(:,:,:)
  logical, external              :: is_in_wavefunction
  PROVIDE H_apply_buffer_allocated

  !$OMP PARALLEL DEFAULT(NONE) SHARED(psi_coef_sorted_bit,N_det_beta_unique,&
      !$OMP N_det_alpha_unique, N_int, psi_det_alpha_unique, psi_det_beta_unique,&
      !$OMP N_det)                                                   &
      !$OMP PRIVATE(i,j,k,l,tmp_det,iproc)
  !$ iproc = omp_get_thread_num()
  allocate (tmp_det(N_int,2,N_det_alpha_unique))
  !$OMP DO SCHEDULE(static,8)
  do j=1,N_det_beta_unique
    l = 1
    do i=1,N_det_alpha_unique
      do k=1,N_int
        tmp_det(k,1,l) = psi_det_alpha_unique(k,i)
        tmp_det(k,2,l) = psi_det_beta_unique (k,j)
      enddo
      if (.not.is_in_wavefunction(tmp_det(1,1,l),N_int)) then
        l = l+1
      endif
    enddo
    call fill_H_apply_buffer_no_selection(l-1, tmp_det, N_int, iproc)
  enddo
  !$OMP END DO
  deallocate(tmp_det)
  !$OMP END PARALLEL
  call copy_H_apply_buffer_to_wf
  SOFT_TOUCH psi_det psi_coef N_det N_det_beta_unique N_det_alpha_unique psi_det_alpha_unique psi_det_beta_unique
end




subroutine get_all_spin_singles_and_doubles(buffer, idx, spindet, Nint, size_buffer, singles, doubles, n_singles, n_doubles)
  use bitmasks
  implicit none
  BEGIN_DOC
  !
  ! Returns the indices of all the single and double excitations in the list of
  ! unique $\alpha$ determinants.
  !
  ! Warning: The buffer is transposed.
  !
  END_DOC
  integer, intent(in)            :: Nint, size_buffer, idx(size_buffer)
  integer(bit_kind), intent(in)  :: buffer(Nint,size_buffer)
  integer(bit_kind), intent(in)  :: spindet(Nint)
  integer, intent(out)           :: singles(size_buffer)
  integer, intent(out)           :: doubles(size_buffer)
  integer, intent(out)           :: n_singles
  integer, intent(out)           :: n_doubles

  select case (Nint)
    case (1)
      call get_all_spin_singles_and_doubles_1(buffer, idx, spindet(1), size_buffer, singles, doubles, n_singles, n_doubles)
    case (2)
      call get_all_spin_singles_and_doubles_2(buffer, idx, spindet, size_buffer, singles, doubles, n_singles, n_doubles)
    case (3)
      call get_all_spin_singles_and_doubles_3(buffer, idx, spindet, size_buffer, singles, doubles, n_singles, n_doubles)
    case (4)
      call get_all_spin_singles_and_doubles_4(buffer, idx, spindet, size_buffer, singles, doubles, n_singles, n_doubles)
      case default
      call get_all_spin_singles_and_doubles_N_int(buffer, idx, spindet, size_buffer, singles, doubles, n_singles, n_doubles)
  end select

end


subroutine get_all_spin_singles(buffer, idx, spindet, Nint, size_buffer, singles, n_singles)
  use bitmasks
  implicit none
  BEGIN_DOC
  !
  ! Returns the indices of all the single excitations in the list of
  ! unique $\alpha$ determinants.
  !
  END_DOC
  integer, intent(in)            :: Nint, size_buffer, idx(size_buffer)
  integer(bit_kind), intent(in)  :: buffer(Nint,size_buffer)
  integer(bit_kind), intent(in)  :: spindet(Nint)
  integer, intent(out)           :: singles(size_buffer)
  integer, intent(out)           :: n_singles

  select case (N_int)
    case (1)
      call get_all_spin_singles_1(buffer, idx, spindet(1), size_buffer, singles, n_singles)
      return
    case (2)
      call get_all_spin_singles_2(buffer, idx, spindet, size_buffer, singles, n_singles)
    case (3)
      call get_all_spin_singles_3(buffer, idx, spindet, size_buffer, singles, n_singles)
    case (4)
      call get_all_spin_singles_4(buffer, idx, spindet, size_buffer, singles, n_singles)
      case default
      call get_all_spin_singles_N_int(buffer, idx, spindet, size_buffer, singles, n_singles)
  end select

end


subroutine get_all_spin_doubles(buffer, idx, spindet, Nint, size_buffer, doubles, n_doubles)
  use bitmasks
  implicit none
  BEGIN_DOC
  !
  ! Returns the indices of all the double excitations in the list of
  ! unique $\alpha$ determinants.
  !
  END_DOC
  integer, intent(in)            :: Nint, size_buffer, idx(size_buffer)
  integer(bit_kind), intent(in)  :: buffer(Nint,size_buffer)
  integer(bit_kind), intent(in)  :: spindet(Nint)
  integer, intent(out)           :: doubles(size_buffer)
  integer, intent(out)           :: n_doubles

  select case (N_int)
    case (1)
      call get_all_spin_doubles_1(buffer, idx, spindet(1), size_buffer, doubles, n_doubles)
    case (2)
      call get_all_spin_doubles_2(buffer, idx, spindet, size_buffer, doubles, n_doubles)
    case (3)
      call get_all_spin_doubles_3(buffer, idx, spindet, size_buffer, doubles, n_doubles)
    case (4)
      call get_all_spin_doubles_4(buffer, idx, spindet, size_buffer, doubles, n_doubles)
      case default
      call get_all_spin_doubles_N_int(buffer, idx, spindet, size_buffer, doubles, n_doubles)
  end select

end





subroutine copy_psi_bilinear_to_psi(psi, isize)
  implicit none
  BEGIN_DOC
  ! Overwrites :c:data:`psi_det` and :c:data:`psi_coef` with the wave function
  ! in bilinear order
  END_DOC
  integer, intent(in)            :: isize
  integer(bit_kind), intent(out) :: psi(N_int,2,isize)
  integer                        :: i,j,k,l
  do k=1,isize
    i = psi_bilinear_matrix_rows(k)
    j = psi_bilinear_matrix_columns(k)
    psi(1:N_int,1,k) = psi_det_alpha_unique(1:N_int,i)
    psi(1:N_int,2,k) = psi_det_beta_unique(1:N_int,j)
  enddo
end


 BEGIN_PROVIDER [ integer*8, singles_alpha_csc_idx, (N_det_alpha_unique+1) ]
&BEGIN_PROVIDER [ integer*8, singles_alpha_csc_size ]
  implicit none
  BEGIN_DOC
  ! singles_alpha_csc_size : Dimension of the :c:data:`singles_alpha_csc` array
  !
  ! singles_alpha_csc_idx  : Index where the single excitations of determinant i start
  END_DOC
  integer                        :: i,j
  integer, allocatable           :: idx0(:), s(:)
  allocate (idx0(N_det_alpha_unique))
  do i=1, N_det_alpha_unique
    idx0(i) = i
  enddo

  !$OMP PARALLEL DEFAULT(NONE)                                       &
      !$OMP   SHARED(N_det_alpha_unique, psi_det_alpha_unique,       &
      !$OMP          idx0, N_int, singles_alpha_csc,                 &
      !$OMP          elec_alpha_num, mo_num, singles_alpha_csc_idx)      &
      !$OMP   PRIVATE(i,s,j)
  allocate (s(elec_alpha_num * (mo_num-elec_alpha_num) ))
  !$OMP DO SCHEDULE(static,64)
  do i=1, N_det_alpha_unique
    call get_all_spin_singles(                                       &
        psi_det_alpha_unique, idx0, psi_det_alpha_unique(1,i), N_int,&
        N_det_alpha_unique, s, j)
    singles_alpha_csc_idx(i+1) = int(j,8)
  enddo
  !$OMP END DO
  deallocate(s)
  !$OMP END PARALLEL
  deallocate(idx0)

  singles_alpha_csc_idx(1) = 1_8
  do i=2, N_det_alpha_unique+1
    singles_alpha_csc_idx(i) = singles_alpha_csc_idx(i) + singles_alpha_csc_idx(i-1)
  enddo
  singles_alpha_csc_size = singles_alpha_csc_idx(N_det_alpha_unique+1)
END_PROVIDER


BEGIN_PROVIDER [ integer, singles_alpha_csc, (singles_alpha_csc_size) ]
  implicit none
  BEGIN_DOC
  ! Indices of all single excitations
  END_DOC
  integer                        :: i, k
  integer, allocatable           :: idx0(:)
  allocate (idx0(N_det_alpha_unique))
  do i=1, N_det_alpha_unique
    idx0(i) = i
  enddo

  !$OMP PARALLEL DO DEFAULT(NONE)                                    &
      !$OMP   SHARED(N_det_alpha_unique, psi_det_alpha_unique,       &
      !$OMP          idx0, N_int, singles_alpha_csc, singles_alpha_csc_idx)&
      !$OMP   PRIVATE(i,k) SCHEDULE(static,1)
  do i=1, N_det_alpha_unique
    call get_all_spin_singles(                                       &
        psi_det_alpha_unique, idx0, psi_det_alpha_unique(1,i), N_int,&
        N_det_alpha_unique, singles_alpha_csc(singles_alpha_csc_idx(i)),&
        k)
  enddo
  !$OMP END PARALLEL DO
  deallocate(idx0)

END_PROVIDER




 BEGIN_PROVIDER [ integer*8, singles_beta_csc_idx, (N_det_beta_unique+1) ]
&BEGIN_PROVIDER [ integer*8, singles_beta_csc_size ]
  implicit none
  BEGIN_DOC
  ! singles_beta_csc_size : Dimension of the :c:data:`singles_beta_csc` array
  !
  ! singles_beta_csc_idx  : Index where the single excitations of determinant i start
  END_DOC
  integer                        :: i,j
  integer, allocatable           :: idx0(:), s(:)
  allocate (idx0(N_det_beta_unique))
  do i=1, N_det_beta_unique
    idx0(i) = i
  enddo

  !$OMP PARALLEL DEFAULT(NONE)                                       &
      !$OMP   SHARED(N_det_beta_unique, psi_det_beta_unique,       &
      !$OMP          idx0, N_int, singles_beta_csc,                 &
      !$OMP          elec_beta_num, mo_num, singles_beta_csc_idx)      &
      !$OMP   PRIVATE(i,s,j)
  allocate (s(elec_beta_num*(mo_num-elec_beta_num)))
  !$OMP DO SCHEDULE(static,1)
  do i=1, N_det_beta_unique
    call get_all_spin_singles(                                       &
        psi_det_beta_unique, idx0, psi_det_beta_unique(1,i), N_int,&
        N_det_beta_unique, s, j)
    singles_beta_csc_idx(i+1) = int(j,8)
  enddo
  !$OMP END DO
  deallocate(s)
  !$OMP END PARALLEL
  deallocate(idx0)

  singles_beta_csc_idx(1) = 1_8
  do i=2, N_det_beta_unique+1
    singles_beta_csc_idx(i) = singles_beta_csc_idx(i) + singles_beta_csc_idx(i-1)
  enddo
  singles_beta_csc_size = singles_beta_csc_idx(N_det_beta_unique+1)
END_PROVIDER


BEGIN_PROVIDER [ integer, singles_beta_csc, (singles_beta_csc_size) ]
  implicit none
  BEGIN_DOC
  ! Indices of all single excitations
  END_DOC
  integer                        :: i, k
  integer, allocatable           :: idx0(:)
  allocate (idx0(N_det_beta_unique))
  do i=1, N_det_beta_unique
    idx0(i) = i
  enddo

  !$OMP PARALLEL DO DEFAULT(NONE)                                    &
      !$OMP   SHARED(N_det_beta_unique, psi_det_beta_unique,       &
      !$OMP          idx0, N_int, singles_beta_csc, singles_beta_csc_idx)&
      !$OMP   PRIVATE(i,k) SCHEDULE(static,64)
  do i=1, N_det_beta_unique
    call get_all_spin_singles(                                       &
        psi_det_beta_unique, idx0, psi_det_beta_unique(1,i), N_int,&
        N_det_beta_unique, singles_beta_csc(singles_beta_csc_idx(i)),&
        k)
  enddo
  !$OMP END PARALLEL DO
  deallocate(idx0)

END_PROVIDER





subroutine get_all_spin_singles_and_doubles_1(buffer, idx, spindet, size_buffer, singles, doubles, n_singles, n_doubles)
  use bitmasks
  implicit none
  BEGIN_DOC
  !
  ! Returns the indices of all the single and double excitations in the list of
  ! unique $\alpha$ determinants.
  !
  ! /!\ : The buffer is transposed !
  !
  END_DOC
  integer, intent(in)            :: size_buffer, idx(size_buffer)
  integer(bit_kind), intent(in)  :: buffer(size_buffer)
  integer(bit_kind), intent(in)  :: spindet
  integer, intent(out)           :: singles(size_buffer)
  integer, intent(out)           :: doubles(size_buffer)
  integer, intent(out)           :: n_singles
  integer, intent(out)           :: n_doubles

  integer                        :: i
  include 'utils/constants.include.F'
  integer                        :: degree
  integer                        :: add_double(0:64) = (/ 0, 0, 0, 0, 1, (0, i=1,60) /)
  integer                        :: add_single(0:64) = (/ 0, 0, 1, 0, 0, (0, i=1,60) /)

  n_singles = 1
  n_doubles = 1
  do i=1,size_buffer
    degree =  popcnt(  xor( spindet, buffer(i) ) )
    doubles(n_doubles) = idx(i)
    singles(n_singles) = idx(i)
    n_doubles = n_doubles+add_double(degree)
    n_singles = n_singles+add_single(degree)
  enddo
  n_singles = n_singles-1
  n_doubles = n_doubles-1

end



subroutine get_all_spin_singles_1(buffer, idx, spindet, size_buffer, singles, n_singles)
  use bitmasks
  implicit none
  BEGIN_DOC
  !
  ! Returns the indices of all the single excitations in the list of
  ! unique $\alpha$ determinants.
  !
  END_DOC
  integer, intent(in)            :: size_buffer, idx(size_buffer)
  integer(bit_kind), intent(in)  :: buffer(size_buffer)
  integer(bit_kind), intent(in)  :: spindet
  integer, intent(out)           :: singles(size_buffer)
  integer, intent(out)           :: n_singles
  integer                        :: i
  integer(bit_kind)              :: v
  integer                        :: degree
  integer                        :: add_single(0:64) = (/ 0, 0, 1, 0, 0, (0, i=1,60) /)
  include 'utils/constants.include.F'

  n_singles = 1
  do i=1,size_buffer
    degree = popcnt(xor( spindet, buffer(i) ))
    singles(n_singles) = idx(i)
    n_singles = n_singles+add_single(degree)
  enddo
  n_singles = n_singles-1

end


subroutine get_all_spin_doubles_1(buffer, idx, spindet, size_buffer, doubles, n_doubles)
  use bitmasks
  implicit none
  BEGIN_DOC
  !
  ! Returns the indices of all the double excitations in the list of
  ! unique $\alpha$ determinants.
  !
  END_DOC
  integer, intent(in)            :: size_buffer, idx(size_buffer)
  integer(bit_kind), intent(in)  :: buffer(size_buffer)
  integer(bit_kind), intent(in)  :: spindet
  integer, intent(out)           :: doubles(size_buffer)
  integer, intent(out)           :: n_doubles
  integer                        :: i
  include 'utils/constants.include.F'
  integer                        :: degree
  integer                        :: add_double(0:64) = (/ 0, 0, 0, 0, 1, (0, i=1,60) /)

  n_doubles = 1
  do i=1,size_buffer
    degree = popcnt(xor( spindet, buffer(i) ))
    doubles(n_doubles) = idx(i)
    n_doubles = n_doubles+add_double(degree)
  enddo
  n_doubles = n_doubles-1

end



BEGIN_TEMPLATE

subroutine get_all_spin_singles_and_doubles_$N_int(buffer, idx, spindet, size_buffer, singles, doubles, n_singles, n_doubles)
  use bitmasks
  implicit none
  BEGIN_DOC
  !
  ! Returns the indices of all the single and double excitations in the list of
  ! unique $\alpha$ determinants.
  !
  ! /!\ : The buffer is transposed !
  !
  END_DOC
  integer, intent(in)            :: size_buffer, idx(size_buffer)
  integer(bit_kind), intent(in)  :: buffer($N_int,size_buffer)
  integer(bit_kind), intent(in)  :: spindet($N_int)
  integer, intent(out)           :: singles(size_buffer)
  integer, intent(out)           :: doubles(size_buffer)
  integer, intent(out)           :: n_singles
  integer, intent(out)           :: n_doubles

  integer                        :: i,k
  integer(bit_kind)              :: xorvec($N_int)
  integer                        :: degree

  n_singles = 1
  n_doubles = 1
  do i=1,size_buffer

    do k=1,$N_int
      xorvec(k) = xor( spindet(k), buffer(k,i) )
    enddo

    degree = 0

    do k=1,$N_int
        degree = degree + popcnt(xorvec(k))
    enddo

    if ( degree == 4 ) then
      doubles(n_doubles) = idx(i)
      n_doubles = n_doubles+1
    else if ( degree == 2 ) then
      singles(n_singles) = idx(i)
      n_singles = n_singles+1
    endif

  enddo
  n_singles = n_singles-1
  n_doubles = n_doubles-1

end


subroutine get_all_spin_singles_$N_int(buffer, idx, spindet, size_buffer, singles, n_singles)
  use bitmasks
  implicit none
  BEGIN_DOC
  !
  ! Returns the indices of all the single excitations in the list of
  ! unique $\alpha$ determinants.
  !
  END_DOC
  integer, intent(in)            :: size_buffer, idx(size_buffer)
  integer(bit_kind), intent(in)  :: buffer($N_int,size_buffer)
  integer(bit_kind), intent(in)  :: spindet($N_int)
  integer, intent(out)           :: singles(size_buffer)
  integer, intent(out)           :: n_singles

  integer                        :: i,k
  include 'utils/constants.include.F'
  integer(bit_kind)              :: xorvec($N_int)
  integer                        :: degree

  n_singles = 1
  do i=1,size_buffer

    do k=1,$N_int
      xorvec(k) = xor( spindet(k), buffer(k,i) )
    enddo

    degree = 0

    do k=1,$N_int
        degree = degree + popcnt(xorvec(k))
    enddo

    if ( degree /= 2 ) then
      cycle
    endif

    singles(n_singles) = idx(i)
    n_singles = n_singles+1

  enddo
  n_singles = n_singles-1

end


subroutine get_all_spin_doubles_$N_int(buffer, idx, spindet, size_buffer, doubles, n_doubles)
  use bitmasks
  implicit none
  BEGIN_DOC
  !
  ! Returns the indices of all the double excitations in the list of
  ! unique $\alpha$ determinants.
  !
  END_DOC
  integer, intent(in)            :: size_buffer, idx(size_buffer)
  integer(bit_kind), intent(in)  :: buffer($N_int,size_buffer)
  integer(bit_kind), intent(in)  :: spindet($N_int)
  integer, intent(out)           :: doubles(size_buffer)
  integer, intent(out)           :: n_doubles

  integer                        :: i,k, degree
  include 'utils/constants.include.F'
  integer(bit_kind)              :: xorvec($N_int)

  n_doubles = 1
  do i=1,size_buffer

    do k=1,$N_int
      xorvec(k) = xor( spindet(k), buffer(k,i) )
    enddo

    degree = 0

    do k=1,$N_int
        degree = degree + popcnt(xorvec(k))
    enddo

    if ( degree /= 4 ) then
      cycle
    endif

    doubles(n_doubles) = idx(i)
    n_doubles = n_doubles+1

  enddo

  n_doubles = n_doubles-1

end

SUBST [ N_int ]
2;;
3;;
4;;
N_int;;

END_TEMPLATE


subroutine wf_of_psi_bilinear_matrix(truncate)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Generate a wave function containing all possible products
  ! of $\alpha$ and $\beta$ determinants
  END_DOC
  logical, intent(in)            :: truncate
  integer                        :: i,j,k
  integer(bit_kind)              :: tmp_det(N_int,2)
  integer                        :: idx
  integer, external              :: get_index_in_psi_det_sorted_bit
  double precision               :: norm(N_states)

  do k=1,N_det
    i = psi_bilinear_matrix_rows(k)
    j = psi_bilinear_matrix_columns(k)
    psi_det(1:N_int,1,k) = psi_det_alpha_unique(1:N_int,i)
    psi_det(1:N_int,2,k) = psi_det_beta_unique (1:N_int,j)
  enddo
  psi_coef(1:N_det,1:N_states) = psi_bilinear_matrix_values(1:N_det,1:N_states)
  TOUCH psi_det psi_coef

  psi_det   = psi_det_sorted
  psi_coef  = psi_coef_sorted
  do while (sum( dabs(psi_coef(N_det,1:N_states)) ) == 0.d0)
    N_det -= 1
  enddo
  SOFT_TOUCH psi_det psi_coef N_det

end

