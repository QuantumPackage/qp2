use bitmasks

BEGIN_PROVIDER [ integer, domain_alpha_det, (N_det) ]
 implicit none
 BEGIN_DOC
 ! Index in psi_det of the k-th determinant of the transposed
 ! bilinear matrix, where the determinants are grouped by alpha
 ! string and sorted by beta string inside each group.
 END_DOC
 integer :: k
 PROVIDE psi_bilinear_matrix_transp_order psi_bilinear_matrix_order
 do k=1,N_det
   domain_alpha_det(k) =                                             &
       psi_bilinear_matrix_order(psi_bilinear_matrix_transp_order(k))
 enddo
END_PROVIDER

double precision function domain_h_phi(det,istate)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! $\langle det | H | \Phi_n \rangle$ for the normalized reference state
  ! n, visiting only the determinants of the reference which can be
  ! connected to det:
  !
  ! * alpha strings doubly excited: the determinants with the same beta
  !   string as det,
  !
  ! * alpha strings at most singly excited: scan of the alpha strings
  !   of the reference, and of the determinants which contain them.
  END_DOC
  integer(bit_kind), intent(in) :: det(N_int,2)
  integer, intent(in)           :: istate
  integer :: ia, ib, k, i, da, db, m
  integer, external :: domain_beta_index
  double precision :: hij

  domain_h_phi = 0.d0

  ! Same beta string, alpha doubly excited
  ib = domain_beta_index(det(1,2))
  if (ib > 0) then
    do k=psi_bilinear_matrix_columns_loc(ib),                         &
          psi_bilinear_matrix_columns_loc(ib+1)-1
      ia = psi_bilinear_matrix_rows(k)
      da = 0
      do m=1,N_int
        da = da + popcnt(ieor(det(m,1),psi_det_alpha_unique(m,ia)))
      enddo
      if (da /= 4) cycle
      i = psi_bilinear_matrix_order(k)
      call i_H_j(det,psi_det(1,1,i),N_int,hij)
      domain_h_phi = domain_h_phi + domain_phi0(i,istate)*hij
    enddo
  endif

  ! Alpha strings at most singly excited
  do ia=1,N_det_alpha_unique
    da = 0
    do m=1,N_int
      da = da + popcnt(ieor(det(m,1),psi_det_alpha_unique(m,ia)))
    enddo
    if (da > 2) cycle
    da = shiftr(da,1)
    do k=psi_bilinear_matrix_transp_rows_loc(ia),                     &
          psi_bilinear_matrix_transp_rows_loc(ia+1)-1
      db = 0
      do m=1,N_int
        db = db + popcnt(ieor(det(m,2),                              &
            psi_det_beta_unique(m,psi_bilinear_matrix_transp_columns(k))))
      enddo
      db = shiftr(db,1)
      if (da+db > 2) cycle
      if (da+db == 0) cycle
      i = domain_alpha_det(k)
      call i_H_j(det,psi_det(1,1,i),N_int,hij)
      domain_h_phi = domain_h_phi + domain_phi0(i,istate)*hij
    enddo
  enddo
end

integer function domain_beta_index(key)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Index of a beta string in psi_det_beta_unique, 0 if it is absent.
  END_DOC
  integer(bit_kind), intent(in) :: key(N_int)
  integer(8), external :: spin_det_search_key
  integer(8) :: ref
  integer :: lo, hi, mid, i
  domain_beta_index = 0
  ref = spin_det_search_key(key,N_int)
  lo = 1
  hi = N_det_beta_unique+1
  do while (lo < hi)
    mid = shiftr(lo+hi,1)
    if (spin_det_search_key(psi_det_beta_unique(1,mid),N_int) < ref) then
      lo = mid+1
    else
      hi = mid
    endif
  enddo
  do i=lo,N_det_beta_unique
    if (spin_det_search_key(psi_det_beta_unique(1,i),N_int) /= ref) exit
    if (all(psi_det_beta_unique(:,i) == key(:))) then
      domain_beta_index = i
      return
    endif
  enddo
end
