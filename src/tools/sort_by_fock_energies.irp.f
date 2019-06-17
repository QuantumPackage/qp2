program sort_by_fock_energies
  BEGIN_DOC
  ! Program that saves the current |MOs| ordered by diagonal element of the Fock operator.
  !
  ! Warning : the Fock operator, and therefore its matrix elements, depends on the occupancy.
  END_DOC
  implicit none
  integer                        :: i,j,k
  integer, allocatable           :: iorder(:)
  double precision, allocatable  :: fock_energies_tmp(:), new_mo_coef(:,:)

  allocate(iorder(mo_num), fock_energies_tmp(mo_num),new_mo_coef(ao_num,mo_num))

  do i = 1, mo_num
    fock_energies_tmp(i) = Fock_matrix_diag_mo(i)
    print*,'fock_energies_tmp(i) = ',fock_energies_tmp(i)
    iorder(i) = i
  enddo

  print*,''
  print*,'Sorting by Fock energies'
  print*,''

  call dsort(fock_energies_tmp,iorder,mo_num)

  do i = 1, mo_num
    k = iorder(i)
    print*,'fock_energies_new(i) = ',fock_energies_tmp(i)
    do j = 1, ao_num
      new_mo_coef(j,i) = mo_coef(j,k)
    enddo
  enddo

  mo_coef = new_mo_coef
  touch mo_coef
  call save_mos
  
end
