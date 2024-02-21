
! ---

program print_j1ecoef_info

  implicit none

  my_grid_becke  = .True.
  PROVIDE tc_grid1_a tc_grid1_r
  my_n_pt_r_grid = tc_grid1_r
  my_n_pt_a_grid = tc_grid1_a
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  if(tc_integ_type .eq. "numeric") then
    my_extra_grid_becke  = .True.
    PROVIDE tc_grid2_a tc_grid2_r
    my_n_pt_r_extra_grid = tc_grid2_r
    my_n_pt_a_extra_grid = tc_grid2_a
    touch my_extra_grid_becke my_n_pt_r_extra_grid my_n_pt_a_extra_grid
  endif

  call print_j1ecoef()
  
end

! ---

subroutine print_j1ecoef()

  implicit none
  integer                       :: i, j, ij
  integer                       :: ierr
  logical                       :: exists
  character(len=10)             :: ni, nj
  double precision, allocatable :: coef_fit2(:)

  PROVIDE ao_l_char_space

  allocate(coef_fit2(ao_num*ao_num))

  if(mpi_master) then
    call ezfio_has_jastrow_j1e_coef_ao2(exists)
  endif 
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    call MPI_BCAST(coef_fit2, ao_num*ao_num, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read j1e_coef_ao2 with MPI'
    endif
  IRP_ENDIF
  if(exists) then
    if(mpi_master) then
      write(6,'(A)') '.. >>>>> [ IO READ: j1e_coef_ao2 ] <<<<< ..'
      call ezfio_get_jastrow_j1e_coef_ao2(coef_fit2)
      IRP_IF MPI
        call MPI_BCAST(coef_fit2, ao_num*ao_num, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        if (ierr /= MPI_SUCCESS) then
          stop 'Unable to read j1e_coef_ao2 with MPI'
        endif
      IRP_ENDIF
    endif
  else

    call get_j1e_coef_fit_ao2(ao_num*ao_num, coef_fit2)
    call ezfio_set_jastrow_j1e_coef_ao2(coef_fit2)

  endif


  do i = 1, ao_num
    write(ni, '(I0)') ao_l(i)+1
    do j = 1, ao_num
      write(nj, '(I0)') ao_l(j)+1
      ij = (i-1)*ao_num + j
      print *, trim(adjustl(ni)) // trim(adjustl(ao_l_char_space(i))), "  " &
             , trim(adjustl(nj)) // trim(adjustl(ao_l_char_space(j))), "  " &
             , dabs(coef_fit2(ij))
    enddo
!    print *, ' '
  enddo


  deallocate(coef_fit2)
 
  return
end

! ---


