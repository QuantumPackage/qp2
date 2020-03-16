integer function zmq_put_psi(zmq_to_qp_run_socket,worker_id)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Put the wave function on the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  character*(256)                :: msg

  integer, external              :: zmq_put_N_states
  integer, external              :: zmq_put_N_det
  integer, external              :: zmq_put_psi_det_size
  integer*8, external            :: zmq_put_psi_det
  integer*8, external            :: zmq_put_psi_coef
  integer*8, external            :: zmq_put_psi_coef_complex

  zmq_put_psi = 0
  if (zmq_put_N_states(zmq_to_qp_run_socket, worker_id) == -1) then
    zmq_put_psi = -1
    return
  endif
  if (zmq_put_N_det(zmq_to_qp_run_socket, worker_id) == -1) then
    zmq_put_psi = -1
    return
  endif
  if (zmq_put_psi_det_size(zmq_to_qp_run_socket, worker_id) == -1) then
    zmq_put_psi = -1
    return
  endif
  if (zmq_put_psi_det(zmq_to_qp_run_socket, worker_id) == -1) then
    zmq_put_psi = -1
    return
  endif
  if (is_complex) then
    if (zmq_put_psi_coef_complex(zmq_to_qp_run_socket, worker_id) == -1) then
      zmq_put_psi = -1
      return
    endif
  else
  if (zmq_put_psi_coef(zmq_to_qp_run_socket, worker_id) == -1) then
    zmq_put_psi = -1
    return
  endif
  endif
end



integer function zmq_get_psi_notouch(zmq_to_qp_run_socket, worker_id)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Get the wave function from the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id

  integer, external              :: zmq_get_N_states
  integer, external              :: zmq_get_N_det
  integer, external              :: zmq_get_psi_det_size
  integer*8, external            :: zmq_get_psi_det
  integer*8, external            :: zmq_get_psi_coef
  integer*8, external            :: zmq_get_psi_coef_complex

  zmq_get_psi_notouch = 0

  if (zmq_get_N_states(zmq_to_qp_run_socket, worker_id) == -1) then
    zmq_get_psi_notouch = -1
    return
  endif
  if (zmq_get_N_det(zmq_to_qp_run_socket, worker_id) == -1) then
    zmq_get_psi_notouch = -1
    return
  endif
  if (zmq_get_psi_det_size(zmq_to_qp_run_socket, worker_id) == -1) then
    zmq_get_psi_notouch = -1
    return
  endif

  if (size(psi_det,kind=8) /= N_int*2_8*psi_det_size*bit_kind) then
    deallocate(psi_det)
    allocate(psi_det(N_int,2,psi_det_size))
  endif

  if (is_complex) then
    if (size(psi_coef_complex,kind=8) /= psi_det_size*N_states) then
      deallocate(psi_coef_complex)
      allocate(psi_coef_complex(psi_det_size,N_states))
    endif
  else
  if (size(psi_coef,kind=8) /= psi_det_size*N_states) then
    deallocate(psi_coef)
    allocate(psi_coef(psi_det_size,N_states))
  endif
  endif

  if (zmq_get_psi_det(zmq_to_qp_run_socket, worker_id) == -1_8) then
    zmq_get_psi_notouch = -1
    return
  endif

  if (is_complex) then
    if (zmq_get_psi_coef_complex(zmq_to_qp_run_socket, worker_id) == -1_8) then
      zmq_get_psi_notouch = -1
      return
    endif
  else
  if (zmq_get_psi_coef(zmq_to_qp_run_socket, worker_id) == -1_8) then
    zmq_get_psi_notouch = -1
    return
  endif
  endif

end


integer function zmq_get_psi(zmq_to_qp_run_socket, worker_id)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Get the wave function from the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  integer, external :: zmq_get_psi_notouch
  zmq_get_psi = zmq_get_psi_notouch(zmq_to_qp_run_socket, worker_id)
  if (is_complex) then
    SOFT_TOUCH psi_det psi_coef_complex psi_det_size N_det N_states
  else
  SOFT_TOUCH psi_det psi_coef psi_det_size N_det N_states
  endif
end





integer function zmq_put_psi_bilinear(zmq_to_qp_run_socket,worker_id)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Put the wave function on the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  character*(256)                :: msg


  zmq_put_psi_bilinear = 0

  integer, external            :: zmq_put_psi
  if (zmq_put_psi(zmq_to_qp_run_socket, worker_id) == -1) then
    zmq_put_psi_bilinear = -1
    return
  endif

  integer*8, external            :: zmq_put_psi_bilinear_matrix_columns
  if (zmq_put_psi_bilinear_matrix_rows(zmq_to_qp_run_socket, worker_id) == -1) then
    zmq_put_psi_bilinear = -1
    return
  endif

  integer*8, external            :: zmq_put_psi_bilinear_matrix_rows
  if (zmq_put_psi_bilinear_matrix_columns(zmq_to_qp_run_socket, worker_id) == -1) then
    zmq_put_psi_bilinear = -1
    return
  endif

  integer*8, external            :: zmq_put_psi_bilinear_matrix_order
  if (zmq_put_psi_bilinear_matrix_order(zmq_to_qp_run_socket, worker_id) == -1) then
    zmq_put_psi_bilinear = -1
    return
  endif
  
  if (is_complex) then
    integer*8, external            :: zmq_put_psi_bilinear_matrix_values_complex
    if (zmq_put_psi_bilinear_matrix_values_complex(zmq_to_qp_run_socket, worker_id) == -1) then
      zmq_put_psi_bilinear = -1
      return
    endif
  else
  integer*8, external            :: zmq_put_psi_bilinear_matrix_values
  if (zmq_put_psi_bilinear_matrix_values(zmq_to_qp_run_socket, worker_id) == -1) then
    zmq_put_psi_bilinear = -1
    return
  endif
  endif

  integer, external            :: zmq_put_N_det_alpha_unique
  if (zmq_put_N_det_alpha_unique(zmq_to_qp_run_socket,worker_id) == -1) then
    zmq_put_psi_bilinear = -1
    return
  endif

  integer, external            :: zmq_put_N_det_beta_unique
  if (zmq_put_N_det_beta_unique(zmq_to_qp_run_socket,worker_id) == -1) then
    zmq_put_psi_bilinear = -1
    return
  endif

  integer*8, external            :: zmq_put_psi_det_alpha_unique
  if (zmq_put_psi_det_alpha_unique(zmq_to_qp_run_socket,worker_id) == -1) then
    zmq_put_psi_bilinear = -1
    return
  endif

  integer*8, external            :: zmq_put_psi_det_beta_unique
  if (zmq_put_psi_det_beta_unique(zmq_to_qp_run_socket,worker_id) == -1) then
    zmq_put_psi_bilinear = -1
    return
  endif

end


integer function zmq_get_psi_bilinear(zmq_to_qp_run_socket, worker_id)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Get the wave function from the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id

  integer, external            :: zmq_get_psi_notouch
  if (zmq_get_psi_notouch(zmq_to_qp_run_socket,1) == -1) then
    zmq_get_psi_bilinear = -1
    return
  endif

  zmq_get_psi_bilinear= 0

  if (is_complex) then
    if (size(psi_bilinear_matrix_values_complex,kind=8) /= N_det*N_states) then
      deallocate(psi_bilinear_matrix_values_complex)
      allocate(psi_bilinear_matrix_values_complex(N_det,N_states))
    endif
  else
  if (size(psi_bilinear_matrix_values,kind=8) /= N_det*N_states) then
    deallocate(psi_bilinear_matrix_values)
    allocate(psi_bilinear_matrix_values(N_det,N_states))
  endif
  endif

  if (size(psi_bilinear_matrix_rows,kind=8) /= N_det) then
    deallocate(psi_bilinear_matrix_rows)
    allocate(psi_bilinear_matrix_rows(N_det))
  endif

  if (size(psi_bilinear_matrix_columns,kind=8) /= N_det) then
    deallocate(psi_bilinear_matrix_columns)
    allocate(psi_bilinear_matrix_columns(N_det))
  endif

  if (size(psi_bilinear_matrix_order,kind=8) /= N_det) then
    deallocate(psi_bilinear_matrix_order)
    allocate(psi_bilinear_matrix_order(N_det))
  endif
  
  if (is_complex) then
    integer*8, external            :: zmq_get_psi_bilinear_matrix_values_complex
    if (zmq_get_psi_bilinear_matrix_values_complex(zmq_to_qp_run_socket, worker_id) == -1_8) then
      zmq_get_psi_bilinear = -1
      return
    endif
  else
  integer*8, external            :: zmq_get_psi_bilinear_matrix_values
  if (zmq_get_psi_bilinear_matrix_values(zmq_to_qp_run_socket, worker_id) == -1_8) then
    zmq_get_psi_bilinear = -1
    return
  endif
  endif

  integer*8, external            :: zmq_get_psi_bilinear_matrix_rows
  if (zmq_get_psi_bilinear_matrix_rows(zmq_to_qp_run_socket, worker_id) == -1_8) then
    zmq_get_psi_bilinear = -1
    return
  endif

  integer*8, external            :: zmq_get_psi_bilinear_matrix_columns
  if (zmq_get_psi_bilinear_matrix_columns(zmq_to_qp_run_socket, worker_id) == -1_8) then
    zmq_get_psi_bilinear = -1
    return
  endif

  integer*8, external            :: zmq_get_psi_bilinear_matrix_order
  if (zmq_get_psi_bilinear_matrix_order(zmq_to_qp_run_socket, worker_id) == -1_8) then
    zmq_get_psi_bilinear = -1
    return
  endif


  integer, external            :: zmq_get_N_det_alpha_unique
  if (zmq_get_N_det_alpha_unique(zmq_to_qp_run_socket,worker_id) == -1) then
    zmq_get_psi_bilinear = -1
    return
  endif

  integer, external            :: zmq_get_N_det_beta_unique
  if (zmq_get_N_det_beta_unique(zmq_to_qp_run_socket,worker_id) == -1) then
    zmq_get_psi_bilinear = -1
    return
  endif

  integer*8, external            :: zmq_get_psi_det_alpha_unique
  if (zmq_get_psi_det_alpha_unique(zmq_to_qp_run_socket,worker_id) == -1) then
    zmq_get_psi_bilinear = -1
    return
  endif

  integer*8, external            :: zmq_get_psi_det_beta_unique
  if (zmq_get_psi_det_beta_unique(zmq_to_qp_run_socket,worker_id) == -1) then
    zmq_get_psi_bilinear = -1
    return
  endif

  if (is_complex) then
    SOFT_TOUCH psi_bilinear_matrix_values_complex psi_bilinear_matrix_rows psi_bilinear_matrix_columns psi_bilinear_matrix_order psi_det psi_coef_complex psi_det_size N_det N_states psi_det_beta_unique psi_det_alpha_unique N_det_beta_unique N_det_alpha_unique
  else
  SOFT_TOUCH psi_bilinear_matrix_values psi_bilinear_matrix_rows psi_bilinear_matrix_columns psi_bilinear_matrix_order psi_det psi_coef psi_det_size N_det N_states psi_det_beta_unique psi_det_alpha_unique N_det_beta_unique N_det_alpha_unique
  endif

end







BEGIN_TEMPLATE

integer function zmq_put_$X(zmq_to_qp_run_socket,worker_id)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Put $X on the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  integer                        :: rc
  character*(256)                :: msg

  zmq_put_$X = 0

  write(msg,'(A,1X,I8,1X,A200)') 'put_data '//trim(zmq_state), worker_id, '$X'
  rc = f77_zmq_send(zmq_to_qp_run_socket,trim(msg),len(trim(msg)),ZMQ_SNDMORE)
  if (rc /= len(trim(msg))) then
    zmq_put_$X = -1
    return
  endif

  rc = f77_zmq_send(zmq_to_qp_run_socket,$X,4,0)
  if (rc /= 4) then
    zmq_put_$X = -1
    return
  endif

  rc = f77_zmq_recv(zmq_to_qp_run_socket,msg,len(msg),0)
  if (msg(1:rc) /= 'put_data_reply ok') then
    zmq_put_$X = -1
    return
  endif

end

integer function zmq_get_$X(zmq_to_qp_run_socket, worker_id)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Get $X from the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  integer                        :: rc
  character*(256)                :: msg

  PROVIDE zmq_state
  zmq_get_$X = 0
  if (mpi_master) then
    write(msg,'(A,1X,I8,1X,A200)') 'get_data '//trim(zmq_state), worker_id, '$X'
    rc = f77_zmq_send(zmq_to_qp_run_socket,trim(msg),len(trim(msg)),0)
    if (rc /= len(trim(msg))) then
      zmq_get_$X = -1
      go to 10
    endif

    rc = f77_zmq_recv(zmq_to_qp_run_socket,msg,len(msg),0)
    if (msg(1:14) /= 'get_data_reply') then
      zmq_get_$X = -1
      go to 10
    endif

    rc = f77_zmq_recv(zmq_to_qp_run_socket,$X,4,0)
    if (rc /= 4) then
      zmq_get_$X = -1
      go to 10
    endif

  endif

 10 continue
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr

    call MPI_BCAST (zmq_get_$X, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to broadcast zmq_get_psi_det'
    endif
    call MPI_BCAST ($X, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to broadcast zmq_get_psi_det'
    endif
  IRP_ENDIF

end

SUBST [ X ]

N_states ;;
N_det ;;
N_det_alpha_unique ;;
N_det_beta_unique ;;
psi_det_size ;;

END_TEMPLATE


BEGIN_TEMPLATE

integer*8 function zmq_put_$X(zmq_to_qp_run_socket,worker_id)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Put $X on the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  integer*8                      :: rc8
  character*(256)                :: msg

  integer*8 :: zmq_put_i8matrix
  integer   :: ni, nj

  if (size($X,kind=8) <= 8388608_8) then
    ni = size($X,kind=4)
    nj = 1
  else
    ni = 8388608_8
    nj = int(size($X,kind=8)/8388608_8,4) + 1
  endif
  rc8 = zmq_put_i8matrix(zmq_to_qp_run_socket, 1, '$X', $X, ni, nj, size($X,kind=8))
  zmq_put_$X = rc8
end

integer*8 function zmq_get_$X(zmq_to_qp_run_socket,worker_id)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Get $X on the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  integer*8                      :: rc8
  character*(256)                :: msg

  integer*8 :: zmq_get_i8matrix
  integer   :: ni, nj

  if (size($X,kind=8) <= 8388608_8) then
    ni = size($X,kind=4)
    nj = 1
  else
    ni = 8388608
    nj = int(size($X,kind=8)/8388608_8,4) + 1
  endif
  rc8 = zmq_get_i8matrix(zmq_to_qp_run_socket, 1, '$X', $X, ni, nj, size($X,kind=8))
  zmq_get_$X = rc8
end

SUBST [ X ]

psi_det ;;
psi_det_alpha_unique ;;
psi_det_beta_unique ;;

END_TEMPLATE

BEGIN_TEMPLATE

integer*8 function zmq_put_$X(zmq_to_qp_run_socket,worker_id)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Put $X on the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  integer*8                      :: rc8
  character*(256)                :: msg

  integer*8 :: zmq_put_imatrix
  integer   :: ni, nj

  if (size($X,kind=8) <= 8388608_8) then
    ni = size($X,kind=4)
    nj = 1
  else
    ni = 8388608_8
    nj = int(size($X,kind=8)/8388608_8,4) + 1
  endif
  rc8 = zmq_put_imatrix(zmq_to_qp_run_socket, 1, '$X', $X, ni, nj, size($X,kind=8))
  zmq_put_$X = rc8
end

integer*8 function zmq_get_$X(zmq_to_qp_run_socket,worker_id)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Get $X on the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  integer*8                      :: rc8
  character*(256)                :: msg

  integer*8 :: zmq_get_imatrix
  integer   :: ni, nj

  if (size($X,kind=8) <= 8388608_8) then
    ni = size($X,kind=4)
    nj = 1
  else
    ni = 8388608
    nj = int(size($X,kind=8)/8388608_8,4) + 1
  endif
  rc8 = zmq_get_imatrix(zmq_to_qp_run_socket, 1, '$X', $X, ni, nj, size($X,kind=8))
  zmq_get_$X = rc8
end

SUBST [ X ]

psi_bilinear_matrix_rows ;;
psi_bilinear_matrix_columns ;;
psi_bilinear_matrix_order ;;

END_TEMPLATE


BEGIN_TEMPLATE

integer*8 function zmq_put_$X(zmq_to_qp_run_socket,worker_id)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Put $X on the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  integer*8                      :: rc8
  character*(256)                :: msg

  zmq_put_$X = 0

  integer*8 :: zmq_put_dmatrix
  integer   :: ni, nj

  if (size($X,kind=8) <= 8388608_8) then
    ni = size($X,kind=4)
    nj = 1
  else
    ni = 8388608
    nj = int(size($X,kind=8)/8388608_8,4) + 1
  endif
  rc8 = zmq_put_dmatrix(zmq_to_qp_run_socket, 1, '$X', $X, ni, nj, size($X,kind=8) )
  zmq_put_$X = rc8
end

integer*8 function zmq_get_$X(zmq_to_qp_run_socket,worker_id)
  use f77_zmq
  implicit none
  BEGIN_DOC
! get $X on the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  integer*8                      :: rc8
  character*(256)                :: msg

  zmq_get_$X = 0_8

  integer*8 :: zmq_get_dmatrix
  integer   :: ni, nj

  if (size($X,kind=8) <= 8388608_8) then
    ni = size($X,kind=4)
    nj = 1
  else
    ni = 8388608
    nj = int(size($X,kind=8)/8388608_8,4) + 1
  endif
  rc8 = zmq_get_dmatrix(zmq_to_qp_run_socket, 1, '$X', $X, ni, nj, size($X,kind=8) )
  zmq_get_$X = rc8
end

SUBST [ X ]

psi_coef ;;
psi_bilinear_matrix_values ;;

END_TEMPLATE

BEGIN_TEMPLATE

integer*8 function zmq_put_$X(zmq_to_qp_run_socket,worker_id)
  use f77_zmq
  implicit none
  BEGIN_DOC
! Put $X on the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  integer*8                      :: rc8
  character*(256)                :: msg

  zmq_put_$X = 0

  integer*8 :: zmq_put_cdmatrix
  integer   :: ni, nj

  if (size($X,kind=8) <= 8388608_8) then
    ni = size($X,kind=4)
    nj = 1
  else
    ni = 8388608
    nj = int(size($X,kind=8)/8388608_8,4) + 1
  endif
  rc8 = zmq_put_cdmatrix(zmq_to_qp_run_socket, 1, '$X', $X, ni, nj, size($X,kind=8) )
  zmq_put_$X = rc8
end

integer*8 function zmq_get_$X(zmq_to_qp_run_socket,worker_id)
  use f77_zmq
  implicit none
  BEGIN_DOC
! get $X on the qp_run scheduler
  END_DOC
  integer(ZMQ_PTR), intent(in)   :: zmq_to_qp_run_socket
  integer, intent(in)            :: worker_id
  integer*8                      :: rc8
  character*(256)                :: msg

  zmq_get_$X = 0_8

  integer*8 :: zmq_get_cdmatrix
  integer   :: ni, nj

  if (size($X,kind=8) <= 8388608_8) then
    ni = size($X,kind=4)
    nj = 1
  else
    ni = 8388608
    nj = int(size($X,kind=8)/8388608_8,4) + 1
  endif
  rc8 = zmq_get_cdmatrix(zmq_to_qp_run_socket, 1, '$X', $X, ni, nj, size($X,kind=8) )
  zmq_get_$X = rc8
end

SUBST [ X ]

psi_coef_complex ;;
psi_bilinear_matrix_values_complex ;;

END_TEMPLATE


!---------------------------------------------------------------------------


