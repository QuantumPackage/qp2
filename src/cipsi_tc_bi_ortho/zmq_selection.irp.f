subroutine ZMQ_selection(N_in, pt2_data)
  use f77_zmq
  use selection_types

  implicit none

  integer(ZMQ_PTR)               :: zmq_to_qp_run_socket , zmq_socket_pull
  integer, intent(in)            :: N_in
  type(selection_buffer)         :: b
  integer                        :: i, l, N
  integer, external              :: omp_get_thread_num
  type(pt2_type), intent(inout)  :: pt2_data

  PROVIDE psi_det psi_coef N_det qp_max_mem N_states pt2_F s2_eig N_det_generators

  N = max(N_in,1)
  N = min(N, (elec_alpha_num * (mo_num-elec_alpha_num))**2)
  if (.True.) then
    PROVIDE pt2_e0_denominator nproc
    PROVIDE psi_bilinear_matrix_columns_loc psi_det_alpha_unique psi_det_beta_unique
    PROVIDE psi_bilinear_matrix_rows psi_det_sorted_tc_order psi_bilinear_matrix_order
    PROVIDE psi_bilinear_matrix_transp_rows_loc psi_bilinear_matrix_transp_columns
    PROVIDE psi_bilinear_matrix_transp_order selection_weight pseudo_sym
    PROVIDE n_act_orb n_inact_orb n_core_orb n_virt_orb n_del_orb seniority_max
    PROVIDE excitation_beta_max  excitation_alpha_max excitation_max

    call new_parallel_job(zmq_to_qp_run_socket,zmq_socket_pull,'selection')

    integer, external              :: zmq_put_psi
    integer, external              :: zmq_put_N_det_generators
    integer, external              :: zmq_put_N_det_selectors
    integer, external              :: zmq_put_dvector

    if (zmq_put_psi(zmq_to_qp_run_socket,1) == -1) then
      stop 'Unable to put psi on ZMQ server'
    endif
    if (zmq_put_N_det_generators(zmq_to_qp_run_socket, 1) == -1) then
      stop 'Unable to put N_det_generators on ZMQ server'
    endif
    if (zmq_put_N_det_selectors(zmq_to_qp_run_socket, 1) == -1) then
      stop 'Unable to put N_det_selectors on ZMQ server'
    endif
    if (zmq_put_dvector(zmq_to_qp_run_socket,1,'energy',pt2_e0_denominator,size(pt2_e0_denominator)) == -1) then
      stop 'Unable to put energy on ZMQ server'
    endif
    if (zmq_put_dvector(zmq_to_qp_run_socket,1,'state_average_weight',state_average_weight,N_states) == -1) then
      stop 'Unable to put state_average_weight on ZMQ server'
    endif
    if (zmq_put_dvector(zmq_to_qp_run_socket,1,'selection_weight',selection_weight,N_states) == -1) then
      stop 'Unable to put selection_weight on ZMQ server'
    endif
    if (zmq_put_dvector(zmq_to_qp_run_socket,1,'threshold_generators',(/threshold_generators/),1) == -1) then
      stop 'Unable to put threshold_generators on ZMQ server'
    endif
    call create_selection_buffer(N, N*2, b)
  endif

  integer, external :: add_task_to_taskserver
  character(len=100000)           :: task
  integer :: j,k,ipos
  ipos=1
  task = ' '

 
 do i= 1, N_det_generators
    do j=1,pt2_F(i)
      write(task(ipos:ipos+30),'(I9,1X,I9,1X,I9,''|'')') j, i, N
      ipos += 30
      if (ipos > 100000-30) then
        if (add_task_to_taskserver(zmq_to_qp_run_socket,trim(task(1:ipos))) == -1) then
          stop 'Unable to add task to task server'
        endif
        ipos=1
      endif
    end do
  enddo
  if (ipos > 1) then
    if (add_task_to_taskserver(zmq_to_qp_run_socket,trim(task(1:ipos))) == -1) then
      stop 'Unable to add task to task server'
    endif
  endif
  N = max(N_in,1)


  ASSERT (associated(b%det))
  ASSERT (associated(b%val))

  integer, external :: zmq_set_running
  if (zmq_set_running(zmq_to_qp_run_socket) == -1) then
    print *,  irp_here, ': Failed in zmq_set_running'
  endif

  integer :: nproc_target
  if (N_det < 3*nproc) then
    nproc_target = N_det/4
  else
    nproc_target = nproc
  endif
  double precision :: mem
  mem = 8.d0 * N_det * (N_int * 2.d0 * 3.d0 +  3.d0 + 5.d0) / (1024.d0**3)
  call write_double(6,mem,'Estimated memory/thread (Gb)')
  if (qp_max_mem > 0) then
    nproc_target = max(1,int(dble(qp_max_mem)/(0.1d0 + mem)))
    nproc_target = min(nproc_target,nproc)
  endif

  f(:) = 1.d0
  if (.not.do_pt2) then
    double precision :: f(N_states), u_dot_u
    do k=1,min(N_det,N_states)
     f(k) = 1.d0 / u_dot_u(psi_selectors_coef(1,k), N_det_selectors)
    enddo
  endif

  !$OMP PARALLEL DEFAULT(shared)  SHARED(b, pt2_data)  PRIVATE(i) NUM_THREADS(nproc_target+1)
  i = omp_get_thread_num()
  if (i==0) then
    call selection_collector(zmq_socket_pull, b, N, pt2_data)
  else
    call selection_slave_inproc(i)
  endif
  !$OMP END PARALLEL

  call end_parallel_job(zmq_to_qp_run_socket, zmq_socket_pull, 'selection')
  if (N_in > 0) then
    if (s2_eig) then
      call make_selection_buffer_s2(b)
    endif
    call fill_H_apply_buffer_no_selection(b%cur,b%det,N_int,0)
  endif
  call delete_selection_buffer(b)

  do k=1,N_states
    pt2_data % pt2(k) = pt2_data % pt2(k) * f(k)
    pt2_data % variance(k) = pt2_data % variance(k) * f(k)
    do l=1,N_states
      pt2_data % overlap(k,l) = pt2_data % overlap(k,l) * dsqrt(f(k)*f(l))
      pt2_data % overlap(l,k) = pt2_data % overlap(l,k) * dsqrt(f(k)*f(l))
    enddo

    pt2_data % rpt2(k) =  &
       pt2_data % pt2(k)/(1.d0 + pt2_data % overlap(k,k))
  enddo

  pt2_overlap(:,:) = pt2_data % overlap(:,:)

  print *, 'Overlap of perturbed states:'
  do l=1,N_states
    print *, pt2_overlap(l,:)
  enddo
  print *, '-------'
  SOFT_TOUCH pt2_overlap
  call update_pt2_and_variance_weights(pt2_data, N_states)

end subroutine


subroutine selection_slave_inproc(i)
  implicit none
  integer, intent(in)            :: i

  call run_selection_slave(1,i,pt2_e0_denominator)
end

subroutine selection_collector(zmq_socket_pull, b, N, pt2_data)
  use f77_zmq
  use selection_types
  use bitmasks
  implicit none


  integer(ZMQ_PTR), intent(in)   :: zmq_socket_pull
  type(selection_buffer), intent(inout) :: b
  integer, intent(in)            :: N
  type(pt2_type), intent(inout)  :: pt2_data
  type(pt2_type)                 :: pt2_data_tmp

  double precision               :: pt2_mwen(N_states)
  double precision               :: variance_mwen(N_states)
  double precision               :: norm2_mwen(N_states)
  integer(ZMQ_PTR),external      :: new_zmq_to_qp_run_socket
  integer(ZMQ_PTR)               :: zmq_to_qp_run_socket

  integer(ZMQ_PTR), external     :: new_zmq_pull_socket

  integer :: msg_size, rc, more
  integer :: acc, i, j, robin, ntask
  double precision, pointer :: val(:)
  integer(bit_kind), pointer :: det(:,:,:)
  integer, allocatable :: task_id(:)
  type(selection_buffer) :: b2




  zmq_to_qp_run_socket = new_zmq_to_qp_run_socket()
  call create_selection_buffer(N, N*2, b2)
  integer :: k
  double precision :: rss
  double precision, external :: memory_of_int
  rss = memory_of_int(N_det_generators)
  call check_mem(rss,irp_here)
  allocate(task_id(N_det_generators))
  more = 1
  pt2_data % pt2(:)      = 0d0
  pt2_data % variance(:) = 0.d0
  pt2_data % overlap(:,:) = 0.d0
  call pt2_alloc(pt2_data_tmp,N_states)
  do while (more == 1)
    call pull_selection_results(zmq_socket_pull, pt2_data_tmp, b2%val(1), b2%det(1,1,1), b2%cur, task_id, ntask)

    call pt2_add(pt2_data, 1.d0, pt2_data_tmp)
    do i=1, b2%cur
      call add_to_selection_buffer(b, b2%det(1,1,i), b2%val(i))
      if (b2%val(i) > b%mini) exit
    end do

    do i=1, ntask
      if(task_id(i) == 0) then
          print *,  "Error in collector"
      endif
      integer, external :: zmq_delete_task
      if (zmq_delete_task(zmq_to_qp_run_socket,zmq_socket_pull,task_id(i),more) == -1) then
        stop 'Unable to delete task'
      endif
    end do
  end do
  call pt2_dealloc(pt2_data_tmp)


  call delete_selection_buffer(b2)
  call sort_selection_buffer(b)
  call end_zmq_to_qp_run_socket(zmq_to_qp_run_socket)
end subroutine

