
! ---

double precision function ao_two_e_integral(i, j, k, l)

  BEGIN_DOC
  !  integral of the AO basis <ik|jl> or (ij|kl)
  !     i(r1) j(r1) 1/r12 k(r2) l(r2)
  END_DOC

  implicit none
  include 'utils/constants.include.F'

  integer, intent(in)            :: i, j, k, l

  integer                        :: p, q, r, s
  integer                        :: num_i,num_j,num_k,num_l,dim1,I_power(3),J_power(3),K_power(3),L_power(3)
  integer                        :: iorder_p(3), iorder_q(3)
  double precision               :: I_center(3), J_center(3), K_center(3), L_center(3)
  double precision               :: integral
  double precision               :: P_new(0:max_dim,3),P_center(3),fact_p,pp
  double precision               :: Q_new(0:max_dim,3),Q_center(3),fact_q,qq

  double precision               :: ao_two_e_integral_schwartz_accel

  double precision               :: ao_two_e_integral_cosgtos


  if(use_cosgtos) then
    !print *, ' use_cosgtos for ao_two_e_integral ?', use_cosgtos

    ao_two_e_integral = ao_two_e_integral_cosgtos(i, j, k, l)

  else

    if (ao_prim_num(i) * ao_prim_num(j) * ao_prim_num(k) * ao_prim_num(l) > 1024 ) then

       ao_two_e_integral = ao_two_e_integral_schwartz_accel(i,j,k,l)

    else

      dim1 = n_pt_max_integrals

      num_i = ao_nucl(i)
      num_j = ao_nucl(j)
      num_k = ao_nucl(k)
      num_l = ao_nucl(l)
      ao_two_e_integral = 0.d0

      if (num_i /= num_j .or. num_k /= num_l .or. num_j /= num_k)then
        do p = 1, 3
          I_power(p) = ao_power(i,p)
          J_power(p) = ao_power(j,p)
          K_power(p) = ao_power(k,p)
          L_power(p) = ao_power(l,p)
          I_center(p) = nucl_coord(num_i,p)
          J_center(p) = nucl_coord(num_j,p)
          K_center(p) = nucl_coord(num_k,p)
          L_center(p) = nucl_coord(num_l,p)
        enddo

        double precision               :: coef1, coef2, coef3, coef4
        double precision               :: p_inv,q_inv
        double precision               :: general_primitive_integral

        do p = 1, ao_prim_num(i)
          coef1 = ao_coef_normalized_ordered_transp(p,i)
          do q = 1, ao_prim_num(j)
            coef2 = coef1*ao_coef_normalized_ordered_transp(q,j)
            call give_explicit_poly_and_gaussian(P_new,P_center,pp,fact_p,iorder_p,&
                ao_expo_ordered_transp(p,i),ao_expo_ordered_transp(q,j),                 &
                I_power,J_power,I_center,J_center,dim1)
            p_inv = 1.d0/pp
            do r = 1, ao_prim_num(k)
              coef3 = coef2*ao_coef_normalized_ordered_transp(r,k)
              do s = 1, ao_prim_num(l)
                coef4 = coef3*ao_coef_normalized_ordered_transp(s,l)
                call give_explicit_poly_and_gaussian(Q_new,Q_center,qq,fact_q,iorder_q,&
                    ao_expo_ordered_transp(r,k),ao_expo_ordered_transp(s,l),             &
                    K_power,L_power,K_center,L_center,dim1)
                q_inv = 1.d0/qq
                integral = general_primitive_integral(dim1,              &
                    P_new,P_center,fact_p,pp,p_inv,iorder_p,             &
                    Q_new,Q_center,fact_q,qq,q_inv,iorder_q)
                ao_two_e_integral = ao_two_e_integral +  coef4 * integral
              enddo ! s
            enddo  ! r
          enddo   ! q
        enddo    ! p

      else

        do p = 1, 3
          I_power(p) = ao_power(i,p)
          J_power(p) = ao_power(j,p)
          K_power(p) = ao_power(k,p)
          L_power(p) = ao_power(l,p)
        enddo
        double  precision              :: ERI

        do p = 1, ao_prim_num(i)
          coef1 = ao_coef_normalized_ordered_transp(p,i)
          do q = 1, ao_prim_num(j)
            coef2 = coef1*ao_coef_normalized_ordered_transp(q,j)
            do r = 1, ao_prim_num(k)
              coef3 = coef2*ao_coef_normalized_ordered_transp(r,k)
              do s = 1, ao_prim_num(l)
                coef4 = coef3*ao_coef_normalized_ordered_transp(s,l)
                integral = ERI(                                          &
                    ao_expo_ordered_transp(p,i),ao_expo_ordered_transp(q,j),ao_expo_ordered_transp(r,k),ao_expo_ordered_transp(s,l),&
                    I_power(1),J_power(1),K_power(1),L_power(1),         &
                    I_power(2),J_power(2),K_power(2),L_power(2),         &
                    I_power(3),J_power(3),K_power(3),L_power(3))
                ao_two_e_integral = ao_two_e_integral + coef4 * integral
              enddo ! s
            enddo  ! r
          enddo   ! q
        enddo    ! p

      endif

    endif

  endif

end

! ---

double precision function ao_two_e_integral_schwartz_accel(i,j,k,l)
  implicit none
  BEGIN_DOC
  !  integral of the AO basis <ik|jl> or (ij|kl)
  !     i(r1) j(r1) 1/r12 k(r2) l(r2)
  END_DOC
  integer,intent(in)             :: i,j,k,l
  integer                        :: p,q,r,s
  double precision               :: I_center(3),J_center(3),K_center(3),L_center(3)
  integer                        :: num_i,num_j,num_k,num_l,dim1,I_power(3),J_power(3),K_power(3),L_power(3)
  double precision               :: integral
  include 'utils/constants.include.F'
  double precision               :: P_new(0:max_dim,3),P_center(3),fact_p,pp
  double precision               :: Q_new(0:max_dim,3),Q_center(3),fact_q,qq
  integer                        :: iorder_p(3), iorder_q(3)
  double precision, allocatable  :: schwartz_kl(:,:)
  double precision               :: schwartz_ij

  dim1 = n_pt_max_integrals

  num_i = ao_nucl(i)
  num_j = ao_nucl(j)
  num_k = ao_nucl(k)
  num_l = ao_nucl(l)
  ao_two_e_integral_schwartz_accel = 0.d0
  double precision               :: thr
  thr = ao_integrals_threshold*ao_integrals_threshold

  allocate(schwartz_kl(0:ao_prim_num(l),0:ao_prim_num(k)))


  if (num_i /= num_j .or. num_k /= num_l .or. num_j /= num_k)then
    do p = 1, 3
      I_power(p) = ao_power(i,p)
      J_power(p) = ao_power(j,p)
      K_power(p) = ao_power(k,p)
      L_power(p) = ao_power(l,p)
      I_center(p) = nucl_coord(num_i,p)
      J_center(p) = nucl_coord(num_j,p)
      K_center(p) = nucl_coord(num_k,p)
      L_center(p) = nucl_coord(num_l,p)
    enddo

    schwartz_kl(0,0) = 0.d0
    do r = 1, ao_prim_num(k)
      coef1 = ao_coef_normalized_ordered_transp(r,k)*ao_coef_normalized_ordered_transp(r,k)
      schwartz_kl(0,r) = 0.d0
      do s = 1, ao_prim_num(l)
        coef2 = coef1 * ao_coef_normalized_ordered_transp(s,l) * ao_coef_normalized_ordered_transp(s,l)
        call give_explicit_poly_and_gaussian(Q_new,Q_center,qq,fact_q,iorder_q,&
            ao_expo_ordered_transp(r,k),ao_expo_ordered_transp(s,l),                 &
            K_power,L_power,K_center,L_center,dim1)
        q_inv = 1.d0/qq
        schwartz_kl(s,r) = general_primitive_integral(dim1,          &
            Q_new,Q_center,fact_q,qq,q_inv,iorder_q,                 &
            Q_new,Q_center,fact_q,qq,q_inv,iorder_q)                 &
            * coef2
        schwartz_kl(0,r) = max(schwartz_kl(0,r),schwartz_kl(s,r))
      enddo
      schwartz_kl(0,0) = max(schwartz_kl(0,r),schwartz_kl(0,0))
    enddo

    do p = 1, ao_prim_num(i)
      double precision               :: coef1
      coef1 = ao_coef_normalized_ordered_transp(p,i)
      do q = 1, ao_prim_num(j)
        double precision               :: coef2
        coef2 = coef1*ao_coef_normalized_ordered_transp(q,j)
        double precision               :: p_inv,q_inv
        call give_explicit_poly_and_gaussian(P_new,P_center,pp,fact_p,iorder_p,&
            ao_expo_ordered_transp(p,i),ao_expo_ordered_transp(q,j),                 &
            I_power,J_power,I_center,J_center,dim1)
        p_inv = 1.d0/pp
        schwartz_ij = general_primitive_integral(dim1,               &
            P_new,P_center,fact_p,pp,p_inv,iorder_p,                 &
            P_new,P_center,fact_p,pp,p_inv,iorder_p) *               &
            coef2*coef2
        if (schwartz_kl(0,0)*schwartz_ij < thr) then
           cycle
        endif
        do r = 1, ao_prim_num(k)
          if (schwartz_kl(0,r)*schwartz_ij < thr) then
             cycle
          endif
          double precision               :: coef3
          coef3 = coef2*ao_coef_normalized_ordered_transp(r,k)
          do s = 1, ao_prim_num(l)
            double precision               :: coef4
            if (schwartz_kl(s,r)*schwartz_ij < thr) then
               cycle
            endif
            coef4 = coef3*ao_coef_normalized_ordered_transp(s,l)
            double precision               :: general_primitive_integral
            call give_explicit_poly_and_gaussian(Q_new,Q_center,qq,fact_q,iorder_q,&
                ao_expo_ordered_transp(r,k),ao_expo_ordered_transp(s,l),             &
                K_power,L_power,K_center,L_center,dim1)
            q_inv = 1.d0/qq
            integral = general_primitive_integral(dim1,              &
                P_new,P_center,fact_p,pp,p_inv,iorder_p,             &
                Q_new,Q_center,fact_q,qq,q_inv,iorder_q)
            ao_two_e_integral_schwartz_accel = ao_two_e_integral_schwartz_accel + coef4 * integral
          enddo ! s
        enddo  ! r
      enddo   ! q
    enddo    ! p

  else

    do p = 1, 3
      I_power(p) = ao_power(i,p)
      J_power(p) = ao_power(j,p)
      K_power(p) = ao_power(k,p)
      L_power(p) = ao_power(l,p)
    enddo
    double  precision              :: ERI

    schwartz_kl(0,0) = 0.d0
    do r = 1, ao_prim_num(k)
      coef1 = ao_coef_normalized_ordered_transp(r,k)*ao_coef_normalized_ordered_transp(r,k)
      schwartz_kl(0,r) = 0.d0
      do s = 1, ao_prim_num(l)
        coef2 = coef1*ao_coef_normalized_ordered_transp(s,l)*ao_coef_normalized_ordered_transp(s,l)
        schwartz_kl(s,r) = ERI(                                      &
            ao_expo_ordered_transp(r,k),ao_expo_ordered_transp(s,l),ao_expo_ordered_transp(r,k),ao_expo_ordered_transp(s,l),&
            K_power(1),L_power(1),K_power(1),L_power(1),             &
            K_power(2),L_power(2),K_power(2),L_power(2),             &
            K_power(3),L_power(3),K_power(3),L_power(3)) * &
            coef2
        schwartz_kl(0,r) = max(schwartz_kl(0,r),schwartz_kl(s,r))
      enddo
      schwartz_kl(0,0) = max(schwartz_kl(0,r),schwartz_kl(0,0))
    enddo

    do p = 1, ao_prim_num(i)
      coef1 = ao_coef_normalized_ordered_transp(p,i)
      do q = 1, ao_prim_num(j)
        coef2 = coef1*ao_coef_normalized_ordered_transp(q,j)
        schwartz_ij = ERI(                                          &
                ao_expo_ordered_transp(p,i),ao_expo_ordered_transp(q,j),ao_expo_ordered_transp(p,i),ao_expo_ordered_transp(q,j),&
                I_power(1),J_power(1),I_power(1),J_power(1),         &
                I_power(2),J_power(2),I_power(2),J_power(2),         &
                I_power(3),J_power(3),I_power(3),J_power(3))*coef2*coef2
        if (schwartz_kl(0,0)*schwartz_ij < thr) then
           cycle
        endif
        do r = 1, ao_prim_num(k)
          if (schwartz_kl(0,r)*schwartz_ij < thr) then
             cycle
          endif
          coef3 = coef2*ao_coef_normalized_ordered_transp(r,k)
          do s = 1, ao_prim_num(l)
            if (schwartz_kl(s,r)*schwartz_ij < thr) then
               cycle
            endif
            coef4 = coef3*ao_coef_normalized_ordered_transp(s,l)
            integral = ERI(                                          &
                ao_expo_ordered_transp(p,i),ao_expo_ordered_transp(q,j),ao_expo_ordered_transp(r,k),ao_expo_ordered_transp(s,l),&
                I_power(1),J_power(1),K_power(1),L_power(1),         &
                I_power(2),J_power(2),K_power(2),L_power(2),         &
                I_power(3),J_power(3),K_power(3),L_power(3))
            ao_two_e_integral_schwartz_accel = ao_two_e_integral_schwartz_accel +  coef4 * integral
          enddo ! s
        enddo  ! r
      enddo   ! q
    enddo    ! p

  endif
  deallocate (schwartz_kl)

end


integer function ao_l4(i,j,k,l)
  implicit none
  BEGIN_DOC
! Computes the product of l values of i,j,k,and l
  END_DOC
  integer, intent(in) :: i,j,k,l
  ao_l4 = ao_l(i)*ao_l(j)*ao_l(k)*ao_l(l)
end



subroutine compute_ao_two_e_integrals(j,k,l,sze,buffer_value)
  implicit none
  use map_module

  BEGIN_DOC
  ! Compute AO 1/r12 integrals for all i and fixed j,k,l
  END_DOC

  include 'utils/constants.include.F'
  integer, intent(in)            :: j,k,l,sze
  real(integral_kind), intent(out) :: buffer_value(sze)
  double precision               :: ao_two_e_integral

  integer                        :: i
  logical, external              :: ao_one_e_integral_zero
  logical, external              :: ao_two_e_integral_zero


  if (ao_one_e_integral_zero(j,l)) then
    buffer_value = 0._integral_kind
    return
  endif

  do i = 1, ao_num
    if (ao_two_e_integral_zero(i,j,k,l)) then
      buffer_value(i) = 0._integral_kind
      cycle
    endif
    !DIR$ FORCEINLINE
    buffer_value(i) = ao_two_e_integral(i,k,j,l)
  enddo

end

BEGIN_PROVIDER [ logical, ao_two_e_integrals_in_map ]
  implicit none
  use f77_zmq
  use map_module
  BEGIN_DOC
  !  Map of Atomic integrals
  !     i(r1) j(r2) 1/r12 k(r1) l(r2)
  END_DOC

  integer                        :: i,j,k,l
  double precision               :: ao_two_e_integral,cpu_1,cpu_2, wall_1, wall_2
  double precision               :: integral, wall_0
  include 'utils/constants.include.F'

  ! For integrals file
  integer(key_kind),allocatable  :: buffer_i(:)
  integer,parameter              :: size_buffer = 1024*64
  real(integral_kind),allocatable :: buffer_value(:)

  integer                        :: n_integrals, rc
  integer                        :: kk, m, j1, i1, lmax
  character*(64)                 :: fmt

  double precision               :: map_mb
  PROVIDE read_ao_two_e_integrals io_ao_two_e_integrals
  if (read_ao_two_e_integrals) then
    print*,'Reading the AO integrals'
    call map_load_from_disk(trim(ezfio_filename)//'/work/ao_ints',ao_integrals_map)
    print*, 'AO integrals provided'
    ao_two_e_integrals_in_map = .True.
  else

    print*, 'Providing the AO integrals'
    call wall_time(wall_0)
    call wall_time(wall_1)
    call cpu_time(cpu_1)

    if (.True.) then
      ! Avoid openMP
      integral = ao_two_e_integral(1,1,1,1)
    endif

    integer(ZMQ_PTR) :: zmq_to_qp_run_socket, zmq_socket_pull
    call new_parallel_job(zmq_to_qp_run_socket,zmq_socket_pull,'ao_integrals')

    character(len=:), allocatable :: task
    allocate(character(len=ao_num*12) :: task)
    write(fmt,*) '(', ao_num, '(I5,X,I5,''|''))'
    do l=1,ao_num
      write(task,fmt) (i,l, i=1,l)
      integer, external :: add_task_to_taskserver
      if (add_task_to_taskserver(zmq_to_qp_run_socket,trim(task)) == -1) then
        stop 'Unable to add task to server'
      endif
    enddo
    deallocate(task)

    integer, external :: zmq_set_running
    if (zmq_set_running(zmq_to_qp_run_socket) == -1) then
      print *,  irp_here, ': Failed in zmq_set_running'
    endif

    PROVIDE nproc
    !$OMP PARALLEL DEFAULT(shared) private(i) num_threads(nproc+1)
        i = omp_get_thread_num()
        if (i==0) then
          call ao_two_e_integrals_in_map_collector(zmq_socket_pull)
        else
          call ao_two_e_integrals_in_map_slave_inproc(i)
        endif
    !$OMP END PARALLEL

    call end_parallel_job(zmq_to_qp_run_socket, zmq_socket_pull, 'ao_integrals')


    print*, 'Sorting the map'
    call map_sort(ao_integrals_map)
    call cpu_time(cpu_2)
    call wall_time(wall_2)
    integer(map_size_kind)         :: get_ao_map_size, ao_map_size
    ao_map_size = get_ao_map_size()

    print*, 'AO integrals provided:'
    print*, ' Size of AO map :         ', map_mb(ao_integrals_map) ,'MB'
    print*, ' Number of AO integrals :', ao_map_size
    print*, ' cpu  time :',cpu_2 - cpu_1, 's'
    print*, ' wall time :',wall_2 - wall_1, 's  ( x ', (cpu_2-cpu_1)/(wall_2-wall_1+tiny(1.d0)), ' )'

    ao_two_e_integrals_in_map = .True.

    if (write_ao_two_e_integrals.and.mpi_master) then
      call ezfio_set_work_empty(.False.)
      call map_save_to_disk(trim(ezfio_filename)//'/work/ao_ints',ao_integrals_map)
      call ezfio_set_ao_two_e_ints_io_ao_two_e_integrals('Read')
    endif

  endif

END_PROVIDER

! ---

BEGIN_PROVIDER [ double precision, ao_two_e_integral_schwartz, (ao_num, ao_num) ]

  BEGIN_DOC
  !  Needed to compute Schwartz inequalities
  END_DOC

  implicit none
  integer          :: i, k
  double precision :: ao_two_e_integral,cpu_1,cpu_2, wall_1, wall_2

  ao_two_e_integral_schwartz(1,1) = ao_two_e_integral(1,1,1,1)
  !$OMP PARALLEL DO PRIVATE(i,k)                                     &
      !$OMP DEFAULT(NONE)                                            &
      !$OMP SHARED (ao_num,ao_two_e_integral_schwartz)              &
      !$OMP SCHEDULE(dynamic)
  do i=1,ao_num
    do k=1,i
      ao_two_e_integral_schwartz(i,k) = dsqrt(ao_two_e_integral(i,i,k,k))
      ao_two_e_integral_schwartz(k,i) = ao_two_e_integral_schwartz(i,k)
    enddo
  enddo
  !$OMP END PARALLEL DO

END_PROVIDER

! ---

double precision function general_primitive_integral(dim,            &
      P_new,P_center,fact_p,p,p_inv,iorder_p,                        &
      Q_new,Q_center,fact_q,q,q_inv,iorder_q)
  implicit none
  BEGIN_DOC
  ! Computes the integral <pq|rs> where p,q,r,s are Gaussian primitives
  END_DOC
  integer,intent(in)             :: dim
  include 'utils/constants.include.F'
  double precision, intent(in)   :: P_new(0:max_dim,3),P_center(3),fact_p,p,p_inv
  double precision, intent(in)   :: Q_new(0:max_dim,3),Q_center(3),fact_q,q,q_inv
  integer, intent(in)            :: iorder_p(3)
  integer, intent(in)            :: iorder_q(3)

  double precision               :: r_cut,gama_r_cut,rho,dist
  double precision               :: dx(0:max_dim),Ix_pol(0:max_dim),dy(0:max_dim),Iy_pol(0:max_dim),dz(0:max_dim),Iz_pol(0:max_dim)
  integer                        :: n_Ix,n_Iy,n_Iz,nx,ny,nz
  double precision               :: bla
  integer                        :: ix,iy,iz,jx,jy,jz,i
  double precision               :: a,b,c,d,e,f,accu,pq,const
  double precision               :: pq_inv, p10_1, p10_2, p01_1, p01_2,pq_inv_2
  integer                        :: n_pt_tmp,n_pt_out, iorder
  double precision               :: d1(0:max_dim),d_poly(0:max_dim),rint,d1_screened(0:max_dim)

  general_primitive_integral = 0.d0

  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: dx,Ix_pol,dy,Iy_pol,dz,Iz_pol
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: d1, d_poly

  ! Gaussian Product
  ! ----------------

  pq = p_inv*0.5d0*q_inv
  pq_inv = 0.5d0/(p+q)
  p10_1 = q*pq  ! 1/(2p)
  p01_1 = p*pq  ! 1/(2q)
  pq_inv_2 = pq_inv+pq_inv
  p10_2 = pq_inv_2 * p10_1*q !0.5d0*q/(pq + p*p)
  p01_2 = pq_inv_2 * p01_1*p !0.5d0*p/(q*q + pq)


  accu = 0.d0
  iorder = iorder_p(1)+iorder_q(1)+iorder_p(1)+iorder_q(1)
  do ix=0,iorder
    Ix_pol(ix) = 0.d0
  enddo
  n_Ix = 0
  do ix = 0, iorder_p(1)
    if (abs(P_new(ix,1)) < thresh) cycle
    a = P_new(ix,1)
    do jx = 0, iorder_q(1)
      d = a*Q_new(jx,1)
      if (abs(d) < thresh) cycle
      !DIR$ FORCEINLINE
      call give_polynom_mult_center_x(P_center(1),Q_center(1),ix,jx,p,q,iorder,pq_inv,pq_inv_2,p10_1,p01_1,p10_2,p01_2,dx,nx)
      !DIR$ FORCEINLINE
      call add_poly_multiply(dx,nx,d,Ix_pol,n_Ix)
    enddo
  enddo
  if (n_Ix == -1) then
    return
  endif
  iorder = iorder_p(2)+iorder_q(2)+iorder_p(2)+iorder_q(2)
  do ix=0, iorder
    Iy_pol(ix) = 0.d0
  enddo
  n_Iy = 0
  do iy = 0, iorder_p(2)
    if (abs(P_new(iy,2)) > thresh) then
      b = P_new(iy,2)
      do jy = 0, iorder_q(2)
        e = b*Q_new(jy,2)
        if (abs(e) < thresh) cycle
        !DIR$ FORCEINLINE
        call   give_polynom_mult_center_x(P_center(2),Q_center(2),iy,jy,p,q,iorder,pq_inv,pq_inv_2,p10_1,p01_1,p10_2,p01_2,dy,ny)
        !DIR$ FORCEINLINE
        call add_poly_multiply(dy,ny,e,Iy_pol,n_Iy)
      enddo
    endif
  enddo
  if (n_Iy == -1) then
    return
  endif

  iorder = iorder_p(3)+iorder_q(3)+iorder_p(3)+iorder_q(3)
  do ix=0,iorder
    Iz_pol(ix) = 0.d0
  enddo
  n_Iz = 0
  do iz = 0, iorder_p(3)
    if (abs(P_new(iz,3)) > thresh) then
      c = P_new(iz,3)
      do jz = 0, iorder_q(3)
        f = c*Q_new(jz,3)
        if (abs(f) < thresh) cycle
        !DIR$ FORCEINLINE
        call   give_polynom_mult_center_x(P_center(3),Q_center(3),iz,jz,p,q,iorder,pq_inv,pq_inv_2,p10_1,p01_1,p10_2,p01_2,dz,nz)
        !DIR$ FORCEINLINE
        call add_poly_multiply(dz,nz,f,Iz_pol,n_Iz)
      enddo
    endif
  enddo
  if (n_Iz == -1) then
    return
  endif

  rho = p*q *pq_inv_2
  dist =  (P_center(1) - Q_center(1))*(P_center(1) - Q_center(1)) +  &
      (P_center(2) - Q_center(2))*(P_center(2) - Q_center(2)) +      &
      (P_center(3) - Q_center(3))*(P_center(3) - Q_center(3))
  const = dist*rho

  n_pt_tmp = n_Ix+n_Iy
  do i=0,n_pt_tmp
    d_poly(i)=0.d0
  enddo

  !DIR$ FORCEINLINE
  call multiply_poly(Ix_pol,n_Ix,Iy_pol,n_Iy,d_poly,n_pt_tmp)
  if (n_pt_tmp == -1) then
    return
  endif
  n_pt_out = n_pt_tmp+n_Iz
  do i=0,n_pt_out
    d1(i)=0.d0
  enddo

  !DIR$ FORCEINLINE
  call multiply_poly(d_poly ,n_pt_tmp ,Iz_pol,n_Iz,d1,n_pt_out)
  double precision               :: rint_sum
  accu = accu + rint_sum(n_pt_out,const,d1)

  general_primitive_integral = fact_p * fact_q * accu *pi_5_2*p_inv*q_inv/dsqrt(p+q)
end


double precision function ERI(alpha,beta,delta,gama,a_x,b_x,c_x,d_x,a_y,b_y,c_y,d_y,a_z,b_z,c_z,d_z)
  implicit none
  BEGIN_DOC
  !  ATOMIC PRIMTIVE two-electron integral between the 4 primitives ::
  !         primitive_1 = x1**(a_x) y1**(a_y) z1**(a_z) exp(-alpha * r1**2)
  !         primitive_2 = x1**(b_x) y1**(b_y) z1**(b_z) exp(- beta * r1**2)
  !         primitive_3 = x2**(c_x) y2**(c_y) z2**(c_z) exp(-delta * r2**2)
  !         primitive_4 = x2**(d_x) y2**(d_y) z2**(d_z) exp(- gama * r2**2)
  END_DOC
  double precision, intent(in)   :: delta,gama,alpha,beta
  integer, intent(in)            :: a_x,b_x,c_x,d_x,a_y,b_y,c_y,d_y,a_z,b_z,c_z,d_z
  integer                        :: a_x_2,b_x_2,c_x_2,d_x_2,a_y_2,b_y_2,c_y_2,d_y_2,a_z_2,b_z_2,c_z_2,d_z_2
  integer                        :: i,j,k,l,n_pt
  integer                        :: n_pt_sup
  double precision               :: p,q,denom,coeff
  double precision               :: I_f
  integer                        :: nx,ny,nz
  include 'utils/constants.include.F'
  nx = a_x+b_x+c_x+d_x
  if(iand(nx,1) == 1) then
    ERI = 0.d0
    return
  endif

  ny = a_y+b_y+c_y+d_y
  if(iand(ny,1) == 1) then
    ERI = 0.d0
    return
  endif

  nz = a_z+b_z+c_z+d_z
  if(iand(nz,1) == 1) then
    ERI = 0.d0
    return
  endif

  ASSERT (alpha >= 0.d0)
  ASSERT (beta >= 0.d0)
  ASSERT (delta >= 0.d0)
  ASSERT (gama >= 0.d0)
  p = alpha + beta
  q = delta + gama
  ASSERT (p+q >= 0.d0)
  n_pt =  shiftl( nx+ny+nz,1 )

  coeff = pi_5_2 / (p * q * dsqrt(p+q))
  if (n_pt == 0) then
    ERI = coeff
    return
  endif

  call integrale_new(I_f,a_x,b_x,c_x,d_x,a_y,b_y,c_y,d_y,a_z,b_z,c_z,d_z,p,q,n_pt)

  ERI = I_f * coeff
end


subroutine integrale_new(I_f,a_x,b_x,c_x,d_x,a_y,b_y,c_y,d_y,a_z,b_z,c_z,d_z,p,q,n_pt)
  BEGIN_DOC
  ! Calculates the integral of the polynomial :
  !
  ! $I_{x_1}(a_x+b_x,c_x+d_x,p,q) \, I_{x_1}(a_y+b_y,c_y+d_y,p,q) \, I_{x_1}(a_z+b_z,c_z+d_z,p,q)$
  ! in $( 0 ; 1)$
  END_DOC


  implicit none
  include 'utils/constants.include.F'
  double precision               :: p,q
  integer                        :: a_x,b_x,c_x,d_x,a_y,b_y,c_y,d_y,a_z,b_z,c_z,d_z
  integer                        :: i, n_pt, j
  double precision               :: I_f, pq_inv, p10_1, p10_2, p01_1, p01_2,rho,pq_inv_2
  integer :: ix,iy,iz, jx,jy,jz, sx,sy,sz

  j = shiftr(n_pt,1)
  ASSERT (n_pt > 1)
  pq_inv = 0.5d0/(p+q)
  pq_inv_2 = pq_inv + pq_inv
  p10_1 = 0.5d0/p
  p01_1 = 0.5d0/q
  p10_2 = 0.5d0 *  q /(p * q + p * p)
  p01_2 = 0.5d0 *  p /(q * q + q * p)
  double precision               :: B00(n_pt_max_integrals)
  double precision               :: B10(n_pt_max_integrals), B01(n_pt_max_integrals)
  double precision               :: t1(n_pt_max_integrals), t2(n_pt_max_integrals)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: t1, t2, B10, B01, B00
  ix = a_x+b_x
  jx = c_x+d_x
  iy = a_y+b_y
  jy = c_y+d_y
  iz = a_z+b_z
  jz = c_z+d_z
  sx = ix+jx
  sy = iy+jy
  sz = iz+jz

  do i = 1,n_pt
    B10(i)  = p10_1 -  gauleg_t2(i,j)* p10_2
    B01(i)  = p01_1 -  gauleg_t2(i,j)* p01_2
    B00(i)  = gauleg_t2(i,j)*pq_inv
  enddo
  if (sx > 0) then
    call I_x1_new(ix,jx,B10,B01,B00,t1,n_pt)
  else
    do i = 1,n_pt
      t1(i) = 1.d0
    enddo
  endif
  if (sy > 0) then
    call I_x1_new(iy,jy,B10,B01,B00,t2,n_pt)
    do i = 1,n_pt
      t1(i) = t1(i)*t2(i)
    enddo
  endif
  if (sz > 0) then
    call I_x1_new(iz,jz,B10,B01,B00,t2,n_pt)
    do i = 1,n_pt
      t1(i) = t1(i)*t2(i)
    enddo
  endif
  I_f= 0.d0
  do i = 1,n_pt
    I_f += gauleg_w(i,j)*t1(i)
  enddo



end

recursive subroutine I_x1_new(a,c,B_10,B_01,B_00,res,n_pt)
  BEGIN_DOC
  !  recursive function involved in the two-electron integral
  END_DOC
  implicit none
  include 'utils/constants.include.F'
  integer, intent(in)            :: a,c,n_pt
  double precision, intent(in)   :: B_10(n_pt_max_integrals),B_01(n_pt_max_integrals),B_00(n_pt_max_integrals)
  double precision, intent(out)  :: res(n_pt_max_integrals)
  double precision               :: res2(n_pt_max_integrals)
  integer                        :: i

  if(c<0)then
    do i=1,n_pt
      res(i) = 0.d0
    enddo
  else if (a==0) then
    call I_x2_new(c,B_10,B_01,B_00,res,n_pt)
  else if (a==1) then
    call I_x2_new(c-1,B_10,B_01,B_00,res,n_pt)
    do i=1,n_pt
      res(i) = c * B_00(i) * res(i)
    enddo
  else
    call I_x1_new(a-2,c,B_10,B_01,B_00,res,n_pt)
    call I_x1_new(a-1,c-1,B_10,B_01,B_00,res2,n_pt)
    do i=1,n_pt
      res(i) = (a-1) * B_10(i) * res(i) &
               + c * B_00(i) * res2(i)
    enddo
  endif
end

recursive subroutine I_x2_new(c,B_10,B_01,B_00,res,n_pt)
  implicit none
  BEGIN_DOC
  !  recursive function involved in the two-electron integral
  END_DOC
  include 'utils/constants.include.F'
  integer, intent(in)            :: c, n_pt
  double precision, intent(in)   :: B_10(n_pt_max_integrals),B_01(n_pt_max_integrals),B_00(n_pt_max_integrals)
  double precision, intent(out)  :: res(n_pt_max_integrals)
  integer                        :: i

  if(c==1)then
    do i=1,n_pt
      res(i) = 0.d0
    enddo
  elseif(c==0) then
    do i=1,n_pt
      res(i) = 1.d0
    enddo
  else
    call I_x1_new(0,c-2,B_10,B_01,B_00,res,n_pt)
    do i=1,n_pt
      res(i) =  (c-1) * B_01(i) * res(i)
    enddo
  endif
end


integer function n_pt_sup(a_x,b_x,c_x,d_x,a_y,b_y,c_y,d_y,a_z,b_z,c_z,d_z)
  implicit none
  BEGIN_DOC
  ! Returns the upper boundary of the degree of the polynomial involved in the
  ! two-electron integral :
  !
  ! $I_x(a_x,b_x,c_x,d_x) \, I_y(a_y,b_y,c_y,d_y) \, I_z(a_z,b_z,c_z,d_z)$
  END_DOC
  integer                        :: a_x,b_x,c_x,d_x,a_y,b_y,c_y,d_y,a_z,b_z,c_z,d_z
  n_pt_sup =  shiftl( a_x+b_x+c_x+d_x + a_y+b_y+c_y+d_y + a_z+b_z+c_z+d_z,1 )
end




subroutine give_polynom_mult_center_x(P_center,Q_center,a_x,d_x,p,q,n_pt_in,pq_inv,pq_inv_2,p10_1,p01_1,p10_2,p01_2,d,n_pt_out)
  implicit none
  BEGIN_DOC
  ! subroutine that returns the explicit polynom in term of the "t"
  ! variable of the following polynomw :
  !
  ! $I_{x_1}(a_x,d_x,p,q) \, I_{x_1}(a_y,d_y,p,q) \ I_{x_1}(a_z,d_z,p,q)$
  END_DOC
  integer, intent(in)            :: n_pt_in
  integer,intent(out)            :: n_pt_out
  integer, intent(in)            :: a_x,d_x
  double precision, intent(in)   :: P_center, Q_center
  double precision, intent(in)   :: p,q,pq_inv,p10_1,p01_1,p10_2,p01_2,pq_inv_2
  include 'utils/constants.include.F'
  double precision,intent(out)   :: d(0:max_dim)
  double precision               :: accu
  accu = 0.d0
  ASSERT (n_pt_in >= 0)
  ! pq_inv = 0.5d0/(p+q)
  ! pq_inv_2 = 1.d0/(p+q)
  ! p10_1 = 0.5d0/p
  ! p01_1 = 0.5d0/q
  ! p10_2 = 0.5d0 *  q /(p * q + p * p)
  ! p01_2 = 0.5d0 *  p /(q * q + q * p)
  double precision               :: B10(0:2), B01(0:2), B00(0:2),C00(0:2),D00(0:2)
  B10(0)  = p10_1
  B10(1)  = 0.d0
  B10(2)  = - p10_2
  ! B10 = p01_1 - t**2 * p10_2
  B01(0)  = p01_1
  B01(1)  = 0.d0
  B01(2)  = - p01_2
  ! B01 = p01_1- t**2 * pq_inv
  B00(0)  = 0.d0
  B00(1)  = 0.d0
  B00(2)  = pq_inv
  ! B00 = t**2 * pq_inv
  do i = 0,n_pt_in
    d(i) = 0.d0
  enddo
  integer                        :: n_pt1,dim,i
  n_pt1 = n_pt_in
  ! C00 = -q/(p+q)*(Px-Qx) * t^2
  C00(0) = 0.d0
  C00(1) = 0.d0
  C00(2) =  -q*(P_center-Q_center) * pq_inv_2
  ! D00 = -p/(p+q)*(Px-Qx) * t^2
  D00(0) = 0.d0
  D00(1) = 0.d0
  D00(2) =  -p*(Q_center-P_center) * pq_inv_2
  !D00(2) =  -p*(Q_center(1)-P_center(1)) /(p+q)
  !DIR$ FORCEINLINE
  call I_x1_pol_mult(a_x,d_x,B10,B01,B00,C00,D00,d,n_pt1,n_pt_in)
  n_pt_out = n_pt1
  if(n_pt1<0)then
    n_pt_out = -1
    do i = 0,n_pt_in
      d(i) = 0.d0
    enddo
    return
  endif

end

subroutine I_x1_pol_mult(a,c,B_10,B_01,B_00,C_00,D_00,d,nd,n_pt_in)
  implicit none
  BEGIN_DOC
  ! Recursive function involved in the two-electron integral
  END_DOC
  integer , intent(in)           :: n_pt_in
  include 'utils/constants.include.F'
  double precision,intent(inout) :: d(0:max_dim)
  integer,intent(inout)          :: nd
  integer, intent(in)            :: a,c
  double precision, intent(in)   :: B_10(0:2),B_01(0:2),B_00(0:2),C_00(0:2),D_00(0:2)
  if( (c>=0).and.(nd>=0) )then

    if (a==1) then
      call I_x1_pol_mult_a1(c,B_10,B_01,B_00,C_00,D_00,d,nd,n_pt_in)
    else if (a==2) then
      call I_x1_pol_mult_a2(c,B_10,B_01,B_00,C_00,D_00,d,nd,n_pt_in)
    else if (a>2) then
      call I_x1_pol_mult_recurs(a,c,B_10,B_01,B_00,C_00,D_00,d,nd,n_pt_in)
    else  ! a == 0

      if( c==0 )then
        nd = 0
        d(0) = 1.d0
        return
      endif

      call I_x2_pol_mult(c,B_10,B_01,B_00,C_00,D_00,d,nd,n_pt_in)
    endif
  else
    nd = -1
  endif
end

recursive subroutine I_x1_pol_mult_recurs(a,c,B_10,B_01,B_00,C_00,D_00,d,nd,n_pt_in)
  implicit none
  BEGIN_DOC
  ! Recursive function involved in the two-electron integral
  END_DOC
  integer , intent(in)           :: n_pt_in
  include 'utils/constants.include.F'
  double precision,intent(inout) :: d(0:max_dim)
  integer,intent(inout)          :: nd
  integer, intent(in)            :: a,c
  double precision, intent(in)   :: B_10(0:2),B_01(0:2),B_00(0:2),C_00(0:2),D_00(0:2)
  double precision               :: X(0:max_dim)
  double precision               :: Y(0:max_dim)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: X,Y
  integer                        :: nx, ix,iy,ny

  ASSERT (a>2)
  !DIR$ LOOP COUNT(8)
  do ix=0,n_pt_in
    X(ix) = 0.d0
  enddo
  nx = 0
  if (a==3) then
    call I_x1_pol_mult_a1(c,B_10,B_01,B_00,C_00,D_00,X,nx,n_pt_in)
  else if (a==4) then
    call I_x1_pol_mult_a2(c,B_10,B_01,B_00,C_00,D_00,X,nx,n_pt_in)
  else
    ASSERT (a>=5)
    call I_x1_pol_mult_recurs(a-2,c,B_10,B_01,B_00,C_00,D_00,X,nx,n_pt_in)
  endif

  !DIR$ LOOP COUNT(8)
  do ix=0,nx
    X(ix) *= dble(a-1)
  enddo

  !DIR$ FORCEINLINE
  call multiply_poly(X,nx,B_10,2,d,nd)

  nx = nd
  !DIR$ LOOP COUNT(8)
  do ix=0,n_pt_in
    X(ix) = 0.d0
  enddo

  if (c>0) then
    if (a==3) then
      call I_x1_pol_mult_a2(c-1,B_10,B_01,B_00,C_00,D_00,X,nx,n_pt_in)
    else
      ASSERT(a >= 4)
      call I_x1_pol_mult_recurs(a-1,c-1,B_10,B_01,B_00,C_00,D_00,X,nx,n_pt_in)
    endif
    if (c>1) then
      !DIR$ LOOP COUNT(8)
      do ix=0,nx
        X(ix) *= c
      enddo
    endif
    !DIR$ FORCEINLINE
    call multiply_poly(X,nx,B_00,2,d,nd)
  endif

  ny=0

  !DIR$ LOOP COUNT(8)
  do ix=0,n_pt_in
    Y(ix) = 0.d0
  enddo
  ASSERT(a > 2)
  if (a==3) then
    call I_x1_pol_mult_a2(c,B_10,B_01,B_00,C_00,D_00,Y,ny,n_pt_in)
  else
    ASSERT(a >= 4)
    call I_x1_pol_mult_recurs(a-1,c,B_10,B_01,B_00,C_00,D_00,Y,ny,n_pt_in)
  endif

  !DIR$ FORCEINLINE
  call multiply_poly(Y,ny,C_00,2,d,nd)

end

recursive subroutine I_x1_pol_mult_a1(c,B_10,B_01,B_00,C_00,D_00,d,nd,n_pt_in)
  implicit none
  BEGIN_DOC
  ! Recursive function involved in the two-electron integral
  END_DOC
  integer , intent(in)           :: n_pt_in
  include 'utils/constants.include.F'
  double precision,intent(inout) :: d(0:max_dim)
  integer,intent(inout)          :: nd
  integer, intent(in)            :: c
  double precision, intent(in)   :: B_10(0:2),B_01(0:2),B_00(0:2),C_00(0:2),D_00(0:2)
  double precision               :: X(0:max_dim)
  double precision               :: Y(0:max_dim)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: X,Y
  integer                        :: nx, ix,iy,ny

  if( (c<0).or.(nd<0) )then
    nd = -1
    return
  endif

  nx = nd
  !DIR$ LOOP COUNT(8)
  do ix=0,n_pt_in
    X(ix) = 0.d0
  enddo
  call I_x2_pol_mult(c-1,B_10,B_01,B_00,C_00,D_00,X,nx,n_pt_in)

  if (c>1) then
    !DIR$ LOOP COUNT(8)
    do ix=0,nx
      X(ix) *= dble(c)
    enddo
  endif

  !DIR$ FORCEINLINE
  call multiply_poly(X,nx,B_00,2,d,nd)

  ny=0

  !DIR$ LOOP COUNT(8)
  do ix=0,n_pt_in
    Y(ix) = 0.d0
  enddo
  call I_x2_pol_mult(c,B_10,B_01,B_00,C_00,D_00,Y,ny,n_pt_in)

  !DIR$ FORCEINLINE
  call multiply_poly(Y,ny,C_00,2,d,nd)

end

recursive subroutine I_x1_pol_mult_a2(c,B_10,B_01,B_00,C_00,D_00,d,nd,n_pt_in)
  implicit none
  BEGIN_DOC
  ! Recursive function involved in the two-electron integral
  END_DOC
  integer , intent(in)           :: n_pt_in
  include 'utils/constants.include.F'
  double precision,intent(inout) :: d(0:max_dim)
  integer,intent(inout)          :: nd
  integer, intent(in)            :: c
  double precision, intent(in)   :: B_10(0:2),B_01(0:2),B_00(0:2),C_00(0:2),D_00(0:2)
  double precision               :: X(0:max_dim)
  double precision               :: Y(0:max_dim)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: X,Y
  integer                        :: nx, ix,iy,ny

  !DIR$ LOOP COUNT(8)
  do ix=0,n_pt_in
    X(ix) = 0.d0
  enddo
  nx = 0
  call I_x2_pol_mult(c,B_10,B_01,B_00,C_00,D_00,X,nx,n_pt_in)

  !DIR$ FORCEINLINE
  call multiply_poly(X,nx,B_10,2,d,nd)

  nx = nd
  !DIR$ LOOP COUNT(8)
  do ix=0,n_pt_in
    X(ix) = 0.d0
  enddo

  !DIR$ FORCEINLINE
  call I_x1_pol_mult_a1(c-1,B_10,B_01,B_00,C_00,D_00,X,nx,n_pt_in)

  if (c>1) then
    !DIR$ LOOP COUNT(8)
    do ix=0,nx
      X(ix) *= dble(c)
    enddo
  endif

  !DIR$ FORCEINLINE
  call multiply_poly(X,nx,B_00,2,d,nd)

  ny=0
  !DIR$ LOOP COUNT(8)
  do ix=0,n_pt_in
    Y(ix) = 0.d0
  enddo
  !DIR$ FORCEINLINE
  call I_x1_pol_mult_a1(c,B_10,B_01,B_00,C_00,D_00,Y,ny,n_pt_in)

  !DIR$ FORCEINLINE
  call multiply_poly(Y,ny,C_00,2,d,nd)

end

recursive subroutine I_x2_pol_mult(c,B_10,B_01,B_00,C_00,D_00,d,nd,dim)
  implicit none
  BEGIN_DOC
  ! Recursive function involved in the two-electron integral
  END_DOC
  integer , intent(in)           :: dim
  include 'utils/constants.include.F'
  double precision               :: d(0:max_dim)
  integer,intent(inout)          :: nd
  integer, intent(in)            :: c
  double precision, intent(in)   :: B_10(0:2),B_01(0:2),B_00(0:2),C_00(0:2),D_00(0:2)
  integer                        :: nx, ix,ny
  double precision               :: X(0:max_dim),Y(0:max_dim)
  !DIR$ ATTRIBUTES ALIGN : $IRP_ALIGN :: X, Y
  integer                        :: i

  select case (c)
    case (0)
      nd = 0
      d(0) = 1.d0
      return

    case (:-1)
      nd = -1
      return

    case (1)
      nd = 2
      d(0) = D_00(0)
      d(1) = D_00(1)
      d(2) = D_00(2)
      return

    case (2)
      nd = 2
      d(0) = B_01(0)
      d(1) = B_01(1)
      d(2) = B_01(2)

      ny = 2
      Y(0) = D_00(0)
      Y(1) = D_00(1)
      Y(2) = D_00(2)

      !DIR$ FORCEINLINE
      call multiply_poly(Y,ny,D_00,2,d,nd)
      return

      case default

      !DIR$ LOOP COUNT(6)
      do ix=0,c+c
        X(ix) = 0.d0
      enddo
      nx = 0
      call I_x2_pol_mult(c-2,B_10,B_01,B_00,C_00,D_00,X,nx,dim)

      !DIR$ LOOP COUNT(6)
      do ix=0,nx
        X(ix) *= dble(c-1)
      enddo

      !DIR$ FORCEINLINE
      call multiply_poly(X,nx,B_01,2,d,nd)

      ny = 0
      !DIR$ LOOP COUNT(6)
      do ix=0,c+c
        Y(ix) = 0.d0
      enddo
      call I_x2_pol_mult(c-1,B_10,B_01,B_00,C_00,D_00,Y,ny,dim)

      !DIR$ FORCEINLINE
      call multiply_poly(Y,ny,D_00,2,d,nd)

  end select
end




subroutine compute_ao_integrals_jl(j,l,n_integrals,buffer_i,buffer_value)
  implicit none
  use map_module
  BEGIN_DOC
  !  Parallel client for AO integrals
  END_DOC

  integer, intent(in)             :: j,l
  integer,intent(out)             :: n_integrals
  integer(key_kind),intent(out)   :: buffer_i(ao_num*ao_num)
  real(integral_kind),intent(out) :: buffer_value(ao_num*ao_num)
  logical, external               :: ao_two_e_integral_zero

  integer                         :: i,k
  double precision                :: ao_two_e_integral,cpu_1,cpu_2, wall_1, wall_2
  double precision                :: integral, wall_0
  double precision                :: thr
  integer                         :: kk, m, j1, i1

  thr = ao_integrals_threshold

  n_integrals = 0

  j1 = j+shiftr(l*l-l,1)
  do k = 1, ao_num           ! r1
    i1 = shiftr(k*k-k,1)
    if (i1 > j1) then
      exit
    endif
    do i = 1, k
      i1 += 1
      if (i1 > j1) then
        exit
      endif
      if (ao_two_e_integral_zero(i,j,k,l)) then
        cycle
      endif
      !DIR$ FORCEINLINE
      integral = ao_two_e_integral(i,k,j,l)  ! i,k : r1    j,l : r2
      if (abs(integral) < thr) then
        cycle
      endif
      n_integrals += 1
      !DIR$ FORCEINLINE
      call two_e_integrals_index(i,j,k,l,buffer_i(n_integrals))
      buffer_value(n_integrals) = integral
    enddo
  enddo

end
