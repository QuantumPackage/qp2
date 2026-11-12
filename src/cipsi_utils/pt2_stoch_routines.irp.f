BEGIN_PROVIDER [ integer, pt2_stoch_istate ]
 implicit none
 BEGIN_DOC
 ! State for stochatsic PT2
 END_DOC
 pt2_stoch_istate = 1
END_PROVIDER

 BEGIN_PROVIDER [ integer, pt2_F, (N_det_generators) ]
&BEGIN_PROVIDER [ integer, pt2_n_tasks_max ]
  implicit none
  BEGIN_DOC
  ! pt2_F(i) is the number of Monte Carlo tasks (combs)
  ! assigned to generator i (i indexes the norm-sorted
  ! generator list):
  ! pt2_n_tasks_max*pt2_min_parallel_tasks for the leading
  ! (largest-norm) generators, pt2_min_parallel_tasks for a
  ! middle range, and 1 for the rest.
  END_DOC
  logical, external :: testTeethBuilding
  integer :: i,j
  pt2_n_tasks_max = elec_alpha_num*elec_alpha_num + elec_alpha_num*elec_beta_num  - n_core_orb*2
  pt2_n_tasks_max = min(pt2_n_tasks_max,1+N_det_generators/10000)
  call write_int(6,pt2_n_tasks_max,'pt2_n_tasks_max')

  pt2_F(:) = max(int(sqrt(float(pt2_n_tasks_max))),1)
  do i=1,pt2_n_0(1+pt2_N_teeth/4)
    pt2_F(i) = pt2_n_tasks_max*pt2_min_parallel_tasks
  enddo
  do i=1+pt2_n_0(pt2_N_teeth-pt2_N_teeth/4), pt2_n_0(pt2_N_teeth-pt2_N_teeth/10)
    pt2_F(i) = pt2_min_parallel_tasks
  enddo
  do i=1+pt2_n_0(pt2_N_teeth-pt2_N_teeth/10), N_det_generators
    pt2_F(i) = 1
  enddo

END_PROVIDER

 BEGIN_PROVIDER [ integer, pt2_N_teeth ]
&BEGIN_PROVIDER [ integer, pt2_minDetInFirstTeeth ]
  implicit none
  BEGIN_DOC
  ! pt2_N_teeth: number of comb teeth of the stochastic PT2
  ! sampling.  1 means "no comb": the caller must then fall back
  ! to the deterministic selection (omp_selection), because a
  ! one-tooth comb samples one generator per comb and costs as
  ! much as the exact pass while adding noise.  This happens when
  ! N_det_generators < 1024, or when no comb can be built (a
  ! warning is printed in that case).
  ! Otherwise: the largest N from 100 down to 2 for which
  ! testTeethBuilding(pt2_minDetInFirstTeeth, N) is true.
  ! pt2_minDetInFirstTeeth: minimum number of determinants
  ! required in the first tooth (1 or min(5, N_det_generators)).
  END_DOC
  logical, external :: testTeethBuilding
  integer :: i
  logical :: found

  ! Never leave pt2_N_teeth to an implicit value: the old loop
  ! `do pt2_N_teeth=100,2,-1' exited with pt2_N_teeth=1 when every
  ! test failed, silently turning a large run into the deterministic
  ! fallback.
  pt2_N_teeth = 1
  if(N_det_generators < 1024) then
    pt2_minDetInFirstTeeth = 1
  else
    pt2_minDetInFirstTeeth = min(5, N_det_generators)
    found = .False.
    do i=100,2,-1
      if(testTeethBuilding(pt2_minDetInFirstTeeth, i)) then
        pt2_N_teeth = i
        found = .True.
        exit
      end if
    end do
    if(.not.found) then
      print *, 'WARNING: no comb can be built: N_det_generators =', N_det_generators
      print *, 'WARNING: falling back to the deterministic PT2 selection (no comb)'
    end if
  end if
  ! Report the effective number of teeth: 0 means that no comb is used.
  ! pt2_N_teeth itself stays 1, so that every array sized pt2_N_teeth+1
  ! remains valid.
  call write_int(6, merge(pt2_N_teeth,0,pt2_N_teeth>=2), &
                 'Number of comb teeth (0: no comb)')
END_PROVIDER


logical function testTeethBuilding(minF, N)
  implicit none
  BEGIN_DOC
! Checks whether the comb can be built with N teeth so that
! the first tooth holds at least minF determinants. Tries
! leading offsets n0 of the norm-sorted generator list and
! returns .true. if the first minF determinants after n0 fit
! in one tooth of width Wt = (1 - u0)/N (u0: cumulative
! weight before n0); .false. if Wt <= 1e-3 or no offset works.
  END_DOC
  integer, intent(in) :: minF, N
  integer :: n0, i
  double precision :: u0, Wt, r

  double precision, allocatable :: tilde_w(:), tilde_cW(:)
  integer, external :: dress_find_sample

  double precision :: rss
  double precision, external :: memory_of_double, memory_of_int

  rss = memory_of_double(2*N_det_generators+1)
  call check_mem(rss,irp_here)

  allocate(tilde_w(N_det_generators), tilde_cW(0:N_det_generators))

  double precision :: norm2
  norm2 = 0.d0
  do i=N_det_generators,1,-1
    tilde_w(i)  = psi_coef_sorted_gen(i,pt2_stoch_istate) * &
                  psi_coef_sorted_gen(i,pt2_stoch_istate)
    norm2 = norm2 + tilde_w(i)
  enddo

  f = 1.d0/norm2
  tilde_w(:) = tilde_w(:) * f

  tilde_cW(0) = -1.d0
  do i=1,N_det_generators
    tilde_cW(i) = tilde_cW(i-1) + tilde_w(i)
  enddo
  tilde_cW(:) = tilde_cW(:) + 1.d0
  deallocate(tilde_w)

  n0 = 0
  testTeethBuilding = .false.
  double precision :: f
  integer :: minFN
  minFN = N_det_generators - minF * N
  f = 1.d0/dble(N)
  do
    u0 = tilde_cW(n0)
    r = tilde_cW(n0 + minF)
    Wt = (1d0 - u0) * f
    if (dabs(Wt) <= 1.d-3) then
      exit
    endif
    if(Wt >= r - u0) then
       testTeethBuilding = .true.
       exit
    end if
    n0 += 1
    if(n0 > minFN) then
      exit
    end if
  end do
  deallocate(tilde_cW)

end function


integer function pt2_find_sample(v, w)
  implicit none
  BEGIN_DOC
! Returns the index where the cumulative weight w crosses v,
! by calling pt2_find_sample_lr on the full range
! 0..N_det_generators.
  END_DOC
  double precision, intent(in) :: v, w(0:N_det_generators)
  integer, external :: pt2_find_sample_lr

  pt2_find_sample = pt2_find_sample_lr(v, w, 0, N_det_generators)
end function


integer function pt2_find_sample_lr(v, w, l_in, r_in)
  implicit none
  BEGIN_DOC
! Binary search on the cumulative weight array w for the
! first index r such that w(r) >= v, within the range
! l_in..r_in. If w is flat there (determinants of zero
! weight), the last index of the flat run is returned.
  END_DOC
  double precision, intent(in) :: v, w(0:N_det_generators)
  integer, intent(in) :: l_in,r_in
  integer :: i,l,r

  l=l_in
  r=r_in

  do while(r-l > 1)
    i = shiftr(r+l,1)
    if(w(i) < v) then
      l = i
    else
      r = i
    end if
  end do
  i = r
  do r=i+1,N_det_generators
    if (w(r) /= w(i)) then
      exit
    endif
  enddo
  pt2_find_sample_lr = r-1
end function


BEGIN_PROVIDER [ integer, pt2_n_tasks ]
 implicit none
 BEGIN_DOC
 ! Number of parallel tasks for the Monte Carlo
 END_DOC
 pt2_n_tasks = N_det_generators
END_PROVIDER

BEGIN_PROVIDER[ double precision, pt2_u, (N_det_generators)]
  implicit none
  BEGIN_DOC
  ! Uniform random numbers in (0,1), one per Monte Carlo
  ! sample (comb) of the stochastic PT2. The random generator
  ! is seeded with (1, 2, ..., m) before being filled.
  END_DOC
  integer, allocatable :: seed(:)
  integer :: m,i
  call random_seed(size=m)
  allocate(seed(m))
  do i=1,m
    seed(i) = i
  enddo
  call random_seed(put=seed)
  deallocate(seed)

  call RANDOM_NUMBER(pt2_u)
 END_PROVIDER

 BEGIN_PROVIDER[ integer, pt2_J, (N_det_generators)]
&BEGIN_PROVIDER[ integer, pt2_R, (N_det_generators)]
  implicit none
  BEGIN_DOC
! pt2_J contains the list of generators after ordering them according to the
! Monte Carlo sampling.
!
! pt2_R(i) is the number of combs drawn when determinant i is computed.
  END_DOC
  integer                :: N_c, N_j
  integer                :: U, t, i
  double precision       :: v
  integer, external      :: pt2_find_sample_lr

  logical, allocatable :: pt2_d(:)
  integer :: m,l,r,k
  integer :: ncache
  integer, allocatable :: ii(:,:)
  double precision :: dt

  ncache = min(N_det_generators,10000)

  double precision :: rss
  double precision, external :: memory_of_double, memory_of_int
  rss = memory_of_int(ncache)*dble(pt2_N_teeth) + memory_of_int(N_det_generators)
  call check_mem(rss,irp_here)

  allocate(ii(pt2_N_teeth,ncache),pt2_d(N_det_generators))

  pt2_R(:) = 0
  pt2_d(:) = .false.
  N_c = 0
  N_j = pt2_n_0(1)
  do i=1,N_j
      pt2_d(i) = .true.
      pt2_J(i) = i
  end do

  U = 0
  do while(N_j < pt2_n_tasks)

    if (N_c+ncache > N_det_generators) then
      ncache = N_det_generators - N_c
    endif

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(dt,v,t,k)
    do k=1, ncache
      dt = pt2_u_0
      do t=1, pt2_N_teeth
        v = dt + pt2_W_T *pt2_u(N_c+k)
        dt = dt + pt2_W_T
        ii(t,k) = pt2_find_sample_lr(v, pt2_cW,pt2_n_0(t),pt2_n_0(t+1))
      end do
    enddo
    !$OMP END PARALLEL DO

    do k=1,ncache
      !ADD_COMB
      N_c = N_c+1
      do t=1, pt2_N_teeth
        i = ii(t,k)
        if(.not. pt2_d(i)) then
          N_j += 1
          pt2_J(N_j) = i
          pt2_d(i) = .true.
        end if
      end do

      pt2_R(N_j) = N_c

      !FILL_TOOTH
      do while(U < N_det_generators)
        U += 1
        if(.not. pt2_d(U)) then
          N_j += 1
          pt2_J(N_j) = U
          pt2_d(U) = .true.
          exit
        end if
      end do
      if (N_j >= pt2_n_tasks) exit
    end do
  enddo

  if(N_det_generators > 1) then
    pt2_R(N_det_generators-1) = 0
    pt2_R(N_det_generators) = N_c
  end if

  deallocate(ii,pt2_d)

END_PROVIDER



 BEGIN_PROVIDER [ double precision, pt2_w, (N_det_generators) ]
&BEGIN_PROVIDER [ double precision, pt2_cW, (0:N_det_generators) ]
&BEGIN_PROVIDER [ double precision, pt2_W_T ]
&BEGIN_PROVIDER [ double precision, pt2_u_0 ]
&BEGIN_PROVIDER [ integer,          pt2_n_0, (pt2_N_teeth+1) ]
   implicit none
   BEGIN_DOC
   ! Weights of the Monte Carlo comb of the stochastic PT2.
   ! pt2_w(i): weight of norm-sorted generator i, rescaled so
   ! that each of the pt2_N_teeth teeth holds the total weight
   ! pt2_W_T. pt2_cW: cumulative sum of pt2_w. pt2_u_0: weight
   ! of the first pt2_n_0(1) generators, always included in the
   ! comb. pt2_n_0: index bounds of the teeth.
   END_DOC
   integer                        :: i, t
   double precision, allocatable  :: tilde_w(:), tilde_cW(:)
   double precision               :: r, tooth_width
   integer, external              :: pt2_find_sample

   double precision               :: rss
   double precision, external     :: memory_of_double, memory_of_int
   rss = memory_of_double(2*N_det_generators+1)
   call check_mem(rss,irp_here)

   if (N_det_generators == 1) then

     pt2_w(1)   = 1.d0
     pt2_cw(1)  = 1.d0
     pt2_u_0    = 1.d0
     pt2_W_T    = 0.d0
     pt2_n_0(1) = 0
     pt2_n_0(2) = 1

   else

     allocate(tilde_w(N_det_generators), tilde_cW(0:N_det_generators))

     tilde_cW(0) = 0d0

     do i=1,N_det_generators
       tilde_w(i)  = psi_coef_sorted_gen(i,pt2_stoch_istate)**2 !+ 1.d-20
     enddo

     double precision               :: norm2
     norm2 = 0.d0
     do i=N_det_generators,1,-1
       norm2 += tilde_w(i)
     enddo

     tilde_w(:) = tilde_w(:) / norm2

     tilde_cW(0) = -1.d0
     do i=1,N_det_generators
       tilde_cW(i) = tilde_cW(i-1) + tilde_w(i)
     enddo
     tilde_cW(:) = tilde_cW(:) + 1.d0

     pt2_n_0(1) = 0
     do
     pt2_u_0 = tilde_cW(pt2_n_0(1))
     r = tilde_cW(pt2_n_0(1) + pt2_minDetInFirstTeeth)
     pt2_W_T = (1d0 - pt2_u_0) / dble(pt2_N_teeth)
     if(pt2_W_T >= r - pt2_u_0) then
       exit
     end if
     pt2_n_0(1) += 1
     if(N_det_generators - pt2_n_0(1) < pt2_minDetInFirstTeeth * pt2_N_teeth) then
       print *, "teeth building failed"
       stop -1
     end if
   end do

   do t=2, pt2_N_teeth
     r = pt2_u_0 + pt2_W_T * dble(t-1)
     pt2_n_0(t) = pt2_find_sample(r, tilde_cW)
   end do
   pt2_n_0(pt2_N_teeth+1) = N_det_generators

   pt2_w(:pt2_n_0(1)) = tilde_w(:pt2_n_0(1))
   do t=1, pt2_N_teeth
     tooth_width = tilde_cW(pt2_n_0(t+1)) - tilde_cW(pt2_n_0(t))
     if (tooth_width == 0.d0) then
       tooth_width = max(1.d-15,sum(tilde_w(pt2_n_0(t):pt2_n_0(t+1))))
     endif
     do i=pt2_n_0(t)+1, pt2_n_0(t+1)
       pt2_w(i) = tilde_w(i) * pt2_W_T / tooth_width
     end do
   end do

   pt2_cW(0) = 0d0
   do i=1,N_det_generators
     pt2_cW(i) = pt2_cW(i-1) + pt2_w(i)
   end do
   pt2_n_0(pt2_N_teeth+1) = N_det_generators

 endif
END_PROVIDER





BEGIN_PROVIDER [ double precision, pt2_overlap, (N_states, N_states) ]
 implicit none
 BEGIN_DOC
 ! Overlap between the perturbed wave functions
 END_DOC
 pt2_overlap(1:N_states,1:N_states) = 0.d0
END_PROVIDER


