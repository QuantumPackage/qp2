use bitmasks

BEGIN_PROVIDER [ integer, N_det_generators_CAS ]
  implicit none
  BEGIN_DOC
  ! Number of generator detetrminants
  END_DOC
  integer                        :: i,k,l
  logical                        :: good
  integer, external :: number_of_holes,number_of_particles
  call write_time(6)
  N_det_generators_CAS = 0
  do i=1,N_det
    good = ( number_of_holes(psi_det_sorted(1,1,i)) ==0).and.(number_of_particles(psi_det_sorted(1,1,i))==0 )
    if (good) then
      N_det_generators_CAS += 1
    endif
  enddo
  N_det_generators_CAS = max(N_det_generators_CAS,1)
  call write_int(6,N_det_generators_CAS,'Number of generators_CAS')
END_PROVIDER

 BEGIN_PROVIDER [ integer(bit_kind), psi_det_generators_CAS, (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ double precision, psi_coef_generators_CAS, (psi_det_size,N_states) ]
&BEGIN_PROVIDER [ integer(bit_kind), psi_det_sorted_gen_CAS, (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ double precision, psi_coef_sorted_gen_CAS, (psi_det_size,N_states) ]
&BEGIN_PROVIDER [ integer, psi_det_sorted_gen_CAS_order, (psi_det_size) ]
  implicit none
  BEGIN_DOC
  ! For Single reference wave functions, the gen_CASerator is the
  ! Hartree-Fock determinant
  END_DOC
  integer                        :: i, k, l, m
  logical                        :: good
  integer, external :: number_of_holes,number_of_particles
  integer, allocatable :: nongen_CAS(:)
  integer :: inongen_CAS

  allocate(nongen_CAS(N_det))

  inongen_CAS = 0
  m=0
  do i=1,N_det
    good = ( number_of_holes(psi_det_sorted(1,1,i)) ==0).and.(number_of_particles(psi_det_sorted(1,1,i))==0 )
    if (good) then
      m = m+1
      psi_det_sorted_gen_CAS_order(i) = m
      do k=1,N_int
        psi_det_generators_CAS(k,1,m) = psi_det_sorted(k,1,i)
        psi_det_generators_CAS(k,2,m) = psi_det_sorted(k,2,i)
      enddo
      psi_coef_generators_CAS(m,:) = psi_coef_sorted(i,:)
    else
      inongen_CAS += 1
      nongen_CAS(inongen_CAS) = i
    endif
  enddo
  ASSERT (m == N_det_generators_CAS)

  psi_det_sorted_gen_CAS(:,:,:N_det_generators_CAS) = psi_det_generators_CAS(:,:,:N_det_generators_CAS)
  psi_coef_sorted_gen_CAS(:N_det_generators_CAS, :) = psi_coef_generators_CAS(:N_det_generators_CAS, :)
  do i=1,inongen_CAS
    psi_det_sorted_gen_CAS_order(nongen_CAS(i)) = N_det_generators_CAS+i
    psi_det_sorted_gen_CAS(:,:,N_det_generators_CAS+i) = psi_det_sorted(:,:,nongen_CAS(i))
    psi_coef_sorted_gen_CAS(N_det_generators_CAS+i, :) = psi_coef_sorted(nongen_CAS(i),:)
  end do

END_PROVIDER

