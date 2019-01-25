use bitmasks

BEGIN_PROVIDER [ integer, N_det_generators ]
  implicit none
  BEGIN_DOC
  ! Number of generator detetrminants
  END_DOC
  integer                        :: i,k,l
  logical                        :: good
  integer, external :: number_of_holes,number_of_particles
  call write_time(6)
  N_det_generators = 0
  do i=1,N_det
    good = ( number_of_holes(psi_det_sorted(1,1,i)) ==0).and.(number_of_particles(psi_det_sorted(1,1,i))==0 )
    if (good) then
      N_det_generators += 1
    endif
  enddo
  N_det_generators = max(N_det_generators,1)
  call write_int(6,N_det_generators,'Number of generators')
END_PROVIDER

 BEGIN_PROVIDER [ integer(bit_kind), psi_det_generators, (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ double precision, psi_coef_generators, (psi_det_size,N_states) ]
&BEGIN_PROVIDER [ integer(bit_kind), psi_det_sorted_gen, (N_int,2,psi_det_size) ]
&BEGIN_PROVIDER [ double precision, psi_coef_sorted_gen, (psi_det_size,N_states) ]
&BEGIN_PROVIDER [ integer, psi_det_sorted_gen_order, (psi_det_size) ]
  implicit none
  BEGIN_DOC
  ! For Single reference wave functions, the generator is the
  ! Hartree-Fock determinant
  END_DOC
  integer                        :: i, k, l, m
  logical                        :: good
  integer, external :: number_of_holes,number_of_particles
  integer, allocatable :: nongen(:)
  integer :: inongen

  allocate(nongen(N_det))

  inongen = 0
  m=0
  do i=1,N_det
    good = ( number_of_holes(psi_det_sorted(1,1,i)) ==0).and.(number_of_particles(psi_det_sorted(1,1,i))==0 )
    if (good) then
      m = m+1
      psi_det_sorted_gen_order(i) = m
      do k=1,N_int
        psi_det_generators(k,1,m) = psi_det_sorted(k,1,i)
        psi_det_generators(k,2,m) = psi_det_sorted(k,2,i)
      enddo
      psi_coef_generators(m,:) = psi_coef_sorted(i,:)
    else
      inongen += 1
      nongen(inongen) = i
    endif
  enddo

  psi_det_sorted_gen(:,:,:N_det_generators) = psi_det_generators(:,:,:N_det_generators)
  psi_coef_sorted_gen(:N_det_generators, :) = psi_coef_generators(:N_det_generators, :)
  do i=1,inongen
    psi_det_sorted_gen_order(nongen(i)) = N_det_generators+i
    psi_det_sorted_gen(:,:,N_det_generators+i) = psi_det_sorted(:,:,nongen(i))
    psi_coef_sorted_gen(N_det_generators+i, :) = psi_coef_sorted(nongen(i),:)
  end do
END_PROVIDER

BEGIN_PROVIDER [ integer, size_select_max]
  implicit none
  BEGIN_DOC
  ! Size of the select_max array
  END_DOC
  size_select_max = 10000
END_PROVIDER

BEGIN_PROVIDER [ double precision, select_max, (size_select_max) ]
  implicit none
  BEGIN_DOC
  ! Memo to skip useless selectors
  END_DOC
  select_max = huge(1.d0)
END_PROVIDER

