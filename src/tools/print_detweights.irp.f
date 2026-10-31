program print_detweights

  implicit none
  BEGIN_DOC
  ! Reads the wave function from the |EZFIO| directory and writes
  ! the excitation degree of each determinant with respect to the
  ! first determinant to the file exc.dat.
  END_DOC

  read_wf = .True.
  touch read_wf

  call print_exc()
  !call main()

end

! ---

subroutine main()

  implicit none
  BEGIN_DOC
  ! Computes the excitation degree of each determinant with
  ! respect to the first determinant, sorts the coefficients of
  ! the first state by decreasing magnitude and writes the
  ! degrees and the sorted coefficients to the file coef.bin.
  END_DOC
  integer                       :: i
  integer                       :: degree
  integer                       :: ios
  integer,          allocatable :: deg(:), ii(:), deg_sorted(:)
  double precision, allocatable :: c(:)

  PROVIDE N_int
  PROVIDE N_det
  PROVIDE psi_det
  PROVIDe psi_coef

  allocate(deg(N_det), ii(N_det), deg_sorted(N_det), c(N_det))

  do i = 1, N_det

    call debug_det(psi_det(1,1,i), N_int)
    call get_excitation_degree(psi_det(1,1,i), psi_det(1,1,1), degree, N_int)

    ii (i) = i
    deg(i) = degree
    c  (i) = dabs(psi_coef(i,1))
  enddo

  call dsort(c, ii, N_det)

  do i = 1, N_det
    deg_sorted(i) = deg(ii(i))
    print *, deg_sorted(i), c(i)
  enddo

  print *, ' saving psi'

  ! Writing output in binary format
  open(unit=10, file="coef.bin", status="replace", action="write", iostat=ios, form="unformatted")

    if(ios /= 0) then
      print *, ' Error opening file!'
      stop
    endif

    write(10) N_det
    write(10) deg_sorted
    write(10) c

  close(10)

  deallocate(deg, ii, deg_sorted, c)

end

! ---

subroutine print_exc()

  implicit none
  BEGIN_DOC
  ! Computes the excitation degree of each determinant with
  ! respect to the first determinant and writes them to the file
  ! exc.dat, preceded by N_det.
  END_DOC

  integer              :: i
  integer, allocatable :: deg(:)

  PROVIDE N_int
  PROVIDE N_det
  PROVIDE psi_det

  allocate(deg(N_det))

  do i = 1, N_det
    call get_excitation_degree(psi_det(1,1,1), psi_det(1,1,i), deg(i), N_int)
  enddo

  open(unit=10, file="exc.dat", action="write")
    write(10,*) N_det
    write(10,*) deg
  close(10)

  deallocate(deg)

end



