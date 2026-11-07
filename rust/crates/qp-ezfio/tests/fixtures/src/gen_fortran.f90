! Reference fixture generator (Fortran producer).
!
! Writes one EZFIO tree with the fixture corpus used by
! `rust/crates/qp-ezfio/tests/golden.rs`. The bytes stored under
! `tests/fixtures/fortran/` are the *decompressed* payloads of the files
! produced here (D6: comparing raw `.gz` bytes is meaningless).
!
! Build and run (see tests/fixtures/README.md for the exact tool versions):
!
!   gfortran -I$HOME/qp2/lib gen_fortran.f90 $HOME/qp2/lib/libezfio_irp.a \
!            -lz -o gen_fortran
!   ./gen_fortran /path/to/empty/dir
!
! The corpus is deliberately written with the library entry points the QP2
! Fortran modules use (`ezfio_write_array_*`, `ezfio_write_*`), not with raw
! `write` statements, so the fixtures capture the real client behaviour.
program gen_fortran
  implicit none
  character*(256) :: dir

  double precision :: reals(6), redge(6), r2d(6), edat(1)
  integer          :: ints(4), i3d(6)
  integer*8        :: i8s(3)
  logical          :: ldat(3)
  character*(8)    :: cdat(3)
  character*(32)   :: labels(3)

  call get_command_argument(1, dir)
  call ezfio_set_file(trim(dir))
  call ezfio_mkdir(trim(dir)//'/oracle')

  ! reals: ordinary values, including both signed zeros.
  reals = (/ 1.d0, -0.5d0, 1.d-3, 1234567.8901234567d0, 0.d0, -0.d0 /)
  call ezfio_write_array_do(trim(dir)//'/oracle', 'reals', 1, (/6/), 6, reals)

  ! redge: exponents that stress the `Ew.d` exponent field (two-digit and
  ! three-digit, the latter drops the `E`) and the smallest normal/subnormal
  ! range.
  redge = (/ 1.d-100, 1.d100, -1.d-5, 3.141592653589793d0, 1.d300, -2.5d-308 /)
  call ezfio_write_array_do(trim(dir)//'/oracle', 'redge', 1, (/6/), 6, redge)

  ! r2d: rank 2, Fortran order (first index fastest).
  r2d = (/ 1.d0, 2.d0, 3.d0, 4.d0, 5.d0, 6.d0 /)
  call ezfio_write_array_do(trim(dir)//'/oracle', 'r2d', 2, (/2, 3/), 6, r2d)

  ! ints: `integer` (the format has one integer width, written with I20).
  ints = (/ 1, -42, 1000000, 0 /)
  call ezfio_write_array_in(trim(dir)//'/oracle', 'ints', 1, (/4/), 4, ints)

  ! i3d: rank 3, all extents > 1 except the first.
  i3d = (/ 1, 2, 3, 4, 5, 6 /)
  call ezfio_write_array_in(trim(dir)//'/oracle', 'i3d', 3, (/1, 2, 3/), 6, i3d)

  ! i8s: `integer*8`, including a value that a 32-bit integer cannot hold.
  i8s = (/ 1_8, -9007199254740993_8, 123456789012345_8 /)
  call ezfio_write_array_i8(trim(dir)//'/oracle', 'i8s', 1, (/3/), 3, i8s)

  ! logicals: `L1` padded to the 24-column element slot.
  ldat = (/.true., .false., .true./)
  call ezfio_write_array_lo(trim(dir)//'/oracle', 'logicals', 1, (/3/), 3, ldat)

  ! chars: blank-padded to len(dat(1)) = 8, including a leading blank.
  cdat = (/ 'abc     ', '  spaced', 'x       ' /)
  call ezfio_write_array_ch(trim(dir)//'/oracle', 'chars', 1, (/3/), 3, cdat)

  ! empty: a zero extent is legal; no element lines follow the header.
  edat = 0.d0
  call ezfio_write_array_do(trim(dir)//'/oracle', 'empty', 2, (/3, 0/), 0, edat)

  ! Scalars.
  call ezfio_write_do(trim(dir)//'/oracle', 'rscalar', 2.5d0)
  call ezfio_write_in(trim(dir)//'/oracle', 'iscalar', 7)
  call ezfio_write_lo(trim(dir)//'/oracle', 'lscalar', .true.)
  call ezfio_write_ch(trim(dir)//'/oracle', 'cscalar', 'hello')

  ! A registry-known key, written the way the `nuclei` module writes it:
  ! `character*(32)`, dims `(nuclei.nucl_num)`. The crate pads to the declared
  ! width under ByteLayout::Fortran, so this fixture is byte-comparable.
  call ezfio_mkdir(trim(dir)//'/nuclei')
  labels = (/ 'H', 'O', 'O' /)
  call ezfio_write_array_ch(trim(dir)//'/nuclei', 'nucl_label', 1, (/3/), 3, labels)

  ! A padded scalar: the Fortran client blanks the file to len(variable), the
  ! readers must trim. Read-only fixture (the crate writes scalars verbatim).
  call ezfio_mkdir(trim(dir)//'/oracle3')
  labels(1) = 'short'
  call ezfio_write_ch(trim(dir)//'/oracle3', 'short_padded', labels(1))
end program
