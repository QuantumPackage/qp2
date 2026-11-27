module sizet
! 64-bit integer kind, the Fortran counterpart of the C size_t on all
! platforms QP3 builds on (x86-64, arm64; gfortran, ifort, nvfortran,
! flang, Cray). Used for array sizes, counts and loop indices that scale
! with system size.
  integer, parameter :: size_t = 8
end module sizet
