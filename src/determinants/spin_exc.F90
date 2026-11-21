! Fortran iso_c_binding interfaces for the three generic C dispatcher
! functions:
!
!   get_all_spin_singles
!   get_all_spin_doubles
!   get_all_spin_singles_and_doubles
!
! Memory layout note:
!   The C functions expect buffer in column-major order:
!     buffer[k + N_int * i]  = word k of determinant i
!   This is identical to Fortran's natural layout for
!     buffer(N_int, size_buffer)
!   so no transposition is needed.
!
! Index convention:
!   idx contains 1-based Fortran indices; the C functions treat them
!   as opaque values and copy them straight into the output arrays,
!   so no adjustment is needed on either side.
!
! Usage example:
!   use spin_excitations_c
!   integer(c_int64_t), allocatable :: buffer(:,:), spindet(:)
!   integer,           allocatable :: idx(:), singles(:), doubles(:)
!   integer :: n_singles, n_doubles, Nint, size_buffer
!   ...
!   call get_all_spin_singles(buffer, idx, spindet, Nint, size_buffer,
!                               singles, n_singles)
!   call get_all_spin_singles_and_doubles(buffer, idx, spindet, Nint,
!                               size_buffer, singles, doubles,
!                               n_singles, n_doubles)

module spin_excitations
  use iso_c_binding
  implicit none


  ! ------------------------------------------------------------------
  ! Raw C interfaces (not intended to be called directly from user code)
  ! ------------------------------------------------------------------
  interface

    subroutine c_get_all_spin_singles_raw(          &
          buffer, idx, spindet,                     &
          N_int, size_buffer,                       &
          singles, n_singles)                       &
        bind(C, name="get_all_spin_singles")
      import :: c_int64_t, c_int
      implicit none
      integer(c_int64_t), intent(in)  :: buffer(*)
      integer(c_int),     intent(in)  :: idx(*)
      integer(c_int64_t), intent(in)  :: spindet(*)
      integer(c_int),     value       :: N_int
      integer(c_int),     value       :: size_buffer
      integer(c_int),     intent(out) :: singles(*)
      integer(c_int),     intent(out) :: n_singles
    end subroutine

    subroutine c_get_all_spin_doubles_raw(          &
          buffer, idx, spindet,                     &
          N_int, size_buffer,                       &
          doubles, n_doubles)                       &
        bind(C, name="get_all_spin_doubles")
      import :: c_int64_t, c_int
      implicit none
      integer(c_int64_t), intent(in)  :: buffer(*)
      integer(c_int),     intent(in)  :: idx(*)
      integer(c_int64_t), intent(in)  :: spindet(*)
      integer(c_int),     value       :: N_int
      integer(c_int),     value       :: size_buffer
      integer(c_int),     intent(out) :: doubles(*)
      integer(c_int),     intent(out) :: n_doubles
    end subroutine

    subroutine c_get_all_spin_singles_and_doubles_raw( &
          buffer, idx, spindet,                        &
          N_int, size_buffer,                          &
          singles, doubles,                            &
          n_singles, n_doubles)                        &
        bind(C, name="get_all_spin_singles_and_doubles")
      import :: c_int64_t, c_int
      implicit none
      integer(c_int64_t), intent(in)  :: buffer(*)
      integer(c_int),     intent(in)  :: idx(*)
      integer(c_int64_t), intent(in)  :: spindet(*)
      integer(c_int),     value       :: N_int
      integer(c_int),     value       :: size_buffer
      integer(c_int),     intent(out) :: singles(*)
      integer(c_int),     intent(out) :: doubles(*)
      integer(c_int),     intent(out) :: n_singles
      integer(c_int),     intent(out) :: n_doubles
    end subroutine

    subroutine c_get_all_spin_singles_raw_1(        &
          buffer, idx, spindet,                     &
          size_buffer,                              &
          singles, n_singles)                       &
        bind(C, name="get_all_spin_singles_1")
      import :: c_int64_t, c_int
      implicit none
      integer(c_int64_t), intent(in)  :: buffer(*)
      integer(c_int),     intent(in)  :: idx(*)
      integer(c_int64_t), intent(in)  :: spindet(1)
      integer(c_int),     value       :: size_buffer
      integer(c_int),     intent(out) :: singles(*)
      integer(c_int),     intent(out) :: n_singles
    end subroutine

    subroutine c_get_all_spin_doubles_raw_1(        &
          buffer, idx, spindet,                     &
          size_buffer,                              &
          doubles, n_doubles)                       &
        bind(C, name="get_all_spin_doubles_1")
      import :: c_int64_t, c_int
      implicit none
      integer(c_int64_t), intent(in)  :: buffer(*)
      integer(c_int),     intent(in)  :: idx(*)
      integer(c_int64_t), intent(in)  :: spindet(1)
      integer(c_int),     value       :: size_buffer
      integer(c_int),     intent(out) :: doubles(*)
      integer(c_int),     intent(out) :: n_doubles
    end subroutine

    subroutine c_get_all_spin_singles_and_doubles_raw_1( &
          buffer, idx, spindet,                        &
          size_buffer,                          &
          singles, doubles,                            &
          n_singles, n_doubles)                        &
        bind(C, name="get_all_spin_singles_and_doubles_1")
      import :: c_int64_t, c_int
      implicit none
      integer(c_int64_t), intent(in)  :: buffer(*)
      integer(c_int),     intent(in)  :: idx(*)
      integer(c_int64_t), intent(in)  :: spindet(1)
      integer(c_int),     value       :: size_buffer
      integer(c_int),     intent(out) :: singles(*)
      integer(c_int),     intent(out) :: doubles(*)
      integer(c_int),     intent(out) :: n_singles
      integer(c_int),     intent(out) :: n_doubles
    end subroutine


    subroutine c_get_all_spin_singles_raw_2(        &
          buffer, idx, spindet,                     &
          size_buffer,                              &
          singles, n_singles)                       &
        bind(C, name="get_all_spin_singles_2")
      import :: c_int64_t, c_int
      implicit none
      integer(c_int64_t), intent(in)  :: buffer(*)
      integer(c_int),     intent(in)  :: idx(*)
      integer(c_int64_t), intent(in)  :: spindet(*)
      integer(c_int),     value       :: size_buffer
      integer(c_int),     intent(out) :: singles(*)
      integer(c_int),     intent(out) :: n_singles
    end subroutine

    subroutine c_get_all_spin_doubles_raw_2(        &
          buffer, idx, spindet,                     &
          size_buffer,                              &
          doubles, n_doubles)                       &
        bind(C, name="get_all_spin_doubles_2")
      import :: c_int64_t, c_int
      implicit none
      integer(c_int64_t), intent(in)  :: buffer(*)
      integer(c_int),     intent(in)  :: idx(*)
      integer(c_int64_t), intent(in)  :: spindet(*)
      integer(c_int),     value       :: size_buffer
      integer(c_int),     intent(out) :: doubles(*)
      integer(c_int),     intent(out) :: n_doubles
    end subroutine

    subroutine c_get_all_spin_singles_and_doubles_raw_2( &
          buffer, idx, spindet,                        &
          size_buffer,                                 &
          singles, doubles,                            &
          n_singles, n_doubles)                        &
        bind(C, name="get_all_spin_singles_and_doubles_2")
      import :: c_int64_t, c_int
      implicit none
      integer(c_int64_t), intent(in)  :: buffer(*)
      integer(c_int),     intent(in)  :: idx(*)
      integer(c_int64_t), intent(in)  :: spindet(*)
      integer(c_int),     value       :: size_buffer
      integer(c_int),     intent(out) :: singles(*)
      integer(c_int),     intent(out) :: doubles(*)
      integer(c_int),     intent(out) :: n_singles
      integer(c_int),     intent(out) :: n_doubles
    end subroutine


    subroutine c_get_all_spin_singles_raw_3(        &
          buffer, idx, spindet,                     &
          size_buffer,                              &
          singles, n_singles)                       &
        bind(C, name="get_all_spin_singles_3")
      import :: c_int64_t, c_int
      implicit none
      integer(c_int64_t), intent(in)  :: buffer(*)
      integer(c_int),     intent(in)  :: idx(*)
      integer(c_int64_t), intent(in)  :: spindet(*)
      integer(c_int),     value       :: size_buffer
      integer(c_int),     intent(out) :: singles(*)
      integer(c_int),     intent(out) :: n_singles
    end subroutine

    subroutine c_get_all_spin_doubles_raw_3(        &
          buffer, idx, spindet,                     &
          size_buffer,                              &
          doubles, n_doubles)                       &
        bind(C, name="get_all_spin_doubles_3")
      import :: c_int64_t, c_int
      implicit none
      integer(c_int64_t), intent(in)  :: buffer(*)
      integer(c_int),     intent(in)  :: idx(*)
      integer(c_int64_t), intent(in)  :: spindet(*)
      integer(c_int),     value       :: size_buffer
      integer(c_int),     intent(out) :: doubles(*)
      integer(c_int),     intent(out) :: n_doubles
    end subroutine

    subroutine c_get_all_spin_singles_and_doubles_raw_3( &
          buffer, idx, spindet,                        &
          size_buffer,                                 &
          singles, doubles,                            &
          n_singles, n_doubles)                        &
        bind(C, name="get_all_spin_singles_and_doubles_3")
      import :: c_int64_t, c_int
      implicit none
      integer(c_int64_t), intent(in)  :: buffer(*)
      integer(c_int),     intent(in)  :: idx(*)
      integer(c_int64_t), intent(in)  :: spindet(*)
      integer(c_int),     value       :: size_buffer
      integer(c_int),     intent(out) :: singles(*)
      integer(c_int),     intent(out) :: doubles(*)
      integer(c_int),     intent(out) :: n_singles
      integer(c_int),     intent(out) :: n_doubles
    end subroutine


    subroutine c_get_all_spin_singles_raw_4(        &
          buffer, idx, spindet,                     &
          size_buffer,                              &
          singles, n_singles)                       &
        bind(C, name="get_all_spin_singles_4")
      import :: c_int64_t, c_int
      implicit none
      integer(c_int64_t), intent(in)  :: buffer(*)
      integer(c_int),     intent(in)  :: idx(*)
      integer(c_int64_t), intent(in)  :: spindet(*)
      integer(c_int),     value       :: size_buffer
      integer(c_int),     intent(out) :: singles(*)
      integer(c_int),     intent(out) :: n_singles
    end subroutine

    subroutine c_get_all_spin_doubles_raw_4(        &
          buffer, idx, spindet,                     &
          size_buffer,                              &
          doubles, n_doubles)                       &
        bind(C, name="get_all_spin_doubles_4")
      import :: c_int64_t, c_int
      implicit none
      integer(c_int64_t), intent(in)  :: buffer(*)
      integer(c_int),     intent(in)  :: idx(*)
      integer(c_int64_t), intent(in)  :: spindet(*)
      integer(c_int),     value       :: size_buffer
      integer(c_int),     intent(out) :: doubles(*)
      integer(c_int),     intent(out) :: n_doubles
    end subroutine

    subroutine c_get_all_spin_singles_and_doubles_raw_4( &
          buffer, idx, spindet,                        &
          size_buffer,                                 &
          singles, doubles,                            &
          n_singles, n_doubles)                        &
        bind(C, name="get_all_spin_singles_and_doubles_4")
      import :: c_int64_t, c_int
      implicit none
      integer(c_int64_t), intent(in)  :: buffer(*)
      integer(c_int),     intent(in)  :: idx(*)
      integer(c_int64_t), intent(in)  :: spindet(*)
      integer(c_int),     value       :: size_buffer
      integer(c_int),     intent(out) :: singles(*)
      integer(c_int),     intent(out) :: doubles(*)
      integer(c_int),     intent(out) :: n_singles
      integer(c_int),     intent(out) :: n_doubles
    end subroutine


  end interface

contains

  ! ------------------------------------------------------------------
  ! get_all_spin_singles
  !
  ! Returns the indices of all single excitations (degree == 2) in
  ! buffer relative to spindet.
  !
  ! Arguments:
  !   buffer(Nint, size_buffer)  [in]  – packed determinant list
  !   idx(size_buffer)           [in]  – caller-supplied indices
  !   spindet(Nint)              [in]  – reference spin-determinant
  !   Nint                       [in]  – number of 64-bit words per det
  !   size_buffer                [in]  – number of determinants
  !   singles(size_buffer)       [out] – indices of single excitations
  !   n_singles                  [out] – number of single excitations
  ! ------------------------------------------------------------------
  subroutine get_all_spin_singles(                 &
        buffer, idx, spindet,                        &
        Nint, size_buffer,                           &
        singles, n_singles)
    implicit none
    integer,            intent(in)  :: Nint, size_buffer
    integer(c_int64_t), intent(in)  :: buffer(Nint, size_buffer)
    integer,            intent(in)  :: idx(size_buffer)
    integer(c_int64_t), intent(in)  :: spindet(Nint)
    integer,            intent(out) :: singles(size_buffer)
    integer,            intent(out) :: n_singles

    call c_get_all_spin_singles_raw(           &
          buffer, idx, spindet,                &
          int(Nint, c_int),                    &
          int(size_buffer, c_int),             &
          singles, n_singles)
  end subroutine

  ! ------------------------------------------------------------------
  ! get_all_spin_doubles
  !
  ! Returns the indices of all double excitations (degree == 4) in
  ! buffer relative to spindet.
  !
  ! Arguments:
  !   buffer(Nint, size_buffer)  [in]  – packed determinant list
  !   idx(size_buffer)           [in]  – caller-supplied indices
  !   spindet(Nint)              [in]  – reference spin-determinant
  !   Nint                       [in]  – number of 64-bit words per det
  !   size_buffer                [in]  – number of determinants
  !   doubles(size_buffer)       [out] – indices of double excitations
  !   n_doubles                  [out] – number of double excitations
  ! ------------------------------------------------------------------
  subroutine get_all_spin_doubles(                 &
        buffer, idx, spindet,                        &
        Nint, size_buffer,                           &
        doubles, n_doubles)
    implicit none
    integer,            intent(in)  :: Nint, size_buffer
    integer(c_int64_t), intent(in)  :: buffer(Nint, size_buffer)
    integer,            intent(in)  :: idx(size_buffer)
    integer(c_int64_t), intent(in)  :: spindet(Nint)
    integer,            intent(out) :: doubles(size_buffer)
    integer,            intent(out) :: n_doubles

    call c_get_all_spin_doubles_raw(           &
          buffer, idx, spindet,                &
          int(Nint, c_int),                    &
          int(size_buffer, c_int),             &
          doubles, n_doubles)
  end subroutine

  ! ------------------------------------------------------------------
  ! get_all_spin_singles_and_doubles
  !
  ! Returns the indices of all single (degree == 2) and double
  ! (degree == 4) excitations in buffer relative to spindet.
  !
  ! Arguments:
  !   buffer(Nint, size_buffer)  [in]  – packed determinant list
  !   idx(size_buffer)           [in]  – caller-supplied indices
  !   spindet(Nint)              [in]  – reference spin-determinant
  !   Nint                       [in]  – number of 64-bit words per det
  !   size_buffer                [in]  – number of determinants
  !   singles(size_buffer)       [out] – indices of single excitations
  !   doubles(size_buffer)       [out] – indices of double excitations
  !   n_singles                  [out] – number of single excitations
  !   n_doubles                  [out] – number of double excitations
  ! ------------------------------------------------------------------
  subroutine get_all_spin_singles_and_doubles(     &
        buffer, idx, spindet,                        &
        Nint, size_buffer,                           &
        singles, doubles,                            &
        n_singles, n_doubles)
    implicit none
    integer,            intent(in)  :: Nint, size_buffer
    integer(c_int64_t), intent(in)  :: buffer(Nint, size_buffer)
    integer,            intent(in)  :: idx(size_buffer)
    integer(c_int64_t), intent(in)  :: spindet(Nint)
    integer,            intent(out) :: singles(size_buffer)
    integer,            intent(out) :: doubles(size_buffer)
    integer,            intent(out) :: n_singles
    integer,            intent(out) :: n_doubles

    call c_get_all_spin_singles_and_doubles_raw( &
          buffer, idx, spindet,                  &
          int(Nint, c_int),                      &
          int(size_buffer, c_int),               &
          singles, doubles,                      &
          n_singles, n_doubles)
  end subroutine

end module spin_excitations

