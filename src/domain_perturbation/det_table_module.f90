module det_table_module
! Sharded hash table determinant -> vector of weights, with one OpenMP
! lock per shard. The weights are the Taylor coefficients (jets) of an
! amplitude, in powers of (E - E_0), from 0 to njet.
  use omp_lib
  implicit none

  integer, parameter :: n_shard = 1024

  type det_shard
    integer(8) :: n = 0_8
    integer(8) :: cap = 0_8
    integer(8) :: nslot = 0_8
    integer(8), allocatable       :: key(:,:)
    double precision, allocatable :: w(:,:)
    integer(8), allocatable       :: slot(:)
  end type

  type det_table
    integer :: nkey = 0
    integer :: njet = 0
    type(det_shard), allocatable          :: sh(:)
    integer(omp_lock_kind), allocatable   :: lock(:)
  end type

contains

  pure integer(8) function det_hash(key,nkey)
    integer, intent(in)    :: nkey
    integer(8), intent(in) :: key(nkey)
    integer(8) :: z
    integer :: i
    z = 1234567891234567_8
    do i=1,nkey
      z = ieor(z,key(i))
      z = z * (-7046029254386353131_8)
      z = ieor(z,shiftr(z,29))
      z = z * (-4658895280553007687_8)
      z = ieor(z,shiftr(z,32))
    end do
    det_hash = z
  end function

  subroutine table_create(t,nkey,njet)
    type(det_table), intent(inout) :: t
    integer, intent(in) :: nkey, njet
    integer :: i
    t%nkey = nkey
    t%njet = njet
    allocate(t%sh(n_shard), t%lock(n_shard))
    do i=1,n_shard
      call omp_init_lock(t%lock(i))
    end do
  end subroutine

  subroutine table_free(t)
    type(det_table), intent(inout) :: t
    integer :: i
    if (.not.allocated(t%sh)) return
    do i=1,n_shard
      call omp_destroy_lock(t%lock(i))
    end do
    deallocate(t%sh, t%lock)
  end subroutine

  ! Moves the contents of a table to another one, with no copy
  subroutine move_alloc_table(a,b)
    type(det_table), intent(inout) :: a, b
    integer :: i
    do i=1,n_shard
      b%sh(i)%n = a%sh(i)%n
      b%sh(i)%cap = a%sh(i)%cap
      b%sh(i)%nslot = a%sh(i)%nslot
      call move_alloc(a%sh(i)%key,b%sh(i)%key)
      call move_alloc(a%sh(i)%w,b%sh(i)%w)
      call move_alloc(a%sh(i)%slot,b%sh(i)%slot)
      a%sh(i)%n = 0_8
      a%sh(i)%cap = 0_8
      a%sh(i)%nslot = 0_8
    end do
  end subroutine

  integer(8) function table_size(t)
    type(det_table), intent(in) :: t
    integer :: i
    table_size = 0_8
    do i=1,n_shard
      table_size = table_size + t%sh(i)%n
    end do
  end function

  subroutine shard_grow(s,nkey,njet)
    type(det_shard), intent(inout) :: s
    integer, intent(in) :: nkey, njet
    integer(8), allocatable       :: k2(:,:)
    double precision, allocatable :: w2(:,:)
    integer(8) :: newcap, i, p, h
    newcap = max(16_8, 2_8*s%cap)
    allocate(k2(nkey,newcap), w2(0:njet,newcap))
    if (s%n > 0_8) then
      k2(:,1:s%n) = s%key(:,1:s%n)
      w2(:,1:s%n) = s%w(:,1:s%n)
    end if
    call move_alloc(k2,s%key)
    call move_alloc(w2,s%w)
    s%cap = newcap
    s%nslot = 2_8*newcap
    if (allocated(s%slot)) deallocate(s%slot)
    allocate(s%slot(0:s%nslot-1))
    s%slot = 0_8
    do i=1,s%n
      h = shiftr(det_hash(s%key(1,i),nkey),10)
      p = iand(h, s%nslot-1_8)
      do while (s%slot(p) /= 0_8)
        p = iand(p+1_8, s%nslot-1_8)
      end do
      s%slot(p) = i
    end do
  end subroutine

  subroutine table_add(t,key,w)
    type(det_table), intent(inout) :: t
    integer(8), intent(in)         :: key(t%nkey)
    double precision, intent(in)   :: w(0:t%njet)
    integer(8) :: h, p, i
    integer    :: is
    h  = det_hash(key,t%nkey)
    is = int(iand(h,int(n_shard-1,8)))+1
    h  = shiftr(h,10)
    call omp_set_lock(t%lock(is))
    associate (s => t%sh(is))
      if (s%n >= s%cap) call shard_grow(s,t%nkey,t%njet)
      p = iand(h, s%nslot-1_8)
      do
        i = s%slot(p)
        if (i == 0_8) then
          s%n = s%n+1_8
          s%key(:,s%n) = key(:)
          s%w(:,s%n) = w(:)
          s%slot(p) = s%n
          exit
        else if (all(s%key(:,i) == key(:))) then
          s%w(:,i) = s%w(:,i) + w(:)
          exit
        end if
        p = iand(p+1_8, s%nslot-1_8)
      end do
    end associate
    call omp_unset_lock(t%lock(is))
  end subroutine

  ! Weight of order 0 of a determinant, 0.d0 if absent. Read only.
  double precision function table_get(t,key)
    type(det_table), intent(in) :: t
    integer(8), intent(in)      :: key(t%nkey)
    integer(8) :: h, p, i
    integer    :: is
    table_get = 0.d0
    h  = det_hash(key,t%nkey)
    is = int(iand(h,int(n_shard-1,8)))+1
    if (t%sh(is)%n == 0_8) return
    h  = shiftr(h,10)
    p = iand(h, t%sh(is)%nslot-1_8)
    do
      i = t%sh(is)%slot(p)
      if (i == 0_8) return
      if (all(t%sh(is)%key(:,i) == key(:))) then
        table_get = t%sh(is)%w(0,i)
        return
      end if
      p = iand(p+1_8, t%sh(is)%nslot-1_8)
    end do
  end function

  ! Mixing of an integer, to obtain independent streams of random numbers
  pure integer(8) function det_mix(x)
    integer(8), intent(in) :: x
    integer(8) :: z
    z = ieor(x,shiftr(x,30)) * (-4658895280553007687_8)
    z = ieor(z,shiftr(z,27)) * (-7723592293110705685_8)
    det_mix = ieor(z,shiftr(z,31))
  end function

  ! Reproducible random numbers, independent of the number of threads
  double precision function det_rnd(state)
    integer(8), intent(inout) :: state
    integer(8) :: z
    state = state + (-7046029254386353131_8)
    z = state
    z = ieor(z,shiftr(z,30)) * (-4658895280553007687_8)
    z = ieor(z,shiftr(z,27)) * (-7723592293110705685_8)
    z = ieor(z,shiftr(z,31))
    det_rnd = dble(shiftr(z,11)) * 2.d0**(-53)
  end function

end module
