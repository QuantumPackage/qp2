use bitmasks

integer, parameter :: hole_     = 1
integer, parameter :: particle_ = 2
integer, parameter :: hole2_    = 3
integer, parameter :: particle2_= 4

BEGIN_PROVIDER [ integer, N_single_exc_bitmasks ]
 implicit none
 BEGIN_DOC
 ! Number of single excitation bitmasks
 END_DOC
 N_single_exc_bitmasks = 1
 !TODO : Read from input!
END_PROVIDER

BEGIN_PROVIDER [ integer(bit_kind), single_exc_bitmask, (N_int, 2, N_single_exc_bitmasks) ]
 implicit none
 BEGIN_DOC
 ! single_exc_bitmask(:,1,i) is the bitmask for holes
 !
 ! single_exc_bitmask(:,2,i) is the bitmask for particles
 !
 ! for a given couple of hole/particle excitations i.
 END_DOC

 single_exc_bitmask(:,hole_,1) = HF_bitmask(:,1)
 single_exc_bitmask(:,particle_,1) = not(HF_bitmask(:,2))
END_PROVIDER


BEGIN_PROVIDER [ integer, N_double_exc_bitmasks ]
 implicit none
 BEGIN_DOC
 ! Number of double excitation bitmasks
 END_DOC
 N_double_exc_bitmasks = 1
END_PROVIDER

BEGIN_PROVIDER [ integer(bit_kind), double_exc_bitmask, (N_int, 4, N_double_exc_bitmasks) ]
 implicit none
 BEGIN_DOC
 ! double_exc_bitmask(:,1,i) is the bitmask for holes of excitation 1
 !
 ! double_exc_bitmask(:,2,i) is the bitmask for particles of excitation 1
 !
 ! double_exc_bitmask(:,3,i) is the bitmask for holes of excitation 2
 !
 ! double_exc_bitmask(:,4,i) is the bitmask for particles of excitation 2
 !
 ! for a given couple of hole/particle excitations i.
 END_DOC

 double_exc_bitmask(:,hole_,1) = HF_bitmask(:,1)
 double_exc_bitmask(:,particle_,1) = not(HF_bitmask(:,2))
 double_exc_bitmask(:,hole2_,1) = HF_bitmask(:,1)
 double_exc_bitmask(:,particle2_,1) = not(HF_bitmask(:,2))

END_PROVIDER

