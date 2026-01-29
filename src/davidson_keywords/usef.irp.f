use bitmasks
use f77_zmq


! ---

BEGIN_PROVIDER [ double precision, threshold_davidson_pt2 ]
 implicit none
 BEGIN_DOC
 ! Threshold of Davidson's algorithm, using PT2 as a guide
 END_DOC
 threshold_davidson_pt2 = threshold_davidson

END_PROVIDER

! ---

