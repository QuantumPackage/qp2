 BEGIN_PROVIDER [ double precision, dressing_column_h, (N_det,N_states) ]
&BEGIN_PROVIDER [ double precision, dressing_column_s, (N_det,N_states) ]
&BEGIN_PROVIDER [ double precision, dressing_delta   , (N_det,N_states) ]
 implicit none
 BEGIN_DOC
 ! Null dressing vectors
 END_DOC
 dressing_column_h(:,:) = 0.d0
 dressing_column_s(:,:) = 0.d0
 dressing_delta   (:,:) = 0.d0
END_PROVIDER

