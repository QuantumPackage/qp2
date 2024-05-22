BEGIN_PROVIDER [ logical, initialize_pt2_E0_denominator ]
 implicit none
 BEGIN_DOC
 ! If true, initialize pt2_E0_denominator
 END_DOC
 initialize_pt2_E0_denominator = .True.
END_PROVIDER

BEGIN_PROVIDER [ double precision, pt2_E0_denominator, (N_states) ]
 implicit none
 BEGIN_DOC
 ! E0 in the denominator of the PT2
 END_DOC
 integer :: i,j

  pt2_E0_denominator = eigval_right_tc_bi_orth

END_PROVIDER

