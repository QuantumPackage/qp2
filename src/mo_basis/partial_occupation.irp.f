BEGIN_PROVIDER [ double precision, partial_occupation, (mo_num) ]
 implicit none
 BEGIN_DOC
 ! Partial occupation numbers for optimal damping algorithm
 END_DOC
 partial_occupation = 0.d0
 integer  :: i
 do i=1,elec_beta_num
   partial_occupation(i) = 2.d0
 enddo
 do i=elec_beta_num+1, elec_alpha_num
   partial_occupation(i) = 1.d0
 enddo
END_PROVIDER


