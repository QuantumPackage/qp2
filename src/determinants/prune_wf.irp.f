BEGIN_PROVIDER [ logical, pruned, (N_det) ]
 implicit none
 BEGIN_DOC
 ! True if determinant is removed by pruning
 END_DOC

 pruned(:) = .False.

 if (pruning == 0.d0) then
   return
 endif

 integer :: i,j,k,ndet_new,nsop_max
 double precision :: thr

 if (s2_eig) then

  nsop_max = max(1,int ( dble(N_occ_pattern) * (1.d0 - pruning) + 0.5d0 ))

  do i=1,N_det
    k = det_to_occ_pattern(i)
    pruned(i) = psi_occ_pattern_sorted_order_reverse(k) > nsop_max
  enddo

 else

  ndet_new = max(1,int( dble(N_det) * (1.d0 - pruning) + 0.5d0 ))
  thr = psi_average_norm_contrib_sorted(ndet_new)
  do i=1, N_det
    pruned(i) = psi_average_norm_contrib(i) < thr
  enddo

 endif

END_PROVIDER
