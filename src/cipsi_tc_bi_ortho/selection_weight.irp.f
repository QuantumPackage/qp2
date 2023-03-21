BEGIN_PROVIDER [ double precision, pt2_match_weight, (N_states) ]
 implicit none
 BEGIN_DOC
 ! Weights adjusted along the selection to make the PT2 contributions
 ! of each state coincide.
 END_DOC
 pt2_match_weight(:) = 1.d0
END_PROVIDER



BEGIN_PROVIDER [ double precision, variance_match_weight, (N_states) ]
 implicit none
 BEGIN_DOC
 ! Weights adjusted along the selection to make the variances
 ! of each state coincide.
 END_DOC
 variance_match_weight(:) = 1.d0
END_PROVIDER



subroutine update_pt2_and_variance_weights(pt2_data, N_st)
  implicit none
  use selection_types
  BEGIN_DOC
! Updates the PT2- and Variance- matching weights.
  END_DOC
  integer, intent(in)          :: N_st
  type(pt2_type), intent(in)   :: pt2_data
  double precision             :: pt2(N_st)
  double precision             :: variance(N_st)

  double precision :: avg, element, dt, x
  integer          :: k
  pt2(:)      = pt2_data % pt2(:)
  variance(:) = pt2_data % variance(:)

  avg = sum(pt2(1:N_st)) / dble(N_st) + 1.d-32 ! Avoid future division by zero

  dt = 8.d0 !* selection_factor
  do k=1,N_st
    element = exp(dt*(pt2(k)/avg - 1.d0))
    element = min(2.0d0 , element)
    element = max(0.5d0 , element)
    pt2_match_weight(k) *= element
  enddo


  avg = sum(variance(1:N_st)) / dble(N_st) + 1.d-32 ! Avoid future division by zero

  do k=1,N_st
    element = exp(dt*(variance(k)/avg -1.d0))
    element = min(2.0d0 , element)
    element = max(0.5d0 , element)
    variance_match_weight(k) *= element
  enddo

  if (N_det < 100) then
    ! For tiny wave functions, weights are 1.d0
    pt2_match_weight(:) = 1.d0
    variance_match_weight(:) = 1.d0
  endif

  threshold_davidson_pt2 = min(1.d-6, &
     max(threshold_davidson, 1.e-1 * PT2_relative_error * minval(abs(pt2(1:N_states)))) )

  SOFT_TOUCH pt2_match_weight variance_match_weight threshold_davidson_pt2
end




BEGIN_PROVIDER [ double precision, selection_weight, (N_states) ]
   implicit none
   BEGIN_DOC
   ! Weights used in the selection criterion
   END_DOC
   select case (weight_selection)

     case (0)
      print *,  'Using input weights in selection'
      selection_weight(1:N_states) = c0_weight(1:N_states) * state_average_weight(1:N_states)

     case (1)
      print *,  'Using 1/c_max^2 weight in selection'
      selection_weight(1:N_states) = c0_weight(1:N_states)

     case (2)
      print *,  'Using pt2-matching weight in selection'
      selection_weight(1:N_states) = c0_weight(1:N_states) * pt2_match_weight(1:N_states)
      print *, '# PT2 weight ', real(pt2_match_weight(:),4)

     case (3)
      print *,  'Using variance-matching weight in selection'
      selection_weight(1:N_states) = c0_weight(1:N_states) * variance_match_weight(1:N_states)
      print *, '# var weight ', real(variance_match_weight(:),4)

     case (4)
      print *,  'Using variance- and pt2-matching weights in selection'
      selection_weight(1:N_states) = c0_weight(1:N_states) * sqrt(variance_match_weight(1:N_states) * pt2_match_weight(1:N_states))
      print *, '# PT2 weight ', real(pt2_match_weight(:),4)
      print *, '# var weight ', real(variance_match_weight(:),4)

     case (5)
      print *,  'Using variance-matching weight in selection'
      selection_weight(1:N_states) = c0_weight(1:N_states) * variance_match_weight(1:N_states)
      print *, '# var weight ', real(variance_match_weight(:),4)

     case (6)
      print *,  'Using CI coefficient-based selection'
      selection_weight(1:N_states) = c0_weight(1:N_states)

     case (7)
      print *,  'Input weights multiplied by variance- and pt2-matching'
      selection_weight(1:N_states) = c0_weight(1:N_states) * sqrt(variance_match_weight(1:N_states) * pt2_match_weight(1:N_states)) * state_average_weight(1:N_states)
      print *, '# PT2 weight ', real(pt2_match_weight(:),4)
      print *, '# var weight ', real(variance_match_weight(:),4)

     case (8)
      print *,  'Input weights multiplied by pt2-matching'
      selection_weight(1:N_states) = c0_weight(1:N_states) * pt2_match_weight(1:N_states) * state_average_weight(1:N_states)
      print *, '# PT2 weight ', real(pt2_match_weight(:),4)

     case (9)
      print *,  'Input weights multiplied by variance-matching'
      selection_weight(1:N_states) = c0_weight(1:N_states) * variance_match_weight(1:N_states) * state_average_weight(1:N_states)
      print *, '# var weight ', real(variance_match_weight(:),4)

    end select
     print *, '# Total weight ', real(selection_weight(:),4)

END_PROVIDER

