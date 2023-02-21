! Algorithm for the trust region

! step_in_trust_region:
! Computes the step in the trust region (delta)
! (automatically sets at the iteration 0 and which evolves during the
! process in function of the evolution of rho). The step is computing by
! constraining its norm with a lagrange multiplier.
! Since the calculation of the step is based on the Newton method, an
! estimation of the gain in energy is given using the Taylors series
! truncated at the second order (criterion_model).
! If (DABS(criterion-criterion_model) < 1d-12) then 
!   must_exit = .True.
! else 
!   must_exit = .False.

! This estimation of the gain in energy is used by
! is_step_cancel_trust_region to say if the step is accepted or cancelled. 

! If the step must be cancelled, the calculation restart from the same
! hessian and gradient and recomputes the step but in a smaller trust
! region and so on until the step is accepted. If the step is accepted
! the hessian and the gradient are recomputed to produce a new step.

! Example:


!  !### Initialization ###
!  delta = 0d0
!  nb_iter = 0 ! Must start at 0 !!!
!  rho = 0.5d0
!  not_converged = .True.
!
!  ! ### TODO ###
!  ! Compute the criterion before the loop
!  call #your_criterion(prev_criterion)
!
!  do while (not_converged)
!    ! ### TODO ## 
!    ! Call your gradient
!    ! Call you hessian
!    call #your_gradient(v_grad) (1D array)
!    call #your_hessian(H) (2D array) 
!
!    ! ### TODO ###
!    ! Diagonalization of the hessian 
!    call diagonalization_hessian(n,H,e_val,w)
!
!    cancel_step = .True. ! To enter in the loop just after 
!    ! Loop to Reduce the trust radius until the criterion decreases and rho >= thresh_rho
!    do while (cancel_step)
!
!      ! Hessian,gradient,Criterion -> x 
!      call trust_region_step_w_expected_e(tmp_n,W,e_val,v_grad,prev_criterion,rho,nb_iter,delta,criterion_model,tmp_x,must_exit)
!
!      if (must_exit) then
!        ! ### Message ###
!        ! if step_in_trust_region sets must_exit on true for numerical reasons
!        print*,'algo_trust1 sends the message : Exit'
!        !### exit ###
!      endif
!
!      !### TODO ###  
!      ! Compute x -> m_x
!      ! Compute m_x -> R
!      ! Apply R and keep the previous MOs...
!      ! Update/touch 
!      ! Compute the new criterion/energy -> criterion
!
!      call #your_routine_1D_to_2D_antisymmetric_array(x,m_x)
!      call #your_routine_2D_antisymmetric_array_to_rotation_matrix(m_x,R)
!      call #your_routine_to_apply_the_rotation_matrix(R,prev_mos)
!
!      TOUCH #your_variables
!      
!      call #your_criterion(criterion)
!
!      ! Criterion -> step accepted or rejected 
!      call trust_region_is_step_cancelled(nb_iter,prev_criterion, criterion, criterion_model,rho,cancel_step)
!
!      ! ### TODO ###
!      !if (cancel_step) then
!      ! Cancel the previous step (mo_coef = prev_mos if you keep them...)
!      !endif
!      #if (cancel_step) then
!        #mo_coef = prev_mos
!      #endif
!
!    enddo
!
!    !call save_mos() !### depend of the time for 1 iteration
!
!    ! To exit the external loop if must_exit = .True.
!    if (must_exit) then
!      !### exit ###
!    endif 
!
!    ! Step accepted, nb iteration + 1
!    nb_iter = nb_iter + 1
!
!    ! ### TODO ###
!    !if (###Conditions###) then
!    ! no_converged = .False.
!    !endif
!    #if (#your_conditions) then
!    #  not_converged = .False.
!    #endif
!
!  enddo



! Variables:

! Input:
! | n              | integer          | m*(m-1)/2                                                                      |
! | m              | integer          | number of mo in the mo_class                                                   |
! | H(n,n)         | double precision | Hessian                                                                        |
! | v_grad(n)      | double precision | Gradient                                                                       |
! | W(n,n)         | double precision | Eigenvectors of the hessian                                                    |
! | e_val(n)       | double precision | Eigenvalues of the hessian                                                     |
! | criterion      | double precision | Actual criterion                                                               |
! | prev_criterion | double precision | Value of the criterion before the first iteration/after the previous iteration |
! | rho            | double precision | Given by is_step_cancel_trus_region                                            |
! |                |                  | Agreement between the real function and the Taylor series (2nd order)          |
! | nb_iter        | integer          | Actual number of iterations                                                    |

! Input/output:
! | delta           | double precision | Radius of the trust region                                                      |

! Output:
! | criterion_model | double precision | Predicted criterion after the rotation                                          |
! | x(n)            | double precision | Step                                                                            |
! | must_exit       | logical          | If the program must exit the loop                                               |


subroutine trust_region_step_w_expected_e(n,H,W,e_val,v_grad,prev_criterion,rho,nb_iter,delta,criterion_model,x,must_exit)

  include 'pi.h'

  BEGIN_DOC
  ! Compute the step and the expected criterion/energy after the step
  END_DOC

  implicit none

  ! in
  integer, intent(in)             :: n, nb_iter
  double precision, intent(in)    :: H(n,n), W(n,n), v_grad(n)
  double precision, intent(in)    :: rho, prev_criterion

  ! inout
  double precision, intent(inout) :: delta, e_val(n)

  ! out
  double precision, intent(out)   :: criterion_model, x(n)
  logical, intent(out)            :: must_exit

  ! internal
  integer :: info

  must_exit = .False.
  
  call trust_region_step(n,nb_iter,v_grad,rho,e_val,W,x,delta)

  call trust_region_expected_e(n,v_grad,H,x,prev_criterion,criterion_model)

  ! exit if DABS(prev_criterion - criterion_model) < 1d-12
  if (DABS(prev_criterion - criterion_model) < thresh_model) then
    print*,''
    print*,'###############################################################################'
    print*,'DABS(prev_criterion - criterion_model) <', thresh_model, 'stop the trust region'
    print*,'###############################################################################'
    print*,''
    must_exit = .True.
  endif

  if (delta < thresh_delta) then
    print*,''
    print*,'##############################################'
    print*,'Delta <', thresh_delta, 'stop the trust region'
    print*,'##############################################'
    print*,''
    must_exit = .True.
  endif

  ! Add after the call to this subroutine, a statement:
  ! "if (must_exit) then
  !   exit
  ! endif"
  ! in order to exit the optimization loop

end subroutine



! Variables:

! Input:
! | nb_iter         | integer          | actual number of iterations                    |
! | prev_criterion  | double precision | criterion before the application of the step x |
! | criterion       | double precision | criterion after the application of the step x  |
! | criterion_model | double precision | predicted criterion after the application of x |

! Output:
! | rho         | double precision | Agreement between the predicted criterion and the real new criterion |
! | cancel_step | logical          | If the step must be cancelled                                        |


subroutine trust_region_is_step_cancelled(nb_iter,prev_criterion, criterion, criterion_model,rho,cancel_step)

  include 'pi.h'

  BEGIN_DOC
  ! Compute if the step should be cancelled
  END_DOC

  implicit none
 
  ! in
  double precision, intent(in)  :: prev_criterion, criterion, criterion_model
  
  ! inout
  integer, intent(inout)        :: nb_iter

  ! out
  logical, intent(out)          :: cancel_step
  double precision, intent(out) :: rho

  ! Computes rho
  call trust_region_rho(prev_criterion,criterion,criterion_model,rho)
  
  if (nb_iter == 0) then
    nb_iter = 1 ! in order to enable the change of delta if the first iteration is cancelled  
  endif

  ! If rho < thresh_rho -> give something in output to cancel the step
  if (rho >= thresh_rho) then !0.1d0) then
     ! The step is accepted
     cancel_step = .False.
  else
     ! The step is rejected
     cancel_step = .True.
     print*, '***********************'
     print*, 'Step cancel : rho <', thresh_rho
     print*, '***********************'
  endif
  
end subroutine
