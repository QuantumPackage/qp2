
BEGIN_PROVIDER [real*8, gradvec_old, (nMonoEx)]
  BEGIN_DOC
  ! calculate the orbital gradient <Psi| H E_pq |Psi> by hand, i.e. for
  ! each determinant I we determine the string E_pq |I> (alpha and beta
  ! separately) and generate <Psi|H E_pq |I>
  ! sum_I c_I <Psi|H E_pq |I> is then the pq component of the orbital
  ! gradient
  ! E_pq = a^+_pa_q + a^+_Pa_Q
  END_DOC
  implicit none
  integer                        :: ii,tt,aa,indx,ihole,ipart,istate
  real*8                         :: res
  
  do indx=1,nMonoEx
    ihole=excit(1,indx)
    ipart=excit(2,indx)
    call calc_grad_elem(ihole,ipart,res)
    gradvec_old(indx)=res
  end do
  
  real*8                         :: norm_grad
  norm_grad=0.d0
  do indx=1,nMonoEx
    norm_grad+=gradvec_old(indx)*gradvec_old(indx)
  end do
  norm_grad=sqrt(norm_grad)
  if (bavard) then
    write(6,*)
    write(6,*) ' Norm of the orbital gradient (via <0|EH|0>) : ', norm_grad
    write(6,*)
  endif
  
  
END_PROVIDER

subroutine calc_grad_elem(ihole,ipart,res)
  BEGIN_DOC
  ! eq 18 of Siegbahn et al, Physica Scripta 1980
  ! we calculate 2 <Psi| H E_pq | Psi>, q=hole, p=particle
  END_DOC
  implicit none
  integer                        :: ihole,ipart,mu,iii,ispin,ierr,nu,istate
  real*8                         :: res
  integer(bit_kind), allocatable :: det_mu(:,:),det_mu_ex(:,:)
  real*8                         :: i_H_psi_array(N_states),phase
  allocate(det_mu(N_int,2))
  allocate(det_mu_ex(N_int,2))
  
  res=0.D0
  
  do mu=1,n_det
    ! get the string of the determinant
    call det_extract(det_mu,mu,N_int)
    do ispin=1,2
      ! do the monoexcitation on it
      call det_copy(det_mu,det_mu_ex,N_int)
      call do_signed_mono_excitation(det_mu,det_mu_ex,nu             &
          ,ihole,ipart,ispin,phase,ierr)
      if (ierr.eq.1) then
        call i_H_psi(det_mu_ex,psi_det,psi_coef,N_int                &
            ,N_det,N_det,N_states,i_H_psi_array)
        do istate=1,N_states
          res+=i_H_psi_array(istate)*psi_coef(mu,istate)*phase
        end do
      end if
    end do
  end do
  
  ! state-averaged gradient
  res*=2.D0/dble(N_states)
  
end subroutine calc_grad_elem

