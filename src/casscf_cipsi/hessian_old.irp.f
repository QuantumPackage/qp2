
use bitmasks
BEGIN_PROVIDER [real*8, hessmat_old, (nMonoEx,nMonoEx)]
  BEGIN_DOC
  ! calculate the orbital hessian 2 <Psi| E_pq H E_rs |Psi>
  ! + <Psi| E_pq E_rs H |Psi> + <Psi| E_rs E_pq H |Psi> by hand,
  ! determinant per determinant, as for the gradient
  !
  ! we assume that we have natural active orbitals
  END_DOC
  implicit none
  integer                        :: indx,ihole,ipart
  integer                        :: jndx,jhole,jpart
  character*3                    :: iexc,jexc
  real*8                         :: res
  
  if (bavard) then
    write(6,*) ' providing Hessian matrix hessmat_old '
    write(6,*) '  nMonoEx = ',nMonoEx
  endif
  
  do indx=1,nMonoEx
    do jndx=1,nMonoEx
      hessmat_old(indx,jndx)=0.D0
    end do
  end do
  
  do indx=1,nMonoEx
    ihole=excit(1,indx)
    ipart=excit(2,indx)
    iexc=excit_class(indx)
    do jndx=indx,nMonoEx
      jhole=excit(1,jndx)
      jpart=excit(2,jndx)
      jexc=excit_class(jndx)
      call calc_hess_elem(ihole,ipart,jhole,jpart,res)
      hessmat_old(indx,jndx)=res
      hessmat_old(jndx,indx)=res
    end do
  end do
  
END_PROVIDER

subroutine calc_hess_elem(ihole,ipart,jhole,jpart,res)
  BEGIN_DOC
  ! eq 19 of Siegbahn et al, Physica Scripta 1980
  ! we calculate  2 <Psi| E_pq H E_rs |Psi>
  !  + <Psi| E_pq E_rs H |Psi> + <Psi| E_rs E_pq H |Psi>
  ! average over all states is performed.
  ! no transition between states.
  END_DOC
  implicit none
  integer                        :: ihole,ipart,ispin,mu,istate
  integer                        :: jhole,jpart,jspin
  integer                        :: mu_pq, mu_pqrs, mu_rs, mu_rspq, nu_rs,nu
  real*8                         :: res
  integer(bit_kind), allocatable :: det_mu(:,:)
  integer(bit_kind), allocatable :: det_nu(:,:)
  integer(bit_kind), allocatable :: det_mu_pq(:,:)
  integer(bit_kind), allocatable :: det_mu_rs(:,:)
  integer(bit_kind), allocatable :: det_nu_rs(:,:)
  integer(bit_kind), allocatable :: det_mu_pqrs(:,:)
  integer(bit_kind), allocatable :: det_mu_rspq(:,:)
  real*8                         :: i_H_psi_array(N_states),phase,phase2,phase3
  real*8                         :: i_H_j_element
  allocate(det_mu(N_int,2))
  allocate(det_nu(N_int,2))
  allocate(det_mu_pq(N_int,2))
  allocate(det_mu_rs(N_int,2))
  allocate(det_nu_rs(N_int,2))
  allocate(det_mu_pqrs(N_int,2))
  allocate(det_mu_rspq(N_int,2))
  integer                        :: mu_pq_possible
  integer                        :: mu_rs_possible
  integer                        :: nu_rs_possible
  integer                        :: mu_pqrs_possible
  integer                        :: mu_rspq_possible
  
  res=0.D0
  
  ! the terms <0|E E H |0>
  do mu=1,n_det
    ! get the string of the determinant
    call det_extract(det_mu,mu,N_int)
    do ispin=1,2
      ! do the monoexcitation pq on it
      call det_copy(det_mu,det_mu_pq,N_int)
      call do_signed_mono_excitation(det_mu,det_mu_pq,mu_pq          &
          ,ihole,ipart,ispin,phase,mu_pq_possible)
      if (mu_pq_possible.eq.1) then
        ! possible, but not necessarily in the list
        ! do the second excitation
        do jspin=1,2
          call det_copy(det_mu_pq,det_mu_pqrs,N_int)
          call do_signed_mono_excitation(det_mu_pq,det_mu_pqrs,mu_pqrs&
              ,jhole,jpart,jspin,phase2,mu_pqrs_possible)
          ! excitation possible
          if (mu_pqrs_possible.eq.1) then
            call i_H_psi(det_mu_pqrs,psi_det,psi_coef,N_int          &
                ,N_det,N_det,N_states,i_H_psi_array)
            do istate=1,N_states
              res+=i_H_psi_array(istate)*psi_coef(mu,istate)*phase*phase2
            end do
          end if
          ! try the de-excitation with opposite sign
          call det_copy(det_mu_pq,det_mu_pqrs,N_int)
          call do_signed_mono_excitation(det_mu_pq,det_mu_pqrs,mu_pqrs&
              ,jpart,jhole,jspin,phase2,mu_pqrs_possible)
          phase2=-phase2
          ! excitation possible
          if (mu_pqrs_possible.eq.1) then
            call i_H_psi(det_mu_pqrs,psi_det,psi_coef,N_int          &
                ,N_det,N_det,N_states,i_H_psi_array)
            do istate=1,N_states
              res+=i_H_psi_array(istate)*psi_coef(mu,istate)*phase*phase2
            end do
          end if
        end do
      end if
      ! exchange the notion of pq and rs
      ! do the monoexcitation rs on the initial determinant
      call det_copy(det_mu,det_mu_rs,N_int)
      call do_signed_mono_excitation(det_mu,det_mu_rs,mu_rs          &
          ,jhole,jpart,ispin,phase2,mu_rs_possible)
      if (mu_rs_possible.eq.1) then
        ! do the second excitation
        do jspin=1,2
          call det_copy(det_mu_rs,det_mu_rspq,N_int)
          call do_signed_mono_excitation(det_mu_rs,det_mu_rspq,mu_rspq&
              ,ihole,ipart,jspin,phase3,mu_rspq_possible)
          ! excitation possible (of course, the result is outside the CAS)
          if (mu_rspq_possible.eq.1) then
            call i_H_psi(det_mu_rspq,psi_det,psi_coef,N_int          &
                ,N_det,N_det,N_states,i_H_psi_array)
            do istate=1,N_states
              res+=i_H_psi_array(istate)*psi_coef(mu,istate)*phase2*phase3
            end do
          end if
          ! we may try the de-excitation, with opposite sign
          call det_copy(det_mu_rs,det_mu_rspq,N_int)
          call do_signed_mono_excitation(det_mu_rs,det_mu_rspq,mu_rspq&
              ,ipart,ihole,jspin,phase3,mu_rspq_possible)
          phase3=-phase3
          ! excitation possible (of course, the result is outside the CAS)
          if (mu_rspq_possible.eq.1) then
            call i_H_psi(det_mu_rspq,psi_det,psi_coef,N_int          &
                ,N_det,N_det,N_states,i_H_psi_array)
            do istate=1,N_states
              res+=i_H_psi_array(istate)*psi_coef(mu,istate)*phase2*phase3
            end do
          end if
        end do
      end if
      !
      ! the operator E H E, we have to do a double loop over the determinants
      !  we still have the determinant mu_pq and the phase in memory
      if (mu_pq_possible.eq.1) then
        do nu=1,N_det
          call det_extract(det_nu,nu,N_int)
          do jspin=1,2
            call det_copy(det_nu,det_nu_rs,N_int)
            call do_signed_mono_excitation(det_nu,det_nu_rs,nu_rs    &
                ,jhole,jpart,jspin,phase2,nu_rs_possible)
            ! excitation possible ?
            if (nu_rs_possible.eq.1) then
              call i_H_j(det_mu_pq,det_nu_rs,N_int,i_H_j_element)
              do istate=1,N_states
                res+=2.D0*i_H_j_element*psi_coef(mu,istate)          &
                    *psi_coef(nu,istate)*phase*phase2
              end do
            end if
          end do
        end do
      end if
    end do
  end do
  
  ! state-averaged Hessian
  res*=1.D0/dble(N_states)
  
end subroutine calc_hess_elem

BEGIN_PROVIDER [real*8, hessmat_peter, (nMonoEx,nMonoEx)]
  BEGIN_DOC
  ! explicit hessian matrix from density matrices and integrals
  ! of course, this will be used for a direct Davidson procedure later
  ! we will not store the matrix in real life
  ! formulas are broken down as functions for the 6 classes of matrix elements
  !
  END_DOC
  implicit none
  integer                        :: i,j,t,u,a,b,indx,jndx,bstart,ustart,indx_shift
  
  real*8                         :: hessmat_itju
  real*8                         :: hessmat_itja
  real*8                         :: hessmat_itua
  real*8                         :: hessmat_iajb
  real*8                         :: hessmat_iatb
  real*8                         :: hessmat_taub
  
  if (bavard) then
    write(6,*) ' providing Hessian matrix hessmat_peter '
    write(6,*) '  nMonoEx = ',nMonoEx
  endif
  provide mo_two_e_integrals_in_map
  
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP SHARED(hessmat_peter,n_core_inact_orb,n_act_orb,n_virt_orb,nMonoEx) &
  !$OMP PRIVATE(i,indx,jndx,j,ustart,t,u,a,bstart,indx_shift)

  !$OMP DO
  ! (DOUBLY OCCUPIED ---> ACT ) 
  do i=1,n_core_inact_orb
    do t=1,n_act_orb
      indx = t + (i-1)*n_act_orb
      jndx=indx
      ! (DOUBLY OCCUPIED ---> ACT )
      do j=i,n_core_inact_orb
        if (i.eq.j) then
          ustart=t
        else
          ustart=1
        end if
        do u=ustart,n_act_orb
          hessmat_peter(jndx,indx)=hessmat_itju(i,t,j,u)
          jndx+=1
        end do
      end do
      ! (DOUBLY OCCUPIED ---> VIRTUAL) 
      do j=1,n_core_inact_orb
        do a=1,n_virt_orb
          hessmat_peter(jndx,indx)=hessmat_itja(i,t,j,a)
          jndx+=1
        end do
      end do
      ! (ACTIVE ---> VIRTUAL) 
      do u=1,n_act_orb
        do a=1,n_virt_orb
          hessmat_peter(jndx,indx)=hessmat_itua(i,t,u,a)
          jndx+=1
        end do
      end do
    end do
  end do
  !$OMP END DO NOWAIT
  
  indx_shift = n_core_inact_orb*n_act_orb
  !$OMP DO
  ! (DOUBLY OCCUPIED ---> VIRTUAL) 
  do a=1,n_virt_orb
    do i=1,n_core_inact_orb
      indx = a + (i-1)*n_virt_orb + indx_shift
      jndx=indx
      ! (DOUBLY OCCUPIED ---> VIRTUAL) 
      do j=i,n_core_inact_orb
        if (i.eq.j) then
          bstart=a
        else
          bstart=1
        end if
        do b=bstart,n_virt_orb
          hessmat_peter(jndx,indx)=hessmat_iajb(i,a,j,b)
          jndx+=1
        end do
      end do
      ! (ACT ---> VIRTUAL) 
      do t=1,n_act_orb
        do b=1,n_virt_orb
          hessmat_peter(jndx,indx)=hessmat_iatb(i,a,t,b)
          jndx+=1
        end do
      end do
    end do
  end do
  !$OMP END DO NOWAIT
  
  indx_shift += n_core_inact_orb*n_virt_orb
  !$OMP DO 
  ! (ACT ---> VIRTUAL) 
  do a=1,n_virt_orb
    do t=1,n_act_orb
      indx = a + (t-1)*n_virt_orb + indx_shift
      jndx=indx
      ! (ACT ---> VIRTUAL) 
      do u=t,n_act_orb
        if (t.eq.u) then
          bstart=a
        else
          bstart=1
        end if
        do b=bstart,n_virt_orb
          hessmat_peter(jndx,indx)=hessmat_taub(t,a,u,b)
          jndx+=1
        end do
      end do
    end do
  end do
  !$OMP END DO 

  !$OMP END PARALLEL
  
  do jndx=1,nMonoEx
    do indx=1,jndx-1
      hessmat_peter(indx,jndx) = hessmat_peter(jndx,indx)
    enddo
  enddo

  
END_PROVIDER

