use bitmasks

BEGIN_PROVIDER [real*8, hessmat, (nMonoEx,nMonoEx)]
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
    write(6,*) ' providing Hessian matrix hessmat '
    write(6,*) '  nMonoEx = ',nMonoEx
  endif
  
  do indx=1,nMonoEx
    do jndx=1,nMonoEx
      hessmat(indx,jndx)=0.D0
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
      hessmat(indx,jndx)=res
      hessmat(jndx,indx)=res
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

BEGIN_PROVIDER [real*8, hessmat2, (nMonoEx,nMonoEx)]
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
    write(6,*) ' providing Hessian matrix hessmat2 '
    write(6,*) '  nMonoEx = ',nMonoEx
  endif
  
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP SHARED(hessmat2,n_core_inact_orb,n_act_orb,n_virt_orb,nMonoEx) &
  !$OMP PRIVATE(i,indx,jndx,j,ustart,t,u,a,bstart,indx_shift)

  !$OMP DO
  do i=1,n_core_inact_orb
    do t=1,n_act_orb
      indx = t + (i-1)*n_act_orb
      jndx=indx
      do j=i,n_core_inact_orb
        if (i.eq.j) then
          ustart=t
        else
          ustart=1
        end if
        do u=ustart,n_act_orb
          hessmat2(jndx,indx)=hessmat_itju(i,t,j,u)
          jndx+=1
        end do
      end do
      do j=1,n_core_inact_orb
        do a=1,n_virt_orb
          hessmat2(jndx,indx)=hessmat_itja(i,t,j,a)
          jndx+=1
        end do
      end do
      do u=1,n_act_orb
        do a=1,n_virt_orb
          hessmat2(jndx,indx)=hessmat_itua(i,t,u,a)
          jndx+=1
        end do
      end do
    end do
  end do
  !$OMP END DO NOWAIT
  
  indx_shift = n_core_inact_orb*n_act_orb
  !$OMP DO
  do a=1,n_virt_orb
    do i=1,n_core_inact_orb
      indx = a + (i-1)*n_virt_orb + indx_shift
      jndx=indx
      do j=i,n_core_inact_orb
        if (i.eq.j) then
          bstart=a
        else
          bstart=1
        end if
        do b=bstart,n_virt_orb
          hessmat2(jndx,indx)=hessmat_iajb(i,a,j,b)
          jndx+=1
        end do
      end do
      do t=1,n_act_orb
        do b=1,n_virt_orb
          hessmat2(jndx,indx)=hessmat_iatb(i,a,t,b)
          jndx+=1
        end do
      end do
    end do
  end do
  !$OMP END DO NOWAIT
  
  indx_shift += n_core_inact_orb*n_virt_orb
  !$OMP DO 
  do a=1,n_virt_orb
    do t=1,n_act_orb
      indx = a + (t-1)*n_virt_orb + indx_shift
      jndx=indx
      do u=t,n_act_orb
        if (t.eq.u) then
          bstart=a
        else
          bstart=1
        end if
        do b=bstart,n_virt_orb
          hessmat2(jndx,indx)=hessmat_taub(t,a,u,b)
          jndx+=1
        end do
      end do
    end do
  end do
  !$OMP END DO 

  !$OMP END PARALLEL
  
  do jndx=1,nMonoEx
    do indx=1,jndx-1
      hessmat2(indx,jndx) = hessmat2(jndx,indx)
    enddo
  enddo

  
END_PROVIDER

real*8 function hessmat_itju(i,t,j,u)
  BEGIN_DOC
  ! the orbital hessian for core/inactive -> active, core/inactive -> active
  ! i, t, j, u are list indices, the corresponding orbitals are ii,tt,jj,uu
  !
  ! we assume natural orbitals
  END_DOC
  implicit none
  integer                        :: i,t,j,u,ii,tt,uu,v,vv,x,xx,y,jj
  real*8                         :: term,t2
  
  ii=list_core_inact(i)
  tt=list_act(t)
  if (i.eq.j) then
    if (t.eq.u) then
      ! diagonal element
      term=occnum(tt)*Fipq(ii,ii)+2.D0*(Fipq(tt,tt)+Fapq(tt,tt))     &
          -2.D0*(Fipq(ii,ii)+Fapq(ii,ii))
      term+=2.D0*(3.D0*bielec_pxxq_no(tt,i,i,tt)-bielec_pqxx_no(tt,tt,i,i))
      term-=2.D0*occnum(tt)*(3.D0*bielec_pxxq_no(tt,i,i,tt)             &
          -bielec_pqxx_no(tt,tt,i,i))
      term-=occnum(tt)*Fipq(tt,tt)
      do v=1,n_act_orb
        vv=list_act(v)
        do x=1,n_act_orb
          xx=list_act(x)
          term+=2.D0*(P0tuvx_no(t,t,v,x)*bielec_pqxx_no(vv,xx,i,i)      &
              +(P0tuvx_no(t,x,v,t)+P0tuvx_no(t,x,t,v))*              &
              bielec_pxxq_no(vv,i,i,xx))
          do y=1,n_act_orb
            term-=2.D0*P0tuvx_no(t,v,x,y)*bielecCI_no(t,v,y,xx)
          end do
        end do
      end do
    else
      ! it/iu, t != u
      uu=list_act(u)
      term=2.D0*(Fipq(tt,uu)+Fapq(tt,uu))
      term+=2.D0*(4.D0*bielec_PxxQ_no(tt,i,j,uu)-bielec_PxxQ_no(uu,i,j,tt) &
          -bielec_PQxx_no(tt,uu,i,j))
      term-=occnum(tt)*Fipq(uu,tt)
      term-=(occnum(tt)+occnum(uu))                                  &
          *(3.D0*bielec_PxxQ_no(tt,i,i,uu)-bielec_PQxx_no(uu,tt,i,i))
      do v=1,n_act_orb
        vv=list_act(v)
        !         term-=D0tu(u,v)*Fipq(tt,vv)   ! published, but inverting t and u seems more correct
        do x=1,n_act_orb
          xx=list_act(x)
          term+=2.D0*(P0tuvx_no(u,t,v,x)*bielec_pqxx_no(vv,xx,i,i)      &
              +(P0tuvx_no(u,x,v,t)+P0tuvx_no(u,x,t,v))               &
              *bielec_pxxq_no(vv,i,i,xx))
          do y=1,n_act_orb
            term-=2.D0*P0tuvx_no(t,v,x,y)*bielecCI_no(u,v,y,xx)
          end do
        end do
      end do
    end if
  else
    ! it/ju
    jj=list_core_inact(j)
    uu=list_act(u)
    if (t.eq.u) then
      term=occnum(tt)*Fipq(ii,jj)
      term-=2.D0*(Fipq(ii,jj)+Fapq(ii,jj))
    else
      term=0.D0
    end if
    term+=2.D0*(4.D0*bielec_PxxQ_no(tt,i,j,uu)-bielec_PxxQ_no(uu,i,j,tt)   &
        -bielec_PQxx_no(tt,uu,i,j))
    term-=(occnum(tt)+occnum(uu))*                                   &
        (4.D0*bielec_PxxQ_no(tt,i,j,uu)-bielec_PxxQ_no(uu,i,j,tt)          &
        -bielec_PQxx_no(uu,tt,i,j))
    do v=1,n_act_orb
      vv=list_act(v)
      do x=1,n_act_orb
        xx=list_act(x)
        term+=2.D0*(P0tuvx_no(u,t,v,x)*bielec_pqxx_no(vv,xx,i,j)        &
            +(P0tuvx_no(u,x,v,t)+P0tuvx_no(u,x,t,v))                 &
            *bielec_pxxq_no(vv,i,j,xx))
      end do
    end do
  end if
  
  term*=2.D0
  hessmat_itju=term
  
end function hessmat_itju

real*8 function hessmat_itja(i,t,j,a)
  BEGIN_DOC
  ! the orbital hessian for core/inactive -> active, core/inactive -> virtual
  END_DOC
  implicit none
  integer                        :: i,t,j,a,ii,tt,jj,aa,v,vv,x,y
  real*8                         :: term
  
  ! it/ja
  ii=list_core_inact(i)
  tt=list_act(t)
  jj=list_core_inact(j)
  aa=list_virt(a)
  term=2.D0*(4.D0*bielec_pxxq_no(aa,j,i,tt)                             &
      -bielec_pqxx_no(aa,tt,i,j) -bielec_pxxq_no(aa,i,j,tt))
  term-=occnum(tt)*(4.D0*bielec_pxxq_no(aa,j,i,tt)                      &
      -bielec_pqxx_no(aa,tt,i,j) -bielec_pxxq_no(aa,i,j,tt))
  if (i.eq.j) then
    term+=2.D0*(Fipq(aa,tt)+Fapq(aa,tt))
    term-=0.5D0*occnum(tt)*Fipq(aa,tt)
    do v=1,n_act_orb
      do x=1,n_act_orb
        do y=1,n_act_orb
          term-=P0tuvx_no(t,v,x,y)*bielecCI_no(x,y,v,aa)
        end do
      end do
    end do
  end if
  term*=2.D0
  hessmat_itja=term
  
end function hessmat_itja

real*8 function hessmat_itua(i,t,u,a)
  BEGIN_DOC
  ! the orbital hessian for core/inactive -> active, active -> virtual
  END_DOC
  implicit none
  integer                        :: i,t,u,a,ii,tt,uu,aa,v,vv,x,xx,u3,t3,v3
  real*8                         :: term
  
  ii=list_core_inact(i)
  tt=list_act(t)
  t3=t+n_core_inact_orb
  uu=list_act(u)
  u3=u+n_core_inact_orb
  aa=list_virt(a)
  if (t.eq.u) then
    term=-occnum(tt)*Fipq(aa,ii)
  else
    term=0.D0
  end if
  term-=occnum(uu)*(bielec_pqxx_no(aa,ii,t3,u3)-4.D0*bielec_pqxx_no(aa,uu,t3,i)&
      +bielec_pxxq_no(aa,t3,u3,ii))
  do v=1,n_act_orb
    vv=list_act(v)
    v3=v+n_core_inact_orb
    do x=1,n_act_orb
      integer                        :: x3
      xx=list_act(x)
      x3=x+n_core_inact_orb
      term-=2.D0*(P0tuvx_no(t,u,v,x)*bielec_pqxx_no(aa,ii,v3,x3)        &
          +(P0tuvx_no(t,v,u,x)+P0tuvx_no(t,v,x,u))                   &
          *bielec_pqxx_no(aa,xx,v3,i))
    end do
  end do
  if (t.eq.u) then
    term+=Fipq(aa,ii)+Fapq(aa,ii)
  end if
  term*=2.D0
  hessmat_itua=term
  
end function hessmat_itua

real*8 function hessmat_iajb(i,a,j,b)
  BEGIN_DOC
  ! the orbital hessian for core/inactive -> virtual, core/inactive -> virtual
  END_DOC
  implicit none
  integer                        :: i,a,j,b,ii,aa,jj,bb
  real*8                         :: term
  
  ii=list_core_inact(i)
  aa=list_virt(a)
  if (i.eq.j) then
    if (a.eq.b) then
      ! ia/ia
      term=2.D0*(Fipq(aa,aa)+Fapq(aa,aa)-Fipq(ii,ii)-Fapq(ii,ii))
      term+=2.D0*(3.D0*bielec_pxxq_no(aa,i,i,aa)-bielec_pqxx_no(aa,aa,i,i))
    else
      bb=list_virt(b)
      ! ia/ib
      term=2.D0*(Fipq(aa,bb)+Fapq(aa,bb))
      term+=2.D0*(3.D0*bielec_pxxq_no(aa,i,i,bb)-bielec_pqxx_no(aa,bb,i,i))
    end if
  else
    ! ia/jb
    jj=list_core_inact(j)
    bb=list_virt(b)
    term=2.D0*(4.D0*bielec_pxxq_no(aa,i,j,bb)-bielec_pqxx_no(aa,bb,i,j)    &
        -bielec_pxxq_no(aa,j,i,bb))
    if (a.eq.b) then
      term-=2.D0*(Fipq(ii,jj)+Fapq(ii,jj))
    end if
  end if
  term*=2.D0
  hessmat_iajb=term
  
end function hessmat_iajb

real*8 function hessmat_iatb(i,a,t,b)
  BEGIN_DOC
  ! the orbital hessian for core/inactive -> virtual, active -> virtual
  END_DOC
  implicit none
  integer                        :: i,a,t,b,ii,aa,tt,bb,v,vv,x,y,v3,t3
  real*8                         :: term
  
  ii=list_core_inact(i)
  aa=list_virt(a)
  tt=list_act(t)
  bb=list_virt(b)
  t3=t+n_core_inact_orb
  term=occnum(tt)*(4.D0*bielec_pxxq_no(aa,i,t3,bb)-bielec_pxxq_no(aa,t3,i,bb)&
      -bielec_pqxx_no(aa,bb,i,t3))
  if (a.eq.b) then
    term-=Fipq(tt,ii)+Fapq(tt,ii)
    term-=0.5D0*occnum(tt)*Fipq(tt,ii)
    do v=1,n_act_orb
      do x=1,n_act_orb
        do y=1,n_act_orb
          term-=P0tuvx_no(t,v,x,y)*bielecCI_no(x,y,v,ii)
        end do
      end do
    end do
  end if
  term*=2.D0
  hessmat_iatb=term
  
end function hessmat_iatb

real*8 function hessmat_taub(t,a,u,b)
  BEGIN_DOC
  ! the orbital hessian for act->virt,act->virt
  END_DOC
  implicit none
  integer                        :: t,a,u,b,tt,aa,uu,bb,v,vv,x,xx,y
  integer                        :: v3,x3
  real*8                         :: term,t1,t2,t3
  
  double precision,allocatable :: P0tuvx_no_t(:,:,:)
  double precision :: bielec_pqxx_no_2(n_act_orb,n_act_orb)
  double precision :: bielec_pxxq_no_2(n_act_orb,n_act_orb)
  tt=list_act(t)
  aa=list_virt(a)
  if (t == u) then
    if (a == b) then
      ! ta/ta
      t1=occnum(tt)*Fipq(aa,aa)
      t2=0.D0
      t3=0.D0
      t1-=occnum(tt)*Fipq(tt,tt)
      do x=1,n_act_orb
        xx=list_act(x)
        x3=x+n_core_inact_orb
        do v=1,n_act_orb
          vv=list_act(v)
          v3=v+n_core_inact_orb
          t2+=P0tuvx_no(t,t,v,x)*bielec_pqxx_no(aa,aa,v3,x3)
        end do
      end do
      do v=1,n_act_orb
        vv=list_act(v)
        v3=v+n_core_inact_orb
        do x=1,n_act_orb
          xx=list_act(x)
          x3=x+n_core_inact_orb
          t2+=(P0tuvx_no(t,x,v,t)+P0tuvx_no(t,x,t,v))*              &
              bielec_pxxq_no(aa,x3,v3,aa)
        end do
      end do
      do y=1,n_act_orb
        do x=1,n_act_orb
          xx=list_act(x)
          do v=1,n_act_orb
            t3-=P0tuvx_no(t,v,x,y)*bielecCI_no(t,v,y,xx)
          end do
        end do
      end do
      term=t1+2.d0*(t2+t3)
    else
      bb=list_virt(b)
      ! ta/tb b/=a
      term=0.5d0*occnum(tt)*Fipq(aa,bb)
      do x=1,n_act_orb
        xx=list_act(x)
          x3=x+n_core_inact_orb
        do v=1,n_act_orb
          vv=list_act(v)
          v3=v+n_core_inact_orb
          term = term + P0tuvx_no(t,t,v,x)*bielec_pqxx_no(aa,bb,v3,x3)
        end do
      end do
      do v=1,n_act_orb
        vv=list_act(v)
        v3=v+n_core_inact_orb
        do x=1,n_act_orb
          xx=list_act(x)
          x3=x+n_core_inact_orb
          term= term + (P0tuvx_no(t,x,v,t)+P0tuvx_no(t,x,t,v))               &
              *bielec_pxxq_no(aa,x3,v3,bb)
        end do
      end do
      term += term
    end if
  else
    ! ta/ub t/=u
    uu=list_act(u)
    bb=list_virt(b)
    allocate(P0tuvx_no_t(n_act_orb,n_act_orb,n_act_orb))
    P0tuvx_no_t(:,:,:) = P0tuvx_no(t,:,:,:)
    do x=1,n_act_orb
      x3=x+n_core_inact_orb
      do v=1,n_act_orb
        v3=v+n_core_inact_orb
        bielec_pqxx_no_2(v,x) = bielec_pqxx_no(aa,bb,v3,x3)
        bielec_pxxq_no_2(v,x) = bielec_pxxq_no(aa,v3,x3,bb)
      end do
    end do
    term=0.D0
    do x=1,n_act_orb
      do v=1,n_act_orb
        term += P0tuvx_no_t(u,v,x)*bielec_pqxx_no_2(v,x)
        term += bielec_pxxq_no_2(x,v) * (P0tuvx_no_t(x,v,u)+P0tuvx_no_t(x,u,v))
      end do
    end do
    term = 6.d0*term
    if (a.eq.b) then
      term-=0.5D0*(occnum(tt)*Fipq(uu,tt)+occnum(uu)*Fipq(tt,uu))
      do v=1,n_act_orb
        do y=1,n_act_orb
          do x=1,n_act_orb
            term-=P0tuvx_no_t(v,x,y)*bielecCI_no(x,y,v,uu)
            term-=P0tuvx_no(u,v,x,y)*bielecCI_no(x,y,v,tt)
          end do
        end do
      end do
    end if
    
  end if
  
  term*=2.D0
  hessmat_taub=term
  
end function hessmat_taub

BEGIN_PROVIDER [real*8, hessdiag, (nMonoEx)]
  BEGIN_DOC
  ! the diagonal of the Hessian, needed for the Davidson procedure
  END_DOC
  implicit none
  integer                        :: i,t,a,indx,indx_shift
  real*8                         :: hessmat_itju,hessmat_iajb,hessmat_taub
  
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP SHARED(hessdiag,n_core_inact_orb,n_act_orb,n_virt_orb,nMonoEx) &
  !$OMP PRIVATE(i,indx,t,a,indx_shift)

  !$OMP DO
  do i=1,n_core_inact_orb
    do t=1,n_act_orb
      indx = t + (i-1)*n_act_orb
      hessdiag(indx)=hessmat_itju(i,t,i,t)
    end do
  end do
  !$OMP END DO NOWAIT
  
  indx_shift = n_core_inact_orb*n_act_orb
  !$OMP DO 
  do a=1,n_virt_orb
    do i=1,n_core_inact_orb
      indx = a + (i-1)*n_virt_orb + indx_shift
      hessdiag(indx)=hessmat_iajb(i,a,i,a)
    end do
  end do
  !$OMP END DO NOWAIT
  
  indx_shift += n_core_inact_orb*n_virt_orb
  !$OMP DO 
  do a=1,n_virt_orb
    do t=1,n_act_orb
      indx = a + (t-1)*n_virt_orb + indx_shift
      hessdiag(indx)=hessmat_taub(t,a,t,a)
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  
END_PROVIDER
