use bitmasks

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
  
  tt=list_act(t)
  aa=list_virt(a)
  if (t == u) then
    if (a == b) then
      ! ta/ta
      t1=occnum(tt)*Fipq(aa,aa)
      t2=0.D0
      t3=0.D0
      t1-=occnum(tt)*Fipq(tt,tt)
      do v=1,n_act_orb
        vv=list_act(v)
        v3=v+n_core_inact_orb
        do x=1,n_act_orb
          xx=list_act(x)
          x3=x+n_core_inact_orb
          t2+=2.D0*(P0tuvx_no(t,t,v,x)*bielec_pqxx_no(aa,aa,v3,x3)      &
              +(P0tuvx_no(t,x,v,t)+P0tuvx_no(t,x,t,v))*              &
              bielec_pxxq_no(aa,x3,v3,aa))
          do y=1,n_act_orb
            t3-=2.D0*P0tuvx_no(t,v,x,y)*bielecCI_no(t,v,y,xx)
          end do
        end do
      end do
      term=t1+t2+t3
    else
      bb=list_virt(b)
      ! ta/tb b/=a
      term=occnum(tt)*Fipq(aa,bb)
      do v=1,n_act_orb
        vv=list_act(v)
        v3=v+n_core_inact_orb
        do x=1,n_act_orb
          xx=list_act(x)
          x3=x+n_core_inact_orb
          term+=2.D0*(P0tuvx_no(t,t,v,x)*bielec_pqxx_no(aa,bb,v3,x3)    &
              +(P0tuvx_no(t,x,v,t)+P0tuvx_no(t,x,t,v))               &
              *bielec_pxxq_no(aa,x3,v3,bb))
        end do
      end do
    end if
  else
    ! ta/ub t/=u
    uu=list_act(u)
    bb=list_virt(b)
    term=0.D0
    do v=1,n_act_orb
      vv=list_act(v)
      v3=v+n_core_inact_orb
      do x=1,n_act_orb
        xx=list_act(x)
        x3=x+n_core_inact_orb
        term+=2.D0*(P0tuvx_no(t,u,v,x)*bielec_pqxx_no(aa,bb,v3,x3)      &
            +(P0tuvx_no(t,x,v,u)+P0tuvx_no(t,x,u,v))                 &
            *bielec_pxxq_no(aa,x3,v3,bb))
      end do
    end do
    if (a.eq.b) then
      term-=0.5D0*(occnum(tt)*Fipq(uu,tt)+occnum(uu)*Fipq(tt,uu))
      do v=1,n_act_orb
        do y=1,n_act_orb
          do x=1,n_act_orb
            term-=P0tuvx_no(t,v,x,y)*bielecCI_no(x,y,v,uu)
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


BEGIN_PROVIDER [double precision, hessmat, (nMonoEx,nMonoEx)]
 implicit none
 integer                        :: i,j,t,u,a,b
 integer                        :: indx,indx_tmp, jndx, jndx_tmp
 integer                        :: ustart,bstart
 real*8                         :: hessmat_itju
 real*8                         :: hessmat_itja
 real*8                         :: hessmat_itua
 real*8                         :: hessmat_iajb
 real*8                         :: hessmat_iatb
 real*8                         :: hessmat_taub
 !       c-a c-v a-v
 !  c-a | X   X  X
 !  c-v |     X  X 
 !  a-v |        X

  provide mo_two_e_integrals_in_map

 hessmat = 0.d0

 !$OMP PARALLEL DEFAULT(NONE) &
 !$OMP SHARED(hessmat,n_c_a_prov,list_idx_c_a,n_core_inact_orb,n_act_orb,mat_idx_c_a) &
 !$OMP PRIVATE(indx_tmp,indx,i,t,j,u,ustart,jndx)

 !$OMP DO
!!!! < Core-active| H |Core-active >
 ! Core-active excitations 
 do indx_tmp = 1, n_c_a_prov
  indx = list_idx_c_a(1,indx_tmp)
  i    = list_idx_c_a(2,indx_tmp)
  t    = list_idx_c_a(3,indx_tmp)
  ! Core-active excitations 
  do j = 1, n_core_inact_orb
   if (i.eq.j) then
     ustart=t
   else
     ustart=1
   end if
   do u=ustart,n_act_orb
    jndx = mat_idx_c_a(j,u)
    hessmat(jndx,indx) = hessmat_itju(i,t,j,u)
    hessmat(indx,jndx) = hessmat(jndx,indx)
   enddo
  enddo
 enddo
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL

 !$OMP PARALLEL DEFAULT(NONE) &
 !$OMP SHARED(hessmat,n_c_a_prov,n_c_v_prov,list_idx_c_a,list_idx_c_v) &
 !$OMP PRIVATE(indx_tmp,jndx_tmp,indx,i,t,j,a,jndx)

 !$OMP DO
!!!! < Core-active| H |Core-VIRTUAL >
 ! Core-active excitations 
 do indx_tmp = 1, n_c_a_prov
  indx = list_idx_c_a(1,indx_tmp)
  i    = list_idx_c_a(2,indx_tmp)
  t    = list_idx_c_a(3,indx_tmp)
  ! Core-VIRTUAL excitations 
  do jndx_tmp = 1, n_c_v_prov
   jndx = list_idx_c_v(1,jndx_tmp)
   j    = list_idx_c_v(2,jndx_tmp)
   a    = list_idx_c_v(3,jndx_tmp)
   hessmat(jndx,indx) = hessmat_itja(i,t,j,a)
   hessmat(indx,jndx) = hessmat(jndx,indx)
  enddo
 enddo
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL

 !$OMP PARALLEL DEFAULT(NONE) &
 !$OMP SHARED(hessmat,n_c_a_prov,n_a_v_prov,list_idx_c_a,list_idx_a_v) &
 !$OMP PRIVATE(indx_tmp,jndx_tmp,indx,i,t,u,a,jndx)

 !$OMP DO
!!!! < Core-active| H |ACTIVE-VIRTUAL >
 ! Core-active excitations 
 do indx_tmp = 1, n_c_a_prov
  indx = list_idx_c_a(1,indx_tmp)
  i    = list_idx_c_a(2,indx_tmp)
  t    = list_idx_c_a(3,indx_tmp)
  ! ACTIVE-VIRTUAL excitations 
  do jndx_tmp = 1, n_a_v_prov
   jndx = list_idx_a_v(1,jndx_tmp)
   u    = list_idx_a_v(2,jndx_tmp)
   a    = list_idx_a_v(3,jndx_tmp)
   hessmat(jndx,indx) = hessmat_itua(i,t,u,a)
   hessmat(indx,jndx) = hessmat(jndx,indx)
  enddo
 enddo

  !$OMP END DO NOWAIT
  !$OMP END PARALLEL


 if(hess_cv_cv)then
 !$OMP PARALLEL DEFAULT(NONE) &
 !$OMP SHARED(hessmat,n_c_v_prov,list_idx_c_v,n_core_inact_orb,n_virt_orb,mat_idx_c_v) &
 !$OMP PRIVATE(indx_tmp,indx,i,a,j,b,bstart,jndx)
  !$OMP DO
!!!!! < Core-VIRTUAL | H |Core-VIRTUAL >
  ! Core-VIRTUAL excitations 
  do indx_tmp = 1, n_c_v_prov
   indx = list_idx_c_v(1,indx_tmp)
   i    = list_idx_c_v(2,indx_tmp)
   a    = list_idx_c_v(3,indx_tmp)
   ! Core-VIRTUAL excitations 
   do j = 1, n_core_inact_orb
    if (i.eq.j) then
      bstart=a
    else
      bstart=1
    end if
    do b=bstart,n_virt_orb
     jndx = mat_idx_c_v(j,b)
     hessmat(jndx,indx) = hessmat_iajb(i,a,j,b)
     hessmat(indx,jndx) = hessmat(jndx,indx)
    enddo
   enddo
  enddo
 
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL
 endif

 !$OMP PARALLEL DEFAULT(NONE) &
 !$OMP SHARED(hessmat,n_c_v_prov,n_a_v_prov,list_idx_c_v,list_idx_a_v) &
 !$OMP PRIVATE(indx_tmp,jndx_tmp,indx,i,a,t,b,jndx)

 !$OMP DO
!!!! < Core-VIRTUAL | H |Active-VIRTUAL >
 ! Core-VIRTUAL excitations 
 do indx_tmp = 1, n_c_v_prov
  indx = list_idx_c_v(1,indx_tmp)
  i    = list_idx_c_v(2,indx_tmp)
  a    = list_idx_c_v(3,indx_tmp)
  ! Active-VIRTUAL excitations 
  do jndx_tmp = 1, n_a_v_prov
   jndx = list_idx_a_v(1,jndx_tmp)
   t    = list_idx_a_v(2,jndx_tmp)
   b    = list_idx_a_v(3,jndx_tmp)
   hessmat(jndx,indx) = hessmat_iatb(i,a,t,b)
   hessmat(indx,jndx) = hessmat(jndx,indx)
  enddo
 enddo
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL


 !$OMP PARALLEL DEFAULT(NONE) &
 !$OMP SHARED(hessmat,n_a_v_prov,list_idx_a_v,n_act_orb,n_virt_orb,mat_idx_a_v) &
 !$OMP PRIVATE(indx_tmp,indx,t,a,u,b,bstart,jndx)

 !$OMP DO
!!!! < Active-VIRTUAL | H |Active-VIRTUAL >
 ! Active-VIRTUAL excitations 
 do indx_tmp = 1, n_a_v_prov
  indx = list_idx_a_v(1,indx_tmp)
  t    = list_idx_a_v(2,indx_tmp)
  a    = list_idx_a_v(3,indx_tmp)
  ! Active-VIRTUAL excitations 
  do u=t,n_act_orb
   if (t.eq.u) then
     bstart=a
   else
     bstart=1
   end if
   do b=bstart,n_virt_orb
    jndx = mat_idx_a_v(u,b)
    hessmat(jndx,indx) = hessmat_taub(t,a,u,b)
    hessmat(indx,jndx) = hessmat(jndx,indx)
   enddo
  enddo
 enddo
  !$OMP END DO NOWAIT
  !$OMP END PARALLEL

END_PROVIDER 
