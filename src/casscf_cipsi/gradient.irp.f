use bitmasks

BEGIN_PROVIDER [ integer, nMonoEx ]
  BEGIN_DOC
  ! Number of single excitations
  END_DOC
  implicit none
  nMonoEx=n_core_inact_orb*n_act_orb+n_core_inact_orb*n_virt_orb+n_act_orb*n_virt_orb
END_PROVIDER

 BEGIN_PROVIDER [integer, n_c_a_prov]
&BEGIN_PROVIDER [integer, n_c_v_prov]
&BEGIN_PROVIDER [integer, n_a_v_prov]
  implicit none
  n_c_a_prov = n_core_inact_orb * n_act_orb
  n_c_v_prov = n_core_inact_orb * n_virt_orb
  n_a_v_prov = n_act_orb * n_virt_orb 
 END_PROVIDER 

 BEGIN_PROVIDER [integer, excit, (2,nMonoEx)]
&BEGIN_PROVIDER [character*3, excit_class, (nMonoEx)]
&BEGIN_PROVIDER [integer, list_idx_c_a, (3,n_c_a_prov) ]
&BEGIN_PROVIDER [integer, list_idx_c_v, (3,n_c_v_prov) ]
&BEGIN_PROVIDER [integer, list_idx_a_v, (3,n_a_v_prov) ]
&BEGIN_PROVIDER [integer, mat_idx_c_a, (n_core_inact_orb,n_act_orb)
&BEGIN_PROVIDER [integer, mat_idx_c_v, (n_core_inact_orb,n_virt_orb)
&BEGIN_PROVIDER [integer, mat_idx_a_v, (n_act_orb,n_virt_orb)
  BEGIN_DOC
  ! a list of the orbitals involved in the excitation
  END_DOC
  
  implicit none
  integer                        :: i,t,a,ii,tt,aa,indx,indx_tmp
  indx=0
  indx_tmp = 0
  do ii=1,n_core_inact_orb
    i=list_core_inact(ii)
    do tt=1,n_act_orb
      t=list_act(tt)
      indx+=1
      excit(1,indx)=i
      excit(2,indx)=t
      excit_class(indx)='c-a'
      indx_tmp += 1
      list_idx_c_a(1,indx_tmp) = indx
      list_idx_c_a(2,indx_tmp) = ii
      list_idx_c_a(3,indx_tmp) = tt
      mat_idx_c_a(ii,tt) = indx
    end do
  end do
  
  indx_tmp = 0
  do ii=1,n_core_inact_orb
    i=list_core_inact(ii)
    do aa=1,n_virt_orb
      a=list_virt(aa)
      indx+=1
      excit(1,indx)=i
      excit(2,indx)=a
      excit_class(indx)='c-v'
      indx_tmp += 1
      list_idx_c_v(1,indx_tmp) = indx
      list_idx_c_v(2,indx_tmp) = ii
      list_idx_c_v(3,indx_tmp) = aa 
      mat_idx_c_v(ii,aa) = indx
    end do
  end do
  
  indx_tmp = 0
  do tt=1,n_act_orb
    t=list_act(tt)
    do aa=1,n_virt_orb
      a=list_virt(aa)
      indx+=1
      excit(1,indx)=t
      excit(2,indx)=a
      excit_class(indx)='a-v'
      indx_tmp += 1
      list_idx_a_v(1,indx_tmp) = indx
      list_idx_a_v(2,indx_tmp) = tt
      list_idx_a_v(3,indx_tmp) = aa
      mat_idx_a_v(tt,aa) = indx
    end do
  end do
  
  if (bavard) then
    write(6,*) ' Filled the table of the Monoexcitations '
    do indx=1,nMonoEx
      write(6,*) ' ex ',indx,' : ',excit(1,indx),' -> '              &
          ,excit(2,indx),'  ',excit_class(indx)
    end do
  end if
  
END_PROVIDER

 BEGIN_PROVIDER [real*8, gradvec2, (nMonoEx)]
&BEGIN_PROVIDER [real*8, norm_grad_vec2]
&BEGIN_PROVIDER [real*8, norm_grad_vec2_tab, (3)]
  BEGIN_DOC
  ! calculate the orbital gradient <Psi| H E_pq |Psi> from density
  ! matrices and integrals; Siegbahn et al, Phys Scr 1980
  ! eqs 14 a,b,c
  END_DOC
  implicit none
  integer                        :: i,t,a,indx
  real*8                         :: gradvec_it,gradvec_ia,gradvec_ta
  
  indx=0
  norm_grad_vec2_tab = 0.d0
  do i=1,n_core_inact_orb
    do t=1,n_act_orb
      indx+=1
      gradvec2(indx)=gradvec_it(i,t)
      norm_grad_vec2_tab(1) += gradvec2(indx)*gradvec2(indx)
    end do
  end do
  
  do i=1,n_core_inact_orb
    do a=1,n_virt_orb
      indx+=1
      gradvec2(indx)=gradvec_ia(i,a)
      norm_grad_vec2_tab(2) += gradvec2(indx)*gradvec2(indx)
    end do
  end do
  
  do t=1,n_act_orb
    do a=1,n_virt_orb
      indx+=1
      gradvec2(indx)=gradvec_ta(t,a)
      norm_grad_vec2_tab(3) += gradvec2(indx)*gradvec2(indx)
    end do
  end do
  
  norm_grad_vec2=0.d0
  do indx=1,nMonoEx
    norm_grad_vec2+=gradvec2(indx)*gradvec2(indx)
  end do
  do i = 1, 3
   norm_grad_vec2_tab(i) = dsqrt(norm_grad_vec2_tab(i))
  enddo
  norm_grad_vec2=sqrt(norm_grad_vec2)
  if(bavard)then
   write(6,*)
   write(6,*) ' Norm of the orbital gradient (via D, P and integrals): ', norm_grad_vec2
   write(6,*)
  endif
  
END_PROVIDER

real*8 function gradvec_it(i,t)
  BEGIN_DOC
  ! the orbital gradient core/inactive -> active
  ! we assume natural orbitals
  END_DOC
  implicit none
  integer                        :: i,t
  
  integer                        :: ii,tt,v,vv,x,y
  integer                        :: x3,y3
  
  ii=list_core_inact(i)
  tt=list_act(t)
  gradvec_it=2.D0*(Fipq(tt,ii)+Fapq(tt,ii))
  gradvec_it-=occnum(tt)*Fipq(ii,tt)
  do v=1,n_act_orb ! active 
    vv=list_act(v)
    do x=1,n_act_orb ! active 
      x3=x+n_core_inact_orb ! list_act(x) 
      do y=1,n_act_orb ! active
        y3=y+n_core_inact_orb ! list_act(y)
        !                Gamma(2)  a a a a      1/r12      i  a  a  a
        gradvec_it-=2.D0*P0tuvx_no(t,v,x,y)*bielec_PQxx_no(ii,vv,x3,y3)
      end do
    end do
  end do
  gradvec_it*=2.D0
end function gradvec_it

real*8 function gradvec_ia(i,a)
  BEGIN_DOC
  ! the orbital gradient core/inactive -> virtual
  END_DOC
  implicit none
  integer                        :: i,a,ii,aa
  
  ii=list_core_inact(i)
  aa=list_virt(a)
  gradvec_ia=2.D0*(Fipq(aa,ii)+Fapq(aa,ii))
  gradvec_ia*=2.D0
  
end function gradvec_ia

real*8 function gradvec_ta(t,a)
  BEGIN_DOC
  ! the orbital gradient active -> virtual
  ! we assume natural orbitals
  END_DOC
  implicit none
  integer                        :: t,a,tt,aa,v,vv,x,y
  
  tt=list_act(t)
  aa=list_virt(a)
  gradvec_ta=0.D0
  gradvec_ta+=occnum(tt)*Fipq(aa,tt)
  do v=1,n_act_orb
    do x=1,n_act_orb
      do y=1,n_act_orb
        gradvec_ta+=2.D0*P0tuvx_no(t,v,x,y)*bielecCI_no(x,y,v,aa)
      end do
    end do
  end do
  gradvec_ta*=2.D0
  
end function gradvec_ta

