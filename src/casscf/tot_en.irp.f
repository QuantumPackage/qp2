 BEGIN_PROVIDER [real*8, etwo]
&BEGIN_PROVIDER [real*8, eone]
&BEGIN_PROVIDER [real*8, eone_bis]
&BEGIN_PROVIDER [real*8, etwo_bis]
&BEGIN_PROVIDER [real*8, etwo_ter]
&BEGIN_PROVIDER [real*8, ecore]
&BEGIN_PROVIDER [real*8, ecore_bis]
   implicit none
   integer                        :: t,u,v,x,i,ii,tt,uu,vv,xx,j,jj,t3,u3,v3,x3
   real*8                         :: e_one_all,e_two_all
   e_one_all=0.D0
   e_two_all=0.D0
   do i=1,n_core_inact_orb
     ii=list_core_inact(i)
     e_one_all+=2.D0*mo_one_e_integrals(ii,ii)
     do j=1,n_core_inact_orb
       jj=list_core_inact(j)
       e_two_all+=2.D0*bielec_PQxx(ii,ii,j,j)-bielec_PQxx(ii,jj,j,i)
     end do
     do t=1,n_act_orb
       tt=list_act(t)
       t3=t+n_core_inact_orb
       do u=1,n_act_orb
         uu=list_act(u)
         u3=u+n_core_inact_orb
         e_two_all+=D0tu(t,u)*(2.D0*bielec_PQxx(tt,uu,i,i)           &
             -bielec_PQxx(tt,ii,i,u3))
       end do
     end do
   end do
   do t=1,n_act_orb
     tt=list_act(t)
     do u=1,n_act_orb
       uu=list_act(u)
       e_one_all+=D0tu(t,u)*mo_one_e_integrals(tt,uu)
       do v=1,n_act_orb
         v3=v+n_core_inact_orb
         do x=1,n_act_orb
           x3=x+n_core_inact_orb
           e_two_all  +=P0tuvx(t,u,v,x)*bielec_PQxx(tt,uu,v3,x3)
         end do
       end do
     end do
   end do
   ecore    =nuclear_repulsion
   ecore_bis=nuclear_repulsion
   do i=1,n_core_inact_orb
     ii=list_core_inact(i)
     ecore    +=2.D0*mo_one_e_integrals(ii,ii)
     ecore_bis+=2.D0*mo_one_e_integrals(ii,ii)
     do j=1,n_core_inact_orb
       jj=list_core_inact(j)
       ecore    +=2.D0*bielec_PQxx(ii,ii,j,j)-bielec_PQxx(ii,jj,j,i)
       ecore_bis+=2.D0*bielec_PxxQ(ii,i,j,jj)-bielec_PxxQ(ii,j,j,ii)
     end do
   end do
   eone    =0.D0
   eone_bis=0.D0
   etwo    =0.D0
   etwo_bis=0.D0
   etwo_ter=0.D0
   do t=1,n_act_orb
     tt=list_act(t)
     t3=t+n_core_inact_orb
     do u=1,n_act_orb
       uu=list_act(u)
       u3=u+n_core_inact_orb
       eone    +=D0tu(t,u)*mo_one_e_integrals(tt,uu)
       eone_bis+=D0tu(t,u)*mo_one_e_integrals(tt,uu)
       do i=1,n_core_inact_orb
         ii=list_core_inact(i)
         eone    +=D0tu(t,u)*(2.D0*bielec_PQxx(tt,uu,i,i)            &
             -bielec_PQxx(tt,ii,i,u3))
         eone_bis+=D0tu(t,u)*(2.D0*bielec_PxxQ(tt,u3,i,ii)           &
             -bielec_PxxQ(tt,i,i,uu))
       end do
       do v=1,n_act_orb
         vv=list_act(v)
         v3=v+n_core_inact_orb
         do x=1,n_act_orb
           xx=list_act(x)
           x3=x+n_core_inact_orb
           real*8                         :: h1,h2,h3
           h1=bielec_PQxx(tt,uu,v3,x3)
           h2=bielec_PxxQ(tt,u3,v3,xx)
           h3=bielecCI(t,u,v,xx)
           etwo    +=P0tuvx(t,u,v,x)*h1
           etwo_bis+=P0tuvx(t,u,v,x)*h2
           etwo_ter+=P0tuvx(t,u,v,x)*h3
           if ((h1.ne.h2).or.(h1.ne.h3)) then
             write(6,9901) t,u,v,x,h1,h2,h3
             9901 format('aie: ',4I4,3E20.12)
           end if
         end do
       end do
     end do
   end do
   
END_PROVIDER
 
 
