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
  
!  if (bavard) then
    write(6,*) ' Filled the table of the Monoexcitations '
    do indx=1,nMonoEx
      write(6,*) ' ex ',indx,' : ',excit(1,indx),' -> '              &
          ,excit(2,indx),'  ',excit_class(indx)
    end do
!  end if
  
END_PROVIDER
