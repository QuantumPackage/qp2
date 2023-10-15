subroutine reorder_orbitals_for_casscf
 implicit none
 BEGIN_DOC
! routine that reorders the orbitals of the CASSCF in terms block of core, active and virtual
 END_DOC
 integer :: i,j,iorb
 integer, allocatable :: iorder(:),array(:)
 allocate(iorder(mo_num),array(mo_num))
 do i = 1, n_core_orb
  iorb = list_core(i)
  array(iorb) = i
 enddo

 do i = 1, n_inact_orb
  iorb = list_inact(i)
  array(iorb) = mo_num + i
 enddo

 do i = 1, n_act_orb
  iorb = list_act(i)
  array(iorb) = 2 * mo_num + i
 enddo

 do i = 1, n_virt_orb
  iorb = list_virt(i)
  array(iorb) = 3 * mo_num + i
 enddo

 do i = 1, mo_num
  iorder(i) = i
 enddo
 call isort(array,iorder,mo_num)
 double precision, allocatable :: mo_coef_new(:,:)
 allocate(mo_coef_new(ao_num,mo_num))
 do i = 1, mo_num
  mo_coef_new(:,i) = mo_coef(:,iorder(i))
 enddo
 mo_coef = mo_coef_new 
 touch mo_coef

 list_core_reverse = 0
 do i = 1, n_core_orb
  list_core(i) = i
  list_core_reverse(i) = i
  mo_class(i) = "Core"
 enddo

 list_inact_reverse = 0
 do i = 1, n_inact_orb
  list_inact(i) = i + n_core_orb
  list_inact_reverse(i+n_core_orb) = i
  mo_class(i+n_core_orb) = "Inactive"
 enddo

 list_act_reverse = 0
 do i = 1, n_act_orb
  list_act(i) = n_core_inact_orb + i
  list_act_reverse(n_core_inact_orb + i) = i
  mo_class(n_core_inact_orb + i) = "Active"
 enddo

 list_virt_reverse = 0
 do i = 1, n_virt_orb
  list_virt(i) = n_core_inact_orb + n_act_orb + i
  list_virt_reverse(n_core_inact_orb + n_act_orb + i) = i
  mo_class(n_core_inact_orb + n_act_orb + i) = "Virtual"
 enddo
 touch list_core_reverse list_core list_inact list_inact_reverse list_act list_act_reverse list_virt list_virt_reverse 

end
