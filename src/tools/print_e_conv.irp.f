program print_e_conv
 implicit none
 BEGIN_DOC
! program that prints in a human readable format the convergence of the CIPSI algorithm.
!
! for all istate, this program produces
!
! * a file "EZFIO.istate.conv" containing the variational and var+PT2 energies as a function of N_det
!
! * for istate > 1, a file EZFIO.istate.delta_e.conv containing the energy difference (both var and var+PT2) with the ground state as a function of N_det
 END_DOC

  provide ezfio_filename
  call routine_e_conv
 end

subroutine routine_e_conv
 implicit none
 BEGIN_DOC
! routine called by :c:func:`print_e_conv`
 END_DOC
 integer :: N_iter_tmp
 integer :: i,istate
 character*(128) :: output
 integer :: i_unit_output,getUnitAndOpen
 character*(128) :: filename

 integer, allocatable :: n_det_tmp(:)
 call ezfio_get_iterations_N_iter(N_iter_tmp)
 print*,'N_iter_tmp = ',N_iter_tmp
 double precision, allocatable :: e(:,:),pt2(:,:)
 allocate(e(N_states, 100),pt2(N_states, 100),n_det_tmp(100))
 call ezfio_get_iterations_energy_iterations(e)
 call ezfio_get_iterations_pt2_iterations(pt2)
 call ezfio_get_iterations_n_det_iterations(n_det_tmp)


  do istate = 1, N_states
   if (istate.lt.10)then
    write (filename, "(I1)")istate
   else
    write (filename, "(I2)")istate
   endif
   print*,filename
   output=trim(ezfio_filename)//'.'//trim(filename)//'.conv'
   output=trim(output)
   print*,'output = ',trim(output)
   i_unit_output = getUnitAndOpen(output,'w')
   write(i_unit_output,*)'# N_det    E_var        E_var + PT2'
   do i = 1, N_iter_tmp
    write(i_unit_output,'(I9,X,3(F16.10,X))')n_det_tmp(i),e(istate,i),e(istate,i) + pt2(istate,i)
   enddo
  enddo

  if(N_states.gt.1)then
   double precision, allocatable :: deltae(:,:),deltae_pt2(:,:)
   allocate(deltae(N_states,100),deltae_pt2(N_states,100))
   do i = 1, N_iter_tmp
    do istate = 1, N_states
     deltae(istate,i) = e(istate,i) - e(1,i)
     deltae_pt2(istate,i) = e(istate,i) + pt2(istate,i) - (e(1,i) + pt2(1,i))
    enddo
   enddo
   do istate = 2, N_states
     if (istate.lt.10)then
      write (filename, "(I1)")istate
     else
      write (filename, "(I2)")istate
     endif
     output=trim(ezfio_filename)//'.'//trim(filename)//'.delta_e.conv'
     print*,'output = ',trim(output)
     i_unit_output = getUnitAndOpen(output,'w')
     write(i_unit_output,*)'# N_det    Delta E_var   Delta (E_var + PT2)'
     do i = 1, N_iter_tmp
      write(i_unit_output,'(I9,X,100(F16.10,X))')n_det_tmp(i),deltae(istate,i),deltae_pt2(istate,i)
     enddo
   enddo
  endif

end
