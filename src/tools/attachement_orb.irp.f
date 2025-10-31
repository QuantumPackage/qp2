program attachement_orb
 implicit none
 read_wf=.True.
 touch read_wf
 call molden_attachment
end

subroutine molden_attachment
  implicit none
  BEGIN_DOC
  ! Produces a Molden file
  END_DOC
  character*(128)                :: output
  integer                        :: i_unit_output,getUnitAndOpen
  integer                        :: i,j,k,l
  double precision, parameter :: a0 = 0.529177249d0

  PROVIDE ezfio_filename

  output=trim(ezfio_filename)//'.attachement.mol'
  print*,'output = ',trim(output)

  i_unit_output = getUnitAndOpen(output,'w')

  write(i_unit_output,'(A)') '[Molden Format]'

  write(i_unit_output,'(A)') '[Atoms] Angs'
  do i = 1, nucl_num
    write(i_unit_output,'(A2,2X,I4,2X,I4,3(2X,F15.10))')                   &
        trim(element_name(int(nucl_charge(i)))),                     &
        i,                                                           &
        int(nucl_charge(i)),                                         &
        nucl_coord(i,1)*a0, nucl_coord(i,2)*a0, nucl_coord(i,3)*a0
  enddo

  write(i_unit_output,'(A)') '[GTO]'

  character*(1)                  :: character_shell
  integer                        :: i_shell,i_prim,i_ao
  integer                        :: iorder(ao_num)
  integer                        :: nsort(ao_num)

  i_shell = 0
  i_prim = 0
  do i=1,nucl_num
    write(i_unit_output,*) i, 0
    do j=1,nucl_num_shell_aos(i)
      i_shell +=1
      i_ao = nucl_list_shell_aos(i,j)
      character_shell = trim(ao_l_char(i_ao))
      write(i_unit_output,*) character_shell, ao_prim_num(i_ao), '1.00'
      do k = 1, ao_prim_num(i_ao)
        i_prim +=1
        write(i_unit_output,'(ES20.10,2X,ES20.10)') ao_expo(i_ao,k), ao_coef(i_ao,k)
      enddo
      l = i_ao
      do while ( ao_l(l) == ao_l(i_ao) )
        nsort(l) = i*10000 + j*100
        l += 1
        if (l > ao_num) exit
      enddo
    enddo
    write(i_unit_output,*)''
  enddo


  do i=1,ao_num
    iorder(i) = i
    ! p
    if      ((ao_power(i,1) == 1 ).and.(ao_power(i,2) == 0 ).and.(ao_power(i,3) == 0 )) then
      nsort(i) += 1
    else if ((ao_power(i,1) == 0 ).and.(ao_power(i,2) == 1 ).and.(ao_power(i,3) == 0 )) then
      nsort(i) += 2
    else if ((ao_power(i,1) == 0 ).and.(ao_power(i,2) == 0 ).and.(ao_power(i,3) == 1 )) then
      nsort(i) += 3
    ! d
    else if ((ao_power(i,1) == 2 ).and.(ao_power(i,2) == 0 ).and.(ao_power(i,3) == 0 )) then
      nsort(i) += 1
    else if ((ao_power(i,1) == 0 ).and.(ao_power(i,2) == 2 ).and.(ao_power(i,3) == 0 )) then
      nsort(i) += 2
    else if ((ao_power(i,1) == 0 ).and.(ao_power(i,2) == 0 ).and.(ao_power(i,3) == 2 )) then
      nsort(i) += 3
    else if ((ao_power(i,1) == 1 ).and.(ao_power(i,2) == 1 ).and.(ao_power(i,3) == 0 )) then
      nsort(i) += 4
    else if ((ao_power(i,1) == 1 ).and.(ao_power(i,2) == 0 ).and.(ao_power(i,3) == 1 )) then
      nsort(i) += 5
    else if ((ao_power(i,1) == 0 ).and.(ao_power(i,2) == 1 ).and.(ao_power(i,3) == 1 )) then
      nsort(i) += 6
    ! f
    else if ((ao_power(i,1) == 3 ).and.(ao_power(i,2) == 0 ).and.(ao_power(i,3) == 0 )) then
      nsort(i) += 1
    else if ((ao_power(i,1) == 0 ).and.(ao_power(i,2) == 3 ).and.(ao_power(i,3) == 0 )) then
      nsort(i) += 2
    else if ((ao_power(i,1) == 0 ).and.(ao_power(i,2) == 0 ).and.(ao_power(i,3) == 3 )) then
      nsort(i) += 3
    else if ((ao_power(i,1) == 1 ).and.(ao_power(i,2) == 2 ).and.(ao_power(i,3) == 0 )) then
      nsort(i) += 4
    else if ((ao_power(i,1) == 2 ).and.(ao_power(i,2) == 1 ).and.(ao_power(i,3) == 0 )) then
      nsort(i) += 5
    else if ((ao_power(i,1) == 2 ).and.(ao_power(i,2) == 0 ).and.(ao_power(i,3) == 1 )) then
      nsort(i) += 6
    else if ((ao_power(i,1) == 1 ).and.(ao_power(i,2) == 0 ).and.(ao_power(i,3) == 2 )) then
      nsort(i) += 7
    else if ((ao_power(i,1) == 0 ).and.(ao_power(i,2) == 1 ).and.(ao_power(i,3) == 2 )) then
      nsort(i) += 8
    else if ((ao_power(i,1) == 0 ).and.(ao_power(i,2) == 2 ).and.(ao_power(i,3) == 1 )) then
      nsort(i) += 9
    else if ((ao_power(i,1) == 1 ).and.(ao_power(i,2) == 1 ).and.(ao_power(i,3) == 1 )) then
      nsort(i) += 10
    ! g
    else if ((ao_power(i,1) == 4 ).and.(ao_power(i,2) == 0 ).and.(ao_power(i,3) == 0 )) then
      nsort(i) += 1
    else if ((ao_power(i,1) == 0 ).and.(ao_power(i,2) == 4 ).and.(ao_power(i,3) == 0 )) then
      nsort(i) += 2
    else if ((ao_power(i,1) == 0 ).and.(ao_power(i,2) == 0 ).and.(ao_power(i,3) == 4 )) then
      nsort(i) += 3
    else if ((ao_power(i,1) == 3 ).and.(ao_power(i,2) == 1 ).and.(ao_power(i,3) == 0 )) then
      nsort(i) += 4
    else if ((ao_power(i,1) == 3 ).and.(ao_power(i,2) == 0 ).and.(ao_power(i,3) == 1 )) then
      nsort(i) += 5
    else if ((ao_power(i,1) == 1 ).and.(ao_power(i,2) == 3 ).and.(ao_power(i,3) == 0 )) then
      nsort(i) += 6
    else if ((ao_power(i,1) == 0 ).and.(ao_power(i,2) == 3 ).and.(ao_power(i,3) == 1 )) then
      nsort(i) += 7
    else if ((ao_power(i,1) == 1 ).and.(ao_power(i,2) == 0 ).and.(ao_power(i,3) == 3 )) then
      nsort(i) += 8
    else if ((ao_power(i,1) == 0 ).and.(ao_power(i,2) == 1 ).and.(ao_power(i,3) == 3 )) then
      nsort(i) += 9
    else if ((ao_power(i,1) == 2 ).and.(ao_power(i,2) == 2 ).and.(ao_power(i,3) == 0 )) then
      nsort(i) += 10
    else if ((ao_power(i,1) == 2 ).and.(ao_power(i,2) == 0 ).and.(ao_power(i,3) == 2 )) then
      nsort(i) += 11
    else if ((ao_power(i,1) == 0 ).and.(ao_power(i,2) == 2 ).and.(ao_power(i,3) == 2 )) then
      nsort(i) += 12
    else if ((ao_power(i,1) == 2 ).and.(ao_power(i,2) == 1 ).and.(ao_power(i,3) == 1 )) then
      nsort(i) += 13
    else if ((ao_power(i,1) == 1 ).and.(ao_power(i,2) == 2 ).and.(ao_power(i,3) == 1 )) then
      nsort(i) += 14
    else if ((ao_power(i,1) == 1 ).and.(ao_power(i,2) == 1 ).and.(ao_power(i,3) == 2 )) then
      nsort(i) += 15
    endif
  enddo

  call isort(nsort,iorder,ao_num)
  write(i_unit_output,'(A)') '[MO]'
  integer :: istate
  istate = 2
   do i=1,n_dettachment(istate)
     write (i_unit_output,*) 'Sym= 1'
     write (i_unit_output,*) 'Ene=', dettachment_numbers_sorted(i,istate)
     write (i_unit_output,*) 'Spin= Alpha'
     write (i_unit_output,*) 'Occup=', dettachment_numbers_sorted(i,istate)
     do j=1,ao_num
       write(i_unit_output, '(I6,2X,ES20.10)') j, dettachment_orbitals(iorder(j),i,istate)
     enddo
   enddo
   do i=1,n_attachment(istate)
     write (i_unit_output,*) 'Sym= 1'
     write (i_unit_output,*) 'Ene=', attachment_numbers_sorted(i,istate)
     write (i_unit_output,*) 'Spin= Alpha'
     write (i_unit_output,*) 'Occup=', attachment_numbers_sorted(i,istate)
     do j=1,ao_num
       write(i_unit_output, '(I6,2X,ES20.10)') j, attachment_orbitals(iorder(j),i,istate)
     enddo
   enddo
  close(i_unit_output)
end

