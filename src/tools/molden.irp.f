program molden
  implicit none
  BEGIN_DOC
  ! Produce a Molden file
  END_DOC
  character*(128)                :: output
  integer                        :: i_unit_output,getUnitAndOpen
  provide ezfio_filename

  integer                        :: i
  print*,trim(ezfio_filename)
  output=trim(ezfio_filename)//'.mol'
  print*,'output = ',trim(output)
  i_unit_output = getUnitAndOpen(output,'w')
  print*,'i_unit_output = ',i_unit_output
  call write_intro_gamess(i_unit_output)
  call write_geometry(i_unit_output)
  call write_Ao_basis(i_unit_output)
  call write_Mo_basis(i_unit_output)


  write(i_unit_output,*)''
  write(i_unit_output,*)''
  write(i_unit_output,*)'          ------------------------'

  close(i_unit_output)
end

subroutine write_intro_gamess(i_unit_output)
  implicit none
  integer, intent(in)            :: i_unit_output
  integer                        :: i,j,k,l

  write(i_unit_output,*)'         *         GAMESS VERSION = 22 FEB 2006 (R5)          *'
  write(i_unit_output,*)'         *             FROM IOWA STATE UNIVERSITY             *'
  write(i_unit_output,*)'         * M.W.SCHMIDT, K.K.BALDRIDGE, J.A.BOATZ, S.T.ELBERT, *'
  write(i_unit_output,*)'         *   M.S.GORDON, J.H.JENSEN, S.KOSEKI, N.MATSUNAGA,   *'
  write(i_unit_output,*)'         *          K.A.NGUYEN, S.J.SU, T.L.WINDUS,           *'
  write(i_unit_output,*)'         *       TOGETHER WITH M.DUPUIS, J.A.MONTGOMERY       *'
  write(i_unit_output,*)'         *         J.COMPUT.CHEM.  14, 1347-1363(1993)        *'
  write(i_unit_output,*)''

end




subroutine write_geometry(i_unit_output)
  implicit none
  integer, intent(in)            :: i_unit_output
  integer                        :: i,j,k,l, getUnitAndOpen



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(i_unit_output,*)'ATOM      ATOMIC                      COORDINATES (BOHR)          '
  write(i_unit_output,*)'          CHARGE         X                   Y                   Z'
  do i = 1, nucl_num
    ! write(i_unit_output,'(A2 I3 X F3.1 X 3(F16.10))') trim(element_name(int(nucl_charge(i)))),i,(nucl_charge(i)), nucl_coord(i,1), nucl_coord(i,2), nucl_coord(i,3)
    write(i_unit_output,'(A2,I1, 9X F5.1 X 3(F16.10 ,4X))') trim(element_name(int(nucl_charge(i)))),i,(nucl_charge(i)), nucl_coord(i,1), nucl_coord(i,2), nucl_coord(i,3)
  enddo
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end


subroutine write_Ao_basis(i_unit_output)
  implicit none
  integer, intent(in)            :: i_unit_output
  integer                        :: i,j,k,l, getUnitAndOpen
  character*(128)                :: character_shell
  integer                        :: i_shell,i_prim,i_ao
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(i_unit_output,*)''
  write(i_unit_output,*)''
  write(i_unit_output,*)'    ATOMIC BASIS SET'
  write(i_unit_output,*)'    ----------------'
  write(i_unit_output,*)'THE CONTRACTED PRIMITIVE FUNCTIONS HAVE BEEN UNNORMALIZED'
  write(i_unit_output,*)'THE CONTRACTED BASIS FUNCTIONS ARE NOW NORMALIZED TO UNITY'
  write(i_unit_output,*)''
  write(i_unit_output,*)'SHELL TYPE  PRIMITIVE        EXPONENT          CONTRACTION COEFFICIENT(S)'
  write(i_unit_output,*)''
  write(i_unit_output,*)''

  i_shell = 0
  i_prim = 0
  do i = 1, Nucl_num
    write(i_unit_output,'(A2,I1)') trim(element_name(int(nucl_charge(i)))),i
    write(i_unit_output,*)' '
    !  write(i_unit_output,*)'Nucl_num_shell_Aos(i) = ',Nucl_num_shell_Aos(i)
    do j = 1, Nucl_num_shell_Aos(i)
      i_shell +=1
      i_ao = Nucl_list_shell_Aos(i,j)
      character_shell = trim(ao_l_char(i_ao))
      !   write(i_unit_output,*),j,i_shell,i_ao!trim(character_shell)
      do k = 1, ao_prim_num(i_ao)
        i_prim +=1
        if(i_prim.lt.100)then
          write(i_unit_output,'(4X,I3,3X,A1,6X,I2,6X,F16.7,2X,F16.12)')i_shell,character_shell,i_prim,ao_expo(i_ao,k),ao_coef(i_ao,k)
        else
          write(i_unit_output,'(4X,I3,3X,A1,5X,I3,6X,F16.7,2X,F16.12)')i_shell,character_shell,i_prim,ao_expo(i_ao,k),ao_coef(i_ao,k)
        endif
      enddo
      write(i_unit_output,*)''
    enddo
  enddo

  write(i_unit_output,*)''
  write(i_unit_output,'(A47,2X,I3)')'TOTAL NUMBER OF BASIS SET SHELLS             =', i_shell
  write(i_unit_output,'(A47,2X,I3)')'NUMBER OF CARTESIAN GAUSSIAN BASIS FUNCTIONS =', ao_num
  ! this is for the new version of molden
  write(i_unit_output,'(A12)')'PP    =NONE'
  write(i_unit_output,*)''


end

subroutine write_Mo_basis(i_unit_output)
  implicit none
  integer, intent(in)            :: i_unit_output
  integer                        :: i,j,k,l, getUnitAndOpen
  integer                        :: i_5,i_mod

  write(i_unit_output,*) '          ----------------------'
  write(i_unit_output,*) '          MCSCF NATURAL ORBITALS'
  write(i_unit_output,*) '          ----------------------'
  write(i_unit_output,*) '                                '

  do j = 1, mo_num
    write(i_unit_output,'(18X,I3)')j
    write(i_unit_output,*)''
    write(i_unit_output,'(18X,F8.5)')-1.d0
    write(i_unit_output,*)''
    do i = 1, ao_num
      !   write(i_unit_output,'(2X,I3, 2X A1,  I3, 2X A4 , F9.6)')i,trim(element_name(int(nucl_charge(ao_nucl(i))))),ao_nucl(i),(ao_l_char_space(i)),mo_coef(i,j)
      ! F12.6 for larger coefficients...
      write(i_unit_output,'(2X,I3, 2X A1,  I3, 2X A4 , F12.6)')i,trim(element_name(int(nucl_charge(ao_nucl(i))))),ao_nucl(i),(ao_l_char_space(i)),mo_coef(i,j)
      !  write(i_unit_output,'(I3, X A1, X I3, X A4 X F16.8)')i,trim(element_name(int(nucl_charge(ao_nucl(i))))),ao_nucl(i),(ao_l_char_space(i))
    enddo
    write(i_unit_output,*)''
  enddo

end
