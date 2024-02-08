
subroutine write_array_6_index_tensor(n_orb,array_tmp,name_file)
 implicit none
 integer, intent(in) :: n_orb
 character*(128),  intent(in) :: name_file 
 double precision, intent(in) :: array_tmp(n_orb,n_orb,n_orb,n_orb,n_orb,n_orb)

 character*(128)                :: output
 integer                        :: i_unit_output,getUnitAndOpen
 PROVIDE ezfio_filename                                                                                                  
 output=trim(ezfio_filename)//'/work/'//trim(name_file)
 i_unit_output = getUnitAndOpen(output,'W')
 write(i_unit_output)array_tmp
 close(unit=i_unit_output)
end

subroutine read_array_6_index_tensor(n_orb,array_tmp,name_file)
 implicit none
 character*(128)                :: output
 integer                        :: i_unit_output,getUnitAndOpen
 integer, intent(in) :: n_orb
 character*(128),  intent(in)  :: name_file 
 double precision, intent(out) :: array_tmp(n_orb,n_orb,n_orb,n_orb,n_orb,n_orb)
 PROVIDE ezfio_filename                                                                                                  
 output=trim(ezfio_filename)//'/work/'//trim(name_file)
 i_unit_output = getUnitAndOpen(output,'R')
 read(i_unit_output)array_tmp
 close(unit=i_unit_output)
end

subroutine read_fcidump_3_tc(array)
 implicit none
 double precision, intent(out) :: array(mo_num, mo_num, mo_num, mo_num, mo_num, mo_num)
 integer :: i,j,k,l,m,n,i_mo, Reason
 double precision :: integral 
 print*,'Reading the THREE-body integrals from a TC FCIDUMP'
 open (unit=15, file="TCDUMP-nosym", status='old',    &
              access='sequential', action='read' )
 read(15,*)i_mo
 if(i_mo.ne.mo_num)then
  print*,'Something went wrong in the read_fcidump_3_tc !'
  print*,'i_mo.ne.mo_num !'
  print*,i_mo,mo_num
  stop
 endif
 do 
  read(15,*,IOSTAT=Reason)integral,i, j, m, k, l, n
  if(Reason > 0)then
   print*,'Something went wrong in the I/O of read_fcidump_3_tc'
   stop
  else if(Reason < 0)then
   exit   
  else
        !                          1 2 3 1 2 3
        !                         <i j m|k l n>
        !                         (ik|jl|mn)
!        integral = integral * 1.d0/3.d0 !!!! For NECI convention 
          array(i,j,m,k,l,n) =  integral * 3.d0
  
   endif
  enddo

end
