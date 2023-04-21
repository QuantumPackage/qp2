BEGIN_PROVIDER [ character*(128), json_filename ]
 implicit none
 BEGIN_DOC
 ! Fortran unit of the JSON file
 END_DOC
 integer, external            :: getUnitAndOpen
 integer                      :: counter
 character*(128)              :: prefix
 logical                      :: exists

 prefix = trim(ezfio_filename)//'/json/'

 exists = .True.
 counter = 0
 do while (exists)
   counter += 1
   write(json_filename, '(A,I5.5,A)') trim(prefix), counter, '.json'
   INQUIRE(FILE=trim(json_filename), EXIST=exists)
 enddo

END_PROVIDER

BEGIN_PROVIDER [ integer, json_unit]
 implicit none
 BEGIN_DOC
 ! Unit file for JSON output
 END_DOC
 integer, external :: getUnitAndOpen
 call ezfio_set_json_empty(.False.)
 json_unit = getUnitAndOpen(json_filename, 'w')
 write(json_unit, '(A)') '{'
END_PROVIDER

subroutine json_close
 write(json_unit, '(A)') '}'
 close(json_unit)
 FREE json_unit
end

