BEGIN_PROVIDER [ character*(1024), trexio_filename ]
 implicit none
 BEGIN_DOC
 ! Name of the TREXIO file
 END_DOC
 character*(1024) :: prefix

 trexio_filename = trexio_file

 if (trexio_file == 'None') then
    prefix = trim(ezfio_work_dir)//trim(ezfio_filename)
    if (backend == 0) then
      trexio_filename = trim(prefix)//'.h5'
    else if (backend == 1) then
      trexio_filename = trim(prefix)
    endif
 endif
END_PROVIDER


