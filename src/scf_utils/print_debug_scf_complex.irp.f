subroutine print_debug_scf_complex
  implicit none
  BEGIN_DOC
! Build the MOs using the extended Huckel model
  END_DOC
  integer                        :: i,j
  
  write(*,'(A)') 'mo_coef_complex'
  write(*,'(A)') '---------------'
  do i=1,ao_num
    write(*,'(200(E24.15))') mo_coef_complex(i,:)
  enddo
  write(*,'(A)') 'scf_density_matrix_ao_alpha_complex'
  write(*,'(A)') '---------------'
  do i=1,ao_num
    write(*,'(200(E24.15))') scf_density_matrix_ao_alpha_complex(i,:)
  enddo
  write(*,'(A)') 'scf_density_matrix_ao_beta_complex'
  write(*,'(A)') '---------------'
  do i=1,ao_num
    write(*,'(200(E24.15))') scf_density_matrix_ao_beta_complex(i,:)
  enddo
  write(*,'(A)') 'ao_two_e_integral_alpha_complex'
  write(*,'(A)') '---------------'
  do i=1,ao_num
    write(*,'(200(E24.15))') ao_two_e_integral_alpha_complex(i,:)
  enddo
  write(*,'(A)') 'ao_two_e_integral_beta_complex'
  write(*,'(A)') '---------------'
  do i=1,ao_num
    write(*,'(200(E24.15))') ao_two_e_integral_beta_complex(i,:)
  enddo
  write(*,'(A)') 'fock_matrix_ao_alpha_complex'
  write(*,'(A)') '---------------'
  do i=1,ao_num
    write(*,'(200(E24.15))') fock_matrix_ao_alpha_complex(i,:)
  enddo
  write(*,'(A)') 'fock_matrix_ao_beta_complex'
  write(*,'(A)') '---------------'
  do i=1,ao_num
    write(*,'(200(E24.15))') fock_matrix_ao_beta_complex(i,:)
  enddo

end
