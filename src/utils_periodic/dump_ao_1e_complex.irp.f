program print_ao_1e_integrals
  call run
end

subroutine run
  use map_module
  implicit none
  
  integer                        :: i,j
  
  write(*,'(A)') 'ao_one_e_integrals_complex'
  write(*,'(A)') '---------------'
  do i=1,ao_num
    write(*,'(200(E24.15))') ao_one_e_integrals_complex(i,:)
  enddo
  write(*,'(A)') 'ao_overlap_complex'
  write(*,'(A)') '---------------'
  do i=1,ao_num
    write(*,'(200(E24.15))') ao_overlap_complex(i,:)
  enddo
  write(*,'(A)') 's_inv_complex'
  write(*,'(A)') '---------------'
  do i=1,ao_num
    write(*,'(200(E24.15))') s_inv_complex(i,:)
  enddo
  write(*,'(A)') 's_half_inv_complex'
  write(*,'(A)') '---------------'
  do i=1,ao_num
    write(*,'(200(E24.15))') s_half_inv_complex(i,:)
  enddo
  write(*,'(A)') 's_half_complex'
  write(*,'(A)') '---------------'
  do i=1,ao_num
    write(*,'(200(E24.15))') s_half_complex(i,:)
  enddo
  write(*,'(A)') 'ao_ortho_canonical_coef_complex'
  write(*,'(A)') '---------------'
  do i=1,ao_num
    write(*,'(200(E24.15))') ao_ortho_canonical_coef_complex(i,:)
  enddo
end
