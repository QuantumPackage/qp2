program fcidump
  implicit none
  BEGIN_DOC
! Produce a regular `FCIDUMP` file from the |MOs| stored in the |EZFIO|
! directory.
!
! To specify an active space, the class of the |MOs| have to set in the
! |EZFIO| directory (see :ref:`qp_set_mo_class`).
!
! The :ref:`fcidump` program supports 3 types of |MO| classes :
!
! * the *core* orbitals which are always doubly occupied in the
!   calculation
!
! * the *deleted* orbitals that are never occupied in the calculation
!
! * the *active* orbitals that are occupied with a varying number of
!   electrons
!
  END_DOC
  if (is_complex) then
    call fcidump_complex
  else
    call fcidump_real
  endif
end

subroutine fcidump_complex
  implicit none
  character*(128) :: output
  integer :: i_unit_output,getUnitAndOpen
  output=trim(ezfio_filename)//'.FCIDUMP'
  i_unit_output = getUnitAndOpen(output,'w')

  integer :: i,j,k,l
  integer :: i1,j1,k1,l1
  integer :: i2,j2,k2,l2,ik2,jl2
  integer :: ki,kj,kk,kl
  integer :: ii,ij,ik,il
  integer*8 :: m
  character*(2), allocatable :: A(:)

  write(i_unit_output,*) '&FCI NORB=', n_act_orb, ', NELEC=', elec_num-n_core_orb*2, &
   ', MS2=', (elec_alpha_num-elec_beta_num), ','
  allocate (A(n_act_orb))
  A = '1,'
  write(i_unit_output,*) 'ORBSYM=', (A(i), i=1,n_act_orb)
  write(i_unit_output,*) 'ISYM=0,'
  write(i_unit_output,*) '/'
  deallocate(A)

  integer(key_kind), allocatable :: keys(:)
  double precision, allocatable  :: values(:)
  integer(cache_map_size_kind)   :: n_elements, n_elements_max
  PROVIDE mo_two_e_integrals_in_map

  complex*16 :: get_two_e_integral_complex, integral

  do kl=1,kpt_num
    do kj=1,kl
      do kk=1,kl
        ki=kconserv(kl,kk,kj)
        if (ki>kl) cycle
        do l1=1,n_act_orb_kpts(kl)
          il=list_act_kpts(l1,kl)
          l = (kl-1)*mo_num_per_kpt + il
          do j1=1,n_act_orb_kpts(kj)
            ij=list_act_kpts(j1,kj)
            j = (kj-1)*mo_num_per_kpt + ij
            if (j>l) exit
            call idx2_tri_int(j,l,jl2)
            do k1=1,n_act_orb_kpts(kk)
              ik=list_act_kpts(k1,kk)
              k = (kk-1)*mo_num_per_kpt + ik
              if (k>l) exit
              do i1=1,n_act_orb_kpts(ki)
                ii=list_act_kpts(i1,ki)
                i = (ki-1)*mo_num_per_kpt + ii
                if ((j==l) .and. (i>k)) exit
                call idx2_tri_int(i,k,ik2)
                if (ik2 > jl2) exit
                integral = get_two_e_integral_complex(i,j,k,l,mo_integrals_map,mo_integrals_map_2)
                if (cdabs(integral) > mo_integrals_threshold) then
                  write(i_unit_output,'(2(E25.15,X),4(I6,X))') dble(integral), dimag(integral),i,k,j,l
                endif
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  do kj=1,kpt_num
    do j1=1,n_act_orb_kpts(kj)
      ij = list_act_kpts(j1,kj)
      j = (kj-1)*mo_num_per_kpt + ij
      do i1=j1,n_act_orb_kpts(kj)
        ii = list_act_kpts(i1,kj)
        i = (kj-1)*mo_num_per_kpt + ii
        integral = mo_one_e_integrals_kpts(ii,ij,kj) + core_fock_operator_complex(i,j)
        if (cdabs(integral) > mo_integrals_threshold) then
          write(i_unit_output,'(2(E25.15,X),4(I6,X))') dble(integral),dimag(integral), i,j,0,0
        endif
      enddo
    enddo
  enddo
  write(i_unit_output,*) core_energy, 0, 0, 0, 0
end
subroutine fcidump_real
  implicit none
  character*(128) :: output
  integer :: i_unit_output,getUnitAndOpen
  output=trim(ezfio_filename)//'.FCIDUMP'
  i_unit_output = getUnitAndOpen(output,'w')

  integer :: i,j,k,l
  integer :: i1,j1,k1,l1
  integer :: i2,j2,k2,l2
  integer*8 :: m
  character*(2), allocatable :: A(:)

  write(i_unit_output,*) '&FCI NORB=', n_act_orb, ', NELEC=', elec_num-n_core_orb*2, &
   ', MS2=', (elec_alpha_num-elec_beta_num), ','
  allocate (A(n_act_orb))
  A = '1,'
  write(i_unit_output,*) 'ORBSYM=', (A(i), i=1,n_act_orb)
  write(i_unit_output,*) 'ISYM=0,'
  write(i_unit_output,*) '/'
  deallocate(A)

  integer(key_kind), allocatable :: keys(:)
  double precision, allocatable  :: values(:)
  integer(cache_map_size_kind)   :: n_elements, n_elements_max
  PROVIDE mo_two_e_integrals_in_map

  double precision :: get_two_e_integral, integral

  do l=1,n_act_orb
   l1 = list_act(l)
   do k=1,n_act_orb
    k1 = list_act(k)
    do j=l,n_act_orb
     j1 = list_act(j)
     do i=k,n_act_orb
      i1 = list_act(i)
       if (i1>=j1) then
          integral = get_two_e_integral(i1,j1,k1,l1,mo_integrals_map)
          if (dabs(integral) > mo_integrals_threshold) then
            write(i_unit_output,*) integral, i,k,j,l
          endif
       end if
     enddo
    enddo
   enddo
  enddo

  do j=1,n_act_orb
   j1 = list_act(j)
   do i=j,n_act_orb
    i1 = list_act(i)
      integral = mo_one_e_integrals(i1,j1) + core_fock_operator(i1,j1)
      if (dabs(integral) > mo_integrals_threshold) then
        write(i_unit_output,*) integral, i,j,0,0
      endif
   enddo
  enddo
  write(i_unit_output,*) core_energy, 0, 0, 0, 0
end
