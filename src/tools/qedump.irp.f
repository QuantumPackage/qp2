program qedump
  read_wf=.True.
  touch read_wf
  call fdump
  call vecdump
  call edump
end
subroutine fdump
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
  character*(128) :: output
  integer :: iunit,getUnitAndOpen
  output=trim(ezfio_filename)//'.FCIDUMP'
  iunit = getUnitAndOpen(output,'w')

  integer :: i,j,k,l
  integer :: i1,j1,k1,l1
  integer :: i2,j2,k2,l2
  integer*8 :: m
  character*(2), allocatable :: A(:)

  write(iunit,*) '&FCI NORB=', n_act_orb, ', NELEC=', elec_num-n_core_orb*2, &
   ', MS2=', (elec_alpha_num-elec_beta_num), ','
  allocate (A(n_act_orb))
  A = '1,'
  write(iunit,*) 'ORBSYM=', (A(i), i=1,n_act_orb)
  write(iunit,*) 'ISYM=0,'
  write(iunit,*) '/'
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
            write(iunit,*) integral, i,k,j,l
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
        write(iunit,*) integral, i,j,0,0
      endif
   enddo
  enddo
  write(iunit,*) core_energy, 0, 0, 0, 0
  close(iunit)
end


subroutine vecdump
  use bitmasks
  implicit none
  BEGIN_DOC
!
  END_DOC
  character*(128) :: output
  integer :: iunit,getUnitAndOpen
  output=trim(ezfio_filename)//'.WFDUMP'
  iunit = getUnitAndOpen(output,'w')

  integer :: i
  character*(2048)                :: output2(2)

  do i=1,N_det
    write(iunit,'(E25.17)') psi_coef(i,1)

    call bitstring_to_str_trim( output2(1), psi_det_sorted(1,1,i), N_int )
    call bitstring_to_str_trim( output2(2), psi_det_sorted(1,2,i), N_int )
    write(iunit,*)  trim(output2(1))
    write(iunit,*)  trim(output2(2))
    write(iunit,*)''
  enddo
  close(iunit)
end

subroutine edump
  implicit none
  BEGIN_DOC
!
  END_DOC
  character*(128) :: output
  integer :: iunit,getUnitAndOpen
  output=trim(ezfio_filename)//'.EDUMP'
  iunit = getUnitAndOpen(output,'w')

  integer, save :: ifirst = 0
  double precision :: Vee, Ven, Vnn, Vecp, T, f
  integer  :: i,j,k
  provide n_states

  double precision :: e(n_states), ept2(n_states), pt2(n_states)
  Vnn = nuclear_repulsion

  call ezfio_get_fci_energy(e)
  call ezfio_get_fci_energy_pt2(ept2)
  write(iunit,*) 'Energy components'
  write(iunit,*) '================='
  write(iunit,*) ''
  do k=1,N_states

    Ven  = 0.d0
    Vecp = 0.d0
    T    = 0.d0

    do j=1,mo_num
      do i=1,mo_num
        f = one_e_dm_mo_alpha(i,j,k) + one_e_dm_mo_beta(i,j,k)
        Ven  = Ven  + f * mo_integrals_n_e(i,j)
        Vecp = Vecp + f * mo_pseudo_integrals(i,j)
        T    = T    + f * mo_kinetic_integrals(i,j)
      enddo
    enddo
    Vee = psi_energy(k) - Ven - Vecp - T
    
    if (ifirst == 0) then
      ifirst = 1
      write(iunit,*) 'Vnn  : Nucleus-Nucleus   potential energy'
      write(iunit,*) 'Ven  : Electron-Nucleus  potential energy'
      write(iunit,*) 'Vee  : Electron-Electron potential energy'
      write(iunit,*) 'Vecp : Potential energy of the pseudo-potentials'
      write(iunit,*) 'T    : Electronic kinetic energy'
      write(iunit,*) ''
    endif

    write(iunit,*) 'State ', k
    write(iunit,*) '---------'
    write(iunit,*) ''
    write(iunit,*) 'Vnn  = ', Vnn
    write(iunit,*) 'Ven  = ', Ven
    write(iunit,*) 'Vee  = ', Vee
    write(iunit,*) 'Vecp = ', Vecp
    write(iunit,*) 'T    = ', T
    write(iunit,*) 'Evar = ',e(k)
    write(iunit,*) 'E+pt2= ',ept2(k)
    write(iunit,*) 'pt2  = ',ept2(k) - e(k)
    write(iunit,*) ''
  enddo

  write(iunit,*) ''
  close(iunit)
end

subroutine bitstring_to_str_trim( output, string, Nint )
  use bitmasks
  implicit none
  BEGIN_DOC
! Transform a bit string to a string for printing
  END_DOC
  character*(*), intent(out)     :: output
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: string(Nint)

  integer                        :: i, j, ibuf
  integer(bit_kind)              :: itemp

  ibuf = 1
  output = ''
!  output(ibuf:ibuf) = '|'
!  ibuf = ibuf+1
  do i=1,Nint
    itemp = 1_bit_kind
    do j=1,bit_kind_size
      if (iand(itemp,string(i)) == itemp) then
        output(ibuf:ibuf) = '+'
      else
        output(ibuf:ibuf) = '-'
      endif
      ibuf = ibuf+1
      itemp = shiftl(itemp,1)
    enddo
  enddo
!  output(ibuf:ibuf) = '|'
end
