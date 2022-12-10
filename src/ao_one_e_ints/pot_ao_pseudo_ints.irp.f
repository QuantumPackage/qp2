BEGIN_PROVIDER [ double precision, ao_pseudo_integrals, (ao_num,ao_num)]
  implicit none
  BEGIN_DOC
  ! Pseudo-potential integrals in the |AO| basis set.
  END_DOC

  if (read_ao_integrals_pseudo) then
    call ezfio_get_ao_one_e_ints_ao_integrals_pseudo(ao_pseudo_integrals)
    print *,  'AO pseudopotential integrals read from disk'
  else

    ao_pseudo_integrals = 0.d0
    if (do_pseudo) then
      if (pseudo_klocmax > 0) then
        ao_pseudo_integrals += ao_pseudo_integrals_local
      endif
      if (pseudo_kmax > 0) then
        ao_pseudo_integrals += ao_pseudo_integrals_non_local
      endif
    endif
  endif

  if (write_ao_integrals_pseudo) then
    call ezfio_set_ao_one_e_ints_ao_integrals_pseudo(ao_pseudo_integrals)
    print *,  'AO pseudopotential integrals written to disk'
  endif

END_PROVIDER

BEGIN_PROVIDER [ double precision, ao_pseudo_integrals_local, (ao_num,ao_num)]
  implicit none
  BEGIN_DOC
  ! Local pseudo-potential
  END_DOC
  include 'utils/constants.include.F'
  double precision               :: alpha, beta, gama, delta
  integer                        :: num_A,num_B
  double precision               :: A_center(3),B_center(3),C_center(3)
  integer                        :: power_A(3),power_B(3)
  integer                        :: i,j,k,l,m
  double precision               :: Vloc, Vpseudo

  double precision               :: wall_1, wall_2, wall_0
  integer                        :: thread_num
  integer, external              :: omp_get_thread_num
  double precision               :: c
  double precision               :: Z

  PROVIDE ao_coef_normalized_ordered_transp
  PROVIDE pseudo_v_k_transp pseudo_n_k_transp pseudo_klocmax pseudo_dz_k_transp

  ao_pseudo_integrals_local = 0.d0

  print*, 'Providing the nuclear electron pseudo integrals (local)'

  ! Dummy iteration for OpenMP
  j=1
  i=1
  l=1
  m=1
  num_A = ao_nucl(j)
  power_A(1:3)= ao_power(j,1:3)
  A_center(1:3) = nucl_coord(num_A,1:3)
  num_B = ao_nucl(i)
  power_B(1:3)= ao_power(i,1:3)
  B_center(1:3) = nucl_coord(num_B,1:3)
  alpha = ao_expo_ordered_transp(l,j)
  beta = ao_expo_ordered_transp(m,i)
  c = 0.d0
  do  k = 1, nucl_num
    Z = nucl_charge(k)

    C_center(1:3) = nucl_coord(k,1:3)

    c = c + Vloc(pseudo_klocmax,                                     &
        pseudo_v_k_transp (1,k),                                     &
        pseudo_n_k_transp (1,k),                                     &
        pseudo_dz_k_transp(1,k),                                     &
        A_center,power_A,alpha,B_center,power_B,beta,C_center)

  enddo


  ao_pseudo_integrals_local = 0.d0
  call wall_time(wall_1)

  thread_num = 0
  !$OMP PARALLEL                                                     &
      !$OMP DEFAULT (NONE)                                           &
      !$OMP PRIVATE (i,j,k,l,m,alpha,beta,A_center,B_center,C_center,power_A,power_B,&
      !$OMP          num_A,num_B,Z,c,                                &
      !$OMP          wall_0,wall_2,thread_num)                       &
      !$OMP SHARED (ao_num,ao_prim_num,ao_expo_ordered_transp,ao_power,ao_nucl,nucl_coord,ao_coef_normalized_ordered_transp,&
      !$OMP         ao_pseudo_integrals_local,nucl_num,nucl_charge,   &
      !$OMP         pseudo_klocmax,pseudo_lmax,pseudo_kmax,pseudo_v_k_transp,pseudo_n_k_transp, pseudo_dz_k_transp,&
      !$OMP         wall_1)

  !$ thread_num = omp_get_thread_num()

  wall_0 = wall_1
  !$OMP DO

  do j = 1, ao_num

    num_A = ao_nucl(j)
    power_A(1:3)= ao_power(j,1:3)
    A_center(1:3) = nucl_coord(num_A,1:3)

    do i = 1, ao_num

      num_B = ao_nucl(i)
      power_B(1:3)= ao_power(i,1:3)
      B_center(1:3) = nucl_coord(num_B,1:3)

      do l=1,ao_prim_num(j)
        alpha = ao_expo_ordered_transp(l,j)

        do m=1,ao_prim_num(i)
          beta = ao_expo_ordered_transp(m,i)
          c = 0.d0

          if (dabs(ao_coef_normalized_ordered_transp(l,j)*ao_coef_normalized_ordered_transp(m,i))&
                < thresh) then
            cycle
          endif
          do  k = 1, nucl_num
            Z = nucl_charge(k)

            C_center(1:3) = nucl_coord(k,1:3)

            c = c + Vloc(pseudo_klocmax, &
                pseudo_v_k_transp (1,k), &
                pseudo_n_k_transp (1,k), &
                pseudo_dz_k_transp(1,k), &
                A_center,power_A,alpha,B_center,power_B,beta,C_center)

          enddo
          ao_pseudo_integrals_local(i,j) = ao_pseudo_integrals_local(i,j) +&
              ao_coef_normalized_ordered_transp(l,j)*ao_coef_normalized_ordered_transp(m,i)*c
        enddo
      enddo
    enddo

    call wall_time(wall_2)
    if (thread_num == 0) then
      if (wall_2 - wall_0 > 1.d0) then
        wall_0 = wall_2
        print*, 100.*float(j)/float(ao_num), '% in ',                &
            wall_2-wall_1, 's'
      endif
    endif
  enddo

 !$OMP END DO
 !$OMP END PARALLEL

 END_PROVIDER


 BEGIN_PROVIDER [ double precision, ao_pseudo_integrals_non_local, (ao_num,ao_num)]
  implicit none
  BEGIN_DOC
  ! Non-local pseudo-potential
  END_DOC
  include 'utils/constants.include.F'
  double precision               :: alpha, beta, gama, delta
  integer                        :: num_A,num_B
  double precision               :: A_center(3),B_center(3),C_center(3)
  integer                        :: power_A(3),power_B(3)
  integer                        :: i,j,k,l,m
  double precision               :: Vloc, Vpseudo
  integer, external              :: omp_get_thread_num

  double precision               :: wall_1, wall_2, wall_0
  integer                        :: thread_num
  double precision               :: c
  double precision               :: Z

  PROVIDE ao_coef_normalized_ordered_transp
  PROVIDE pseudo_lmax pseudo_kmax pseudo_v_kl_transp pseudo_n_kl_transp pseudo_dz_kl_transp
  ao_pseudo_integrals_non_local = 0.d0

  print*, 'Providing the nuclear electron pseudo integrals (non-local)'

  call wall_time(wall_1)
  thread_num = 0

  !$OMP PARALLEL                                                     &
      !$OMP DEFAULT (NONE)                                           &
      !$OMP PRIVATE (i,j,k,l,m,alpha,beta,A_center,B_center,C_center,power_A,power_B,&
      !$OMP          num_A,num_B,Z,c,                                &
      !$OMP          wall_0,wall_2,thread_num)                       &
      !$OMP SHARED (ao_num,ao_prim_num,ao_expo_ordered_transp,ao_power,ao_nucl,nucl_coord,ao_coef_normalized_ordered_transp,&
      !$OMP         ao_pseudo_integrals_non_local,nucl_num,nucl_charge,&
      !$OMP         pseudo_klocmax,pseudo_lmax,pseudo_kmax,pseudo_n_kl_transp, pseudo_v_kl_transp, pseudo_dz_kl_transp,&
      !$OMP         wall_1)

  !$ thread_num = omp_get_thread_num()

  wall_0 = wall_1
  !$OMP DO SCHEDULE (guided)
!
  do j = 1, ao_num

    num_A = ao_nucl(j)
    power_A(1:3)= ao_power(j,1:3)
    A_center(1:3) = nucl_coord(num_A,1:3)

    do i = 1, ao_num

      num_B = ao_nucl(i)
      power_B(1:3)= ao_power(i,1:3)
      B_center(1:3) = nucl_coord(num_B,1:3)

      do l=1,ao_prim_num(j)
        alpha = ao_expo_ordered_transp(l,j)

        do m=1,ao_prim_num(i)
          beta = ao_expo_ordered_transp(m,i)
          c = 0.d0

          if (dabs(ao_coef_normalized_ordered_transp(l,j)*ao_coef_normalized_ordered_transp(m,i))&
                < thresh) then
            cycle
          endif

          do  k = 1, nucl_num
            Z = nucl_charge(k)

            C_center(1:3) = nucl_coord(k,1:3)

            c = c + Vpseudo(pseudo_lmax,pseudo_kmax, &
                    pseudo_v_kl_transp(1,0,k),  &
                    pseudo_n_kl_transp(1,0,k),  &
                    pseudo_dz_kl_transp(1,0,k), &
                    A_center,power_A,alpha,B_center,power_B,beta,C_center)
          enddo
          ao_pseudo_integrals_non_local(i,j) = ao_pseudo_integrals_non_local(i,j) +&
              ao_coef_normalized_ordered_transp(l,j)*ao_coef_normalized_ordered_transp(m,i)*c
        enddo
      enddo
    enddo

    call wall_time(wall_2)
    if (thread_num == 0) then
      if (wall_2 - wall_0 > 1.d0) then
        wall_0 = wall_2
        print*, 100.*float(j)/float(ao_num), '% in ',                &
            wall_2-wall_1, 's'
      endif
    endif
  enddo

  !$OMP END DO

  !$OMP END PARALLEL


END_PROVIDER

 BEGIN_PROVIDER [ double precision, pseudo_v_k_transp, (pseudo_klocmax,nucl_num) ]
&BEGIN_PROVIDER [ integer         , pseudo_n_k_transp, (pseudo_klocmax,nucl_num) ]
&BEGIN_PROVIDER [ double precision, pseudo_dz_k_transp, (pseudo_klocmax,nucl_num)]
 implicit none
 BEGIN_DOC
 ! Transposed arrays for pseudopotentials
 END_DOC

 integer :: i,j
 do j=1,nucl_num
   do i=1,pseudo_klocmax
     pseudo_v_k_transp (i,j) = pseudo_v_k (j,i)
     pseudo_n_k_transp (i,j) = pseudo_n_k (j,i)
     pseudo_dz_k_transp(i,j) = pseudo_dz_k(j,i)
   enddo
 enddo
END_PROVIDER

 BEGIN_PROVIDER [ double precision, pseudo_v_kl_transp, (pseudo_kmax,0:pseudo_lmax,nucl_num) ]
&BEGIN_PROVIDER [ integer         , pseudo_n_kl_transp, (pseudo_kmax,0:pseudo_lmax,nucl_num) ]
&BEGIN_PROVIDER [ double precision, pseudo_dz_kl_transp, (pseudo_kmax,0:pseudo_lmax,nucl_num)]
 implicit none
 BEGIN_DOC
 ! Transposed arrays for pseudopotentials
 END_DOC

 integer :: i,j,l
 do j=1,nucl_num
   do l=0,pseudo_lmax
     do i=1,pseudo_kmax
       pseudo_v_kl_transp (i,l,j) = pseudo_v_kl (j,i,l)
       pseudo_n_kl_transp (i,l,j) = pseudo_n_kl (j,i,l)
       pseudo_dz_kl_transp(i,l,j) = pseudo_dz_kl(j,i,l)
     enddo
   enddo
 enddo
END_PROVIDER

