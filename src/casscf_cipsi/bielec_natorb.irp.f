 BEGIN_PROVIDER [real*8, bielec_PQxx_no_array, (mo_num, mo_num,n_core_inact_act_orb,n_core_inact_act_orb)]
  BEGIN_DOC
  ! WARNING !!! Old version !!! NOT USED ANYMORE IN THE PROGRAM !!! TOO BIG TO BE STORED ON LARGE SYSTEMS !!! 
  ! 
  ! Replaced by the Cholesky-based function bielec_PQxx_no
  !
  ! integral (pq|xx) in the basis of natural MOs
  ! indices are unshifted orbital numbers
  ! 
  END_DOC
  implicit none
  integer                        :: i,j,k,l,t,u,p,q
  double precision, allocatable  :: f(:,:,:), d(:,:,:)

  print*,''
  print*,'Providing bielec_PQxx_no_array, WARNING IT CAN BE A VERY BIG ARRAY WHEN MO_NUM IS LARGE !!!'
  print*,''


  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(j,k,l,p,d,f) &
  !$OMP SHARED(n_core_inact_act_orb,mo_num,n_act_orb,n_core_inact_orb, &
  !$OMP   bielec_PQxx_no_array,bielec_PQxx_array,list_act,natorbsCI)

  allocate (f(n_act_orb,mo_num,n_core_inact_act_orb), &
      d(n_act_orb,mo_num,n_core_inact_act_orb))

  !$OMP DO
  do l=1,n_core_inact_act_orb
    bielec_PQxx_no_array(:,:,:,l) = bielec_PQxx_array(:,:,:,l)

    do k=1,n_core_inact_act_orb
      do j=1,mo_num
        do p=1,n_act_orb
          f(p,j,k)=bielec_PQxx_no_array(list_act(p),j,k,l)
        end do
      end do
    end do
    call dgemm('T','N',n_act_orb,mo_num*n_core_inact_act_orb,n_act_orb,1.d0,  &
          natorbsCI, size(natorbsCI,1),                              &
          f, n_act_orb,                                              &
          0.d0,                                                      &
          d, n_act_orb)
    do k=1,n_core_inact_act_orb
      do j=1,mo_num
        do p=1,n_act_orb
          bielec_PQxx_no_array(list_act(p),j,k,l)=d(p,j,k)
        end do
      end do

      do j=1,mo_num
        do p=1,n_act_orb
          f(p,j,k)=bielec_PQxx_no_array(j,list_act(p),k,l)
        end do
      end do
    end do
    call dgemm('T','N',n_act_orb,mo_num*n_core_inact_act_orb,n_act_orb,1.d0,  &
          natorbsCI, n_act_orb,                                      &
          f, n_act_orb,                                              &
          0.d0,                                                      &
          d, n_act_orb)
    do k=1,n_core_inact_act_orb
      do p=1,n_act_orb
        do j=1,mo_num
          bielec_PQxx_no_array(j,list_act(p),k,l)=d(p,j,k)
        end do
      end do
    end do
  end do
  !$OMP END DO NOWAIT

  deallocate (f,d)

  allocate (f(mo_num,mo_num,n_act_orb),d(mo_num,mo_num,n_act_orb))

  !$OMP DO 
  do l=1,n_core_inact_act_orb

    do p=1,n_act_orb
      do k=1,mo_num
        do j=1,mo_num
          f(j,k,p) = bielec_PQxx_no_array(j,k,n_core_inact_orb+p,l)
        end do
      end do
    end do
    call dgemm('N','N',mo_num*mo_num,n_act_orb,n_act_orb,1.d0,       &
          f, mo_num*mo_num,                                          &
          natorbsCI, n_act_orb,                                      &
          0.d0,                                                      &
          d, mo_num*mo_num)
    do p=1,n_act_orb
      do k=1,mo_num
        do j=1,mo_num
          bielec_PQxx_no_array(j,k,n_core_inact_orb+p,l)=d(j,k,p)
        end do
      end do
    end do
  end do
  !$OMP END DO NOWAIT

  !$OMP BARRIER 

  !$OMP DO
  do l=1,n_core_inact_act_orb
    do p=1,n_act_orb
      do k=1,mo_num
        do j=1,mo_num
          f(j,k,p) = bielec_PQxx_no_array(j,k,l,n_core_inact_orb+p)
        end do
      end do
    end do
    call dgemm('N','N',mo_num*mo_num,n_act_orb,n_act_orb,1.d0,       &
          f, mo_num*mo_num,                                          &
          natorbsCI, n_act_orb,                                      &
          0.d0,                                                      &
          d, mo_num*mo_num)
    do p=1,n_act_orb
      do k=1,mo_num
        do j=1,mo_num
          bielec_PQxx_no_array(j,k,l,n_core_inact_orb+p)=d(j,k,p)
        end do
      end do
    end do
  end do
  !$OMP END DO

  deallocate (f,d)
  !$OMP END PARALLEL

END_PROVIDER



BEGIN_PROVIDER [real*8, bielec_PxxQ_no_array, (mo_num,n_core_inact_act_orb,n_core_inact_act_orb, mo_num)]
  BEGIN_DOC
  ! WARNING !!! Old version !!! NOT USED ANYMORE IN THE PROGRAM !!! TOO BIG TO BE STORED ON LARGE SYSTEMS !!! 
  ! 
  ! Replaced by the Cholesky-based function bielec_PxxQ_no
  !
  ! integral (px|xq) in the basis of natural MOs
  ! indices are unshifted orbital numbers
  END_DOC
  implicit none
  integer                        :: i,j,k,l,t,u,p,q
  double precision, allocatable  :: f(:,:,:), d(:,:,:)

  print*,''
  print*,'Providing bielec_PxxQ_no_array, WARNING IT CAN BE A VERY BIG ARRAY WHEN MO_NUM IS LARGE !!!'
  print*,''

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(j,k,l,p,d,f) &
  !$OMP SHARED(n_core_inact_act_orb,mo_num,n_act_orb,n_core_inact_orb, &
  !$OMP   bielec_PxxQ_no_array,bielec_PxxQ_array,list_act,natorbsCI)


  allocate (f(n_act_orb,n_core_inact_act_orb,n_core_inact_act_orb), &
      d(n_act_orb,n_core_inact_act_orb,n_core_inact_act_orb))

  !$OMP DO
  do j=1,mo_num
    bielec_PxxQ_no_array(:,:,:,j) = bielec_PxxQ_array(:,:,:,j)
    do l=1,n_core_inact_act_orb
      do k=1,n_core_inact_act_orb
        do p=1,n_act_orb
            f(p,k,l) = bielec_PxxQ_no_array(list_act(p),k,l,j)
        end do
      end do
    end do
    call dgemm('T','N',n_act_orb,n_core_inact_act_orb**2,n_act_orb,1.d0,  &
          natorbsCI, size(natorbsCI,1),                              &
          f, n_act_orb,                                              &
          0.d0,                                                      &
          d, n_act_orb)
    do l=1,n_core_inact_act_orb
      do k=1,n_core_inact_act_orb
        do p=1,n_act_orb
          bielec_PxxQ_no_array(list_act(p),k,l,j)=d(p,k,l)
        end do
      end do
    end do
  end do
  !$OMP END DO NOWAIT

  deallocate (f,d)

  allocate (f(n_act_orb,mo_num,n_core_inact_act_orb), &
    d(n_act_orb,mo_num,n_core_inact_act_orb))

  !$OMP DO
  do k=1,mo_num
    do l=1,n_core_inact_act_orb
      do j=1,mo_num
        do p=1,n_act_orb
          f(p,j,l) = bielec_PxxQ_no_array(j,n_core_inact_orb+p,l,k)
        end do
      end do
    end do
    call dgemm('T','N',n_act_orb,mo_num*n_core_inact_act_orb,n_act_orb,1.d0,  &
          natorbsCI, size(natorbsCI,1),                              &
          f, n_act_orb,                                              &
          0.d0,                                                      &
          d, n_act_orb)
    do l=1,n_core_inact_act_orb
      do j=1,mo_num
        do p=1,n_act_orb
          bielec_PxxQ_no_array(j,n_core_inact_orb+p,l,k)=d(p,j,l)
        end do
      end do
    end do
  end do
  !$OMP END DO NOWAIT

  deallocate(f,d)

  allocate(f(mo_num,n_core_inact_act_orb,n_act_orb), &
    d(mo_num,n_core_inact_act_orb,n_act_orb) )

  !$OMP DO
  do k=1,mo_num
    do p=1,n_act_orb
      do l=1,n_core_inact_act_orb
        do j=1,mo_num
          f(j,l,p) = bielec_PxxQ_no_array(j,l,n_core_inact_orb+p,k)
        end do
      end do
    end do
    call dgemm('N','N',mo_num*n_core_inact_act_orb,n_act_orb,n_act_orb,1.d0,  &
          f, mo_num*n_core_inact_act_orb,                                       &
          natorbsCI, size(natorbsCI,1),                              &
          0.d0,                                                      &
          d, mo_num*n_core_inact_act_orb)
    do p=1,n_act_orb
      do l=1,n_core_inact_act_orb
        do j=1,mo_num
          bielec_PxxQ_no_array(j,l,n_core_inact_orb+p,k)=d(j,l,p)
        end do
      end do
    end do
  end do
  !$OMP END DO NOWAIT

  !$OMP BARRIER

  !$OMP DO 
  do l=1,n_core_inact_act_orb
    do p=1,n_act_orb
      do k=1,n_core_inact_act_orb
        do j=1,mo_num
          f(j,k,p) = bielec_PxxQ_no_array(j,k,l,n_core_inact_orb+p)
        end do
      end do
    end do
    call dgemm('N','N',mo_num*n_core_inact_act_orb,n_act_orb,n_act_orb,1.d0,  &
          f, mo_num*n_core_inact_act_orb,                                       &
          natorbsCI, size(natorbsCI,1),                              &
          0.d0,                                                      &
          d, mo_num*n_core_inact_act_orb)
    do p=1,n_act_orb
      do k=1,n_core_inact_act_orb
        do j=1,mo_num
          bielec_PxxQ_no_array(j,k,l,n_core_inact_orb+p)=d(j,k,p)
        end do
      end do
    end do
  end do
  !$OMP END DO NOWAIT
  deallocate(f,d)
  !$OMP END PARALLEL

END_PROVIDER


BEGIN_PROVIDER [real*8, bielecCI_no, (n_act_orb,n_act_orb,n_act_orb, mo_num)]
  BEGIN_DOC
  ! integrals (tu|vp) in the basis of natural MOs
  ! index p runs over the whole basis, t,u,v only over the active orbitals
  ! 
  ! This array can be stored anyway. Ex: 50 active orbitals, 1500 MOs ==> 8x50^3x1500 = 1.5 Gb
  END_DOC
  implicit none
  integer                        :: i,j,k,l,t,u,p,q
  double precision, allocatable  :: f(:,:,:), d(:,:,:)
 
  double precision :: wall0, wall1
  call wall_time(wall0)
  print*,'Providing bielecCI_no'
  
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(j,k,l,p,d,f) &
  !$OMP SHARED(n_core_inact_act_orb,mo_num,n_act_orb,n_core_inact_orb, &
  !$OMP   bielecCI_no,bielecCI,list_act,natorbsCI)

  allocate (f(n_act_orb,n_act_orb,mo_num), &
      d(n_act_orb,n_act_orb,mo_num))

  !$OMP DO
  do l=1,mo_num
    bielecCI_no(:,:,:,l) = bielecCI(:,:,:,l)
    do k=1,n_act_orb
      do j=1,n_act_orb
        do p=1,n_act_orb
          f(p,j,k)=bielecCI_no(p,j,k,l)
        end do
      end do
    end do
    call dgemm('T','N',n_act_orb,n_act_orb*n_act_orb,n_act_orb,1.d0,  &
          natorbsCI, size(natorbsCI,1),                              &
          f, n_act_orb,                                              &
          0.d0,                                                      &
          d, n_act_orb)
    do k=1,n_act_orb
      do j=1,n_act_orb
        do p=1,n_act_orb
          bielecCI_no(p,j,k,l)=d(p,j,k)
        end do
      end do

      do j=1,n_act_orb
        do p=1,n_act_orb
          f(p,j,k)=bielecCI_no(j,p,k,l)
        end do
      end do
    end do
    call dgemm('T','N',n_act_orb,n_act_orb*n_act_orb,n_act_orb,1.d0,  &
          natorbsCI, n_act_orb,                                      &
          f, n_act_orb,                                              &
          0.d0,                                                      &
          d, n_act_orb)
    do k=1,n_act_orb
      do p=1,n_act_orb
        do j=1,n_act_orb
          bielecCI_no(j,p,k,l)=d(p,j,k)
        end do
      end do
    end do

    do p=1,n_act_orb
      do k=1,n_act_orb
        do j=1,n_act_orb
          f(j,k,p)=bielecCI_no(j,k,p,l)
        end do
      end do
    end do
    call dgemm('N','N',n_act_orb*n_act_orb,n_act_orb,n_act_orb,1.d0,       &
          f, n_act_orb*n_act_orb,                                          &
          natorbsCI, n_act_orb,                                      &
          0.d0,                                                      &
          d, n_act_orb*n_act_orb)

    do p=1,n_act_orb
      do k=1,n_act_orb
        do j=1,n_act_orb
          bielecCI_no(j,k,p,l)=d(j,k,p)
        end do
      end do 
    end do  
  end do  
  !$OMP END DO

  !$OMP DO
  do l=1,n_act_orb
    do p=1,n_act_orb
      do k=1,n_act_orb
        do j=1,n_act_orb
          f(j,k,p)=bielecCI_no(j,k,l,list_act(p))
        end do
      end do
    end do
    call dgemm('N','N',n_act_orb*n_act_orb,n_act_orb,n_act_orb,1.d0,       &
          f, n_act_orb*n_act_orb,                                          &
          natorbsCI, n_act_orb,                                      &
          0.d0,                                                      &
          d, n_act_orb*n_act_orb)

    do p=1,n_act_orb
      do k=1,n_act_orb
        do j=1,n_act_orb
          bielecCI_no(j,k,l,list_act(p))=d(j,k,p)
        end do
      end do 
    end do  
  end do  
  !$OMP END DO
  
  deallocate(d,f)
  !$OMP END PARALLEL
  call wall_time(wall1)
  print*,'Time to provide bielecCI_no = ',wall1-wall0


END_PROVIDER

