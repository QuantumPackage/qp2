subroutine u_0_S2_u_0_complex(e_0,u_0,n,keys_tmp,Nint,N_st,sze_8)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes e_0 = <u_0|S2|u_0>/<u_0|u_0>
  !
  ! n : number of determinants
  !
  END_DOC
  integer, intent(in)            :: n,Nint, N_st, sze_8
  double precision, intent(out)  :: e_0(N_st)
  complex*16, intent(in)   :: u_0(sze_8,N_st)
  integer(bit_kind),intent(in)   :: keys_tmp(Nint,2,n)

  complex*16, allocatable  :: v_0(:,:)
  double precision               :: u_dot_u_complex
  complex*16                     :: u_dot_v_complex
  integer :: i,j
  allocate (v_0(sze_8,N_st))

  call s2_u_0_nstates_complex(v_0,u_0,n,keys_tmp,Nint,N_st,sze_8)
  do i=1,N_st
    e_0(i) = dble(u_dot_v_complex(u_0(1,i),v_0(1,i),n))/u_dot_u_complex(u_0(1,i),n)
  enddo
end



subroutine S2_u_0_complex(v_0,u_0,n,keys_tmp,Nint)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes v_0 = S^2|u_0>
  !
  ! n : number of determinants
  !
  END_DOC
  integer, intent(in)            :: n,Nint
  complex*16, intent(out)  :: v_0(n)
  complex*16, intent(in)   :: u_0(n)
  integer(bit_kind),intent(in)   :: keys_tmp(Nint,2,n)
  call s2_u_0_nstates_complex(v_0,u_0,n,keys_tmp,Nint,1,n)
end

subroutine S2_u_0_nstates_complex(v_0,u_0,n,keys_tmp,Nint,N_st,sze_8)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes v_0  = S^2|u_0>
  !
  ! n : number of determinants
  !
  END_DOC
  integer, intent(in)            :: N_st,n,Nint, sze_8
  complex*16, intent(out)  :: v_0(sze_8,N_st)
  complex*16, intent(in)   :: u_0(sze_8,N_st)
  integer(bit_kind),intent(in)   :: keys_tmp(Nint,2,n)
  double precision               :: s2_tmp
  complex*16, allocatable  :: vt(:,:)
  integer                        :: i,j,k,l, jj,ii
  integer                        :: i0, j0

  integer, allocatable           :: shortcut(:,:), sort_idx(:,:)
  integer(bit_kind), allocatable :: sorted(:,:,:), version(:,:,:)
  integer(bit_kind)              :: sorted_i(Nint)

  integer                        :: sh, sh2, ni, exa, ext, org_i, org_j, endi, istate


  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  ASSERT (n>0)
  PROVIDE ref_bitmask_energy

  allocate (shortcut(0:n+1,2), sort_idx(n,2), sorted(Nint,n,2), version(Nint,n,2))
  v_0 = (0.d0,0.d0)

  call sort_dets_ab_v(keys_tmp, sorted(1,1,1), sort_idx(1,1), shortcut(0,1), version(1,1,1), n, Nint)
  call sort_dets_ba_v(keys_tmp, sorted(1,1,2), sort_idx(1,2), shortcut(0,2), version(1,1,2), n, Nint)

  !$OMP PARALLEL DEFAULT(NONE)                                       &
      !$OMP PRIVATE(i,s2_tmp,j,k,jj,vt,ii,sh,sh2,ni,exa,ext,org_i,org_j,endi,sorted_i,istate)&
      !$OMP SHARED(n,u_0,keys_tmp,Nint,v_0,sorted,shortcut,sort_idx,version,N_st,sze_8)
  allocate(vt(sze_8,N_st))
  vt = (0.d0,0.d0)

  do sh=1,shortcut(0,1)
    !$OMP DO SCHEDULE(static,1)
    do sh2=sh,shortcut(0,1)
      exa = 0
      do ni=1,Nint
        exa = exa + popcnt(xor(version(ni,sh,1), version(ni,sh2,1)))
      end do
      if(exa > 2) then
        cycle
      end if

      do i=shortcut(sh,1),shortcut(sh+1,1)-1
        org_i = sort_idx(i,1)
        if(sh==sh2) then
          endi = i-1
        else
          endi = shortcut(sh2+1,1)-1
        end if
        do ni=1,Nint
          sorted_i(ni) = sorted(ni,i,1)
        enddo

        do j=shortcut(sh2,1),endi
          org_j = sort_idx(j,1)
          ext = exa
          do ni=1,Nint
            ext = ext + popcnt(xor(sorted_i(ni), sorted(ni,j,1)))
          end do
          if(ext <= 4) then
            call get_s2(keys_tmp(1,1,org_i),keys_tmp(1,1,org_j),Nint,s2_tmp)
            do istate=1,N_st
              vt (org_i,istate) = vt (org_i,istate) + s2_tmp*u_0(org_j,istate)
              vt (org_j,istate) = vt (org_j,istate) + s2_tmp*u_0(org_i,istate)
            enddo
          endif
        enddo
      enddo
    enddo
    !$OMP END DO NOWAIT
  enddo

  do sh=1,shortcut(0,2)
    !$OMP DO
    do i=shortcut(sh,2),shortcut(sh+1,2)-1
      org_i = sort_idx(i,2)
      do j=shortcut(sh,2),i-1
        org_j = sort_idx(j,2)
        ext = 0
        do ni=1,Nint
          ext = ext + popcnt(xor(sorted(ni,i,2), sorted(ni,j,2)))
        end do
        if(ext == 4) then
          call get_s2(keys_tmp(1,1,org_i),keys_tmp(1,1,org_j),Nint,s2_tmp)
          do istate=1,N_st
            vt (org_i,istate) = vt (org_i,istate) + s2_tmp*u_0(org_j,istate)
            vt (org_j,istate) = vt (org_j,istate) + s2_tmp*u_0(org_i,istate)
          enddo
        end if
      end do
    end do
    !$OMP END DO NOWAIT
  enddo
  !$OMP BARRIER

  do istate=1,N_st
    do i=n,1,-1
      !$OMP ATOMIC
      v_0(i,istate) = v_0(i,istate) + vt(i,istate)
    enddo
  enddo

  deallocate(vt)
  !$OMP END PARALLEL

  do i=1,n
    call get_s2(keys_tmp(1,1,i),keys_tmp(1,1,i),Nint,s2_tmp)
    do istate=1,N_st
      v_0(i,istate) += s2_tmp * u_0(i,istate)
    enddo
  enddo

  deallocate (shortcut, sort_idx, sorted, version)
end







subroutine get_uJ_s2_uI_complex(psi_keys_tmp,psi_coefs_tmp,n,nmax_coefs,nmax_keys,s2,nstates)
  !todo: modify/implement for complex
  print*,irp_here,' not implemented for complex'
  stop -1
!  implicit none
!  use bitmasks
!  integer, intent(in)            :: n,nmax_coefs,nmax_keys,nstates
!  integer(bit_kind), intent(in)  :: psi_keys_tmp(N_int,2,nmax_keys)
!  complex*16, intent(in)   :: psi_coefs_tmp(nmax_coefs,nstates)
!  complex*16, intent(out)  :: s2(nstates,nstates)
!  double precision               :: s2_tmp
!  complex*16                     :: accu
!  integer                        :: i,j,l,jj,ll,kk
!  integer, allocatable           :: idx(:)
!  BEGIN_DOC
!  ! returns the matrix elements of S^2 "s2(i,j)" between the "nstates" states
!  ! psi_coefs_tmp(:,i) and psi_coefs_tmp(:,j)
!  END_DOC
!  s2 = (0.d0,0.d0)
!  do ll = 1, nstates
!    do jj = 1, nstates
!      accu = (0.d0,0.d0)
!      !$OMP PARALLEL DEFAULT(NONE)                                   &
!          !$OMP PRIVATE (i,j,kk,idx,s2_tmp)                          &
!          !$OMP SHARED (ll,jj,psi_keys_tmp,psi_coefs_tmp,N_int,n,nstates)&
!          !$OMP REDUCTION(+:accu)
!      allocate(idx(0:n))
!      !$OMP DO SCHEDULE(dynamic)
!      do i = n,1,-1   ! Better OMP scheduling
!        call get_s2(psi_keys_tmp(1,1,i),psi_keys_tmp(1,1,i),N_int,s2_tmp)
!        accu += dconjg(psi_coefs_tmp(i,ll)) * s2_tmp * psi_coefs_tmp(i,jj)
!        call filter_connected(psi_keys_tmp,psi_keys_tmp(1,1,i),N_int,i-1,idx)
!        do kk=1,idx(0)
!          j = idx(kk)
!          call get_s2(psi_keys_tmp(1,1,i),psi_keys_tmp(1,1,j),N_int,s2_tmp)
!          accu += dconjg(psi_coefs_tmp(i,ll)) * s2_tmp * psi_coefs_tmp(j,jj) + psi_coefs_tmp(i,jj) * s2_tmp * psi_coefs_tmp(j,ll)
!        enddo
!      enddo
!      !$OMP END DO
!      deallocate(idx)
!      !$OMP END PARALLEL
!      s2(ll,jj) += accu
!    enddo
!  enddo
!  do i = 1, nstates
!    do j =i+1,nstates
!      accu = 0.5d0 * (s2(i,j) + s2(j,i))
!      s2(i,j) = accu
!      s2(j,i) = accu
!    enddo
!  enddo
end


subroutine i_S2_psi_minilist_complex(key,keys,idx_key,N_minilist,coef,Nint,Ndet,Ndet_max,Nstate,i_S2_psi_array)
  !todo: modify/implement for complex
  print*,irp_here,' not implemented for complex'
  stop -1
!  use bitmasks
!  implicit none
!  integer, intent(in)            :: Nint, Ndet,Ndet_max,Nstate,idx_key(Ndet), N_minilist
!  integer(bit_kind), intent(in)  :: keys(Nint,2,Ndet)
!  integer(bit_kind), intent(in)  :: key(Nint,2)
!  double precision, intent(in)   :: coef(Ndet_max,Nstate)
!  double precision, intent(out)  :: i_S2_psi_array(Nstate)
!
!  integer                        :: i, ii,j, i_in_key, i_in_coef
!  double precision               :: phase
!  integer                        :: exc(0:2,2,2)
!  double precision               :: s2ij
!  integer                        :: idx(0:Ndet)
!  BEGIN_DOC
!! Computes $\langle i|S^2|\Psi \rangle = \sum_J c_J \langle i|S^2|J \rangle$.
!!
!! Uses filter_connected_i_H_psi0 to get all the $|J\rangle$ to which $|i\rangle$
!! is connected. The $|J\rangle$ are searched in short pre-computed lists.
!  END_DOC
!
!  ASSERT (Nint > 0)
!  ASSERT (N_int == Nint)
!  ASSERT (Nstate > 0)
!  ASSERT (Ndet > 0)
!  ASSERT (Ndet_max >= Ndet)
!  i_S2_psi_array = 0.d0
!
!  call filter_connected_i_H_psi0(keys,key,Nint,N_minilist,idx)
!  if (Nstate == 1) then
!
!    do ii=1,idx(0)
!      i_in_key = idx(ii)
!      i_in_coef = idx_key(idx(ii))
!      !DIR$ FORCEINLINE
!      call get_s2(keys(1,1,i_in_key),key,Nint,s2ij)
!      ! TODO : Cache misses
!      i_S2_psi_array(1) = i_S2_psi_array(1) + coef(i_in_coef,1)*s2ij
!    enddo
!
!  else
!
!    do ii=1,idx(0)
!      i_in_key = idx(ii)
!      i_in_coef = idx_key(idx(ii))
!      !DIR$ FORCEINLINE
!      call get_s2(keys(1,1,i_in_key),key,Nint,s2ij)
!      do j = 1, Nstate
!        i_S2_psi_array(j) = i_S2_psi_array(j) + coef(i_in_coef,j)*s2ij
!      enddo
!    enddo
!
!  endif
!
end
