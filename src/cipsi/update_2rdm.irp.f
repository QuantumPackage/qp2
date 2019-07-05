use bitmasks

subroutine give_2rdm_pert_contrib(det,coef,psi_det_connection,psi_coef_connection,n_det_connection,nkeys,keys,values,sze_buff)
 implicit none
 integer, intent(in) :: n_det_connection,sze_buff
 double precision, intent(in) :: coef(N_states)
 integer(bit_kind), intent(in) :: det(N_int,2)
 integer(bit_kind), intent(in) :: psi_det_connection(N_int,2,n_det_connection) 
 double precision, intent(in)  :: psi_coef_connection(n_det_connection, N_states)
 integer,           intent(inout) :: keys(4,sze_buff),nkeys
 double precision,  intent(inout) :: values(sze_buff)
 integer :: i,j
 integer                        :: exc(0:2,2,2)
 integer                        :: degree
 double precision               :: phase, contrib
 do i = 1, n_det_connection
  call get_excitation(det,psi_det_connection(1,1,i),exc,degree,phase,N_int) 
  if(degree.gt.2)cycle
  contrib = 0.d0
  do j = 1, N_states
   contrib += state_average_weight(j) * psi_coef_connection(i,j) * phase * coef(j)
  enddo
  ! case of single excitations 
  if(degree == 1)then
   if (nkeys+ 2 * elec_alpha_num .ge. sze_buff)then
    call update_rdms(nkeys,keys,values,sze_buff)
    nkeys = 0
   endif
   call update_buffer_single_exc_rdm(det,psi_det_connection(1,1,i),exc,phase,contrib,nkeys,keys,values,sze_buff)
  else 
  ! case of double excitations 
   if (nkeys+ 4  .ge. sze_buff)then
    call update_rdms(nkeys,keys,values,sze_buff)
    nkeys = 0
   endif
   call update_buffer_double_exc_rdm(exc,phase,contrib,nkeys,keys,values,sze_buff)
  endif
 enddo

end

subroutine update_buffer_single_exc_rdm(det1,det2,exc,phase,contrib,nkeys,keys,values,sze_buff)
 implicit none
 integer, intent(in) :: sze_buff
 integer(bit_kind), intent(in) :: det1(N_int,2)
 integer(bit_kind), intent(in) :: det2(N_int,2)
 integer,intent(in)             :: exc(0:2,2,2)
 double precision,intent(in)    :: phase, contrib
 integer, intent(inout)         :: nkeys, keys(4,sze_buff)
 double precision, intent(inout):: values(sze_buff)
 


end

subroutine update_buffer_double_exc_rdm(exc,phase,contrib,nkeys,keys,values,sze_buff)
 implicit none
 integer, intent(in) :: sze_buff
 integer,intent(in)             :: exc(0:2,2,2)
 double precision,intent(in)    :: phase, contrib
 integer, intent(inout)         :: nkeys, keys(4,sze_buff)
 double precision, intent(inout):: values(sze_buff)


end


subroutine update_rdms(nkeys,keys,values,sze_buff)
 implicit none
 integer, intent(in) :: nkeys, keys(4,sze_buff),sze_buff
 double precision, intent(in) :: values(sze_buff)
 integer :: i,h1,h2,p1,p2
 do i = 1, nkeys
  h1 = keys(1,i)
  h2 = keys(2,i)
  p1 = keys(3,i)
  p2 = keys(4,i)
  pert_2rdm_provider(h1,h2,p1,p2) += values(i)
 enddo
end
