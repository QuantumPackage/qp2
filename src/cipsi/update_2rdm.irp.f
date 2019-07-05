use bitmasks

subroutine give_2rdm_pert_contrib(det,coef,psi_det_connection,psi_coef_connection,n_det_connection,nkeys,keys,values,sze_buff)
 implicit none
 integer, intent(in) :: n_det_connection,nkeys
 double precision, intent(in) :: coef(N_states)
 integer(bit_kind), intent(in) :: det(N_int,2)
 integer(bit_kind), intent(in) :: psi_det_connection(N_int,2,n_det_connection) 
 double precision, intent(in)  :: psi_coef_connection(n_det_connection, N_states)
 integer,           intent(inout) :: keys(4,sze_buff),sze_buff
 double precision,  intent(inout) :: values(sze_buff)

end
