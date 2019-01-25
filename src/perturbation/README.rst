============
perturbation
============


All subroutines in ``*.irp.f`` starting with `pt2_` in the current directory are
perturbation computed using the routine `i_H_psi`. Other cases are not allowed.
The arguments of the `pt2_` are always:

.. code-block:: fortran

  subroutine pt2_...(                                          &
      psi_ref,                                                 &
      psi_ref_coefs,                                           &
      E_refs,                                                  &
      det_pert,                                                &
      c_pert,                                                  &
      e_2_pert,                                                &
      H_pert_diag,                                             &
      Nint,                                                    &
      Ndet,                                                    &
      N_st )


  integer          , intent(in)  :: Nint,Ndet,N_st
  integer(bit_kind), intent(in)  :: psi_ref(Nint,2,Ndet)
  double precision , intent(in)  :: psi_ref_coefs(Ndet,N_st)
  double precision , intent(in)  :: E_refs(N_st)
  integer(bit_kind), intent(in)  :: det_pert(Nint,2)
  double precision , intent(out) :: c_pert(N_st),e_2_pert(N_st),H_pert_diag


`psi_ref`
  bitstring of the determinants present in the various `N_st` states

`psi_ref_coefs`
  coefficients of the determinants on the various `N_st` states

`E_refs`
  Energy of the various `N_st` states

`det_pert`
  Perturber determinant

`c_pert`
  Perturbative coefficients for the various states

`e_2_pert`
  Perturbative energetic contribution for the various states

`H_pert_diag`
  Diagonal |H| matrix element of the perturber

`Nint`
  Should be equal to `N_int`

`Ndet`
  Number of determinants `i` in |Psi| on which we apply <det_pert | |H| | `i`>

`N_st`
  Number of states



