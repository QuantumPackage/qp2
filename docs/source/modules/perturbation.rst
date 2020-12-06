.. _module_perturbation: 
 
.. program:: perturbation 
 
.. default-role:: option 
 
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



 
 
 
EZFIO parameters 
---------------- 
 
.. option:: do_pt2
 
    If `True`, compute the |PT2| contribution
 
    Default: True
 
.. option:: pt2_max
 
    The selection process stops when the largest |PT2| (for all the state) is lower than `pt2_max` in absolute value
 
    Default: 0.0001
 
.. option:: variance_max
 
    The selection process stops when the largest variance (for all the state) is lower than `variance_max` in absolute value
 
    Default: 0.0
 
.. option:: pt2_relative_error
 
    Stop stochastic |PT2| when the relative error is smaller than `pT2_relative_error`
 
    Default: 0.002
 
.. option:: correlation_energy_ratio_max
 
    The selection process stops at a fixed correlation ratio (useful for getting same accuracy between molecules). Defined as :math:`(E_{CI}-E_{HF})/(E_{CI}+E_{PT2} - E_{HF})`.
 
    Default: 1.00
 
.. option:: h0_type
 
    Type of denominator in PT2. [EN | SOP | HF]
 
    Default: EN
 
 
Providers 
--------- 
 
.. c:var:: max_exc_pert


    File : :file:`perturbation/exc_max.irp.f`

    .. code:: fortran

        integer	:: max_exc_pert	




 
.. c:var:: selection_criterion


    File : :file:`perturbation/selection.irp.f`

    .. code:: fortran

        double precision	:: selection_criterion	
        double precision	:: selection_criterion_min	
        double precision	:: selection_criterion_factor	


    Threshold to select determinants. Set by selection routines.


 
.. c:var:: selection_criterion_factor


    File : :file:`perturbation/selection.irp.f`

    .. code:: fortran

        double precision	:: selection_criterion	
        double precision	:: selection_criterion_min	
        double precision	:: selection_criterion_factor	


    Threshold to select determinants. Set by selection routines.


 
.. c:var:: selection_criterion_min


    File : :file:`perturbation/selection.irp.f`

    .. code:: fortran

        double precision	:: selection_criterion	
        double precision	:: selection_criterion_min	
        double precision	:: selection_criterion_factor	


    Threshold to select determinants. Set by selection routines.


 
.. c:var:: var_pt2_ratio


    File : :file:`perturbation/var_pt2_ratio_provider.irp.f`

    .. code:: fortran

        double precision	:: var_pt2_ratio	


    The selection process stops when the energy ratio variational/(variational+PT2)
    is equal to var_pt2_ratio

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`correlation_energy_ratio_max`


 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: fill_h_apply_buffer_selection:


    File : :file:`perturbation/selection.irp.f`

    .. code:: fortran

        subroutine fill_H_apply_buffer_selection(n_selected,det_buffer,e_2_pert_buffer,coef_pert_buffer, N_st,Nint,iproc,select_max_out)


    Fill the H_apply buffer with determiants for the selection

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`h_apply_buffer_allocated`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`selection_criterion`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`omp_set_lock`
       * :c:func:`omp_unset_lock`
       * :c:func:`resize_h_apply_buffer`

 
.. c:function:: perturb_buffer_by_mono_dummy:


    File : :file:`perturbation/perturbation.irp.f_shell_13`

    .. code:: fortran

        subroutine perturb_buffer_by_mono_dummy(i_generator,buffer,buffer_size,e_2_pert_buffer,coef_pert_buffer,sum_e_2_pert,sum_norm_pert,sum_H_pert_diag,N_st,Nint,key_mask,fock_diag_tmp,electronic_energy)


    Apply pertubration ``dummy`` to the buffer of determinants generated in the H_apply
    routine.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_det_generators`
       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`psi_det_generators`
       * :c:data:`psi_selectors`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`create_minilist`
       * :c:func:`create_minilist_find_previous`
       * :c:func:`pt2_dummy`

 
.. c:function:: perturb_buffer_by_mono_epstein_nesbet:


    File : :file:`perturbation/perturbation.irp.f_shell_13`

    .. code:: fortran

        subroutine perturb_buffer_by_mono_epstein_nesbet(i_generator,buffer,buffer_size,e_2_pert_buffer,coef_pert_buffer,sum_e_2_pert,sum_norm_pert,sum_H_pert_diag,N_st,Nint,key_mask,fock_diag_tmp,electronic_energy)


    Apply pertubration ``epstein_nesbet`` to the buffer of determinants generated in the H_apply
    routine.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_det_generators`
       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`psi_det_generators`
       * :c:data:`psi_selectors`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`create_minilist`
       * :c:func:`create_minilist_find_previous`
       * :c:func:`pt2_epstein_nesbet`

 
.. c:function:: perturb_buffer_by_mono_epstein_nesbet_2x2:


    File : :file:`perturbation/perturbation.irp.f_shell_13`

    .. code:: fortran

        subroutine perturb_buffer_by_mono_epstein_nesbet_2x2(i_generator,buffer,buffer_size,e_2_pert_buffer,coef_pert_buffer,sum_e_2_pert,sum_norm_pert,sum_H_pert_diag,N_st,Nint,key_mask,fock_diag_tmp,electronic_energy)


    Apply pertubration ``epstein_nesbet_2x2`` to the buffer of determinants generated in the H_apply
    routine.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_det_generators`
       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`psi_det_generators`
       * :c:data:`psi_selectors`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`create_minilist`
       * :c:func:`create_minilist_find_previous`
       * :c:func:`pt2_epstein_nesbet_2x2`

 
.. c:function:: perturb_buffer_by_mono_epstein_nesbet_2x2_no_ci_diag:


    File : :file:`perturbation/perturbation.irp.f_shell_13`

    .. code:: fortran

        subroutine perturb_buffer_by_mono_epstein_nesbet_2x2_no_ci_diag(i_generator,buffer,buffer_size,e_2_pert_buffer,coef_pert_buffer,sum_e_2_pert,sum_norm_pert,sum_H_pert_diag,N_st,Nint,key_mask,fock_diag_tmp,electronic_energy)


    Apply pertubration ``epstein_nesbet_2x2_no_ci_diag`` to the buffer of determinants generated in the H_apply
    routine.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_det_generators`
       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`psi_det_generators`
       * :c:data:`psi_selectors`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`create_minilist`
       * :c:func:`create_minilist_find_previous`
       * :c:func:`pt2_epstein_nesbet_2x2_no_ci_diag`

 
.. c:function:: perturb_buffer_by_mono_moller_plesset:


    File : :file:`perturbation/perturbation.irp.f_shell_13`

    .. code:: fortran

        subroutine perturb_buffer_by_mono_moller_plesset(i_generator,buffer,buffer_size,e_2_pert_buffer,coef_pert_buffer,sum_e_2_pert,sum_norm_pert,sum_H_pert_diag,N_st,Nint,key_mask,fock_diag_tmp,electronic_energy)


    Apply pertubration ``moller_plesset`` to the buffer of determinants generated in the H_apply
    routine.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_det_generators`
       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`psi_det_generators`
       * :c:data:`psi_selectors`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`create_minilist`
       * :c:func:`create_minilist_find_previous`
       * :c:func:`pt2_moller_plesset`

 
.. c:function:: perturb_buffer_by_mono_qdpt:


    File : :file:`perturbation/perturbation.irp.f_shell_13`

    .. code:: fortran

        subroutine perturb_buffer_by_mono_qdpt(i_generator,buffer,buffer_size,e_2_pert_buffer,coef_pert_buffer,sum_e_2_pert,sum_norm_pert,sum_H_pert_diag,N_st,Nint,key_mask,fock_diag_tmp,electronic_energy)


    Apply pertubration ``qdpt`` to the buffer of determinants generated in the H_apply
    routine.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_det`
       * :c:data:`n_det_generators`
       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`psi_det_generators`
       * :c:data:`psi_selectors`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`create_minilist`
       * :c:func:`create_minilist_find_previous`
       * :c:func:`pt2_qdpt`

 
.. c:function:: perturb_buffer_dummy:


    File : :file:`perturbation/perturbation.irp.f_shell_13`

    .. code:: fortran

        subroutine perturb_buffer_dummy(i_generator,buffer,buffer_size,e_2_pert_buffer,coef_pert_buffer,sum_e_2_pert,sum_norm_pert,sum_H_pert_diag,N_st,Nint,key_mask,fock_diag_tmp,electronic_energy)


    Apply pertubration ``dummy`` to the buffer of determinants generated in the H_apply
    routine.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_det_generators`
       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`psi_det_generators`
       * :c:data:`psi_selectors`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`create_microlist`
       * :c:func:`create_minilist`
       * :c:func:`create_minilist_find_previous`
       * :c:func:`getmobiles`
       * :c:func:`pt2_dummy`

 
.. c:function:: perturb_buffer_epstein_nesbet:


    File : :file:`perturbation/perturbation.irp.f_shell_13`

    .. code:: fortran

        subroutine perturb_buffer_epstein_nesbet(i_generator,buffer,buffer_size,e_2_pert_buffer,coef_pert_buffer,sum_e_2_pert,sum_norm_pert,sum_H_pert_diag,N_st,Nint,key_mask,fock_diag_tmp,electronic_energy)


    Apply pertubration ``epstein_nesbet`` to the buffer of determinants generated in the H_apply
    routine.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_det_generators`
       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`psi_det_generators`
       * :c:data:`psi_selectors`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`create_microlist`
       * :c:func:`create_minilist`
       * :c:func:`create_minilist_find_previous`
       * :c:func:`getmobiles`
       * :c:func:`pt2_epstein_nesbet`

 
.. c:function:: perturb_buffer_epstein_nesbet_2x2:


    File : :file:`perturbation/perturbation.irp.f_shell_13`

    .. code:: fortran

        subroutine perturb_buffer_epstein_nesbet_2x2(i_generator,buffer,buffer_size,e_2_pert_buffer,coef_pert_buffer,sum_e_2_pert,sum_norm_pert,sum_H_pert_diag,N_st,Nint,key_mask,fock_diag_tmp,electronic_energy)


    Apply pertubration ``epstein_nesbet_2x2`` to the buffer of determinants generated in the H_apply
    routine.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_det_generators`
       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`psi_det_generators`
       * :c:data:`psi_selectors`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`create_microlist`
       * :c:func:`create_minilist`
       * :c:func:`create_minilist_find_previous`
       * :c:func:`getmobiles`
       * :c:func:`pt2_epstein_nesbet_2x2`

 
.. c:function:: perturb_buffer_epstein_nesbet_2x2_no_ci_diag:


    File : :file:`perturbation/perturbation.irp.f_shell_13`

    .. code:: fortran

        subroutine perturb_buffer_epstein_nesbet_2x2_no_ci_diag(i_generator,buffer,buffer_size,e_2_pert_buffer,coef_pert_buffer,sum_e_2_pert,sum_norm_pert,sum_H_pert_diag,N_st,Nint,key_mask,fock_diag_tmp,electronic_energy)


    Apply pertubration ``epstein_nesbet_2x2_no_ci_diag`` to the buffer of determinants generated in the H_apply
    routine.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_det_generators`
       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`psi_det_generators`
       * :c:data:`psi_selectors`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`create_microlist`
       * :c:func:`create_minilist`
       * :c:func:`create_minilist_find_previous`
       * :c:func:`getmobiles`
       * :c:func:`pt2_epstein_nesbet_2x2_no_ci_diag`

 
.. c:function:: perturb_buffer_moller_plesset:


    File : :file:`perturbation/perturbation.irp.f_shell_13`

    .. code:: fortran

        subroutine perturb_buffer_moller_plesset(i_generator,buffer,buffer_size,e_2_pert_buffer,coef_pert_buffer,sum_e_2_pert,sum_norm_pert,sum_H_pert_diag,N_st,Nint,key_mask,fock_diag_tmp,electronic_energy)


    Apply pertubration ``moller_plesset`` to the buffer of determinants generated in the H_apply
    routine.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_det_generators`
       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`psi_det_generators`
       * :c:data:`psi_selectors`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`create_microlist`
       * :c:func:`create_minilist`
       * :c:func:`create_minilist_find_previous`
       * :c:func:`getmobiles`
       * :c:func:`pt2_moller_plesset`

 
.. c:function:: perturb_buffer_qdpt:


    File : :file:`perturbation/perturbation.irp.f_shell_13`

    .. code:: fortran

        subroutine perturb_buffer_qdpt(i_generator,buffer,buffer_size,e_2_pert_buffer,coef_pert_buffer,sum_e_2_pert,sum_norm_pert,sum_H_pert_diag,N_st,Nint,key_mask,fock_diag_tmp,electronic_energy)


    Apply pertubration ``qdpt`` to the buffer of determinants generated in the H_apply
    routine.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_det_generators`
       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`psi_det_generators`
       * :c:data:`psi_selectors`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`create_microlist`
       * :c:func:`create_minilist`
       * :c:func:`create_minilist_find_previous`
       * :c:func:`getmobiles`
       * :c:func:`pt2_qdpt`

 
.. c:function:: pt2_dummy:


    File : :file:`perturbation/pt2_equations.irp.f_template_305`

    .. code:: fortran

        subroutine pt2_dummy (electronic_energy,det_ref,det_pert,fock_diag_tmp,c_pert,e_2_pert,H_pert_diag,Nint,ndet,N_st,minilist,idx_minilist,N_minilist)


    Dummy perturbation to add all connected determinants.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_det_selectors`
       * :c:data:`psi_selectors`
       * :c:data:`psi_selectors_size`
       * :c:data:`selection_criterion`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`perturb_buffer_by_mono_dummy`
       * :c:func:`perturb_buffer_dummy`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_h_psi_minilist`

 
.. c:function:: pt2_epstein_nesbet:


    File : :file:`perturbation/pt2_equations.irp.f_template_305`

    .. code:: fortran

        subroutine pt2_epstein_nesbet (electronic_energy,det_ref,det_pert,fock_diag_tmp,c_pert,e_2_pert,H_pert_diag,Nint,ndet,N_st,minilist,idx_minilist,N_minilist)


    Compute the standard Epstein-Nesbet perturbative first order coefficient and
    second order energetic contribution for the various N_st states.
    
    `c_pert(i)` = $\frac{\langle i|H|\alpha \rangle}{ E_n - \langle \alpha|H|\alpha \rangle }$.
    
    `e_2_pert(i)` = $\frac{\langle i|H|\alpha \rangle^2}{ E_n - \langle \alpha|H|\alpha \rangle }$.
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`psi_selectors`
       * :c:data:`psi_selectors_size`
       * :c:data:`selection_criterion`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`perturb_buffer_by_mono_epstein_nesbet`
       * :c:func:`perturb_buffer_epstein_nesbet`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_h_psi_minilist`

 
.. c:function:: pt2_epstein_nesbet_2x2:


    File : :file:`perturbation/pt2_equations.irp.f_template_305`

    .. code:: fortran

        subroutine pt2_epstein_nesbet_2x2 (electronic_energy,det_ref,det_pert,fock_diag_tmp,c_pert,e_2_pert,H_pert_diag,Nint,ndet,N_st,minilist,idx_minilist,N_minilist)


    Computes the Epstein-Nesbet 2x2 diagonalization coefficient and energetic contribution
    for the various N_st states.
    
    `e_2_pert(i)` = $\frac{1}{2} ( \langle \alpha|H|\alpha \rangle -  E_n) - \sqrt{ (\langle \alpha|H|\alpha \rangle -  E_n)^2 + 4 \langle i|H|\alpha \rangle^2 }$.
    
    `c_pert(i)` = `e_2_pert(i)` $\times \frac{1}{ \langle i|H|\alpha \rangle}$.
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`psi_selectors`
       * :c:data:`psi_selectors_size`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`perturb_buffer_by_mono_epstein_nesbet_2x2`
       * :c:func:`perturb_buffer_epstein_nesbet_2x2`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_h_psi`

 
.. c:function:: pt2_epstein_nesbet_2x2_no_ci_diag:


    File : :file:`perturbation/pt2_equations.irp.f_template_305`

    .. code:: fortran

        subroutine pt2_epstein_nesbet_2x2_no_ci_diag(electronic_energy,det_ref,det_pert,fock_diag_tmp,c_pert,e_2_pert,H_pert_diag,Nint,ndet,N_st,minilist,idx_minilist,N_minilist)


    compute the Epstein-Nesbet 2x2 diagonalization coefficient and energetic contribution
    
    for the various N_st states.
    
    e_2_pert(i) = 0.5 * (( <det_pert|H|det_pert> -  E(i) )  - sqrt( ( <det_pert|H|det_pert> -  E(i)) ^2 + 4 <psi(i)|H|det_pert>^2  )
    
    c_pert(i) = e_2_pert(i)/ <psi(i)|H|det_pert>
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`psi_energy`
       * :c:data:`psi_selectors`
       * :c:data:`psi_selectors_size`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`perturb_buffer_by_mono_epstein_nesbet_2x2_no_ci_diag`
       * :c:func:`perturb_buffer_epstein_nesbet_2x2_no_ci_diag`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`i_h_psi`

 
.. c:function:: pt2_moller_plesset:


    File : :file:`perturbation/pt2_equations.irp.f_template_305`

    .. code:: fortran

        subroutine pt2_moller_plesset (electronic_energy,det_ref,det_pert,fock_diag_tmp,c_pert,e_2_pert,H_pert_diag,Nint,ndet,N_st,minilist,idx_minilist,N_minilist)


    Computes the standard Moller-Plesset perturbative first order coefficient and second
    order energetic contribution for the various N_st states.
    
    `c_pert(i)` = $\frac{\langle i|H|\alpha \rangle}{\text{difference of orbital energies}}$.
    
    `e_2_pert(i)` = $\frac{\langle i|H|\alpha \rangle^2}{\text{difference of orbital energies}}$.
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_mo`
       * :c:data:`mo_num`
       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`psi_selectors`
       * :c:data:`psi_selectors_size`
       * :c:data:`ref_bitmask`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`perturb_buffer_by_mono_moller_plesset`
       * :c:func:`perturb_buffer_moller_plesset`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`decode_exc`
       * :c:func:`get_excitation`
       * :c:func:`i_h_psi_minilist`

 
.. c:function:: pt2_qdpt:


    File : :file:`perturbation/pt2_equations.irp.f_template_305`

    .. code:: fortran

        subroutine pt2_qdpt (electronic_energy,det_ref,det_pert,fock_diag_tmp,c_pert,e_2_pert,H_pert_diag,Nint,ndet,N_st,minilist,idx_minilist,N_minilist)


    Computes the QDPT first order coefficient and second order energetic contribution
    for the various N_st states.
    
    `c_pert(i)` = $\frac{\langle i|H|\alpha \rangle}{\langle i|H|i \rangle - \langle \alpha|H|\alpha \rangle}$.
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_num`
       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`psi_selectors`
       * :c:data:`psi_selectors_size`
       * :c:data:`selection_criterion`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`perturb_buffer_by_mono_qdpt`
       * :c:func:`perturb_buffer_qdpt`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`get_excitation_degree`
       * :c:func:`i_h_j`
       * :c:func:`i_h_psi_minilist`

 
.. c:function:: remove_small_contributions:


    File : :file:`perturbation/selection.irp.f`

    .. code:: fortran

        subroutine remove_small_contributions


    Remove determinants with small contributions. N_states is assumed to be
    provided.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det`
       * :c:data:`n_det_generators`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_coef`
       * :c:data:`psi_det_sorted`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`
       * :c:data:`psi_det_sorted`
       * :c:data:`selection_criterion`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`diagonalize_ci`
       * :c:func:`i_h_psi`
       * :c:func:`write_int`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`ci_electronic_energy`
       * :c:data:`ci_electronic_energy`
       * :c:data:`ci_energy`
       * :c:data:`ci_electronic_energy`
       * :c:data:`n_det`
       * :c:data:`psi_coef`
       * :c:data:`psi_det`
       * :c:data:`psi_energy`
       * :c:data:`psi_energy`

