.. _module_selectors_utils: 
 
.. program:: selectors_utils 
 
.. default-role:: option 
 
===============
selectors_utils
===============

Quantities derived from the selectors.

From the selectors of the wave function, this module computes:

* the excitation degree of each selector with respect to the |HF|
  reference, :c:data:`exc_degree_per_selectors`, and the list
  :c:data:`double_index_selectors` and the number
  :c:data:`n_double_selectors` of the double excitations among them
* the correlation energy per selector,
  :c:data:`e_corr_per_selectors`: for a double excitation
  :math:`D_i`, it is
  :math:`\langle D_i | H | \text{HF} \rangle c(D_i)/c(\text{HF})`
* the second order estimate of the correlation energy,
  :c:data:`e_corr_second_order`, and its variant restricted to the
  double excitations, :c:data:`e_corr_double_only`
* the diagonal elements of the Hamiltonian matrix for each selector,
  :c:data:`psi_selectors_diag_h_mat`
 
 
 
Providers 
--------- 
 
.. c:var:: double coef_hf_selector


    File : :file:`selectors_utils/e_corr_selectors.irp.f`

    .. code:: fortran

        double precision	:: coef_hf_selector	
        double precision	:: inv_selectors_coef_hf	
        double precision	:: inv_selectors_coef_hf_squared	
        double precision, allocatable	:: e_corr_per_selectors	(N_det_selectors)
        double precision, allocatable	:: i_h_hf_per_selectors	(N_det_selectors)
        double precision, allocatable	:: delta_e_per_selector	(N_det_selectors)
        double precision	:: e_corr_double_only	
        double precision	:: e_corr_second_order	


    Correlation energy per determinant with respect to the Hartree-Fock determinant
    for the all the double excitations in the selectors determinants.
    
    E_corr_per_selectors(i) = :math:`\langle D_i | H | \text{HF}\rangle  c(D_i)/c(HF)` if :math:`| D_i \rangle` is a double excitation.
    
    E_corr_per_selectors(i) = -1000.d0 if it is not a double excitation
    
    coef_hf_selector = coefficient of the Hartree Fock determinant in the selectors determinants

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`exc_degree_per_selectors`
       * :c:data:`mo_integrals_map`
       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`psi_selectors`
       * :c:data:`ref_bitmask`
       * :c:data:`ref_bitmask_energy`


 
.. c:var:: double delta_e_per_selector


    File : :file:`selectors_utils/e_corr_selectors.irp.f`

    .. code:: fortran

        double precision	:: coef_hf_selector	
        double precision	:: inv_selectors_coef_hf	
        double precision	:: inv_selectors_coef_hf_squared	
        double precision, allocatable	:: e_corr_per_selectors	(N_det_selectors)
        double precision, allocatable	:: i_h_hf_per_selectors	(N_det_selectors)
        double precision, allocatable	:: delta_e_per_selector	(N_det_selectors)
        double precision	:: e_corr_double_only	
        double precision	:: e_corr_second_order	


    Correlation energy per determinant with respect to the Hartree-Fock determinant
    for the all the double excitations in the selectors determinants.
    
    E_corr_per_selectors(i) = :math:`\langle D_i | H | \text{HF}\rangle  c(D_i)/c(HF)` if :math:`| D_i \rangle` is a double excitation.
    
    E_corr_per_selectors(i) = -1000.d0 if it is not a double excitation
    
    coef_hf_selector = coefficient of the Hartree Fock determinant in the selectors determinants

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`exc_degree_per_selectors`
       * :c:data:`mo_integrals_map`
       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`psi_selectors`
       * :c:data:`ref_bitmask`
       * :c:data:`ref_bitmask_energy`


 
.. c:var:: int double_index_selectors


    File : :file:`selectors_utils/e_corr_selectors.irp.f`

    .. code:: fortran

        integer, allocatable	:: exc_degree_per_selectors	(N_det_selectors)
        integer, allocatable	:: double_index_selectors	(N_det_selectors)
        integer	:: n_double_selectors	


    Degree of excitation respect to Hartree Fock for the wave function
    for the all the selectors determinants.
    
    double_index_selectors = list of the index of the double excitations
    
    n_double_selectors = number of double excitations in the selectors determinants

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`psi_selectors`
       * :c:data:`ref_bitmask`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`coef_hf_selector`

 
.. c:var:: double e_corr_double_only


    File : :file:`selectors_utils/e_corr_selectors.irp.f`

    .. code:: fortran

        double precision	:: coef_hf_selector	
        double precision	:: inv_selectors_coef_hf	
        double precision	:: inv_selectors_coef_hf_squared	
        double precision, allocatable	:: e_corr_per_selectors	(N_det_selectors)
        double precision, allocatable	:: i_h_hf_per_selectors	(N_det_selectors)
        double precision, allocatable	:: delta_e_per_selector	(N_det_selectors)
        double precision	:: e_corr_double_only	
        double precision	:: e_corr_second_order	


    Correlation energy per determinant with respect to the Hartree-Fock determinant
    for the all the double excitations in the selectors determinants.
    
    E_corr_per_selectors(i) = :math:`\langle D_i | H | \text{HF}\rangle  c(D_i)/c(HF)` if :math:`| D_i \rangle` is a double excitation.
    
    E_corr_per_selectors(i) = -1000.d0 if it is not a double excitation
    
    coef_hf_selector = coefficient of the Hartree Fock determinant in the selectors determinants

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`exc_degree_per_selectors`
       * :c:data:`mo_integrals_map`
       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`psi_selectors`
       * :c:data:`ref_bitmask`
       * :c:data:`ref_bitmask_energy`


 
.. c:var:: double e_corr_per_selectors


    File : :file:`selectors_utils/e_corr_selectors.irp.f`

    .. code:: fortran

        double precision	:: coef_hf_selector	
        double precision	:: inv_selectors_coef_hf	
        double precision	:: inv_selectors_coef_hf_squared	
        double precision, allocatable	:: e_corr_per_selectors	(N_det_selectors)
        double precision, allocatable	:: i_h_hf_per_selectors	(N_det_selectors)
        double precision, allocatable	:: delta_e_per_selector	(N_det_selectors)
        double precision	:: e_corr_double_only	
        double precision	:: e_corr_second_order	


    Correlation energy per determinant with respect to the Hartree-Fock determinant
    for the all the double excitations in the selectors determinants.
    
    E_corr_per_selectors(i) = :math:`\langle D_i | H | \text{HF}\rangle  c(D_i)/c(HF)` if :math:`| D_i \rangle` is a double excitation.
    
    E_corr_per_selectors(i) = -1000.d0 if it is not a double excitation
    
    coef_hf_selector = coefficient of the Hartree Fock determinant in the selectors determinants

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`exc_degree_per_selectors`
       * :c:data:`mo_integrals_map`
       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`psi_selectors`
       * :c:data:`ref_bitmask`
       * :c:data:`ref_bitmask_energy`


 
.. c:var:: double e_corr_second_order


    File : :file:`selectors_utils/e_corr_selectors.irp.f`

    .. code:: fortran

        double precision	:: coef_hf_selector	
        double precision	:: inv_selectors_coef_hf	
        double precision	:: inv_selectors_coef_hf_squared	
        double precision, allocatable	:: e_corr_per_selectors	(N_det_selectors)
        double precision, allocatable	:: i_h_hf_per_selectors	(N_det_selectors)
        double precision, allocatable	:: delta_e_per_selector	(N_det_selectors)
        double precision	:: e_corr_double_only	
        double precision	:: e_corr_second_order	


    Correlation energy per determinant with respect to the Hartree-Fock determinant
    for the all the double excitations in the selectors determinants.
    
    E_corr_per_selectors(i) = :math:`\langle D_i | H | \text{HF}\rangle  c(D_i)/c(HF)` if :math:`| D_i \rangle` is a double excitation.
    
    E_corr_per_selectors(i) = -1000.d0 if it is not a double excitation
    
    coef_hf_selector = coefficient of the Hartree Fock determinant in the selectors determinants

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`exc_degree_per_selectors`
       * :c:data:`mo_integrals_map`
       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`psi_selectors`
       * :c:data:`ref_bitmask`
       * :c:data:`ref_bitmask_energy`


 
.. c:var:: int exc_degree_per_selectors


    File : :file:`selectors_utils/e_corr_selectors.irp.f`

    .. code:: fortran

        integer, allocatable	:: exc_degree_per_selectors	(N_det_selectors)
        integer, allocatable	:: double_index_selectors	(N_det_selectors)
        integer	:: n_double_selectors	


    Degree of excitation respect to Hartree Fock for the wave function
    for the all the selectors determinants.
    
    double_index_selectors = list of the index of the double excitations
    
    n_double_selectors = number of double excitations in the selectors determinants

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`psi_selectors`
       * :c:data:`ref_bitmask`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`coef_hf_selector`

 
.. c:var:: double i_h_hf_per_selectors


    File : :file:`selectors_utils/e_corr_selectors.irp.f`

    .. code:: fortran

        double precision	:: coef_hf_selector	
        double precision	:: inv_selectors_coef_hf	
        double precision	:: inv_selectors_coef_hf_squared	
        double precision, allocatable	:: e_corr_per_selectors	(N_det_selectors)
        double precision, allocatable	:: i_h_hf_per_selectors	(N_det_selectors)
        double precision, allocatable	:: delta_e_per_selector	(N_det_selectors)
        double precision	:: e_corr_double_only	
        double precision	:: e_corr_second_order	


    Correlation energy per determinant with respect to the Hartree-Fock determinant
    for the all the double excitations in the selectors determinants.
    
    E_corr_per_selectors(i) = :math:`\langle D_i | H | \text{HF}\rangle  c(D_i)/c(HF)` if :math:`| D_i \rangle` is a double excitation.
    
    E_corr_per_selectors(i) = -1000.d0 if it is not a double excitation
    
    coef_hf_selector = coefficient of the Hartree Fock determinant in the selectors determinants

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`exc_degree_per_selectors`
       * :c:data:`mo_integrals_map`
       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`psi_selectors`
       * :c:data:`ref_bitmask`
       * :c:data:`ref_bitmask_energy`


 
.. c:var:: double inv_selectors_coef_hf


    File : :file:`selectors_utils/e_corr_selectors.irp.f`

    .. code:: fortran

        double precision	:: coef_hf_selector	
        double precision	:: inv_selectors_coef_hf	
        double precision	:: inv_selectors_coef_hf_squared	
        double precision, allocatable	:: e_corr_per_selectors	(N_det_selectors)
        double precision, allocatable	:: i_h_hf_per_selectors	(N_det_selectors)
        double precision, allocatable	:: delta_e_per_selector	(N_det_selectors)
        double precision	:: e_corr_double_only	
        double precision	:: e_corr_second_order	


    Correlation energy per determinant with respect to the Hartree-Fock determinant
    for the all the double excitations in the selectors determinants.
    
    E_corr_per_selectors(i) = :math:`\langle D_i | H | \text{HF}\rangle  c(D_i)/c(HF)` if :math:`| D_i \rangle` is a double excitation.
    
    E_corr_per_selectors(i) = -1000.d0 if it is not a double excitation
    
    coef_hf_selector = coefficient of the Hartree Fock determinant in the selectors determinants

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`exc_degree_per_selectors`
       * :c:data:`mo_integrals_map`
       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`psi_selectors`
       * :c:data:`ref_bitmask`
       * :c:data:`ref_bitmask_energy`


 
.. c:var:: double inv_selectors_coef_hf_squared


    File : :file:`selectors_utils/e_corr_selectors.irp.f`

    .. code:: fortran

        double precision	:: coef_hf_selector	
        double precision	:: inv_selectors_coef_hf	
        double precision	:: inv_selectors_coef_hf_squared	
        double precision, allocatable	:: e_corr_per_selectors	(N_det_selectors)
        double precision, allocatable	:: i_h_hf_per_selectors	(N_det_selectors)
        double precision, allocatable	:: delta_e_per_selector	(N_det_selectors)
        double precision	:: e_corr_double_only	
        double precision	:: e_corr_second_order	


    Correlation energy per determinant with respect to the Hartree-Fock determinant
    for the all the double excitations in the selectors determinants.
    
    E_corr_per_selectors(i) = :math:`\langle D_i | H | \text{HF}\rangle  c(D_i)/c(HF)` if :math:`| D_i \rangle` is a double excitation.
    
    E_corr_per_selectors(i) = -1000.d0 if it is not a double excitation
    
    coef_hf_selector = coefficient of the Hartree Fock determinant in the selectors determinants

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_mo_integrals`
       * :c:data:`big_array_coulomb_integrals`
       * :c:data:`exc_degree_per_selectors`
       * :c:data:`mo_integrals_map`
       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`psi_selectors`
       * :c:data:`ref_bitmask`
       * :c:data:`ref_bitmask_energy`


 
.. c:var:: int n_double_selectors


    File : :file:`selectors_utils/e_corr_selectors.irp.f`

    .. code:: fortran

        integer, allocatable	:: exc_degree_per_selectors	(N_det_selectors)
        integer, allocatable	:: double_index_selectors	(N_det_selectors)
        integer	:: n_double_selectors	


    Degree of excitation respect to Hartree Fock for the wave function
    for the all the selectors determinants.
    
    double_index_selectors = list of the index of the double excitations
    
    n_double_selectors = number of double excitations in the selectors determinants

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`psi_selectors`
       * :c:data:`ref_bitmask`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`coef_hf_selector`

 
.. c:var:: double psi_selectors_coef_transp


    File : :file:`selectors_utils/selectors.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_selectors_coef_transp	(N_states,psi_selectors_size)


    Transposed psi_selectors

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_selectors`
       * :c:data:`n_states`
       * :c:data:`psi_selectors`
       * :c:data:`psi_selectors_size`


 
.. c:var:: double psi_selectors_diag_h_mat


    File : :file:`selectors_utils/selectors.irp.f`

    .. code:: fortran

        double precision, allocatable	:: psi_selectors_diag_h_mat	(psi_selectors_size)


    Diagonal elements of the H matrix for each selectors

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_num`
       * :c:data:`n_det_selectors`
       * :c:data:`n_int`
       * :c:data:`psi_selectors`
       * :c:data:`psi_selectors_size`
       * :c:data:`ref_bitmask`
       * :c:data:`ref_bitmask_energy`


 
.. c:var:: int psi_selectors_size


    File : :file:`selectors_utils/selectors.irp.f`

    .. code:: fortran

        integer	:: psi_selectors_size	



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`psi_det_size`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`psi_selectors`
       * :c:data:`psi_selectors_coef_transp`
       * :c:data:`psi_selectors_diag_h_mat`

 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: void zmq_get_n_det_generators(void)


    File : :file:`selectors_utils/zmq.irp.f_template_102`

    .. code:: fortran

        integer function zmq_get_N_det_generators(zmq_to_qp_run_socket, worker_id)


    Get N_det_generators from the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`n_det_generators`
       * :c:data:`zmq_state`

 
.. c:function:: void zmq_get_n_det_selectors(void)


    File : :file:`selectors_utils/zmq.irp.f_template_102`

    .. code:: fortran

        integer function zmq_get_N_det_selectors(zmq_to_qp_run_socket, worker_id)


    Get N_det_selectors from the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mpi_master`
       * :c:data:`n_det_selectors`
       * :c:data:`zmq_state`

 
.. c:function:: void zmq_put_n_det_generators(void)


    File : :file:`selectors_utils/zmq.irp.f_template_102`

    .. code:: fortran

        integer function zmq_put_N_det_generators(zmq_to_qp_run_socket,worker_id)


    Put N_det_generators on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_generators`
       * :c:data:`zmq_state`

 
.. c:function:: void zmq_put_n_det_selectors(void)


    File : :file:`selectors_utils/zmq.irp.f_template_102`

    .. code:: fortran

        integer function zmq_put_N_det_selectors(zmq_to_qp_run_socket,worker_id)


    Put N_det_selectors on the qp_run scheduler

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_selectors`
       * :c:data:`zmq_state`

