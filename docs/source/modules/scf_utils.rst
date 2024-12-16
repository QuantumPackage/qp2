.. _module_scf_utils: 
 
.. program:: scf_utils 
 
.. default-role:: option 
 
=========
scf_utils
=========



The scf_utils module is an abstract module which contains the basics to perform *Restricted* SCF calculations (the
spatial part of the |MOs| is common for alpha and beta spinorbitals) based on a single-determinant wave function.

This module does not produce any executable *and must not do*, but instead it contains everything one needs to perform an orbital optimization based on an Fock matrix.
The ``scf_utils`` module is meant to be included in the :file:`NEED` of the various single determinant SCF procedures, such as ``hartree_fock`` or ``kohn_sham``, where a specific definition of the Fock matrix is given (see :file:`hartree_fock fock_matrix_hf.irp.f` for an example).

All SCF programs perform the following actions:


#. Compute/Read all the one- and two-electron integrals, and store them in memory

#. Check in the |EZFIO| database if there is a set of |MOs|. If there is, it
   will read them as initial guess. Otherwise, it will create a guess.
#. Perform the |SCF| iterations based on the definition of the Fock matrix


The main keywords/options are:

* :option:`scf_utils thresh_scf`
* :option:`scf_utils level_shift`

At each iteration, the |MOs| are saved in the |EZFIO| database. Hence, if the calculation
crashes for any unexpected reason, the calculation can be restarted by running again
the |SCF| with the same |EZFIO| database.

The `DIIS`_ algorithm is implemented, as well as the `level-shifting`_ method.
If the |SCF| does not converge, try again with a higher value of :option:`level_shift`.

To start a calculation from scratch, the simplest way is to remove the
``mo_basis`` directory from the |EZFIO| database, and run the |SCF| again.

.. _DIIS: https://en.wikipedia.org/w/index.php?title=DIIS
.. _level-shifting: https://doi.org/10.1002/qua.560070407

 
 
 
EZFIO parameters 
---------------- 
 
.. option:: max_dim_diis
 
    Maximum size of the DIIS extrapolation procedure
 
    Default: 15
 
.. option:: threshold_diis
 
    Threshold on the convergence of the DIIS error vector during a Hartree-Fock calculation. If 0. is chosen, the square root of thresh_scf will be used.
 
    Default: 0.
 
.. option:: thresh_scf
 
    Threshold on the convergence of the Hartree Fock energy.
 
    Default: 1.e-10
 
.. option:: n_it_scf_max
 
    Maximum number of SCF iterations
 
    Default: 500
 
.. option:: level_shift
 
    Energy shift on the virtual MOs to improve SCF convergence
 
    Default: 0.
 
.. option:: scf_algorithm
 
    Type of SCF algorithm used. Possible choices are [ Simple | DIIS]
 
    Default: DIIS
 
.. option:: mo_guess_type
 
    Initial MO guess. Can be [ Huckel | HCore ]
 
    Default: Huckel
 
.. option:: energy
 
    Calculated HF energy
 
 
.. option:: do_mom
 
    If true, this will run a MOM calculation. The overlap will be computed at each step with respect to the initial MOs. After an initial Hartree-Fock calculation, the guess can be created by swapping molecular orbitals through the qp run swap_mos command.
 
    Default: False
 
.. option:: frozen_orb_scf
 
    If true, leave untouched all the orbitals defined as core and optimize all the orbitals defined as active with qp_set_mo_class
 
    Default: False
 
.. option:: no_oa_or_av_opt
 
    If true, you set to zero all Fock elements between the orbital set to active and all the other orbitals
 
    Default: False
 
 
Providers 
--------- 
 
.. c:var:: all_shells_closed


    File : :file:`scf_utils/scf_density_matrix_ao.irp.f`

    .. code:: fortran

        logical	:: all_shells_closed	


    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_ao`
       * :c:data:`fock_matrix_mo`
       * :c:data:`scf_density_matrix_ao`

 
.. c:var:: eigenvalues_fock_matrix_ao


    File : :file:`scf_utils/diis.irp.f`

    .. code:: fortran

        double precision, allocatable	:: eigenvalues_fock_matrix_ao	(AO_num)
        double precision, allocatable	:: eigenvectors_fock_matrix_ao	(AO_num,AO_num)


    Eigenvalues and eigenvectors of the Fock matrix over the AO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`fock_matrix_ao`
       * :c:data:`s_half_inv`


 
.. c:var:: eigenvectors_fock_matrix_ao


    File : :file:`scf_utils/diis.irp.f`

    .. code:: fortran

        double precision, allocatable	:: eigenvalues_fock_matrix_ao	(AO_num)
        double precision, allocatable	:: eigenvectors_fock_matrix_ao	(AO_num,AO_num)


    Eigenvalues and eigenvectors of the Fock matrix over the AO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`fock_matrix_ao`
       * :c:data:`s_half_inv`


 
.. c:var:: eigenvectors_fock_matrix_mo


    File : :file:`scf_utils/diagonalize_fock.irp.f`

    .. code:: fortran

        double precision, allocatable	:: eigenvectors_fock_matrix_mo	(ao_num,mo_num)


    Eigenvectors of the Fock matrix in the |MO| basis obtained with level shift.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`fock_matrix_mo`
       * :c:data:`frozen_orb_scf`
       * :c:data:`level_shift`
       * :c:data:`list_act`
       * :c:data:`list_core`
       * :c:data:`list_inact`
       * :c:data:`list_virt`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_orb`
       * :c:data:`n_inact_orb`
       * :c:data:`n_virt_orb`
       * :c:data:`no_oa_or_av_opt`


 
.. c:function:: extrapolate_fock_matrix:


    File : :file:`scf_utils/roothaan_hall_scf.irp.f`

    .. code:: fortran

        subroutine extrapolate_Fock_matrix(               &
           error_matrix_DIIS,Fock_matrix_DIIS,    &
           Fock_matrix_AO_,size_Fock_matrix_AO,   &
           iteration_SCF,dim_DIIS                 &
           )


    Compute the extrapolated Fock matrix using the DIIS procedure

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`max_dim_diis`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`roothaan_hall_scf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgecon`
       * :c:func:`dgesv`
       * :c:func:`dgetrf`

 
.. c:var:: fock_matrix_ao


    File : :file:`scf_utils/fock_matrix.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_matrix_ao	(ao_num,ao_num)


    Fock matrix in AO basis set

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_shells_closed`
       * :c:data:`ao_num`
       * :c:data:`fock_matrix_ao_alpha`
       * :c:data:`fock_matrix_mo`
       * :c:data:`frozen_orb_scf`
       * :c:data:`level_shift`
       * :c:data:`mo_num`
       * :c:data:`s_mo_coef`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`eigenvalues_fock_matrix_ao`
       * :c:data:`fps_spf_matrix_ao`

 
.. c:var:: fock_matrix_diag_mo


    File : :file:`scf_utils/fock_matrix.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_matrix_mo	(mo_num,mo_num)
        double precision, allocatable	:: fock_matrix_diag_mo	(mo_num)


    Fock matrix on the MO basis.
    For open shells, the ROHF Fock Matrix is ::
    
          |  Rcc  |  F^b  |  Fcv  |
          |-----------------------|
          |  F^b  |  Roo  |  F^a  |
          |-----------------------|
          |  Fcv  |  F^a  |  Rvv  |
    
    C: Core, O: Open, V: Virtual
    
    Rcc = Acc Fcc^a + Bcc Fcc^b
    Roo = Aoo Foo^a + Boo Foo^b
    Rvv = Avv Fvv^a + Bvv Fvv^b
    Fcv = (F^a + F^b)/2
    
    F^a: Fock matrix alpha (MO), F^b: Fock matrix beta (MO)
    A,B: Coupling parameters
    
    J. Chem. Phys. 133, 141102 (2010), https://doi.org/10.1063/1.3503173
    Coupling parameters from J. Chem. Phys. 125, 204110 (2006); https://doi.org/10.1063/1.2393223.
          cc   oo   vv
     A  -0.5  0.5  1.5
     B   1.5  0.5 -0.5
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_shells_closed`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`fock_matrix_mo_alpha`
       * :c:data:`fock_matrix_mo_beta`
       * :c:data:`frozen_orb_scf`
       * :c:data:`list_act`
       * :c:data:`list_core`
       * :c:data:`list_inact`
       * :c:data:`list_virt`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_orb`
       * :c:data:`n_inact_orb`
       * :c:data:`n_virt_orb`
       * :c:data:`no_oa_or_av_opt`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`fock_matrix_ao`

 
.. c:var:: fock_matrix_mo


    File : :file:`scf_utils/fock_matrix.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_matrix_mo	(mo_num,mo_num)
        double precision, allocatable	:: fock_matrix_diag_mo	(mo_num)


    Fock matrix on the MO basis.
    For open shells, the ROHF Fock Matrix is ::
    
          |  Rcc  |  F^b  |  Fcv  |
          |-----------------------|
          |  F^b  |  Roo  |  F^a  |
          |-----------------------|
          |  Fcv  |  F^a  |  Rvv  |
    
    C: Core, O: Open, V: Virtual
    
    Rcc = Acc Fcc^a + Bcc Fcc^b
    Roo = Aoo Foo^a + Boo Foo^b
    Rvv = Avv Fvv^a + Bvv Fvv^b
    Fcv = (F^a + F^b)/2
    
    F^a: Fock matrix alpha (MO), F^b: Fock matrix beta (MO)
    A,B: Coupling parameters
    
    J. Chem. Phys. 133, 141102 (2010), https://doi.org/10.1063/1.3503173
    Coupling parameters from J. Chem. Phys. 125, 204110 (2006); https://doi.org/10.1063/1.2393223.
          cc   oo   vv
     A  -0.5  0.5  1.5
     B   1.5  0.5 -0.5
    

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_shells_closed`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`fock_matrix_mo_alpha`
       * :c:data:`fock_matrix_mo_beta`
       * :c:data:`frozen_orb_scf`
       * :c:data:`list_act`
       * :c:data:`list_core`
       * :c:data:`list_inact`
       * :c:data:`list_virt`
       * :c:data:`mo_num`
       * :c:data:`n_act_orb`
       * :c:data:`n_core_orb`
       * :c:data:`n_inact_orb`
       * :c:data:`n_virt_orb`
       * :c:data:`no_oa_or_av_opt`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`fock_matrix_ao`

 
.. c:var:: fock_matrix_mo_alpha


    File : :file:`scf_utils/fock_matrix.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_matrix_mo_alpha	(mo_num,mo_num)


    Fock matrix on the MO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`fock_matrix_ao_alpha`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_mo`

 
.. c:var:: fock_matrix_mo_beta


    File : :file:`scf_utils/fock_matrix.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fock_matrix_mo_beta	(mo_num,mo_num)


    Fock matrix on the MO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`fock_matrix_ao_alpha`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_mo`

 
.. c:var:: fps_spf_matrix_ao


    File : :file:`scf_utils/diis.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fps_spf_matrix_ao	(AO_num,AO_num)


    Commutator FPS - SPF

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`fock_matrix_ao`
       * :c:data:`scf_density_matrix_ao`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`fps_spf_matrix_mo`

 
.. c:var:: fps_spf_matrix_mo


    File : :file:`scf_utils/diis.irp.f`

    .. code:: fortran

        double precision, allocatable	:: fps_spf_matrix_mo	(mo_num,mo_num)


    Commutator FPS - SPF in MO basis

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`fps_spf_matrix_ao`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`


 
.. c:var:: scf_density_matrix_ao


    File : :file:`scf_utils/scf_density_matrix_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: scf_density_matrix_ao	(ao_num,ao_num)


    Sum of :math:`\alpha`  and :math:`\beta`  density matrices

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`all_shells_closed`
       * :c:data:`ao_num`
       * :c:data:`scf_density_matrix_ao_alpha`
       * :c:data:`scf_density_matrix_ao_beta`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_integral_alpha_chol`
       * :c:data:`fps_spf_matrix_ao`

 
.. c:var:: scf_density_matrix_ao_alpha


    File : :file:`scf_utils/scf_density_matrix_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: scf_density_matrix_ao_alpha	(ao_num,ao_num)


    :math:`C.C^t`  over :math:`\alpha`  MOs

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`elec_alpha_num`
       * :c:data:`mo_coef`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_integral_alpha`
       * :c:data:`ao_two_e_integral_alpha_chol`
       * :c:data:`hf_energy`
       * :c:data:`hf_kinetic_energy`
       * :c:data:`mcscf_fock_alpha_ao`
       * :c:data:`scf_density_matrix_ao`
       * :c:data:`scf_energy`

 
.. c:var:: scf_density_matrix_ao_beta


    File : :file:`scf_utils/scf_density_matrix_ao.irp.f`

    .. code:: fortran

        double precision, allocatable	:: scf_density_matrix_ao_beta	(ao_num,ao_num)


    :math:`C.C^t`  over :math:`\beta`  MOs

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_coef`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`ao_two_e_integral_alpha`
       * :c:data:`ao_two_e_integral_alpha_chol`
       * :c:data:`hf_energy`
       * :c:data:`hf_kinetic_energy`
       * :c:data:`mcscf_fock_alpha_ao`
       * :c:data:`scf_density_matrix_ao`
       * :c:data:`scf_energy`

 
.. c:var:: scf_energy


    File : :file:`scf_utils/fock_matrix.irp.f`

    .. code:: fortran

        double precision	:: scf_energy	


    Hartree-Fock energy

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_one_e_integrals`
       * :c:data:`extra_e_contrib_density`
       * :c:data:`fock_matrix_ao_alpha`
       * :c:data:`nuclear_repulsion`
       * :c:data:`scf_density_matrix_ao_alpha`
       * :c:data:`scf_density_matrix_ao_beta`


 
.. c:var:: threshold_diis_nonzero


    File : :file:`scf_utils/diis.irp.f`

    .. code:: fortran

        double precision	:: threshold_diis_nonzero	


    If threshold_DIIS is zero, choose sqrt(thresh_scf)

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`thresh_scf`
       * :c:data:`threshold_diis`


 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: damping_scf:


    File : :file:`scf_utils/damping_scf.irp.f`

    .. code:: fortran

        subroutine damping_SCF



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`fock_matrix_ao`
       * :c:data:`fock_matrix_mo`
       * :c:data:`frozen_orb_scf`
       * :c:data:`mo_coef`
       * :c:data:`mo_label`
       * :c:data:`n_it_scf_max`
       * :c:data:`scf_density_matrix_ao_alpha`
       * :c:data:`scf_density_matrix_ao_beta`
       * :c:data:`scf_energy`
       * :c:data:`thresh_scf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`ezfio_set_hartree_fock_energy`
       * :c:func:`initialize_mo_coef_begin_iteration`
       * :c:func:`mo_as_eigvectors_of_mo_matrix`
       * :c:func:`reorder_core_orb`
       * :c:func:`save_mos`
       * :c:func:`write_double`
       * :c:func:`write_time`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`scf_density_matrix_ao_alpha`
       * :c:data:`scf_density_matrix_ao_beta`
       * :c:data:`mo_coef`

 
.. c:function:: huckel_guess:


    File : :file:`scf_utils/huckel.irp.f`

    .. code:: fortran

        subroutine huckel_guess


    Build the MOs using the extended Huckel model

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_one_e_integrals`
       * :c:data:`ao_overlap`
       * :c:data:`ao_two_e_integral_alpha`
       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`fock_matrix_ao_alpha`
       * :c:data:`fock_matrix_ao_alpha`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`create_guess`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`orthonormalize_mos`
       * :c:func:`restore_symmetry`
       * :c:func:`save_mos`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_ao_alpha`
       * :c:data:`fock_matrix_ao_alpha`
       * :c:data:`mo_coef`

 
.. c:function:: reorder_mo_max_overlap:


    File : :file:`scf_utils/reorder_mo_max_overlap.irp.f`

    .. code:: fortran

        subroutine reorder_mo_max_overlap


    routines that compute the projection of each MO of the current `mo_coef` on the space spanned by the occupied orbitals of `mo_coef_begin_iteration`

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_num`
       * :c:data:`ao_overlap`
       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_coef_begin_iteration`
       * :c:data:`mo_num`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`roothaan_hall_scf`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`dsort`

 
.. c:function:: roothaan_hall_scf:


    File : :file:`scf_utils/roothaan_hall_scf.irp.f`

    .. code:: fortran

        subroutine Roothaan_Hall_SCF


    Roothaan-Hall algorithm for SCF Hartree-Fock calculation

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_md5`
       * :c:data:`ao_num`
       * :c:data:`do_mom`
       * :c:data:`eigenvectors_fock_matrix_mo`
       * :c:data:`fock_matrix_ao`
       * :c:data:`fock_matrix_ao_alpha`
       * :c:data:`fock_matrix_ao_alpha`
       * :c:data:`fock_matrix_mo`
       * :c:data:`fps_spf_matrix_ao`
       * :c:data:`fps_spf_matrix_mo`
       * :c:data:`frozen_orb_scf`
       * :c:data:`json_int_fmt`
       * :c:data:`json_unit`
       * :c:data:`level_shift`
       * :c:data:`max_dim_diis`
       * :c:data:`mo_coef`
       * :c:data:`mo_label`
       * :c:data:`mo_num`
       * :c:data:`mo_occ`
       * :c:data:`n_it_scf_max`
       * :c:data:`scf_algorithm`
       * :c:data:`scf_energy`
       * :c:data:`thresh_scf`
       * :c:data:`threshold_diis_nonzero`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`run`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`dgemm`
       * :c:func:`extrapolate_fock_matrix`
       * :c:func:`initialize_mo_coef_begin_iteration`
       * :c:func:`lock_io`
       * :c:func:`mo_as_eigvectors_of_mo_matrix`
       * :c:func:`orthonormalize_mos`
       * :c:func:`pivoted_cholesky`
       * :c:func:`reorder_core_orb`
       * :c:func:`reorder_mo_max_overlap`
       * :c:func:`restore_symmetry`
       * :c:func:`save_mos`
       * :c:func:`sleep`
       * :c:func:`unlock_io`
       * :c:func:`write_double`
       * :c:func:`write_time`

    Touches:

    .. hlist::
       :columns: 3

       * :c:data:`fock_matrix_ao_alpha`
       * :c:data:`fock_matrix_ao_alpha`
       * :c:data:`mo_coef`
       * :c:data:`level_shift`
       * :c:data:`mo_coef`

