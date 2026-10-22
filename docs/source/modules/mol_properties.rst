.. _module_mol_properties: 
 
.. program:: mol_properties 
 
.. default-role:: option 
 
mol_properties
==============

No README.rst in this module; the providers and
subroutines below are documented from the irpf90 man pages. 
 
 
EZFIO parameters 
---------------- 
 
.. option:: dft_grid_mo_dipole
 
    If true, you compute the MO x/y/z dipole with DFT's grid
 
    Default: false
 
.. option:: print_all_transitions
 
    If true, print the transition between all the states
 
    Default: false
 
.. option:: calc_dipole_moment
 
    If true, the electric dipole moment will be computed
 
    Default: false
 
.. option:: calc_tr_dipole_moment
 
    If true and N_states > 1, the transition electric dipole moment will be computed
 
    Default: false
 
.. option:: calc_osc_str
 
    If true and N_states > 1, the oscillator strength will be computed
 
    Default: false
 
.. option:: calc_energy_components
 
    If true, the components of the energy (1e, 2e, kinetic) will be computed
 
    Default: false
 
 
Programs 
-------- 
 
 * :ref:`properties` 
 
Providers 
--------- 
 
.. c:var:: ci_energy_no_diag


    File : :file:`ci_energy_no_diag.irp.f`

    .. code:: fortran

        double precision, allocatable	:: ci_energy_no_diag	(N_states)


    CI energy from density matrices and integrals
    Avoid the rediagonalization for ci_energy

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`n_states`
       * :c:data:`nuclear_repulsion`
       * :c:data:`psi_energy`


 
.. c:var:: mo_deriv_1_x


    File : :file:`mo_deriv_1.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_deriv_1_x	(mo_num,mo_num)
        double precision, allocatable	:: mo_deriv_1_y	(mo_num,mo_num)
        double precision, allocatable	:: mo_deriv_1_z	(mo_num,mo_num)


    array of the integrals of MO_i * d/dx MO_j
    array of the integrals of MO_i * d/dy MO_j
    array of the integrals of MO_i * d/dz MO_j

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_deriv_1_x`
       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`multi_s_deriv_1`

 
.. c:var:: mo_deriv_1_y


    File : :file:`mo_deriv_1.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_deriv_1_x	(mo_num,mo_num)
        double precision, allocatable	:: mo_deriv_1_y	(mo_num,mo_num)
        double precision, allocatable	:: mo_deriv_1_z	(mo_num,mo_num)


    array of the integrals of MO_i * d/dx MO_j
    array of the integrals of MO_i * d/dy MO_j
    array of the integrals of MO_i * d/dz MO_j

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_deriv_1_x`
       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`multi_s_deriv_1`

 
.. c:var:: mo_deriv_1_z


    File : :file:`mo_deriv_1.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_deriv_1_x	(mo_num,mo_num)
        double precision, allocatable	:: mo_deriv_1_y	(mo_num,mo_num)
        double precision, allocatable	:: mo_deriv_1_z	(mo_num,mo_num)


    array of the integrals of MO_i * d/dx MO_j
    array of the integrals of MO_i * d/dy MO_j
    array of the integrals of MO_i * d/dz MO_j

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ao_deriv_1_x`
       * :c:data:`ao_num`
       * :c:data:`mo_coef`
       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`multi_s_deriv_1`

 
.. c:var:: mo_prop_dipole_x


    File : :file:`mo_dipole_properties.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_prop_dipole_x	(mo_num,mo_num)
        double precision, allocatable	:: mo_prop_dipole_y	(mo_num,mo_num)
        double precision, allocatable	:: mo_prop_dipole_z	(mo_num,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`dft_grid_mo_dipole`
       * :c:data:`mo_dipole_x`
       * :c:data:`mo_dipole_x_becke`
       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`multi_s_dipole_moment`

 
.. c:var:: mo_prop_dipole_y


    File : :file:`mo_dipole_properties.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_prop_dipole_x	(mo_num,mo_num)
        double precision, allocatable	:: mo_prop_dipole_y	(mo_num,mo_num)
        double precision, allocatable	:: mo_prop_dipole_z	(mo_num,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`dft_grid_mo_dipole`
       * :c:data:`mo_dipole_x`
       * :c:data:`mo_dipole_x_becke`
       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`multi_s_dipole_moment`

 
.. c:var:: mo_prop_dipole_z


    File : :file:`mo_dipole_properties.irp.f`

    .. code:: fortran

        double precision, allocatable	:: mo_prop_dipole_x	(mo_num,mo_num)
        double precision, allocatable	:: mo_prop_dipole_y	(mo_num,mo_num)
        double precision, allocatable	:: mo_prop_dipole_z	(mo_num,mo_num)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`dft_grid_mo_dipole`
       * :c:data:`mo_dipole_x`
       * :c:data:`mo_dipole_x_becke`
       * :c:data:`mo_num`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`multi_s_dipole_moment`

 
.. c:var:: multi_s_deriv_1


    File : :file:`multi_s_deriv_1.irp.f`

    .. code:: fortran

        double precision, allocatable	:: multi_s_deriv_1	(N_states,N_states)
        double precision, allocatable	:: multi_s_x_deriv_1	(N_states,N_states)
        double precision, allocatable	:: multi_s_y_deriv_1	(N_states,N_states)
        double precision, allocatable	:: multi_s_z_deriv_1	(N_states,N_states)


    Providers for :
    <Psi_m|v_x|Psi_n>
    <Psi_m|v_y|Psi_n>
    <Psi_m|v_z|Psi_n>
    ||v|| = sqrt(v_x^2 + v_y^2 + v_z^2)
    v_x = d/dx
    Cf. multi_s_dipole_moment for the equations

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_deriv_1_x`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`one_e_tr_dm_mo`


 
.. c:var:: multi_s_dipole_moment


    File : :file:`multi_s_dipole_moment.irp.f`

    .. code:: fortran

        double precision, allocatable	:: multi_s_dipole_moment	(N_states,N_states)
        double precision, allocatable	:: multi_s_x_dipole_moment	(N_states,N_states)
        double precision, allocatable	:: multi_s_y_dipole_moment	(N_states,N_states)
        double precision, allocatable	:: multi_s_z_dipole_moment	(N_states,N_states)


    Providers for :
    <\Psi_m|\mu_x|\Psi_n>
    <\Psi_m|\mu_y|\Psi_n>
    <\Psi_m|\mu_z|\Psi_n>
    ||\mu|| = \sqrt{\mu_x^2 + \mu_y^2 + \mu_z^2}
    
    <\Psi_n|x| \Psi_m > = \sum_p \gamma_{pp}^{nm} \bra{\phi_p} x \ket{\phi_p}
      + \sum_{pq, p \neq q} \gamma_{pq}^{nm} \bra{\phi_p} x \ket{\phi_q}
    \Psi: wf
    n,m indexes for the states
    p,q: general spatial MOs
    gamma^{nm}: density matrix \bra{\Psi^n} a^{\dagger}_a a_i \ket{\Psi^m}

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_num`
       * :c:data:`mo_prop_dipole_x`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`one_e_tr_dm_mo`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`qp_max_mem`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`multi_s_x_dipole_moment_eigenvec`

 
.. c:var:: multi_s_x_deriv_1


    File : :file:`multi_s_deriv_1.irp.f`

    .. code:: fortran

        double precision, allocatable	:: multi_s_deriv_1	(N_states,N_states)
        double precision, allocatable	:: multi_s_x_deriv_1	(N_states,N_states)
        double precision, allocatable	:: multi_s_y_deriv_1	(N_states,N_states)
        double precision, allocatable	:: multi_s_z_deriv_1	(N_states,N_states)


    Providers for :
    <Psi_m|v_x|Psi_n>
    <Psi_m|v_y|Psi_n>
    <Psi_m|v_z|Psi_n>
    ||v|| = sqrt(v_x^2 + v_y^2 + v_z^2)
    v_x = d/dx
    Cf. multi_s_dipole_moment for the equations

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_deriv_1_x`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`one_e_tr_dm_mo`


 
.. c:var:: multi_s_x_dipole_moment


    File : :file:`multi_s_dipole_moment.irp.f`

    .. code:: fortran

        double precision, allocatable	:: multi_s_dipole_moment	(N_states,N_states)
        double precision, allocatable	:: multi_s_x_dipole_moment	(N_states,N_states)
        double precision, allocatable	:: multi_s_y_dipole_moment	(N_states,N_states)
        double precision, allocatable	:: multi_s_z_dipole_moment	(N_states,N_states)


    Providers for :
    <\Psi_m|\mu_x|\Psi_n>
    <\Psi_m|\mu_y|\Psi_n>
    <\Psi_m|\mu_z|\Psi_n>
    ||\mu|| = \sqrt{\mu_x^2 + \mu_y^2 + \mu_z^2}
    
    <\Psi_n|x| \Psi_m > = \sum_p \gamma_{pp}^{nm} \bra{\phi_p} x \ket{\phi_p}
      + \sum_{pq, p \neq q} \gamma_{pq}^{nm} \bra{\phi_p} x \ket{\phi_q}
    \Psi: wf
    n,m indexes for the states
    p,q: general spatial MOs
    gamma^{nm}: density matrix \bra{\Psi^n} a^{\dagger}_a a_i \ket{\Psi^m}

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_num`
       * :c:data:`mo_prop_dipole_x`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`one_e_tr_dm_mo`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`qp_max_mem`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`multi_s_x_dipole_moment_eigenvec`

 
.. c:var:: multi_s_x_dipole_moment_eigenval


    File : :file:`multi_s_dipole_moment.irp.f`

    .. code:: fortran

        double precision, allocatable	:: multi_s_x_dipole_moment_eigenvec	(N_states,N_states)
        double precision, allocatable	:: multi_s_y_dipole_moment_eigenvec	(N_states,N_states)
        double precision, allocatable	:: multi_s_z_dipole_moment_eigenvec	(N_states,N_states)
        double precision, allocatable	:: multi_s_x_dipole_moment_eigenval	(N_states)
        double precision, allocatable	:: multi_s_y_dipole_moment_eigenval	(N_states)
        double precision, allocatable	:: multi_s_z_dipole_moment_eigenval	(N_states)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`multi_s_dipole_moment`
       * :c:data:`multi_s_dipole_moment`
       * :c:data:`multi_s_dipole_moment`
       * :c:data:`multi_s_dipole_moment`
       * :c:data:`n_states`


 
.. c:var:: multi_s_x_dipole_moment_eigenvec


    File : :file:`multi_s_dipole_moment.irp.f`

    .. code:: fortran

        double precision, allocatable	:: multi_s_x_dipole_moment_eigenvec	(N_states,N_states)
        double precision, allocatable	:: multi_s_y_dipole_moment_eigenvec	(N_states,N_states)
        double precision, allocatable	:: multi_s_z_dipole_moment_eigenvec	(N_states,N_states)
        double precision, allocatable	:: multi_s_x_dipole_moment_eigenval	(N_states)
        double precision, allocatable	:: multi_s_y_dipole_moment_eigenval	(N_states)
        double precision, allocatable	:: multi_s_z_dipole_moment_eigenval	(N_states)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`multi_s_dipole_moment`
       * :c:data:`multi_s_dipole_moment`
       * :c:data:`multi_s_dipole_moment`
       * :c:data:`multi_s_dipole_moment`
       * :c:data:`n_states`


 
.. c:var:: multi_s_y_deriv_1


    File : :file:`multi_s_deriv_1.irp.f`

    .. code:: fortran

        double precision, allocatable	:: multi_s_deriv_1	(N_states,N_states)
        double precision, allocatable	:: multi_s_x_deriv_1	(N_states,N_states)
        double precision, allocatable	:: multi_s_y_deriv_1	(N_states,N_states)
        double precision, allocatable	:: multi_s_z_deriv_1	(N_states,N_states)


    Providers for :
    <Psi_m|v_x|Psi_n>
    <Psi_m|v_y|Psi_n>
    <Psi_m|v_z|Psi_n>
    ||v|| = sqrt(v_x^2 + v_y^2 + v_z^2)
    v_x = d/dx
    Cf. multi_s_dipole_moment for the equations

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_deriv_1_x`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`one_e_tr_dm_mo`


 
.. c:var:: multi_s_y_dipole_moment


    File : :file:`multi_s_dipole_moment.irp.f`

    .. code:: fortran

        double precision, allocatable	:: multi_s_dipole_moment	(N_states,N_states)
        double precision, allocatable	:: multi_s_x_dipole_moment	(N_states,N_states)
        double precision, allocatable	:: multi_s_y_dipole_moment	(N_states,N_states)
        double precision, allocatable	:: multi_s_z_dipole_moment	(N_states,N_states)


    Providers for :
    <\Psi_m|\mu_x|\Psi_n>
    <\Psi_m|\mu_y|\Psi_n>
    <\Psi_m|\mu_z|\Psi_n>
    ||\mu|| = \sqrt{\mu_x^2 + \mu_y^2 + \mu_z^2}
    
    <\Psi_n|x| \Psi_m > = \sum_p \gamma_{pp}^{nm} \bra{\phi_p} x \ket{\phi_p}
      + \sum_{pq, p \neq q} \gamma_{pq}^{nm} \bra{\phi_p} x \ket{\phi_q}
    \Psi: wf
    n,m indexes for the states
    p,q: general spatial MOs
    gamma^{nm}: density matrix \bra{\Psi^n} a^{\dagger}_a a_i \ket{\Psi^m}

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_num`
       * :c:data:`mo_prop_dipole_x`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`one_e_tr_dm_mo`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`qp_max_mem`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`multi_s_x_dipole_moment_eigenvec`

 
.. c:var:: multi_s_y_dipole_moment_eigenval


    File : :file:`multi_s_dipole_moment.irp.f`

    .. code:: fortran

        double precision, allocatable	:: multi_s_x_dipole_moment_eigenvec	(N_states,N_states)
        double precision, allocatable	:: multi_s_y_dipole_moment_eigenvec	(N_states,N_states)
        double precision, allocatable	:: multi_s_z_dipole_moment_eigenvec	(N_states,N_states)
        double precision, allocatable	:: multi_s_x_dipole_moment_eigenval	(N_states)
        double precision, allocatable	:: multi_s_y_dipole_moment_eigenval	(N_states)
        double precision, allocatable	:: multi_s_z_dipole_moment_eigenval	(N_states)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`multi_s_dipole_moment`
       * :c:data:`multi_s_dipole_moment`
       * :c:data:`multi_s_dipole_moment`
       * :c:data:`multi_s_dipole_moment`
       * :c:data:`n_states`


 
.. c:var:: multi_s_y_dipole_moment_eigenvec


    File : :file:`multi_s_dipole_moment.irp.f`

    .. code:: fortran

        double precision, allocatable	:: multi_s_x_dipole_moment_eigenvec	(N_states,N_states)
        double precision, allocatable	:: multi_s_y_dipole_moment_eigenvec	(N_states,N_states)
        double precision, allocatable	:: multi_s_z_dipole_moment_eigenvec	(N_states,N_states)
        double precision, allocatable	:: multi_s_x_dipole_moment_eigenval	(N_states)
        double precision, allocatable	:: multi_s_y_dipole_moment_eigenval	(N_states)
        double precision, allocatable	:: multi_s_z_dipole_moment_eigenval	(N_states)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`multi_s_dipole_moment`
       * :c:data:`multi_s_dipole_moment`
       * :c:data:`multi_s_dipole_moment`
       * :c:data:`multi_s_dipole_moment`
       * :c:data:`n_states`


 
.. c:var:: multi_s_z_deriv_1


    File : :file:`multi_s_deriv_1.irp.f`

    .. code:: fortran

        double precision, allocatable	:: multi_s_deriv_1	(N_states,N_states)
        double precision, allocatable	:: multi_s_x_deriv_1	(N_states,N_states)
        double precision, allocatable	:: multi_s_y_deriv_1	(N_states,N_states)
        double precision, allocatable	:: multi_s_z_deriv_1	(N_states,N_states)


    Providers for :
    <Psi_m|v_x|Psi_n>
    <Psi_m|v_y|Psi_n>
    <Psi_m|v_z|Psi_n>
    ||v|| = sqrt(v_x^2 + v_y^2 + v_z^2)
    v_x = d/dx
    Cf. multi_s_dipole_moment for the equations

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_deriv_1_x`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`one_e_tr_dm_mo`


 
.. c:var:: multi_s_z_dipole_moment


    File : :file:`multi_s_dipole_moment.irp.f`

    .. code:: fortran

        double precision, allocatable	:: multi_s_dipole_moment	(N_states,N_states)
        double precision, allocatable	:: multi_s_x_dipole_moment	(N_states,N_states)
        double precision, allocatable	:: multi_s_y_dipole_moment	(N_states,N_states)
        double precision, allocatable	:: multi_s_z_dipole_moment	(N_states,N_states)


    Providers for :
    <\Psi_m|\mu_x|\Psi_n>
    <\Psi_m|\mu_y|\Psi_n>
    <\Psi_m|\mu_z|\Psi_n>
    ||\mu|| = \sqrt{\mu_x^2 + \mu_y^2 + \mu_z^2}
    
    <\Psi_n|x| \Psi_m > = \sum_p \gamma_{pp}^{nm} \bra{\phi_p} x \ket{\phi_p}
      + \sum_{pq, p \neq q} \gamma_{pq}^{nm} \bra{\phi_p} x \ket{\phi_q}
    \Psi: wf
    n,m indexes for the states
    p,q: general spatial MOs
    gamma^{nm}: density matrix \bra{\Psi^n} a^{\dagger}_a a_i \ket{\Psi^m}

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`elec_alpha_num`
       * :c:data:`elec_beta_num`
       * :c:data:`mo_num`
       * :c:data:`mo_prop_dipole_x`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_coord`
       * :c:data:`nucl_num`
       * :c:data:`one_e_tr_dm_mo`
       * :c:data:`psi_bilinear_matrix_transp_values`
       * :c:data:`psi_bilinear_matrix_values`
       * :c:data:`psi_det_alpha_unique`
       * :c:data:`psi_det_beta_unique`
       * :c:data:`qp_max_mem`

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`multi_s_x_dipole_moment_eigenvec`

 
.. c:var:: multi_s_z_dipole_moment_eigenval


    File : :file:`multi_s_dipole_moment.irp.f`

    .. code:: fortran

        double precision, allocatable	:: multi_s_x_dipole_moment_eigenvec	(N_states,N_states)
        double precision, allocatable	:: multi_s_y_dipole_moment_eigenvec	(N_states,N_states)
        double precision, allocatable	:: multi_s_z_dipole_moment_eigenvec	(N_states,N_states)
        double precision, allocatable	:: multi_s_x_dipole_moment_eigenval	(N_states)
        double precision, allocatable	:: multi_s_y_dipole_moment_eigenval	(N_states)
        double precision, allocatable	:: multi_s_z_dipole_moment_eigenval	(N_states)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`multi_s_dipole_moment`
       * :c:data:`multi_s_dipole_moment`
       * :c:data:`multi_s_dipole_moment`
       * :c:data:`multi_s_dipole_moment`
       * :c:data:`n_states`


 
.. c:var:: multi_s_z_dipole_moment_eigenvec


    File : :file:`multi_s_dipole_moment.irp.f`

    .. code:: fortran

        double precision, allocatable	:: multi_s_x_dipole_moment_eigenvec	(N_states,N_states)
        double precision, allocatable	:: multi_s_y_dipole_moment_eigenvec	(N_states,N_states)
        double precision, allocatable	:: multi_s_z_dipole_moment_eigenvec	(N_states,N_states)
        double precision, allocatable	:: multi_s_x_dipole_moment_eigenval	(N_states)
        double precision, allocatable	:: multi_s_y_dipole_moment_eigenval	(N_states)
        double precision, allocatable	:: multi_s_z_dipole_moment_eigenval	(N_states)



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`multi_s_dipole_moment`
       * :c:data:`multi_s_dipole_moment`
       * :c:data:`multi_s_dipole_moment`
       * :c:data:`multi_s_dipole_moment`
       * :c:data:`n_states`


 
 
Subroutines / functions 
----------------------- 
 
.. c:macro:: print_dipole_moment:


    File : :file:`print_properties.irp.f`

    .. code:: fortran

        subroutine print_dipole_moment


    To print the dipole moment ||<\Psi_i|\mu|\Psi_i>|| and its x,y,z components

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ha_to_ev`
       * :c:data:`multi_s_dipole_moment`
       * :c:data:`n_states`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`print_mol_properties`

 
.. c:macro:: print_energy_components:


    File : :file:`print_e_components.irp.f`

    .. code:: fortran

        subroutine print_energy_components()


    Prints the different components of the energy.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_integrals_n_e`
       * :c:data:`mo_kinetic_integrals`
       * :c:data:`mo_num`
       * :c:data:`n_states`
       * :c:data:`nuclear_repulsion`
       * :c:data:`one_e_dm_mo_alpha`
       * :c:data:`psi_energy`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`print_mol_properties`

 
.. c:macro:: print_mol_properties:


    File : :file:`print_mol_properties.irp.f`

    .. code:: fortran

        subroutine print_mol_properties()


    Run the propertie calculations

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`calc_dipole_moment`
       * :c:data:`calc_energy_components`
       * :c:data:`calc_osc_str`
       * :c:data:`calc_tr_dipole_moment`
       * :c:data:`n_states`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`mol_properties`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`print_dipole_moment`
       * :c:func:`print_energy_components`
       * :c:func:`print_oscillator_strength`
       * :c:func:`print_transition_dipole_moment`

 
.. c:macro:: print_oscillator_strength:


    File : :file:`print_properties.irp.f`

    .. code:: fortran

        subroutine print_oscillator_strength


    https://doi.org/10.1016/j.cplett.2004.03.126
    Oscillator strength in:
    - length gauge, f^l_{ij} = 2/3 (E_i - E_j) <\Psi_i|r|\Psi_j> <\Psi_j|r|\Psi_i>
    - velocity gauge, f^v_{ij} = 2/3 (E_i - E_j)^(-1) <\Psi_i|v|\Psi_j> <\Psi_j|v|\Psi_i>
    - mixed gauge, f^m_{ij} = -2i/3 <\Psi_i|r|\Psi_j> <\Psi_j|v|\Psi_i>

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ci_energy_no_diag`
       * :c:data:`ha_to_ev`
       * :c:data:`multi_s_deriv_1`
       * :c:data:`multi_s_dipole_moment`
       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`print_all_transitions`
       * :c:data:`s2_values`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`print_mol_properties`

 
.. c:macro:: print_transition_dipole_moment:


    File : :file:`print_properties.irp.f`

    .. code:: fortran

        subroutine print_transition_dipole_moment


    To print the transition dipole moment ||<\Psi_i|\mu|\Psi_j>|| and its components along x, y and z

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ci_energy_no_diag`
       * :c:data:`ha_to_ev`
       * :c:data:`multi_s_dipole_moment`
       * :c:data:`n_det`
       * :c:data:`n_states`
       * :c:data:`print_all_transitions`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`print_mol_properties`

