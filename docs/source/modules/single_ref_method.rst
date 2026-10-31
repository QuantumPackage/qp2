.. _module_single_ref_method: 
 
.. program:: single_ref_method 
 
.. default-role:: option 
 
=================
single_ref_method
=================

Definition of the generators for single reference methods.

Including this module in the :file:`NEED` file of a program makes
the only generator the reference (|HF|) determinant:
:c:data:`n_det_generators` is 1 and :c:data:`psi_det_generators` is
the reference determinant.
 
 
 
Providers 
--------- 
 
.. c:var:: int n_det_generators


    File : :file:`single_ref_method/generators.irp.f`

    .. code:: fortran

        integer	:: n_det_generators	


    For Single reference wave functions, the number of generators is 1 : the
    Hartree-Fock determinant

    Needed by:

    .. hlist::
       :columns: 3

       * :c:data:`n_det_selectors`

 
.. c:var:: double psi_coef_generators


    File : :file:`single_ref_method/generators.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_generators	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_coef_generators	(psi_det_size,N_states)


    For Single reference wave functions, the generator is the
    Hartree-Fock determinant

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`hf_bitmask`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`


 
.. c:var:: int psi_det_generators


    File : :file:`single_ref_method/generators.irp.f`

    .. code:: fortran

        integer(bit_kind), allocatable	:: psi_det_generators	(N_int,2,psi_det_size)
        double precision, allocatable	:: psi_coef_generators	(psi_det_size,N_states)


    For Single reference wave functions, the generator is the
    Hartree-Fock determinant

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`hf_bitmask`
       * :c:data:`n_det`
       * :c:data:`n_int`
       * :c:data:`n_states`
       * :c:data:`psi_det`
       * :c:data:`psi_det_size`


 
.. c:var:: double select_max


    File : :file:`single_ref_method/generators.irp.f`

    .. code:: fortran

        double precision, allocatable	:: select_max	(1)


    Memo to skip useless selectors


 
.. c:var:: int size_select_max


    File : :file:`single_ref_method/generators.irp.f`

    .. code:: fortran

        integer	:: size_select_max	


    Size of select_max


