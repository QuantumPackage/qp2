.. _module_utils_periodic: 
 
.. program:: utils_periodic 
 
.. default-role:: option 
 
utils_periodic
==============

No README.rst in this module; the providers and
subroutines below are documented from the irpf90 man pages. 
 
 
Programs 
-------- 
 
 * :ref:`import_integrals_ao_periodic` 
 * :ref:`import_mo_coef_periodic` 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: void run(void)


    File : :file:`import_mo_coef_periodic.irp.f`

    .. code:: fortran

        subroutine run


    Reads C.qp line by line: a line "i j re im" is stored
    in mo_coef(i,j) and mo_coef_imag(i,j) (real and
    imaginary parts). Sets mo_label to "None" and writes
    the MO basis to the |EZFIO| database via save_mos.
    TODO(doc): is mo_coef_imag stored in the database?
    save_mos writes mo_coef but not mo_coef_imag.

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_coef`
       * :c:data:`mo_coef_imag`
       * :c:data:`mo_label`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`import_integrals_ao_periodic`
       * :c:func:`import_mo_coef_periodic`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`save_mos`

