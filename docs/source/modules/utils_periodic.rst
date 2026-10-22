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
 
.. c:macro:: run:


    File : :file:`import_mo_coef_periodic.irp.f`

    .. code:: fortran

        subroutine run



    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`mo_coef`
       * :c:data:`mo_coef_imag`
       * :c:data:`mo_label`

    Called by:

    .. hlist::
       :columns: 3

       * :c:func:`import_mo_coef_periodic`
       * :c:func:`print_integrals`

    Calls:

    .. hlist::
       :columns: 3

       * :c:func:`save_mos`

