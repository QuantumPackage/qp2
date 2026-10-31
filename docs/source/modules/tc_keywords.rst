.. _module_tc_keywords: 
 
.. program:: tc_keywords 
 
.. default-role:: option 
 
tc_keywords
===========

No README.rst in this module; the providers and
subroutines below are documented from the irpf90 man pages. 
 
 
Providers 
--------- 
 
.. c:var:: j1b_coeff


    File : :file:`j1b_pen.irp.f`

    .. code:: fortran

        double precision, allocatable	:: j1b_coeff	(nucl_num)


    coefficients of the 1-body Jastrow

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`mpi_master`
       * :c:data:`nucl_num`


 
.. c:var:: j1b_pen


    File : :file:`j1b_pen.irp.f`

    .. code:: fortran

        double precision, allocatable	:: j1b_pen	(nucl_num)


    exponents of the 1-body Jastrow

    Needs:

    .. hlist::
       :columns: 3

       * :c:data:`ezfio_filename`
       * :c:data:`mpi_master`
       * :c:data:`nucl_charge`
       * :c:data:`nucl_num`


 
 
Subroutines / functions 
----------------------- 
 
.. c:function:: tc_keywords:


    File : :file:`tc_keywords.irp.f`

    .. code:: fortran

        subroutine tc_keywords


    TODO : Put the documentation of the program here

