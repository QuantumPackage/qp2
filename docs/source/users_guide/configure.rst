=========
configure
=========


.. program:: configure

Program that can either configure the compilations options or download/install
external dependencies (see the installation description). 

Usage
-----

.. code:: bash

  ./configure [-h | -c <file> | -i <package>] 

.. option:: -c <file>, --config <file>

   Define a configuration file, in :file`${QP_ROOT}/config/`

.. option::  -h, --help

   Print the help message

.. option::  -i <package>, --install <package>

   Try to install <package>. Use at your own risk.

Example
-------

.. code:: bash

     ./configure 
     ./configure -c config/gfortran.cfg

