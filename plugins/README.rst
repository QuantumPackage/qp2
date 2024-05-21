==============================
Tutorial for creating a plugin
==============================

Introduction: what is a plugin, and what tutorial will be about ?
=================================================================

The |QP| is split into two kinds of routines/global variables (i.e. *providers*): 
 1)  the **core modules** locatedin qp2/src/, which contains all the bulk of a quantum chemistry software (integrals, matrix elements between Slater determinants, linear algebra routines, DFT stuffs etc..)
 2) the **plugins** which are external routines/*providers* connected to the qp2/src/ routines/*providers*.
 
More precisely, a **plugin** of the |QP| is a directory where you can create routines, 
providers and executables that use all the global variables/functions/routines already created 
in the modules of qp2/src or in other plugins. 

Instead of giving a theoretical lecture on what is a plugin, 
we will go through a series of examples that allow you to do the following thing: 

1)   print out **one- and two-electron integrals** on the AO/MO basis, creates two providers which manipulate these objects, print out these providers, 

2)  browse the **Slater determinants stored** in the |EZFIO| wave function and compute their matrix elements, 

3) build the **Hamiltonian matrix** and **diagonalize** it either with **Lapack or Davidson**,

4)  print out the **one- and two-electron rdms**, 

5)   obtain the **AOs** and **MOs** on the **DFT grid**, together with the **density**,

How the tutorial will be done
-----------------------------

This tuto is as follows: 

 1)  you **READ THIS FILE UNTIL THE END** in order to get the big picture and vocabulary, 

 2) you go to the directory :file:`qp2/plugins/tuto_plugins/` and you will find detailed tutorials for each of the 5 examples. 

Creating a plugin: the basic
----------------------------

The first thing to do is to be in the QPSH mode: you execute the qp2/bin/qpsh script that essentially loads all 
the environement variables and allows for the completion of command lines in bash (that is an AMAZING feature :) 

Then, you need to known **where** you want to create your plugin, and what is the **name** of the plugin. 

.. important::

  The plugins are **NECESSARILY** located in qp2/plugins/, and from there you can create any structures of directories.


Ex: If you want to create a plugin named "my_fancy_plugin" in the directory plugins/plugins_test/, 
this goes with the command 

.. code:: bash

    qp plugins create -n my_fancy_plugin -r plugins_test/

Then, to create the plugin of your dreams, the two questions you need to answer are the following: 

1) What do I **need** to compute what I want, which means what are the **objects** that I need ?

   There are two kind of objects:

   + the *routines/functions*: 

     Ex: Linear algebra routines, integration routines etc ...

   + the global variables which are called the *providers*: 

     Ex: one-electron integrals, Slater determinants, density matrices etc ...

2) **Where do I find** these objects ? 

   The objects (routines/functions/providers) are necessarily created in other *modules/plugins*.  

.. seealso::

   The routine :c:func:`lapack_diagd` (which diagonalises a real hermitian matrix) is located in the file 
   :file:`qp2/src/utils/linear_algebra.irp.f` 
   therefore it "belongs" to the module :ref:`module_utils`

   The routine :c:func:`ao_to_mo` (which converts a given matrix A from the AO basis to the MO basis) is located in the file 
   :file:`qp2/src/mo_one_e_ints/ao_to_mo.irp.f`
   therefore it "belongs" to the module :ref:`module_mo_one_e_ints`

   The provider :c:data:`ao_one_e_integrals` (which is the integrals of one-body part of H on the AO basis) is located in the file 
   :file:`qp2/src/ao_one_e_ints/ao_one_e_ints.irp.f` 
   therefore it belongs to the module :ref:`module_ao_one_e_ints`

   The provider :c:data:`one_e_dm_mo_beta_average` (which is the state average beta density matrix on the MO basis) is located in the file 
   :file:`qp2/src/determinants/density_matrix.irp.f`
   therefore it belongs to the module :ref:`module_determinants`

To import all the variables that you need, you just need to write the name of the plugins in the :file:`NEED` file .

To import all the variables/routines of the module :ref:`module_utils`, :ref:`module_determinants` and :ref:`module_mo_one_e_ints`, the  :file:`NEED` file you will need is simply the following:

.. code:: bash

 cat NEED

 utils
 determinants
 mo_one_e_ints


.. important::

  There are **many** routines/providers in the core modules of QP. 

  Nevertheless, as everything is coded with the |IRPF90|, you can use the following amazing tools: :command:`irpman`

  :command:`irpman` can be used in command line in bash to obtain all the info on a routine or variable ! 


Example: execute the following command line : 

.. code:: bash

  irpman ao_one_e_integrals

Then all the information you need on :c:data:`ao_one_e_integrals` will appear on the screen. 
This includes

 - **where** the provider is created, (*i.e.* the actual file where the provider is designed)
 - the **type** of the provider (*i.e.* a logical, integer etc ...)
 - the **dimension** if it is an array, 
 - what other *providers* are **needed** to build this provider, 
 - what other *providers* **need** this provider. 


