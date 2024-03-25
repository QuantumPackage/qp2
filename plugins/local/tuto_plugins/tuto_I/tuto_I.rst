=============================================
Tuto I: One- and two-e integrals (20 minutes)
=============================================

Requirements 
------------
1) You know how to create an |EZFIO| file and run calculations with |QP| (check the tuto: `<https://quantumpackage.github.io/qp2/post/hartree-fock/>`_),

2) You have an |EZFIO| file with MOs created (with the :ref:`scf` executable for instance). As we are going to print out some integrals, don't take a too large system/basis (Ex: H2, cc-pVDZ is ok :)

3) You made an qp set_file YOUR_EZFIO_FILE_FOR_H2 in order to work on that ezfio folder. 

4) You have READ the :file:`qp2/plugins/README.rst` file to HAVE THE **VOCABULARY**. 

Our goals:
----------
We want to create a plugin to do the following things: 
 1) print out one- and two-electron integrals on the AO/MO basis,  

 2) creates two providers which manipulate these objects, 

 3) print out these providers. 

I) Getting started: creating the plugin
---------------------------------------
We will go step-by-step through these plugins. 

We will create a plugin named "plugin_I", and its location will be in "tuto_plugins". 
Therefore to create the plugin, we do: 

.. code:: bash

 qp plugins create -n plugin_I -r tuto_plugins

Then do an "ls" in qp2/plugins/tuto_plugins/ and you will find a directory called "plugin_I". 

In that directory you will find: 

1) a :file:`NEED` file that will eventually contain all the other modules/plugins needed by our "plugin_I",

2) a :file:`README.rst` file that you can and **SHOULD** modify in order to **DOCUMENT** what is doing the plugin,

3) a :file:`plugin_I.irp.f` file that is a program to be compiled and just printing "Hello world" 

II) Specifying the dependencies
-------------------------------
The next step is to know what are the other modules/plugins that we need to do what we want. 
We need here 

a) the one-electron integrals on the AO basis, which are computed in :file:`qp2/src/ao_one_e_ints/`

b) the one-electron integrals on the MO basis, which are computed in :file:`qp2/src/mo_one_e_ints/`

c) the two-electron integrals on the AO basis, which are computed in :file:`qp2/src/ao_two_e_ints/`

d) the two-electron integrals on the MO basis, which are computed in :file:`qp2/src/mo_two_e_ints/`

Therefore, we will need the following four modules:

 a) ao_one_e_ints
 b) mo_one_e_ints
 c) ao_two_e_ints
 d) mo_two_e_ints

You can then create the following "NEED" file by executing the following command 

.. code:: bash

 cat <<EOF > NEED
 ao_one_e_ints
 mo_one_e_ints 
 ao_two_e_ints
 mo_two_e_ints
 EOF

II) Installing the plugin 
-------------------------
Now that we have specified the various depenencies we need now to INSTALL the plugin, which means to create the equivalent of a Makefile for the compilation. 

To do it we simply do 

.. code:: bash

 qp plugins install plugin_I


III) Compiling the void plugin 
------------------------------
It is customary to compile first your "void" plugin, void in the sense that it does not contain anything else than the program printing "Hello world". 

To do so, just go in the plugin and execute the following command

.. code:: bash

  ninja 

It does a lot of stuffs, but it must conclude with something like 

.. code:: bash

 make: Leaving directory 'SOME_PATH_TOWARD_YOUR_QP2_DIRECTORY/qp2/ocaml'


Since that it has compiled, an executable "plugin_I" has been created. 

Also, if you make "ls" in the "plugin_I" you will notice that many symbolink links have been created, and among which the four modules that you included in the NEED file. 

All the other modules (Ex::ref:`module_ao_basis`, :ref:`module_utils`) are here because they are need by some of the  four modules that you need. 
The variables that we need are 

:data:`ao_one_e_integrals`

:data:`mo_one_e_integrals`

You can check them with 

.. code:: bash

 irpman ao_one_e_integrals


.. code:: bash

 irpman mo_one_e_integrals

in order to get some information on where they are created, and many more information. 
We will now create an executable such that it prints out the integrals. 


IV) Printing out the one-electron integrals
--------------------------------------------
We will now create a program that will print out the one-electron integrals on the AO and MO basis.

You can then copy the file :file:`qp2/plugins/tuto_plugins/tuto_I/print_one_e_h.irp.f` in your plugin. 

In this file you will see that we simply browse the two arrays :data:`ao_one_e_integrals` and :data:`mo_one_e_integrals`, which are the providers and we browse them until either :data:`ao_num` or :data:`mo_num` which are also providers representing the number of AOs or MOs. 


.. seealso::

 You can check these variables with :command:`irpman` ! 

If you recompile using |ninja| as before, and another executable has been created "print_one_e_h". 
Then, you can run the program on the ezfio file by doing 

.. code:: bash

 qp run print_one_e_h 

and will print out the data you need :) 

By the way, as the file :file:`plugin_I.irp.f` contains nothing but a "Hello world" print, you can simply remove it if you want. 

V) Printing out the two-electron integrals
------------------------------------------
We will now create a file that prints out the two-electron integrals in the AO and MO basis.
These can be accessed with the following subroutines :

1- :c:func:`get_ao_two_e_integral` for the AO basis 

2- :c:func:`get_two_e_integral` for the MO basis 


.. seealso::

 check them with irpman !

To print the two-electron integrals, you can copy the file :file:`qp2/plugins/tuto_plugins/tuto_I/print_two_e_h.irp.f` in your plugin and recompile with |ninja|. 
Then just run the program 

.. code:: bash

 qp run print_two_e_h

and it will print all the things you want :)

VI) Creating new providers and a program to print them
------------------------------------------------------
We will now create new providers that manipulates the objects that we just printed.
As an example, we will compute the trace of the one electron integrals in the AO and MO basis. 
In the file :file:`qp2/plugins/tuto_plugins/tuto_I/traces_one_e.irp.f` you will find the several new providers among which 

 1- :c:data:`trace_mo_one_e_ints` : simply the sum of the diagonal matrix element of the one-electron integrals

 2- :c:data:`trace_ao_one_e_ints` : the corresponding trace on the AO basis 
   .. math::                                                                                                              

      \text{Tr}({\bf h}{\bf S}^{-1}) =  \sum_{m,n} S^{-1}_{mn} h_{mn}


 3- :c:data:`trace_ao_one_e_ints_from_mo` : the trace on the AO basis with the integrals obtained first from the MO basis 
   .. math::                                                                                                              

      \text{Tr}({\bf \tilde{h}}{\bf S}^{-1}) = \text{Tr}\big({\bf SC h}({\bf SC }^T){\bf S}^{-1}\big) 

Just copy the :file:`qp2/plugins/tuto_plugins/tuto_I/traces_one_e.irp.f` in your plugin and recompile. 

.. seealso::

 Once it has compiled, check your new providers with :command:`irpman` !

As explained in the files :file:`qp2/plugins/tuto_plugins/tuto_I/traces_one_e.irp.f` and :file:`qp2/plugins/tuto_plugins/tuto_I/print_traces_on_e.irp.f`, :c:data:`trace_mo_one_e_ints` is equal to :c:data:`trace_ao_one_e_ints` only if the number of AO basis functions is equal to the number of MO basis functions, which means if you work with cartesian functions. 


.. seealso::

 You can check with :command:`qp create_ezfio -h` for the option to create an |EZFIO| with cartesian basis functions

In the file :file:`qp2/plugins/tuto_plugins/tuto_I/print_traces_on_e.irp.f` you will find an example of executable that prints out the various providers. 
Copy these two files in your plugin and recompile to execute it.

Execute the program print_traces_on_e and check for the results with 

.. code:: bash

 qp run print_traces_on_e

The code in  :file:`qp2/plugins/tuto_plugins/tuto_I/print_traces_on_e.irp.f` should be easy to read, I let the reader interpret it.
