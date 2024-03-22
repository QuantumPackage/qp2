=====================================================================
Tutorial for plugin I: One-e integrals (duration: 20 minutes at most)
=====================================================================

Requirements 
------------
 a) You know how to create an EZFIO file and run calculations with QP 
 (check the tuto: `<https://quantumpackage.github.io/qp2/post/hartree-fock/>`),
 b) You have an EZFIO file with MOs created (with the 'scf' executable for instance). 
    As we are going to print out some integrals, don't take a too large system/basis (Ex: H2, cc-pVDZ is ok :)
 c) You made an qp set_file YOUR_EZFIO_FILE_FOR_H2 in order to work on that ezfio folder, 
 d) You have READ the ../README.rst file to HAVE THE VOCABULARY.  

Our goals:
----------
We want to create a plugin to do the following things: 
 a) print out one- and two-electron integrals on the AO/MO basis,  
 b) creates two providers which manipulate these objects, 
 c) print out these providers, 

I) Starting: creating the plugin
--------------------------------
We will go step-by-step through these plugins. 

The name of the plugin will be "plugin_I", and its location is in "tuto_plugins". 
Therefore to create the plugin, we do: 

qp plugins create -n plugin_I -r tuto_plugins

Then do an "ls" in qp2/plugins/tuto_plugins/ and you will find a directory called "plugin_I". 
In that directory you will find: 
 i)   a "NEED" file that will eventually contain all the other modules/plugins needed by our "plugin_I"
 ii)  a "README.rst" file that you can AND SHOULD modify in order to document what is doing the plugin. 
 iii) a "plugin_I.irp.f" file that is a program to be compiled and just printing "Hello world" 

II) Specifying the dependencies
-------------------------------
The next step is to know what are the other modules/plugins that we need to do what we want. 
We need here 
 a) the one-electron integrals on the AO basis, which are computed in qp2/src/ao_one_e_ints/
 b) the one-electron integrals on the MO basis, which are computed in qp2/src/mo_one_e_ints/
 c) the two-electron integrals on the AO basis, which are computed in qp2/src/ao_two_e_ints/
 d) the two-electron integrals on the MO basis, which are computed in qp2/src/mo_two_e_ints/

Therefore, we will need the following four modules:
a) ao_one_e_ints
b) mo_one_e_ints
c) ao_two_e_ints
d) mo_two_e_ints

You can then create the following "NEED" file by executing the following command 
$ cat <<EOF > NEED
ao_one_e_ints
mo_one_e_ints 
ao_two_e_ints
mo_two_e_ints
EOF

II) Installing the plugin 
-------------------------
Now that we have specified the various depenencies we need now to INSTALL the plugin, which means to create the equivalent of a Makefile for the compilation. 
To do it we simply do 
$ qp plugins install plugin_I

III) Compiling the void plugin 
------------------------------
It is customary to compile first your "void" plugin, void in the sense that it does not contain anything else than the program printing "Hello world". 
To do so, just go in the plugin and execute the following command
$ ninja 
It does a lot of stuffs, but it must conclude with something like 
"
make: Leaving directory 'SOME_PATH_TOWARD_YOUR_QP2_DIRECTORY/qp2/ocaml'
"

Since that it has compiled, an executable "plugin_I" has been created. 
Also, if you make "ls" in the "plugin_I" you will notice that many symbolink links have been created, and among which the four modules that you included in the NEED file. 
All the other modules (Ex:"ao_basis", "utils") are here because they are need by some of the  four modules that you need. 
The variables that we need are 
ao_one_e_integrals 
mo_one_e_integrals
You can check them with 
irpman ao_one_e_integrals
irpman mo_one_e_integrals
in order to get some information on where they are created, and many more information. 
We will modify the executable such that it prints out the integrals. 


IV) Printing out the one-electron integrals
--------------------------------------------
We will create a program that will print out the one-electron integrals on the AO and MO basis.
You can then copy the file "print_one_e_h.irp.f" located in "plugins/tuto_plugins/tuto_I" in your plugin. 
In the file you will see that we simply browse the two arrays "ao_one_e_integrals" and "mo_one_e_integrals", which are global variables (providers) and we browse them until either "ao_num" or "mo_num" which are also providers representing the number of AOs or MOs. 
You can check these variables with irpman ! 
If you recompile using "ninja" as before, and another executable has been created "print_one_e_h". 
Then, you can run the program on the ezfio file by doing 
qp run print_one_e_h 
and will print out the data you need :) 

By the way, as the file "plugin_I.irp.f" contains nothing but a "Hello world" print, you can simply remove it if you want. 
V) Printing out the two-electron integrals
------------------------------------------
We will now create a file that prints out the two-electron integrals in the AO and MO basis.
These can be accessed with the following subroutines :
+) get_ao_two_e_integral for the AO basis 
+) get_two_e_integral for the MO basis 
check them with irpman !
To print the two-electron integrals, you can copy the file "print_two_e_h.irp.f" in your plugin and recompile. 
Then just run the program 
qp run print_two_e_h
and it will print all the things you want :)

VI) Creating new providers and a program to print them
------------------------------------------------------
We will now create new providers that manipulates the objects that we just printed.
As an example, we will compute the trace of the one electron integrals in the AO and MO basis. 
In the file "traces_one_e.irp.f" you will find the several new providers among which 
 a) trace_mo_one_e_ints         : simply the sum of the diagonal matrix element of the one-electron integrals
 b) trace_ao_one_e_ints         : the corresponding trace on the AO basis : Sum(m,n) S^{-1}_{mn} h_{mn}
 c) trace_ao_one_e_ints_from_mo : the trace on the AO basis with the integrals obtained first from the MO basis 
As explained in these files, "trace_mo_one_e_ints" is equal to "trace_ao_one_e_ints" only if the number of AO basis functions is equal to the number of MO basis functions, which means if you work with cartesian functions. 
(You can check with "qp create_ezfio -h" for the option to create an EZFIO with cartesian basis functions)

In the file "print_traces_on_e.irp.f" you will find an example of executable that prints out the various providers. 
Copy these two files in your plugin and recompile to execute it.

Execute the program print_traces_on_e and check for the results !
