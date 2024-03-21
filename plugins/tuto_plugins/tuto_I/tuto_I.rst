======================================
Tutorial for plugin I: One-e integrals
======================================

!!! Requirements: 
 a) you know how to create an EZFIO file and run calculations with QP 
 (check the tuto: `<https://quantumpackage.github.io/qp2/post/hartree-fock/>`),
 b) you have an EZFIO file in the sto-3g from the file H2.xyz in plugins/tuto_plugins,  
    and you have run an HF calculation giving an energy of -1.116759 a.u.,
 c) you made an qp set_file YOUR_EZFIO_FILE_FOR_H2 in order to be, 
 d) you have READ the ../README.rst file to HAVE THE VOCABULARY.  

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
Therefore to create the plugin, we do 

$ qp plugins create -n plugin_I -r tuto_plugins
Then to an "ls" in qp2/plugins/tuto_plugins/ 
and you will find a directory called "plugin_I". 
In that directory you will find: 
 i)  a "NEED" file that will eventually contain all the other modules/plugins needed by our "plugin_I"
 ii) a "README.rst" file that you can AND SHOULD modify in order to document what is doing the plugin. 
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
irpman ao_one_e_integral
irpman mo_one_e_integral
in order to get some information on where they are created, and many more information. 
We will modify the executable such that it prints out the integrals. 


IV) Printing out the one-electron integrals
--------------------------------------------
We will create a program that will print out the one-electron integrals on the AO and MO basis.
You can then copy the file "print_one_e_h.irp.f" in your plugin. 
In the file you will see that we simply browse the two arrays "ao_one_e_integrals" and "mo_one_e_integrals", which are global variables (providers) and we browse them until either "ao_num" or "mo_num" which are also providers representing the number of AOs or MOs. 
You can check these variables with irpman ! 
If you recompile using "ninja" as before, and another executable has been created "print_one_e_h". 
Then, you can run the program on the ezfio file by doing 
qp run print_one_e_h 
and will print out the data you need :) 

