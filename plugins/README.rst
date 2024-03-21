==============================
Tutorial for creating a plugin
==============================

Introduction: what is a plugin, and what this tuto will be about ?
============================================================
The QP is split into two kinds of routines/global variables (i.e. providers): 
 i)  the core modules locatedin qp2/src/, which contains all the bulk of a quantum chemistry software (integrals, matrix elements between Slater determinants, linear algebra routines, DFT stuffs etc..)
 ii) the plugins which are external stuffs connected to the qp2/src/ stuffs.
 
More precisely, a plugin of the QP is a directory where you can create routines, 
providers and executables that use all the global variables/functions/routines already created 
in the modules ofqp2/src or in other plugins. 

Instead of giving a theoretical lecture on what is a plugin, 
we will go through a series of examples that allow you to do the following thing: 
 I)   print out one- and two-electron integrals on the AO/MO basis,  
      creates two providers which manipulate these objects, 
      print out these providers, 
 II)  browse the Slater determinants stored in the EZFIO wave function and compute their matrix elements, 
 III) build the Hamiltonian matrix and diagonalize it either with Lapck or Davidson,
 IV)  print out the one- and two-electron rdms, 
 V)   obtain the AOs and MOs on the DFT grid, together with the density,

This tuto is as follows: 
 i)  you READ THIS FILE UNTIL THE END in order to get the big picture and vocabulary, 
 ii) you go to the directory qp2/plugins/tuto_plugins/ and you will find detailed tuto there for each of the 5 examples. 

Creating a plugin: the basic
----------------------------
The first thing to do is to be in the QPSH mode: you execute the qp2/bin/qpsh script that essentially loads all 
the environement variables and allows for the completion of command lines in bash (that is an AMAZING feature :) 

Then, you need to known where you want to create your plugin, and what is the name of the plugin. 
!!!! WARINING: The plugins are NECESSARILY located in qp2/plugins/ !!!!
Ex: If you want to create a plugin named "my_fancy_plugin" in the directory plugins/plugins_test/, 
this goes with the command 
qp plugins create -n my_fancy_plugin -r plugins_test/

Then, to create plugin of your dreams, the two questions you need to answer are the following: 
a) What do I need to compute what I want, which means what are the objects that I need ?
   There are two kind of objects:
   + the routines/functions 
     Ex: Linear algebra routines, integration routines etc ...
   + the global variables which are called the PROVIDERS
     Ex: one-electron integrals, Slater determinants, density matrices etc ...
b) Where do I find these objects ? 
   The objects (routines/functions/providers) are necessarily created in other modules/plugins 
     Ex: the routine "lapack_diagd" (which diagonalises a real hermitian matrix) is located in the file 
         qp2/src/utils/linear_algebra.irp.f
         therefore it "belongs" to the module "utils" 
       : the routine "ao_to_mo" (which converts a given matrix A from the AO basis to the MO basis) is located in the file
         qp2/src/mo_one_e_ints/ao_to_mo.irp.f
         therefore it "belongs" to the module "mo_one_e_ints"
       : the provider "ao_one_e_integrals" (which is the integrals of one-body part of H on the AO basis) is located in the file 
         qp2/src/mo_one_e_ints/ao_to_mo.irp.f 
         therefore it belongs to the module "mo_one_e_ints" 
       : the provider "one_e_dm_mo_beta_average" (which is the state average beta density matrix on the MO basis) is located in the file 
         qp2/src/determinants/density_matrix.irp.f 
         therefore it belongs to the module "determinants"

To import all the variables that you need, you just need to write the name of the plugins in the file "NEED"
Ex: to import all the variables/routines of the module "utils", "determinants" and "mo_one_e_ints" you will have the following NEED file:
utils
determinants
mo_one_e_ints

TIPS 
----
There are many many routines/providers in the core modules of QP. Nevertheless, as everything is coded with the IRPF90, you can use the following amazing tools: irpman
irpman can be used in command line in bash to obtain all the info on a routine or variable ! 
Ex: execute the following command line : 
irpman ao_one_e_integrals
Then it appears all the information you want on ao_one_e_integrals, including where it is created, the type, dimension if it is an array, what providers it needs to be built, and what providers need this provider. 


