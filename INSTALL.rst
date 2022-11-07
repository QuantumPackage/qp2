============
Installation
============

|qp| can be downloaded on GitHub as an `archive
<https://github.com/QuantumPackage/qp2/releases>`_ or as a `git
repository <https://github.com/QuantumPackage/qp2>`_.

.. code:: bash

  git clone https://github.com/QuantumPackage/qp2

QP has a long list of dependencies:

- GNU Linux OS
- Fortran compiler : GNU Fortran, Intel Fortran or IBM XL Fortran
- `GNU make`_
- `Autoconf`_
- `Python`_ > 3.7
- |IRPF90| : Fortran code generator
- |EZFIO| : Easy Fortran Input/Output library generator
- |BLAS| and |LAPACK|
- `Zlib`_
- `GNU Patch`_
- |ZeroMQ| : networking library
- `GMP <https://gmplib.org/>`_ : Gnu Multiple Precision Arithmetic Library
- |OCaml| compiler with |OPAM| package manager
- |Ninja| : a parallel build system
- |pkg-config| : a tool which returns information about installed libraries

so we don't recommend the installation of all the dependencies by hand.

Preferred installation method: Guix
===================================

|qp| is very easy to install on systems that have
access to `Guix <https://guix.gnu.org/en>`_. To install Guix
on your system, follow the `Guix installation instructions
<https://guix.gnu.org/manual/en/html_node/Binary-Installation.html>`_.

One Guix is installed, all the dependencies required for compiling
and running |qp| can be installed using:

.. code:: bash

  guix package -L guix -m guix/manifest.scm 

in the :file:`qp2` root directory.

Installation on a x86_64 system without network access
======================================================

Fetch the :file:`qp2-environment.tar.gz` archive containing all the
dependencies `here <http://localhost>`_.  Unpack the archive in the
:file:`guix` directory:

.. code:: bash

   $ cd qp2/guix
   $ wget qp2-environment.tar.gz
   $ tar -zxf qp2-environment.tar.gz
   $ ls
   gnu  manifest.scm  qp2-environment.tar.gz  qp2.scm

That's it!

Compilation
===========

Before anything, go into your :file:`quantum_package` directory and run

.. code:: bash

   ./configure


This script will create the :file:`quantum_package.rc` bash script, which
sets all the environment variables required for the normal operation of
|qp|. It will also initialize the git submodules that are
required, and tell you which external dependencies are missing and need to be
installed. The required dependencies are located in the
`external/qp2-dependencies` directory, such that once |qp| is configured the
internet connection is not needed any more.

When all dependencies have been installed, (the :command:`configure` will
inform you what is missing) source the :file:`quantum_package.rc` in order to
load all environment variables and compile |QP|.

Now all the requirements are met, you can compile the programs using

.. code:: bash

   make



When all the dependencies have been installed, go into the :file:`config`
directory, and copy the configuration file that corresponds to your
architecture. Modify it if needed, and run :command:`configure` with
:option:`configure -c`.

.. code:: bash

   cp ./config/gfortran.example config/gfortran_avx.cfg
   ./configure -c config/gfortran_avx.cfg


.. note::

   The ``popcnt`` instruction accelerates *a lot* the programs, so the
   SSE4.2, AVX or AVX2 instruction sets should be enabled in the
   configuration file if possible.


Help for installing external dependencies
=========================================

Using the :command:`configure` executable
-----------------------------------------

The :command:`configure` executable can help you in installing the minimal dependencies you will need to compile the |QP|.
The command is to be used as follows:

.. code:: bash

   ./configure -i <package>

The following packages are supported by the :command:`configure` installer:

* ninja
* zeromq
* f77zmq
* gmp
* ocaml  (:math:`\approx` 5 minutes)
* docopt
* resultsFile
* bats
* zlib

Example:

.. code:: bash

   ./configure -i ninja



If the :command:`configure` executable fails to install a specific dependency
-----------------------------------------------------------------------------

If the :command:`configure` executable does not succeed in installing a specific
dependency, you should try to install the dependency on your system by yourself.

Before doing anything below, try to install the packages with your package manager
(:command:`apt`, :command:`yum`, etc).


Ninja
-----

*Ninja* is a build system (like GNU make), with a focus on speed.

* Download the latest binary version of Ninja
  here : `<https://github.com/ninja-build/ninja/releases/latest>`_

* Unzip the ninja-linux.zip file, and move the ninja binary into
  the :file:`${QP_ROOT}/bin` directory.



IRPF90
------

*IRPF90* is a Fortran code generator for programming using the Implicit Reference
to Parameters (IRP) method.

If you have *pip* for Python2, you can do

.. code:: bash

   python3 -m pip install --user irpf90

Otherwise,

* Download the latest version of IRPF90
  here : `<https://gitlab.com/scemama/irpf90/-/archive/v1.7.2/irpf90-v1.7.2.tar.gz>`_ and move
  the downloaded archive in the :file:`${QP_ROOT}/external` directory

* Extract the archive and go into the :file:`irpf90-*` directory to run
  :command:`make`

.. note::

    The :envvar:`IRPF90_PATH` variable may need to be updated in the configuration
    file :file:`${QP_ROOT}/etc/irpf90.rc`.



ZeroMQ and its Fortran binding
------------------------------

*ZeroMQ* is a high-performance asynchronous messaging library.

* Download the latest stable version of ZeroMQ
  here : `<https://github.com/zeromq/libzmq/releases/latest>`_ and move the
  downloaded archive in the :file:`${QP_ROOT}/external` directory

* Extract the archive, go into the :file:`zeromq-*` directory and run
  the following commands

.. code:: bash

   ./configure --prefix="${QP_ROOT}" --without-libsodium
   make
   make install


* Download the Fortran binding
  here : `<https://github.com/zeromq/f77_zmq/releases/latest>`_ and move
  the downloaded archive in the :file:`${QP_ROOT}/external` directory

* Extract the archive, go into the :file:`f77_zmq-*` directory and run
  the following commands

.. code:: bash

   export ZMQ_H=${QP_ROOT}/include/zmq.h
   make
   cp libf77zmq.a ${QP_ROOT}/lib
   cp libf77zmq.so ${QP_ROOT}/lib


* Copy the :file:`f77_zmq_free.h` file in the ``ZMQ`` module as follows:

.. code:: bash

   cp f77_zmq_free.h ${QP_ROOT}/src/zmq/f77_zmq.h


Zlib
----

*Zlib* is the compression library used by *gzip*.

* Download the latest version of Zlib here:
  `<https://www.zlib.net/zlib-1.2.11.tar.gz>`_
  and move it in the :file:`${QP_ROOT}/external` directory

* Extract the archive, go into the :file:`zlib-*` directory and run
  the following commands


.. code:: bash

   ./configure --prefix=${QP_ROOT}
   make
   make install

With Debian or Ubuntu, you can use

.. code:: bash

   sudo apt install zlib1g-dev

GMP
---

GMP is the GNU Multiple Precision Arithmetic Library.

* Download the latest version of gmp here:
  `<ftp://ftp.gnu.org/gnu/gmp/gmp-6.1.2.tar.bz2>`_
  and move it in the :file:`${QP_ROOT}/external` directory

* Extract the archive, go into the :file:`gmp-*` directory and run
  the following commands

.. code:: bash

   ./configure --prefix=${QP_ROOT}
   make
   make install

With Debian or Ubuntu, you can use

.. code:: bash

   sudo apt install libgmp-dev




OCaml
-----

*OCaml* is a general purpose programming language with an emphasis on expressiveness and safety.

* The following packages are required (Debian or Ubuntu):

  .. code:: bash

    sudo apt install libncurses5-dev pkg-config libgmp3-dev m4


* Download the installer of the OPAM package manager here :
  `<https://raw.githubusercontent.com/ocaml/opam/master/shell/install.sh>`_
  and move it in the :file:`${QP_ROOT}/external` directory

* If you use OCaml only with |qp|, you can install the OPAM directory
  containing the compiler and all the installed libraries in the
  :file:`${QP_ROOT}/external` directory as

  .. code:: bash

     export OPAMROOT=${QP_ROOT}/external/opam


* Run the installer

  .. code:: bash

     echo ${QP_ROOT}/bin
     ${QP_ROOT}/external/opam_installer.sh --no-backup --fresh

  The :command:`opam` command can be installed in the :file:`${QP_ROOT}/bin`
  directory. To do this, take the output of ``echo ${QP_ROOT}/bin`` and
  use it as an answer to where :command:`opam` should be installed.


* Install the OCaml compiler

  .. code:: bash

      opam init --comp=4.11.1
      eval `${QP_ROOT}/bin/opam env`

  If the installation fails because of bwrap, you can initialize opam using:

  .. code:: bash

      opam init --disable-sandboxing --comp=4.11.1
      eval `${QP_ROOT}/bin/opam env`

* Install the required external OCaml libraries

  .. code:: bash

      opam install ocamlbuild cryptokit zmq sexplib ppx_sexp_conv ppx_deriving getopt


Docopt
------

*Docopt* is a Python package defining a command-line interface description language.

If you have *pip* for Python3, you can do

.. code:: bash

   python3 -m pip install --user docopt

Otherwise,

* Download the archive here : `<https://github.com/docopt/docopt/releases/tag/0.6.2>`_

* Extract the archive

* Copy :file:`docopt-0.6.2/docopt.py` in the :file:`${QP_ROOT}/scripts` directory


resultsFile
-----------

*resultsFile* is a Python package to extract data from output files of quantum chemistry
codes.

If you have *pip* for Python3, you can do

.. code:: bash

   python3 -m pip install --user resultsFile




Installation of dependencies via a Conda environment
====================================================

.. code:: bash

  conda env create -f qp2.yml
￼￼


