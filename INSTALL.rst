============
Installation
============

|qp| can be downloaded on GitHub as an `archive
<https://github.com/QuantumPackage/qp2/releases>`_ or as a `git
repository <https://github.com/QuantumPackage/qp2>`_.

.. code:: bash

  git clone https://github.com/QuantumPackage/qp2


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


.. note::

   No Conda environment is needed. The easiest way to install the external
   dependencies is to let :command:`configure` do it for you, see
   "Help for installing external dependencies" below.


Requirements
============

- Linux OS
- Fortran compiler : GNU Fortran, Intel Fortran, NVIDIA Fortran or LLVM flang
- `GNU make`_
- `Python`_ > 3.7
- |IRPF90| : Fortran code generator
- |EZFIO| : Easy Fortran Input/Output library generator
- |BLAS| and |LAPACK|
- cmake
- TREXIO : library to import and export TREXIO format files
- `Zlib`_
- |ZeroMQ| : networking library
- `GMP <https://gmplib.org/>`_ : Gnu Multiple Precision Arithmetic Library
- Rust toolchain (:command:`cargo`, edition 2024) : https://rustup.rs
- |Ninja| : a parallel build system
- |pkg-config| : a tool which returns information about installed libraries


When all the dependencies have been installed, pick in the :file:`config`
directory the configuration file that corresponds to your machine, copy
it to a new file, modify the copy, and run :command:`configure` with the
:file:`-c` option. The configuration files, and the choice for each
machine, are detailed on the :ref:`configure` page.

.. code:: bash

   cp config/gfortran.cfg config/gfortran_mybox.cfg
   ./configure -c config/gfortran_mybox.cfg


Help for installing external dependencies
=========================================

Using the :command:`configure` executable
-----------------------------------------

The :command:`configure` executable can help you in installing the minimal dependencies you will need to compile the |QP|.
The command is to be used as follows:

.. code:: bash

   ./configure -i <package>

The following packages are supported by the :command:`configure` installer:

* cmake
* zlib
* ninja
* zeromq
* f77zmq
* gmp
* docopt
* resultsFile
* bats
* trexio
* armadillo
* openorb
* simple_gpu
* bse (adds the Basis Set Exchange CLI)

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

If you have *pip* for Python3, you can do

.. code:: bash

   python3 -m pip install --user irpf90

Otherwise,

* Download the latest version of IRPF90 (v2.2.0)
  here : `<https://gitlab.com/scemama/irpf90/-/archive/v2.2.0/irpf90-v2.2.0.tar.gz>`_ and move
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




Rust
----

*Rust* is the language of the command-line tools of |qp| (:command:`qp
create_ezfio`, :command:`qp edit`, :command:`qp run`, :command:`qp
set_mo_class`, :command:`qp tunnel`, :command:`qp print_basis` and
:command:`element_create_db`).  The Fortran/IRPF90 modules are built with
the Fortran compiler, as before; only the tools need Rust.

* Install `rustup <https://rustup.rs>`_:

  .. code:: bash

     curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

  The version of the compiler is pinned by
  :file:`${QP_ROOT}/rust/rust-toolchain.toml`; :command:`cargo` downloads
  the pinned toolchain automatically on the first build.  To install it
  explicitly:

  .. code:: bash

     rustup toolchain install 1.92.0

  The toolchain needs network access the first time it is used.  On a
  machine without network access, install the toolchain before copying
  |qp|.

* Check that :command:`cargo` is visible in the :envvar:`PATH`.  The
  :command:`configure` script stops with an error message if it is not.


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




