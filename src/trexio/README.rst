======
trexio
======

Reading and writing of TREXIO files
(see https://github.com/trex-coe/trexio).

The :ref:`export_trexio` program exports the wave function in TREXIO
format, to the file named by :c:data:`trexio_file`, with the back-end
selected by :c:data:`backend` (0 for HDF5, 1 for text). The
`export_*` providers gather the basis, the |MOs|, the one- and
two-electron integrals (including the Cholesky decompositions), the
determinants and the density matrices of the wave function for the
export.

The program :ref:`import_trexio_integrals` imports the one- and
two-electron integrals from a TREXIO file into the |EZFIO| database,
and :ref:`import_trexio_determinants` imports determinants from a
TREXIO file.
