===========
ezfio_files
===========

Access to the |EZFIO| database and to the basic execution
information.

The provider :c:data:`ezfio_filename` is the name of the |EZFIO|
directory: it is read from the :envvar:`QPACKAGE_INPUT` environment
variable, or from the first argument of the command line.

The other providers and subroutines of this module are:

* the process and thread counts: :c:data:`nproc`,
  :c:data:`nproc_max`, :c:data:`nthreads_pt2`,
  :c:data:`nthreads_davidson`
* the file locking: :c:func:`lock_io`, :c:func:`unlock_io`,
  :c:data:`file_lock`, and :c:func:`getunitandopen` to open a file
  with a free unit number
* the output helpers: :c:func:`write_time`, :c:func:`wall_time`,
  :c:data:`output_wall_time_0`, :c:data:`output_cpu_time_0`, and
  :c:func:`write_bool`, :c:func:`write_int`, :c:func:`write_double`
* the clean stopping of the calculation: :c:data:`qp_stop`,
  :c:data:`qp_stop_variable`, :c:data:`qp_stop_filename` and
  :c:data:`qp_kill_filename`
