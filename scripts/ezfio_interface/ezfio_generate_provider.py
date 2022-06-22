#!/usr/bin/env python3

__author__ = "Applencourt PEP8"
__date__ = "jeudi 26 mars 2015, 12:49:35 (UTC+0100)"

"""
Creates the provider of a variable that has to be
fetched from the EZFIO file.
"""

import sys

def gen_ezfio_provider_disk_access(name_ref, name):
    data = """
  BEGIN_PROVIDER [ logical, read_{name} ]
 &BEGIN_PROVIDER [ logical, write_{name} ]

  BEGIN_DOC
  ! One level of abstraction for {name}
  END_DOC

   if ({name_ref}.EQ.'Read') then
     read_{name} =  .True.
     write_{name} = .False.
   else if  ({name_ref}.EQ.'Write') then
     read_{name} = .False.
     write_{name} =  .True.
   else if ({name_ref}.EQ.'None') then
     read_{name} = .False.
     write_{name} = .False.
   else
     print *, '{name_ref} has a bad type'
     stop 1
   endif

 END_PROVIDER
"""
    return data.format(name=name,name_ref=name_ref)

class EZFIO_Provider(object):

    data = """
BEGIN_PROVIDER [ %(type)s, %(name)s %(size)s ]
  implicit none
  BEGIN_DOC
! %(doc)s
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  if (mpi_master) then
    %(test_null_size)s
    call ezfio_has_%(ezfio_dir)s_%(ezfio_name)s(has)
    if (has) then
      write(6,'(A)') '.. >>>>> [ IO READ: %(name)s ] <<<<< ..'
      call ezfio_get_%(ezfio_dir)s_%(ezfio_name)s(%(name)s)
    else
      print *, '%(ezfio_dir)s/%(ezfio_name)s not found in EZFIO file'
      stop 1
    endif
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( %(name)s, %(size_mpi)s, %(type_mpi)s, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read %(name)s with MPI'
    endif
  IRP_ENDIF
%(write)s
END_PROVIDER
""".strip()

    write_correspondance = {"integer": "write_int",
                            "logical": "write_bool",
                            "double precision": "write_double"}

    mpi_correspondance   = {"integer": "MPI_INTEGER",
                            "integer*8": "MPI_INTEGER8",
                            "character*(32)": "MPI_CHARACTER",
                            "character*(64)": "MPI_CHARACTER",
                            "character*(256)": "MPI_CHARACTER",
                            "logical": "MPI_LOGICAL",
                            "double precision": "MPI_DOUBLE_PRECISION"}

    def __init__(self):
        self.values = "type doc name ezfio_dir ezfio_name write output".split()
        for v in self.values:
            exec("self.{0} = None".format(v))

    def __repr__(self):
        self.set_write()
        self.set_test_null_size()
        for v in self.values:
            if not v:
                msg = "Error : %s is not set in EZFIO.cfg" % (v)
                print(msg, file=sys.stderr)
                sys.exit(1)
        if "size" not in self.__dict__:
            self.__dict__["size"] = ""

        return self.data % self.__dict__

    def set_test_null_size(self):
        if "size" not in self.__dict__:
            self.__dict__["size"] = ""
        if self.size != "":
            self.test_null_size = "if (size(%s) == 0) return\n" % ( self.name )
        else:
            self.test_null_size = ""

    def set_write(self):
        output = self.output
        name = self.name
        l_write = ["",
                   "  call write_time(%(output)s)",
                   ""]

        self.write = "\n".join(l_write) % locals()
        self.type_mpi = self.mpi_correspondance[self.type]
        if "size" in self.__dict__:
            return
        else:
            if self.type in self.write_correspondance:
                write = self.write_correspondance[self.type]

                l_write = ["",
                          "  call write_time(%(output)s)",
                          "  call %(write)s(%(output)s, %(name)s, &",
                          "     '%(name)s')",
                          ""]

    def set_type(self, t):
        self.type = t.lower()

    def set_doc(self, t):
        self.doc = t.strip().replace('\n', '\n! ')

    def set_name(self, t):
        self.name = t

    def set_ezfio_dir(self, t):
        self.ezfio_dir = t.lower()

    def set_ezfio_name(self, t):
        self.ezfio_name = t.lower()

    def set_output(self, t):
        self.output = t

    def set_size(self, t):
        self.size_mpi = t.replace(',',')*(').replace('0:','1+')
        if (self.type == "character*(32)"):
            self.size_mpi += "*32"
        if t != "1":
            self.size = ", " + t
        else:
            self.size = ""

def test_module():
    T = EZFIO_Provider()
    T.set_type("double precision")
    T.set_name("thresh_SCF")
    T.set_doc("Threshold on the convergence of the Hartree Fock energy")
    T.set_ezfio_dir("Hartree_Fock")
    T.set_ezfio_name("thresh_SCF")
    T.set_output("output_Hartree_Fock")
    print(T)

if __name__ == '__main__':
    test_module()
