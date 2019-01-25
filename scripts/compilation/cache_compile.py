#!/usr/bin/env python2
"""
Save the .o from a .f90
and is the .o is asked a second time, retur it
Take in argv command like:
     ifort -g  -openmp -I IRPF90_temp/Ezfio_files/ -c IRPF90_temp/Integrals_Monoelec/kin_ao_ints.irp.module.F90 -o IRPF90_temp/Integrals_Monoelec/kin_ao_ints.irp.module.o
"""

import os
import sys
import hashlib
import re
import shutil
import subprocess

r = re.compile(ur'-c\s+(\S+\.[fF]90)\s+-o\s+(\S+\.o)')
p = re.compile(ur'-I IRPF90_temp/\S*\s+')
mod = re.compile(ur'module\s+(?P<mod>\S+).+end\s?module\s+(?P=mod)?',
                 re.MULTILINE | re.IGNORECASE)

tmpdir_root = os.environ.get("TMPDIR", failobj="/dev/shm")
TMPDIR = os.path.join(tmpdir_root, os.environ["USER"], "qp_compiler")


def return_filename_to_cache(command):
    """
    For a irp_command:
        ifort -g  -openmp -I IRPF90_temp/Ezfio_files/ -c IRPF90_temp/Integrals_Monoelec/kin_ao_ints.irp.module.F90 -o IRPF90_temp/Integrals_Monoelec/kin_ao_ints.irp.module.o

    Return the *.F90 and the *.o
    """
    command_clean = p.sub('', command)
    match = r.search(command_clean)

    input = match.group(1)
    output = match.group(2)

    return (input, output)


def get_hash_key(command, input_data):
    """
    Return the hash of command + input_data
    """
    m = hashlib.md5()
    m.update(command)
    m.update(input_data)

    # Md5 Key containing command + content of Fread
    return m.hexdigest()


def run_and_save_the_data(command, path_output, path_key, is_mod):

    # Compile the file -> .o
    process = subprocess.Popen(command, shell=True)

    if process.wait() != 0:
        sys.exit(1)
    elif not is_mod:
        try:
            shutil.copyfile(path_output, path_key)
        except:
            pass


def cache_utility(command):
    # Create temp directory

    try:
        os.makedirs(TMPDIR)
    except OSError:
        pass

    # Get the filename of the input.f.90
    # and the otput .o

    try:
        (path_input, path_output) = return_filename_to_cache(command)
    except:
        # Canot parse the arg of command
        raise OSError

    try:
        with open(path_input, 'r') as f:
            input_data = f.read()

        # Get the hash
        key = get_hash_key(command, input_data)
        path_key = os.path.join(TMPDIR, key)

        # Try to return the content of the .o file
        try:
            shutil.copyfile(path_key, path_output)
        except IOError:
            is_mod = mod.search(input_data.replace('\n', ' '))
            run_and_save_the_data(command, path_output, path_key, is_mod)
    except:
        raise

if __name__ == '__main__':

    line = sys.argv[1:]
    command = " ".join(line)

    try:
        cache_utility(command)
    except OSError:
        process = subprocess.Popen(command, shell=True)
