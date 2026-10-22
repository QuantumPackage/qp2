#!/usr/bin/env python3


import os
import sys
import configparser

from module_handler import get_binaries


def generate_modules(abs_module, entities, local_entities=None):
  """Generates the doc for modules. local_entities is the entity dict read
     from the module's own IRPF90_man (built modules); when it is None the
     module is provider-only and its entities are looked up in the global
     `entities` dict under its module name."""
  MODULE = os.path.split(abs_module)[-1]
  module = MODULE.lower()
  if module == "dummy":
    return

  readme_file = os.path.join(abs_module, 'README.rst')
  if os.path.exists(readme_file):
    with open(readme_file, 'r') as f:
      readme = f.read()
  else:
    readme = ("%s\n%s\n\nNo README.rst in this module; the providers and\n"
              "subroutines below are documented from the irpf90 man "
              "pages." % (module, "="*len(module)))
  rst = [
    ".. _module_%s:"%(module), "", 
    ".. program:: %s"%(module), "", 
    ".. default-role:: option", "", 
    readme, "", 
  ]

  EZFIO = os.path.join(abs_module,'EZFIO.cfg')
  if os.path.exists(EZFIO):
    rst += ["", "EZFIO parameters", "----------------", ""]
    config_file = configparser.ConfigParser()
    with open(EZFIO, 'r') as f:
#        config_file.readfp(f)
        config_file.read_file(f)
        for section in config_file.sections():
            doc = config_file.get(section, "doc")
            doc = "    " + doc.replace("\n", "\n\n    ")+"\n"
            try:
                default = config_file.get(section, "default")
                default = "    " + "Default: %s\n"%default
            except:
                default = ""
            rst += [".. option:: %s\n"%(section), doc, default]

  src = local_entities if local_entities is not None \
      else {k: e for k, e in entities.items()
            if e.get("module", "").lower() == module.lower()}
  providers = []
  subroutines = {}
  for k in sorted(src.keys()):
    e = src[k]
    if e.get("module", "").lower() == module.lower():
        if "/" not in e["file"] and e["file"] != "ezfio_interface.irp.f":
            if e["type"] == 's':
                subroutines[e["name"]] = e
            elif e["type"] == 'p':
                providers.append(e)

  progs = get_binaries(abs_module)
  binaries = [os.path.basename(p[0]) for p in progs]
  prog_names = set(p[1].lower() for p in progs)

  if progs:
    rst += ["", "Programs", "--------", ""]
    for binary, prog in progs:
        try:
            subroutines[prog.lower()]
        except KeyError:
            print("Error: The program %s (binary %s) in %s has no man page; \
you probably did not run ninja at the root."%
                (prog, os.path.basename(binary), abs_module))
            sys.exit(1)
        # The program page is anchored by the binary (file) name, which is
        # unique. The program statement name need not be: two tools
        # programs both declare 'program print_hamiltonian'.
        rst += [" * :ref:`%s`"%(os.path.basename(binary))]

  if providers:
    rst += ["", "Providers", "---------", ""]
    for p in providers:
        rst += [p["rst"]]

  if subroutines:
    rst += [ "", "Subroutines / functions", "-----------------------", "" ]
    for p in sorted(subroutines.keys()):
        p = subroutines[p]
        if p["name"] in binaries or p["name"] in prog_names:
           continue
        rst += [p["rst"]]

  rst_file = os.path.join('modules', module+".rst")
  with open(rst_file,'w') as f: 
      f.write(" \n".join(rst))

  for binary, prog in progs:
    p = subroutines.get(prog.lower())
    if p is None:
      continue
    b = os.path.basename(binary)
    rst = [".. _%s:"%(b), "", 
           ".. program:: %s"%(b), "", 
           "="*len(b), b, "="*len(b), "", ""]
    rst += [line[3:] for line in p["rst"].splitlines()[8:]]
    rst_file = os.path.join('programs', b+".rst")
    with open(rst_file,'w') as f: 
        f.write(" \n".join(rst))



def generate_providers(abs_module):
  """ Reads the IRPF90_man pages and returns a dict of dicts describing the
      providers.
  """
  MODULE = os.path.split(abs_module)[-1]
  module = MODULE.lower()
  if module == "dummy":
    return

  files    = {}
  entities = {}
  mandir = os.path.join(abs_module, 'IRPF90_man') 
  if not os.path.exists(mandir):
    return {}

  for f in os.listdir(mandir):
        if f.endswith('.rst'):
          continue
        filename = os.path.join(mandir, f)
        if f not in files:
            files[f] = 0
            name = f.split('.')[0] 
            with open(os.path.join(mandir, name+".rst"), 'r') as g:
              rst = g.read()
            with open(filename, 'r') as f:
                state = 0
                entity = {"decl": [], "doc": [] ,
                    "name": name , "module": module, "rst":rst}
                text=f.read()
                text_old = None
                while text_old != text:
                    text_old = text
                    text = text.replace("$"," :math:`",1).replace("$","` ",1)
                for line in text.splitlines():
                    line = line.rstrip()
                    if line.startswith(".SH Declaration"):
                        state = 1
                        continue
                    elif line.startswith(".nf"): continue
                    elif line.startswith(".ni"): continue
                    elif line.startswith(".P"): continue
                    if line.startswith(".SH Description"):
                        state = 2
                        continue
                    elif line.startswith(".SH File"):
                        state = 3
                        continue
                    if line.startswith(".SH Need"):
                        state = 0
                        continue
                    if line.startswith(".SH Instability"):
                        state = 0
                        continue
                    if line.startswith(".SH Call"):
                        state = 0
                        continue

                    if state == 1:
                        entity["decl"] += [line]
                        # Multi-line declarations (Fortran continuation '&')
                        # would otherwise let a continuation line override
                        # the type decided by the first line.
                        if "type" not in entity:
                            if line.startswith("subroutine") \
                            or line.startswith("function ") \
                            or " function " in line:
                                entity["type"] = 's'
                            else:
                                entity["type"] = 'p'
                    elif state == 2:
                        if line.startswith(".br"):
                          line = "\n\n"
                        entity["doc"] += [line]
                    elif state == 3:
                        if line.startswith(".br"):
                            continue
                        entity["file"] = line.split("/")[-1]
                        try:
                            entity["module"] = line.split("/")[-2]
                        except: pass
                        break

            entities[entity["name"]] = entity

  return entities


def generate_index(entities):

  rst_file = os.path.join('programmers_guide','index_providers.rst')

  with open(rst_file,'w') as f: 
    rst = [ "Index of Providers", 
            "------------------", 
            "",
            ".. hlist::",
            ""]

    for e in sorted(entities.keys()):
        e = entities[e]
        # Entities declared in the auto-generated ezfio_interface.irp.f are
        # not documented here: they are the EZFIO parameters of their module
        # and are rendered from EZFIO.cfg in the module page.
        if e.get("file") == "ezfio_interface.irp.f":
            continue
        if e["type"] == 'p':
            rst.append("    * :c:data:`%s`" % (e["name"]))

    rst += [ "","",
             "Index of Subroutines/Functions", 
             "------------------------------", 
             "",
             ".. hlist::",
             "" ]

    for e in sorted(entities.keys()):
        e = entities[e]
        if e.get("file") == "ezfio_interface.irp.f":
            continue
        if e["type"] == 's':
            rst.append("    * :c:func:`%s`" % (e["name"]))

    rst.append("\n")
    f.write(" \n".join(rst))



def main():

  if "QP_ROOT" in os.environ:
    QP_ROOT=os.environ["QP_ROOT"]
  else:
    QP_ROOT="../../"

  SRC = os.path.join(QP_ROOT, "src")

  entities = {}
  local_entities = {}
  for abs_module in os.listdir(SRC):
    if os.path.islink(os.path.join(SRC,abs_module)):
      continue
    abs_module = os.path.join(SRC,abs_module)
    # Collect every entity from every IRPF90_man: provider-only modules
    # have no main, are never built, and only get their man pages inside
    # the IRPF90_man of the programs that NEED them; the owning module of
    # each page is recovered from the .SH File section of the .l file.
    if os.path.isdir( os.path.join(abs_module, "IRPF90_man") ):
      read_entities = generate_providers(abs_module)
      if read_entities:
        local_entities[os.path.basename(abs_module).lower()] = read_entities
        for k in read_entities:
            entities[k] = read_entities[k]

  owned = set(e["module"].lower() for e in entities.values()
              if isinstance(e.get("module"), str))

  for abs_module in os.listdir(SRC):
    abs_module = os.path.join(SRC,abs_module)
    if os.path.islink(abs_module):
      continue
    base = os.path.basename(abs_module).lower()
    if os.path.exists( os.path.join(abs_module, "README.rst") ) \
       or base in owned:
        generate_modules(abs_module,entities,local_entities.get(base))

  generate_index(entities)

if __name__ == '__main__':
  main()

