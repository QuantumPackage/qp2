#!/usr/bin/env python3


import os
import sys
import configparser

from module_handler import get_binaries


def generate_modules(abs_module, entities):
  """Generates the doc for modules"""
  MODULE = os.path.split(abs_module)[-1]
  module = MODULE.lower()
  if module == "dummy":
    return

  with open(os.path.join(abs_module, 'README.rst'), 'r') as f:
    readme = f.read()
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
        config_file.readfp(f)
        for section in config_file.sections():
            doc = config_file.get(section, "doc")
            doc = "    " + doc.replace("\n", "\n\n    ")+"\n"
            try:
                default = config_file.get(section, "default")
                default = "    " + "Default: %s\n"%default
            except:
                default = ""
            rst += [".. option:: %s\n"%(section), doc, default]

  providers = []
  subroutines = {}
  for k in sorted(entities.keys()):
    e = entities[k]
    if e["module"].lower() == module.lower():
        if "/" not in e["file"] and e["file"] != "ezfio_interface.irp.f":
            if e["type"] == 's':
                subroutines[e["name"]] = e
            elif e["type"] == 'p':
                providers.append(e)

  binaries = [os.path.basename(f) for f in get_binaries(abs_module)]

  if binaries:
    rst += ["", "Programs", "--------", ""]
    for b in binaries:
        try:
            b = subroutines[b]
        except KeyError:
            print("Error: The program %s in %s does not have the same name as the file, \
or you did not run ninja at the root."%
                (b, abs_module))
            sys.exit(1)
        rst += [" * :ref:`%s`"%(b["name"])]

  if providers:
    rst += ["", "Providers", "---------", ""]
    for p in providers:
        rst += [p["rst"]]

  if subroutines:
    rst += [ "", "Subroutines / functions", "-----------------------", "" ]
    for p in sorted(subroutines.keys()):
        p = subroutines[p]
        if p["name"] in binaries:
           continue
        rst += [p["rst"]]

  rst_file = os.path.join('modules', module+".rst")
  with open(rst_file,'w') as f: 
      f.write(" \n".join(rst))

  for b in subroutines:
    if b not in binaries:
      continue
    p = subroutines[b]
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
  for abs_module in os.listdir(SRC):
    if os.path.islink(os.path.join(SRC,abs_module)):
      continue
    abs_module = os.path.join(SRC,abs_module)
    if os.path.exists( os.path.join(abs_module, "README.rst") ):
      read_entities = generate_providers(abs_module)
      if read_entities:
        for k in read_entities:
            entities[k] = read_entities[k]

  for abs_module in os.listdir(SRC):
    abs_module = os.path.join(SRC,abs_module)
    if os.path.islink(os.path.join(SRC,abs_module)):
      continue
    if os.path.exists( os.path.join(abs_module, "README.rst") ):
        generate_modules(abs_module,entities)

  generate_index(entities)

if __name__ == '__main__':
  main()

