#!/usr/bin/env python3

import os
import json

def fix_json(s):
   """Properly termitates an incomplete JSON file"""

   s = s.replace(' ','')
   s = s.replace('\n','')
   s = s.replace('\t','')
   s = s.replace(",{}",'')
   tmp = [ c for c in s if c in "[]{}" ]
   tmp = "".join(tmp)
   tmp_old = ""
   while tmp != tmp_old:
      tmp_old = tmp
      tmp = tmp.replace("{}","")
      tmp = tmp.replace("[]","")
   while s[-1] in [ ',', '\n', ' ', '\t' ]:
      s = s[:-1]
   tmp = [ c for c in tmp ]
   tmp.reverse()
   for c in tmp:
      if c == '[': s += "]"
      elif c == '{': s += "}"
   return s


def load(filename):
  """Loads a JSON file after calling the fix_json function."""
  with open(filename,'r') as f:
    data = f.read()
  new_data = fix_json(data)
  return json.loads(new_data)


def load_all(ezfio_filename):
  """Loads all JSON files of an EZFIO."""
  d = {}
  prefix = ezfio_filename+'/json/'
  for filename in [ x for x in os.listdir(prefix) if x.endswith(".json")]:
    d[filename] = load(prefix+filename)
  return d


def load_last(ezfio_filename):
  """Loads last JSON file of an EZFIO."""
  d = {}
  prefix = ezfio_filename+'/json/'
  l = [ x for x in os.listdir(prefix) if x.endswith(".json")]
  l.sort()
  filename = l[-1]
  print(filename)
  return load(prefix+filename)


def fix(ezfio_filename):
  """Fixes all JSON files in an EZFIO."""
  d = load_all(ezfio_filename)
  prefix = ezfio_filename+'/json/'
  for filename in d.keys():
        with open(prefix+filename, 'w') as json_file:
            json.dump(d[filename], json_file)



