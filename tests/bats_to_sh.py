#!/usr/bin/env python3

import sys

with open(sys.argv[1],'r') as f:
  raw_data = f.read()

print("set -x")

output = []
inside = False
level = 0
for i in raw_data:
  new_i = i
  if i == "@":
    inside = True
  elif i == "{" and inside and level == 0:
    new_i = "\nfunction _run_test() {\n setup\n"
  elif i == "}" and inside and level == 1:
    inside = False
    new_i = "}\n_run_test || exit 1"
  if i == "{":
    level += 1
  elif i == "}":
    level -= 1
  output.append(new_i)

print("".join(output).replace("@test ",
"""[[ -z $BATS_TEST_NUMBER ]] && BATS_TEST_NUMBER=0 || ((++BATS_TEST_NUMBER)) ;
export BATS_TEST_DESCRIPTION=""").replace("skip","return"))



