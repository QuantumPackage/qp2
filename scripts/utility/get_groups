#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import urllib
import sys
from bs4 import BeautifulSoup

address="http://gernot-katzers-spice-pages.com/character_tables/%s.html?fmt=simple"

def clean_up(text):
    soup = BeautifulSoup(text, "lxml")
    pre = soup.pre
    group = pre.b.get_text()
    sop = {}
    irred = {}
    irred_count = 0
    sop_count = 0
    for span in pre.find_all('span'):
        cls = span.get('class')
        if cls == ['sop']:
          a = span.decode_contents() 
          if a not in sop:
            sop[a] = sop_count
            sop_count += 1
        elif cls == ['irred']:
          a = span.decode_contents() 
          if a not in irred:
            irred[a] = irred_count
            irred_count += 1
    table = [ [] for j in sop ]
    data = pre.get_text().splitlines()
    def f(x):
      y = x.split()
      if len(y) == 0:
        return False
      else:
        return y[0] in irred
    data = filter(f,data)[:len(irred)]
    for line in data:
       s = line.replace('*','').split()
       l = irred[s[0]]
       data[l] = map(float,s[1:len(irred)+1])

    d = {}
    e = {}
    for k in irred:
      d[irred[k]] = k
    for k in sop:
      e[sop[k]] = k
    n = len(irred)
    print "Group\t", group, "\nn\t", n
    print "\n     \tIrred     \tOperation"
    for i in range(n):
      print "%4d \t %s \t %s"%(i+1, d[i].ljust(10), e[i].ljust(10))

    print "\nTable\n     ",
    for j in range(n):
        print "%8s "%(str(j+1).center(8)),
    for i in range(n):
      print "\n%4d "%(i+1), 
      for j in range(n):
        print "%8.5f "%(data[i][j]), 
    print "\n"

def main():
    for group in sys.argv[1:]:
      f = urllib.urlopen(address%(group))
      html = f.read().split('\n',1)[1]
      text = clean_up(html)

if __name__ == "__main__":
    main()
