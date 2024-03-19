#!/usr/bin/env python3

import zmq
import sys, os

def main():
  context = zmq.Context()
  socket = context.socket(zmq.SUB)
  socket.connect("tcp://127.0.0.1:41280")
  socket.setsockopt(zmq.SUBSCRIBE, "")
  while True:
    print socket.recv()

if __name__ == '__main__':
  main()
