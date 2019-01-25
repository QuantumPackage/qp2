#!/usr/bin/python

import zmq
import sys, os

def main():
  context = zmq.Context()
  socket = context.socket(zmq.REQ)
  socket.connect(os.environ["QP_RUN_ADDRESS"])

  def send(msg,expected):
    print "Send  : ", msg
    socket.send(msg)
    reply = socket.recv()
    print "Reply   : ", ':'+reply+':'
    if (reply != expected):
      print "Expected: ", ':'+expected+':'
    print ""
    assert (reply == expected)
       

  send("new_job ao_integrals tcp://130.120.229.139:12345 inproc://ao_integrals",
       "ok")
  send("new_job ao_integrals tcp://130.120.229.139:12345 inproc://ao_integrals",
       "error A job is already running")

#  send("connect","error Message not understood : connect")

  send("connect  tcp","connect_reply ao_integrals 1 tcp://130.120.229.139:12345")
  send("connect  inproc","connect_reply ao_integrals 2 inproc://ao_integrals")
  send("disconnect ao_integrals 3","error Queuing_system.ml:68:2 : disconnect ao_integrals 3")
  send("disconnect ao_integrals 2","disconnect_reply ao_integrals")
  send("connect  inproc","connect_reply ao_integrals 3 inproc://ao_integrals")

  send("add_task ao_integrals triangle 3", "ok")
  send("add_task ao_integrals range 4 7", "ok")

  for i in range(8,11):
     send("add_task ao_integrals %d %d"%(i,i+10), "ok")

  send("get_task ao_integrals 3", "get_task_reply 10 10 20")
  send("get_task ao_integrals 3", "get_task_reply 9 9 19")
  send("get_task ao_integrals 3", "get_task_reply 8 8 18")

  send("task_done ao_integrals 3 10", "ok")
  send("task_done ao_integrals 3 9", "ok")
  send("task_done ao_integrals 3 8", "ok")
  send("del_task ao_integrals 10", "del_task_reply more 10")
  send("del_task ao_integrals 9", "del_task_reply more 9")
  send("del_task ao_integrals 8", "del_task_reply more 8")
  send("del_task ao_integrals 10", "error Task 10 is already deleted : del_task ao_integrals 10")

  send("get_task ao_integrals 1", "get_task_reply 7 4")
  send("get_task ao_integrals 3", "get_task_reply 6 5")
  send("get_task ao_integrals 1", "get_task_reply 5 6")
  send("get_task ao_integrals 3", "get_task_reply 4 7")
  send("get_task ao_integrals 3", "get_task_reply 3 1 3")
  send("get_task ao_integrals 1", "get_task_reply 2 2 3")
  send("get_task ao_integrals 1", "get_task_reply 1 3 3")

  send("task_done ao_integrals 1 1", "ok")
  send("task_done ao_integrals 1 2", "ok")
  send("task_done ao_integrals 3 3", "ok")
  send("task_done ao_integrals 3 4", "ok")
  send("task_done ao_integrals 1 5", "ok")
  send("task_done ao_integrals 1 6", "error Queuing_system.ml:81:30 : task_done ao_integrals 1 6")
  send("task_done ao_integrals 3 6", "ok")
  send("task_done ao_integrals 1 7", "ok")

  send("del_task ao_integrals 1", "del_task_reply more 1")
  send("del_task ao_integrals 2", "del_task_reply more 2")
  send("del_task ao_integrals 3", "del_task_reply more 3")
  send("del_task ao_integrals 4", "del_task_reply more 4")
  send("del_task ao_integrals 5", "del_task_reply more 5")
  send("del_task ao_integrals 6", "del_task_reply more 6")
  send("del_task ao_integrals 7", "del_task_reply done 7")

  send("end_job ao_integrals","ok")
  send("end_job ao_integrals","error No job is running")
  send("terminate","ok")

if __name__ == '__main__':
  main()
