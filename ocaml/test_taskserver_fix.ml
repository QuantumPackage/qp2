(* test_taskserver_fix.ml — deterministic regression tests for GitHub
   issue #393 (qp_run hangs at 0% CPU after the calculation completes
   when the TaskServer PUB socket cannot be bound).

   Usage: ./test_taskserver_fix <case>
     1: TaskServer.stop against a dead server: must fail with a clear
        error after ~30 s (a dead server must not hang the launcher).
     2: legacy TaskServer.run with the PUB port pre-occupied: the REP
        must still answer; after the 60 s PUB bind retries fail,
        run() must return (no Zmq.Context.terminate deadlock) and
        pub_error() must be set.
     3: TaskServer.run with pre-bound REP+PUB sockets (the qp_run
        path): clean start/stop, pub_error() = None.

   Every case prints "PASS <case> ..." and exits 0 on success.
   The test harness wraps each case in a `timeout` as a backstop. *)

let port_base = 43_210

let fail msg =
  Printf.eprintf "FAIL: %s\n%!" msg;
  exit 1

let expect_ok port =
  let ctx = Zmq.Context.create () in
  let sock = Zmq.Socket.create ctx Zmq.Socket.req in
  Zmq.Socket.set_linger_period sock 1000;
  Zmq.Socket.connect sock (Printf.sprintf "tcp://localhost:%d" port);
  (try
     Printf.sprintf "new_job scf tcp://localhost:%d inproc://scf" port
     |> Zmq.Socket.send sock;
     let reply = Zmq.Socket.recv sock in
     if reply <> "ok" then fail (Printf.sprintf "expected ok, got %s" reply)
   with other ->
     fail (Printf.sprintf "REQ failed: %s" (Printexc.to_string other));
   (* fail exits, so the cleanup below runs on every non-failing path
      (OCaml has no try...finally). *)
   Zmq.Socket.close sock;
   Zmq.Context.terminate ctx)

let () =
  if Array.length Sys.argv <> 2 then
    failwith "usage: test_taskserver_fix <1|2|3>"
  else
    let case = int_of_string Sys.argv.(1) in
    let port = port_base + 10 * case in
    match case with
    | 1 ->
        let t0 = Unix.time () in
        (try
           TaskServer.stop ~port:port;
           fail "stop() returned against a dead server (expected a failure)"
         with
         | Failure msg ->
             let dt = Unix.time () -. t0 in
             if dt < 20. || dt > 60. then
               fail (Printf.sprintf "stop() gave up after %.0f s (expected ~30 s)" dt)
             else
               (Printf.printf "PASS 1: stop() failed cleanly after %.0f s: %s\n%!" dt msg;
                exit 0))
    | 2 ->
        (* The harness occupies port+1 (PUB) externally for 300 s.
           None/None: legacy path, the server creates and binds its
           sockets itself (the PUB bind must fail after 60 s). *)
        let thread =
          Thread.create
            (fun () -> TaskServer.run
                         ~rep_socket:None ~pub_socket:None ~port:port
                         ~address:(Printf.sprintf "tcp://localhost:%d" port))
            ()
        in
        Thread.delay 3.;
        (match TaskServer.pub_error () with
         | Some _ -> fail "pub_error set before the 60 s PUB retries finished"
         | None -> ());
        expect_ok port;  (* REP must answer mid-failure *)
        let t_stop = Unix.time () in
        TaskServer.stop ~port:port;
        Thread.join thread;
        let dt = Unix.time () -. t_stop in
        let () =
          if dt > 120. then
            fail (Printf.sprintf "run() took %.0f s after stop (deadlock?)" dt)
        in
        (match TaskServer.pub_error () with
         | None -> fail "pub_error() is None but the PUB port was occupied"
         | Some msg ->
             Printf.printf "PASS 2: run() returned after %.0f s; pub_error: %s\n%!" dt msg;
             exit 0)
    | 3 ->
        let ctx = Zmq.Context.create () in
        let rep = Zmq.Socket.create ctx Zmq.Socket.rep in
        let pub = Zmq.Socket.create ctx Zmq.Socket.pub in
        Zmq.Socket.set_linger_period rep 1_000_000;
        Zmq.Socket.bind rep (Printf.sprintf "tcp://*:%d" port);
        Zmq.Socket.bind pub (Printf.sprintf "tcp://*:%d" (port + 1));
        let thread =
          Thread.create (fun () ->
              TaskServer.run ~port:port
                            ~rep_socket:(Some rep) ~pub_socket:(Some pub)
                            ~address:(Printf.sprintf "tcp://localhost:%d" port)) ()
        in
        Thread.delay 2.;
        expect_ok port;
        let t_stop = Unix.time () in
        TaskServer.stop ~port:port;
        Thread.join thread;
        Zmq.Context.terminate ctx;
        let dt = Unix.time () -. t_stop in
        let () =
          if dt > 10. then fail (Printf.sprintf "case 3: stop/join took %.0f s" dt)
        in
        (match TaskServer.pub_error () with
         | Some msg -> fail (Printf.sprintf "case 3: pub_error set: %s" msg)
         | None ->
             Printf.printf "PASS 3: pre-bound run/stop clean in %.1f s\n%!" dt;
             exit 0)
    | _ -> failwith "unknown case"
