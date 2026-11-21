let () =
  TaskServer.run ~rep_socket:None ~pub_socket:None ~port:12345
                 ~address:(Printf.sprintf "tcp://%s:%d" (Lazy.force TaskServer.ip_address) 12345)
