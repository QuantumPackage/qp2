include String

(** Split a string on a given character *)
let split ?(on=' ') str =
  split_on_char on str


(** Strip blanks on the left of a string *)
let ltrim s =
    let rec do_work s l =
        match s.[0] with
        | '\n'
        | ' ' -> do_work (sub s 1 (l-1)) (l-1)
        | _ -> s
    in
    let l =
        length s
    in
    if (l > 0) then
        do_work s l
    else
        s

(** Strip blanks on the right of a string *)
let rtrim s =
    let rec do_work s l =
        let newl =
            l-1
        in
        match s.[newl] with
        | '\n'
        | ' ' -> do_work (sub s 0 (newl)) (newl)
        | _ -> s
    in
    let l =
        length s
    in
    if (l > 0) then
        do_work s l
    else
        s


(** Strip blanks on the right and left of a string *)
let strip = String.trim


(** Split a string in two pieces when a character is found the 1st time from the left *)
let lsplit2_exn ?(on=' ') s =
    let length =
        String.length s
    in
    let rec do_work i =
        if (i = length) then
           begin
              raise Not_found
           end
        else if (s.[i] = on) then
           ( String.sub s 0 i,
             String.sub s (i+1) (length-i-1) )
        else
           do_work (i+1)
    in
    do_work 0


(** Split a string in two pieces when a character is found the 1st time from the right *)
let rsplit2_exn ?(on=' ') s =
    let length =
        String.length s 
    in
    let rec do_work i =
        if (i = -1) then
           begin
              raise Not_found
           end
        else if (s.[i] = on) then
           ( String.sub s 0 i,
             String.sub s (i+1) (length-i-1) )
        else
           do_work (i-1)
    in
    do_work (length-1)


let lsplit2 ?(on=' ') s =
  try
    Some (lsplit2_exn ~on s)
  with
  | Not_found -> None


let rsplit2 ?(on=' ') s =
  try
    Some (rsplit2_exn ~on s)
  with
  | Not_found -> None


let to_list s =
  Array.init (String.length s) (fun i -> s.[i])
  |> Array.to_list


let of_list l =
  let a = Array.of_list l in
  String.init (Array.length a) (fun i -> a.(i))

let rev s =
  to_list s
  |> List.rev
  |> of_list

let fold ~init ~f s =
    to_list s
    |> List.fold_left f init


let is_prefix ~prefix s =
  let len =
    String.length prefix
  in
  if len > String.length s then
    false
  else
    prefix = String.sub s 0 len


let of_char c =
  String.make 1 c

let tr ~target ~replacement s =
  String.map (fun c -> if c = target then replacement else c) s


let substr_index ?(pos=0) ~pattern s =
  try
    let regexp = 
      Str.regexp pattern
    in
    Some (Str.search_forward regexp s pos)
  with Not_found -> None


let substr_replace_all ~pattern ~with_ s =
  let regexp =
    Str.regexp pattern
  in
  Str.global_replace regexp with_ s

