module Id = struct
  type t = int

  let of_int x =
    assert (x>0); x

  let to_int x = x

  let of_string x =
    int_of_string x
    |> of_int

  let to_string x =
    string_of_int x

  let increment x = x + 1
  let decrement x = x - 1

  let compare = compare
end

module Task = struct
  include Id
end

module Client = struct
  include Id
end

