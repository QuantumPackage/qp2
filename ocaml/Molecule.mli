exception MultiplicityError of string

type t = {
  nuclei     : Atom.t list;
  elec_alpha : Qptypes.Elec_alpha_number.t;
  elec_beta  : Qptypes.Elec_beta_number.t;
} [@@deriving sexp]

(** Returns the charge of the molecule *)
val get_charge       : t -> Charge.t

(** Returns the multiplicity of the molecule *)
val get_multiplicity : t -> Multiplicity.t

(** Returns the number of nuclei *)
val get_nucl_num : t -> Qptypes.Nucl_number.t

(** The name of the molecule *)
val name : t -> string

(** Conversion for printing *)
val to_string : t -> string
val to_xyz    : t -> string


(** Creates a molecule from an xyz file *)
val of_xyz_file :
  ?charge:Charge.t ->
  ?multiplicity:Multiplicity.t ->
  ?units:Units.units -> string -> t

(** Creates a molecule from a zmt file *)
val of_zmt_file :
  ?charge:Charge.t ->
  ?multiplicity:Multiplicity.t ->
  ?units:Units.units -> string -> t

(** Creates a molecule from a file (xyz or zmt) *)
val of_file :
  ?charge:Charge.t ->
  ?multiplicity:Multiplicity.t ->
  ?units:Units.units -> string -> t

(** Creates a molecule from an xyz file in a string *)
val of_xyz_string :
  ?charge:Charge.t ->
  ?multiplicity:Multiplicity.t ->
  ?units:Units.units -> string -> t

(** Creates the distance matrix between all the atoms *)
val distance_matrix :
  t -> (float array) array

(** Computes the MD5 hash *)
val to_md5 : t -> Qptypes.MD5.t
