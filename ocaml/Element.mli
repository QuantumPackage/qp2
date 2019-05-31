exception ElementError of string

type t = X

|H                                                 |He
|Li|Be                              |B |C |N |O |F |Ne
|Na|Mg                              |Al|Si|P |S |Cl|Ar
|K |Ca|Sc|Ti|V |Cr|Mn|Fe|Co|Ni|Cu|Zn|Ga|Ge|As|Se|Br|Kr
|Rb|Sr|Y |Zr|Nb|Mo|Tc|Ru|Rh|Pd|Ag|Cd|In|Sn|Sb|Te|I |Xe
|Cs|Ba|La|Hf|Ta|W |Re|Os|Ir|Pt|Au|Hg|Tl|Pb|Bi|Po|At|Rn
|Fr|Ra|Ac|Rf|Db|Sg|Bh|Hs|Mt|Ds|Rg|Cn|Nh|Fl|Mc|Lv|Ts|Og

      |Ce|Pr|Nd|Pm|Sm|Eu|Gd|Tb|Dy|Ho|Er|Tm|Yb|Lu
      |Th|Pa|U |Np|Pu|Am|Cm|Bk|Cf|Es|Fm|Md|No|Lr
        
[@@deriving sexp]

(** String conversion functions *)
val of_string : string -> t
val to_string : t -> string
val to_long_string : t -> string

(** Properties *)
val to_charge : t -> Charge.t
val of_charge : Charge.t -> t
val covalent_radius : t -> Qptypes.Positive_float.t
val vdw_radius : t -> Qptypes.Positive_float.t option
val mass : t -> Qptypes.Positive_float.t
