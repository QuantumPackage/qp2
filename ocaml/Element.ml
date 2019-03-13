open Sexplib.Std
open Qptypes

exception ElementError of string

type t =
|X
|H                                                 |He
|Li|Be                              |B |C |N |O |F |Ne
|Na|Mg                              |Al|Si|P |S |Cl|Ar
|K |Ca|Sc|Ti|V |Cr|Mn|Fe|Co|Ni|Cu|Zn|Ga|Ge|As|Se|Br|Kr
|Rb|Sr|Y |Zr|Nb|Mo|Tc|Ru|Rh|Pd|Ag|Cd|In|Sn|Sb|Te|I |Xe
                           |Pt
[@@deriving sexp]

let of_string x =
  match (String.capitalize_ascii (String.lowercase_ascii x)) with
|  "X"   |  "Dummy"       ->  X
|  "H"   |  "Hydrogen"    ->  H
|  "He"  |  "Helium"      ->  He
|  "Li"  |  "Lithium"     ->  Li
|  "Be"  |  "Beryllium"   ->  Be
|  "B"   |  "Boron"       ->  B
|  "C"   |  "Carbon"      ->  C
|  "N"   |  "Nitrogen"    ->  N
|  "O"   |  "Oxygen"      ->  O
|  "F"   |  "Fluorine"    ->  F
|  "Ne"  |  "Neon"        ->  Ne
|  "Na"  |  "Sodium"      ->  Na
|  "Mg"  |  "Magnesium"   ->  Mg
|  "Al"  |  "Aluminum"    ->  Al
|  "Si"  |  "Silicon"     ->  Si
|  "P"   |  "Phosphorus"  ->  P
|  "S"   |  "Sulfur"      ->  S
|  "Cl"  |  "Chlorine"    ->  Cl
|  "Ar"  |  "Argon"       ->  Ar
|  "K"   |  "Potassium"   ->  K
|  "Ca"  |  "Calcium"     ->  Ca
|  "Sc"  |  "Scandium"    ->  Sc
|  "Ti"  |  "Titanium"    ->  Ti
|  "V"   |  "Vanadium"    ->  V
|  "Cr"  |  "Chromium"    ->  Cr
|  "Mn"  |  "Manganese"   ->  Mn
|  "Fe"  |  "Iron"        ->  Fe
|  "Co"  |  "Cobalt"      ->  Co
|  "Ni"  |  "Nickel"      ->  Ni
|  "Cu"  |  "Copper"      ->  Cu
|  "Zn"  |  "Zinc"        ->  Zn
|  "Ga"  |  "Gallium"     ->  Ga
|  "Ge"  |  "Germanium"   ->  Ge
|  "As"  |  "Arsenic"     ->  As
|  "Se"  |  "Selenium"    ->  Se
|  "Br"  |  "Bromine"     ->  Br
|  "Kr"  |  "Krypton"     ->  Kr
|  "Rb"  |  "Rubidium"    ->  Rb
|  "Sr"  |  "Strontium"   ->  Sr
|  "Y"   |  "Yttrium"     ->  Y
|  "Zr"  |  "Zirconium"   ->  Zr
|  "Nb"  |  "Niobium"     ->  Nb
|  "Mo"  |  "Molybdenum"  ->  Mo
|  "Tc"  |  "Technetium"  ->  Tc
|  "Ru"  |  "Ruthenium"   ->  Ru
|  "Rh"  |  "Rhodium"     ->  Rh
|  "Pd"  |  "Palladium"   ->  Pd
|  "Ag"  |  "Silver"      ->  Ag
|  "Cd"  |  "Cadmium"     ->  Cd
|  "In"  |  "Indium"      ->  In
|  "Sn"  |  "Tin"         ->  Sn
|  "Sb"  |  "Antimony"    ->  Sb
|  "Te"  |  "Tellurium"   ->  Te
|  "I"   |  "Iodine"      ->  I
|  "Xe"  |  "Xenon"       ->  Xe
|  "Pt"  |  "Platinum"    ->  Pt
| x -> raise (ElementError ("Element "^x^" unknown"))


let to_string = function
| X   -> "X"
| H   -> "H"
| He  -> "He"
| Li  -> "Li"
| Be  -> "Be"
| B   -> "B"
| C   -> "C"
| N   -> "N"
| O   -> "O"
| F   -> "F"
| Ne  -> "Ne"
| Na  -> "Na"
| Mg  -> "Mg"
| Al  -> "Al"
| Si  -> "Si"
| P   -> "P"
| S   -> "S"
| Cl  -> "Cl"
| Ar  -> "Ar"
| K   -> "K"
| Ca  -> "Ca"
| Sc  -> "Sc"
| Ti  -> "Ti"
| V   -> "V"
| Cr  -> "Cr"
| Mn  -> "Mn"
| Fe  -> "Fe"
| Co  -> "Co"
| Ni  -> "Ni"
| Cu  -> "Cu"
| Zn  -> "Zn"
| Ga  -> "Ga"
| Ge  -> "Ge"
| As  -> "As"
| Se  -> "Se"
| Br  -> "Br"
| Kr  -> "Kr"
| Rb  -> "Rb"
| Sr  -> "Sr"
| Y   -> "Y"
| Zr  -> "Zr"
| Nb  -> "Nb"
| Mo  -> "Mo"
| Tc  -> "Tc"
| Ru  -> "Ru"
| Rh  -> "Rh"
| Pd  -> "Pd"
| Ag  -> "Ag"
| Cd  -> "Cd"
| In  -> "In"
| Sn  -> "Sn"
| Sb  -> "Sb"
| Te  -> "Te"
| I   -> "I"
| Xe  -> "Xe"
| Pt  -> "Pt"


let to_long_string = function
| X   -> "Dummy"
| H   -> "Hydrogen"
| He  -> "Helium"
| Li  -> "Lithium"
| Be  -> "Beryllium"
| B   -> "Boron"
| C   -> "Carbon"
| N   -> "Nitrogen"
| O   -> "Oxygen"
| F   -> "Fluorine"
| Ne  -> "Neon"
| Na  -> "Sodium"
| Mg  -> "Magnesium"
| Al  -> "Aluminum"
| Si  -> "Silicon"
| P   -> "Phosphorus"
| S   -> "Sulfur"
| Cl  -> "Chlorine"
| Ar  -> "Argon"
| K   -> "Potassium"
| Ca  -> "Calcium"
| Sc  -> "Scandium"
| Ti  -> "Titanium"
| V   -> "Vanadium"
| Cr  -> "Chromium"
| Mn  -> "Manganese"
| Fe  -> "Iron"
| Co  -> "Cobalt"
| Ni  -> "Nickel"
| Cu  -> "Copper"
| Zn  -> "Zinc"
| Ga  -> "Gallium"
| Ge  -> "Germanium"
| As  -> "Arsenic"
| Se  -> "Selenium"
| Br  -> "Bromine"
| Kr  -> "Krypton"
| Rb  ->  "Rubidium"
| Sr  ->  "Strontium"
| Y   ->  "Yttrium"
| Zr  ->  "Zirconium"
| Nb  ->  "Niobium"
| Mo  ->  "Molybdenum"
| Tc  ->  "Technetium"
| Ru  ->  "Ruthenium"
| Rh  ->  "Rhodium"
| Pd  ->  "Palladium"
| Ag  ->  "Silver"
| Cd  ->  "Cadmium"
| In  ->  "Indium"
| Sn  ->  "Tin"
| Sb  ->  "Antimony"
| Te  ->  "Tellurium"
| I   ->  "Iodine"
| Xe  ->  "Xenon"
| Pt  ->  "Platinum"


let to_charge c =
  let result = match c with
  | X   -> 0
  | H   -> 1
  | He  -> 2
  | Li  -> 3
  | Be  -> 4
  | B   -> 5
  | C   -> 6
  | N   -> 7
  | O   -> 8
  | F   -> 9
  | Ne  -> 10
  | Na  -> 11
  | Mg  -> 12
  | Al  -> 13
  | Si  -> 14
  | P   -> 15
  | S   -> 16
  | Cl  -> 17
  | Ar  -> 18
  | K   -> 19
  | Ca  -> 20
  | Sc  -> 21
  | Ti  -> 22
  | V   -> 23
  | Cr  -> 24
  | Mn  -> 25
  | Fe  -> 26
  | Co  -> 27
  | Ni  -> 28
  | Cu  -> 29
  | Zn  -> 30
  | Ga  -> 31
  | Ge  -> 32
  | As  -> 33
  | Se  -> 34
  | Br  -> 35
  | Kr  -> 36
  | Rb  -> 37
  | Sr  -> 38
  | Y   -> 39
  | Zr  -> 40
  | Nb  -> 41
  | Mo  -> 42
  | Tc  -> 43
  | Ru  -> 44
  | Rh  -> 45
  | Pd  -> 46
  | Ag  -> 47
  | Cd  -> 48
  | In  -> 49
  | Sn  -> 50
  | Sb  -> 51
  | Te  -> 52
  | I   -> 53
  | Xe  -> 54
  | Pt  -> 78
  in Charge.of_int result


let of_charge c = match (Charge.to_int c) with
|  0   -> X
|  1   -> H
|  2   -> He
|  3   -> Li
|  4   -> Be
|  5   -> B
|  6   -> C
|  7   -> N
|  8   -> O
|  9   -> F
|  10  -> Ne
|  11  -> Na
|  12  -> Mg
|  13  -> Al
|  14  -> Si
|  15  -> P
|  16  -> S
|  17  -> Cl
|  18  -> Ar
|  19  -> K
|  20  -> Ca
|  21  -> Sc
|  22  -> Ti
|  23  -> V
|  24  -> Cr
|  25  -> Mn
|  26  -> Fe
|  27  -> Co
|  28  -> Ni
|  29  -> Cu
|  30  -> Zn
|  31  -> Ga
|  32  -> Ge
|  33  -> As
|  34  -> Se
|  35  -> Br
|  36  -> Kr
|  37  -> Rb
|  38  -> Sr
|  39  -> Y
|  40  -> Zr
|  41  -> Nb
|  42  -> Mo
|  43  -> Tc
|  44  -> Ru
|  45  -> Rh
|  46  -> Pd
|  47  -> Ag
|  48  -> Cd
|  49  -> In
|  50  -> Sn
|  51  -> Sb
|  52  -> Te
|  53  -> I
|  54  -> Xe
|  78  -> Pt
| x -> raise (ElementError ("Element of charge "^(string_of_int x)^" unknown"))


let covalent_radius x =
  let result = function
  | X   -> 0.
  | H   -> 0.37
  | He  -> 0.70
  | Li  -> 1.23
  | Be  -> 0.89
  | B   -> 0.90
  | C   -> 0.85
  | N   -> 0.74
  | O   -> 0.74
  | F   -> 0.72
  | Ne  -> 0.70
  | Na  -> 1.00
  | Mg  -> 1.36
  | Al  -> 1.25
  | Si  -> 1.17
  | P   -> 1.10
  | S   -> 1.10
  | Cl  -> 0.99
  | Ar  -> 0.70
  | K   -> 2.03
  | Ca  -> 1.74
  | Sc  -> 1.44
  | Ti  -> 1.32
  | V   -> 1.22
  | Cr  -> 0.00
  | Mn  -> 1.16
  | Fe  -> 0.00
  | Co  -> 1.15
  | Ni  -> 1.17
  | Cu  -> 1.25
  | Zn  -> 1.25
  | Ga  -> 1.20
  | Ge  -> 1.21
  | As  -> 1.16
  | Se  -> 0.70
  | Br  -> 1.24
  | Kr  -> 1.91
  | Rb  -> 2.20
  | Sr  -> 1.95
  | Y   -> 1.90
  | Zr  -> 1.75
  | Nb  -> 1.64
  | Mo  -> 1.54
  | Tc  -> 1.47
  | Ru  -> 1.46
  | Rh  -> 1.42
  | Pd  -> 1.39
  | Ag  -> 1.45
  | Cd  -> 1.44
  | In  -> 1.42
  | Sn  -> 1.39
  | Sb  -> 1.39
  | Te  -> 1.38
  | I   -> 1.39
  | Xe  -> 1.40
  | Pt  -> 1.30
  in
  Units.angstrom_to_bohr *. (result x)
  |> Positive_float.of_float

let vdw_radius x =
  let result = function
  | X    -> 0.
  | H    -> 1.20
  | He  ->  1.70
  | Li  ->  1.70
  | Be  ->  1.70
  | B   ->  1.70
  | C   ->  1.70
  | N   ->  1.55
  | O   ->  1.52
  | F   ->  1.47
  | Ne  ->  1.70
  | Na  ->  1.70
  | Mg  ->  1.70
  | Al  ->  1.94
  | Si  ->  2.10
  | P   ->  1.80
  | S   ->  1.80
  | Cl  ->  1.75
  | Ar  ->  1.70
  | K   ->  1.70
  | Ca  ->  1.70
  | Sc  ->  1.70
  | Ti  ->  1.70
  | V   ->  1.98
  | Cr  ->  1.94
  | Mn  ->  1.93
  | Fe  ->  1.93
  | Co  ->  1.92
  | Ni  ->  1.70
  | Cu  ->  1.70
  | Zn  ->  1.70
  | Ga  ->  2.02
  | Ge  ->  1.70
  | As  ->  1.96
  | Se  ->  1.70
  | Br  ->  2.10
  | Kr  ->  1.70
  | Rb  -> 3.03
  | Sr  -> 2.49
  | Y   -> 0.
  | Zr  -> 0.
  | Nb  -> 0.
  | Mo  -> 0.
  | Tc  -> 0.
  | Ru  -> 0.
  | Rh  -> 0.
  | Pd  -> 1.63
  | Ag  -> 1.72
  | Cd  -> 1.58
  | In  -> 1.93
  | Sn  -> 2.17
  | Sb  -> 2.06
  | Te  -> 2.06
  | I   -> 1.98
  | Xe  -> 2.16
  | Pt  -> 1.75
  in
  Units.angstrom_to_bohr *. (result x)
  |> Positive_float.of_float

let mass x =
  let result = function
  | X   ->    0.
  | H   ->    1.0079
  | He  ->    4.00260
  | Li  ->    6.941
  | Be  ->    9.01218
  | B   ->   10.81
  | C   ->   12.011
  | N   ->   14.0067
  | O   ->   15.9994
  | F   ->   18.998403
  | Ne  ->   20.179
  | Na  ->   22.98977
  | Mg  ->   24.305
  | Al  ->   26.98154
  | Si  ->   28.0855
  | P   ->   30.97376
  | S   ->   32.06
  | Cl  ->   35.453
  | Ar  ->   39.948
  | K   ->   39.0983
  | Ca  ->   40.08
  | Sc  ->   44.9559
  | Ti  ->   47.90
  | V   ->   50.9415
  | Cr  ->   51.996
  | Mn  ->   54.9380
  | Fe  ->   55.9332
  | Co  ->   58.9332
  | Ni  ->   58.70
  | Cu  ->   63.546
  | Zn  ->   65.38
  | Ga  ->   69.72
  | Ge  ->   72.59
  | As  ->   74.9216
  | Se  ->   78.96
  | Br  ->   79.904
  | Kr  ->   83.80
  | Rb  -> 85.4678
  | Sr  -> 87.62
  | Y   -> 88.90584
  | Zr  -> 91.224
  | Nb  -> 92.90637
  | Mo  -> 95.95
  | Tc  -> 98.
  | Ru  -> 101.07
  | Rh  -> 102.90550
  | Pd  -> 106.42
  | Ag  -> 107.8682
  | Cd  -> 112.414
  | In  -> 114.818
  | Sn  -> 118.710
  | Sb  -> 121.760
  | Te  -> 127.60
  | I   -> 126.90447
  | Xe  -> 131.293
  | Pt  -> 195.084
  in
  result x
  |> Positive_float.of_float


