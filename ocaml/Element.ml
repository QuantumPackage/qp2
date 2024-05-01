open Sexplib.Std
open Qptypes

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

let of_string x =
  match (String.capitalize_ascii (String.lowercase_ascii x)) with
|  "X"   |   "Ghost"        ->  X
|  "H"   |   "Hydrogen"     ->  H
|  "He"  |   "Helium"       ->  He
|  "Li"  |   "Lithium"      ->  Li
|  "Be"  |   "Beryllium"    ->  Be
|  "B"   |   "Boron"        ->  B
|  "C"   |   "Carbon"       ->  C
|  "N"   |   "Nitrogen"     ->  N
|  "O"   |   "Oxygen"       ->  O
|  "F"   |   "Fluorine"     ->  F
|  "Ne"  |   "Neon"         ->  Ne
|  "Na"  |   "Sodium"       ->  Na
|  "Mg"  |   "Magnesium"    ->  Mg
|  "Al"  |   "Aluminum"     ->  Al
|  "Si"  |   "Silicon"      ->  Si
|  "P"   |   "Phosphorus"   ->  P
|  "S"   |   "Sulfur"       ->  S
|  "Cl"  |   "Chlorine"     ->  Cl
|  "Ar"  |   "Argon"        ->  Ar
|  "K"   |   "Potassium"    ->  K
|  "Ca"  |   "Calcium"      ->  Ca
|  "Sc"  |   "Scandium"     ->  Sc
|  "Ti"  |   "Titanium"     ->  Ti
|  "V"   |   "Vanadium"     ->  V
|  "Cr"  |   "Chromium"     ->  Cr
|  "Mn"  |   "Manganese"    ->  Mn
|  "Fe"  |   "Iron"         ->  Fe
|  "Co"  |   "Cobalt"       ->  Co
|  "Ni"  |   "Nickel"       ->  Ni
|  "Cu"  |   "Copper"       ->  Cu
|  "Zn"  |   "Zinc"         ->  Zn
|  "Ga"  |   "Gallium"      ->  Ga
|  "Ge"  |   "Germanium"    ->  Ge
|  "As"  |   "Arsenic"      ->  As
|  "Se"  |   "Selenium"     ->  Se
|  "Br"  |   "Bromine"      ->  Br
|  "Kr"  |   "Krypton"      ->  Kr
|  "Rb"  |   "Rubidium"     ->  Rb
|  "Sr"  |   "Strontium"    ->  Sr
|  "Y"   |   "Yttrium"      ->  Y
|  "Zr"  |   "Zirconium"    ->  Zr
|  "Nb"  |   "Niobium"      ->  Nb
|  "Mo"  |   "Molybdenum"   ->  Mo
|  "Tc"  |   "Technetium"   ->  Tc
|  "Ru"  |   "Ruthenium"    ->  Ru
|  "Rh"  |   "Rhodium"      ->  Rh
|  "Pd"  |   "Palladium"    ->  Pd
|  "Ag"  |   "Silver"       ->  Ag
|  "Cd"  |   "Cadmium"      ->  Cd
|  "In"  |   "Indium"       ->  In
|  "Sn"  |   "Tin"          ->  Sn
|  "Sb"  |   "Antimony"     ->  Sb
|  "Te"  |   "Tellurium"    ->  Te
|  "I"   |   "Iodine"       ->  I
|  "Xe"  |   "Xenon"        ->  Xe
|  "Cs"  |   "Cesium"       ->  Cs
|  "Ba"  |   "Barium"       ->  Ba
|  "La"  |   "Lanthanum"    ->  La
|  "Ce"  |   "Cerium"       ->  Ce
|  "Pr"  |   "Praseodymium"  ->  Pr
|  "Nd"  |   "Neodymium"    ->  Nd
|  "Pm"  |   "Promethium"   ->  Pm
|  "Sm"  |   "Samarium"     ->  Sm
|  "Eu"  |   "Europium"     ->  Eu
|  "Gd"  |   "Gadolinium"   ->  Gd
|  "Tb"  |   "Terbium"      ->  Tb
|  "Dy"  |   "Dysprosium"   ->  Dy
|  "Ho"  |   "Holmium"      ->  Ho
|  "Er"  |   "Erbium"       ->  Er
|  "Tm"  |   "Thulium"      ->  Tm
|  "Yb"  |   "Ytterbium"    ->  Yb
|  "Lu"  |   "Lutetium"     ->  Lu
|  "Hf"  |   "Hafnium"      ->  Hf
|  "Ta"  |   "Tantalum"     ->  Ta
|  "W"   |   "Tungsten"     ->  W
|  "Re"  |   "Rhenium"      ->  Re
|  "Os"  |   "Osmium"       ->  Os
|  "Ir"  |   "Iridium"      ->  Ir
|  "Pt"  |   "Platinum"     ->  Pt
|  "Au"  |   "Gold"         ->  Au
|  "Hg"  |   "Mercury"      ->  Hg
|  "Tl"  |   "Thallium"     ->  Tl
|  "Pb"  |   "Lead"         ->  Pb
|  "Bi"  |   "Bismuth"      ->  Bi
|  "Po"  |   "Polonium"     ->  Po
|  "At"  |   "Astatine"     ->  At
|  "Rn"  |   "Radon"        ->  Rn
|  "Fr"  |   "Francium"     ->  Fr
|  "Ra"  |   "Radium"       ->  Ra
|  "Ac"  |   "Actinium"     ->  Ac
|  "Th"  |   "Thorium"      ->  Th
|  "Pa"  |   "Protactinium" ->  Pa
|  "U"   |   "Uranium"      ->  U
|  "Np"  |   "Neptunium"    ->  Np
|  "Pu"  |   "Plutonium"    ->  Pu
|  "Am"  |   "Americium"    ->  Am
|  "Cm"  |   "Curium"       ->  Cm
|  "Bk"  |   "Berkelium"    ->  Bk
|  "Cf"  |   "Californium"  ->  Cf
|  "Es"  |   "Einsteinium"  ->  Es
|  "Fm"  |   "Fermium"      ->  Fm
|  "Md"  |   "Mendelevium"  ->  Md
|  "No"  |   "Nobelium"     ->  No
|  "Lr"  |   "Lawrencium"   ->  Lr
|  "Rf"  |   "Rutherfordium"->  Rf
|  "Db"  |   "Dubnium"      ->  Db
|  "Sg"  |   "Seaborgium"   ->  Sg
|  "Bh"  |   "Bohrium"      ->  Bh
|  "Hs"  |   "Hassium"      ->  Hs
|  "Mt"  |   "Meitnerium"   ->  Mt
|  "Ds"  |   "Darmstadtium" ->  Ds
|  "Rg"  |   "Roentgenium"  ->  Rg
|  "Cn"  |   "Copernicium"  ->  Cn
|  "Nh"  |   "Nihonium"     ->  Nh
|  "Fl"  |   "Flerovium"    ->  Fl
|  "Mc"  |   "Moscovium"    ->  Mc
|  "Lv"  |   "Livermorium"  ->  Lv
|  "Ts"  |   "Tennessine"   ->  Ts
|  "Og"  |   "Oganesson"    ->  Og
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
| Cs  -> "Cs"
| Ba  -> "Ba"
| La  -> "La"
| Hf  -> "Hf"
| Ta  -> "Ta"
| W   -> "W"
| Re  -> "Re"
| Os  -> "Os"
| Ir  -> "Ir"
| Pt  -> "Pt"
| Au  -> "Au"
| Hg  -> "Hg"
| Tl  -> "Tl"
| Pb  -> "Pb"
| Bi  -> "Bi"
| Po  -> "Po"
| At  -> "At"
| Rn  -> "Rn"
| Fr  -> "Fr"
| Ra  -> "Ra"
| Ac  -> "Ac"
| Rf  -> "Rf"
| Db  -> "Db"
| Sg  -> "Sg"
| Bh  -> "Bh"
| Hs  -> "Hs"
| Mt  -> "Mt"
| Ds  -> "Ds"
| Rg  -> "Rg"
| Cn  -> "Cn"
| Nh  -> "Nh"
| Fl  -> "Fl"
| Mc  -> "Mc"
| Lv  -> "Lv"
| Ts  -> "Ts"
| Og  -> "Og"
| Ce  -> "Ce"
| Pr  -> "Pr"
| Nd  -> "Nd"
| Pm  -> "Pm"
| Sm  -> "Sm"
| Eu  -> "Eu"
| Gd  -> "Gd"
| Tb  -> "Tb"
| Dy  -> "Dy"
| Ho  -> "Ho"
| Er  -> "Er"
| Tm  -> "Tm"
| Yb  -> "Yb"
| Lu  -> "Lu"
| Th  -> "Th"
| Pa  -> "Pa"
| U   -> "U"
| Np  -> "Np"
| Pu  -> "Pu"
| Am  -> "Am"
| Cm  -> "Cm"
| Bk  -> "Bk"
| Cf  -> "Cf"
| Es  -> "Es"
| Fm  -> "Fm"
| Md  -> "Md"
| No  -> "No"
| Lr  -> "Lr"


let to_long_string = function
| X   -> "Ghost"
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
| Rb  -> "Rubidium"
| Sr  -> "Strontium"
| Y   -> "Yttrium"
| Zr  -> "Zirconium"
| Nb  -> "Niobium"
| Mo  -> "Molybdenum"
| Tc  -> "Technetium"
| Ru  -> "Ruthenium"
| Rh  -> "Rhodium"
| Pd  -> "Palladium"
| Ag  -> "Silver"
| Cd  -> "Cadmium"
| In  -> "Indium"
| Sn  -> "Tin"
| Sb  -> "Antimony"
| Te  -> "Tellurium"
| I   -> "Iodine"
| Xe  -> "Xenon"
| Cs  -> "Cesium"
| Ba  -> "Barium"
| La  -> "Lanthanum"
| Ce  -> "Cerium"
| Pr  -> "Praseodymium"
| Nd  -> "Neodymium"
| Pm  -> "Promethium"
| Sm  -> "Samarium"
| Eu  -> "Europium"
| Gd  -> "Gadolinium"
| Tb  -> "Terbium"
| Dy  -> "Dysprosium"
| Ho  -> "Holmium"
| Er  -> "Erbium"
| Tm  -> "Thulium"
| Yb  -> "Ytterbium"
| Lu  -> "Lutetium"
| Hf  -> "Hafnium"
| Ta  -> "Tantalum"
| W   -> "Tungsten"
| Re  -> "Rhenium"
| Os  -> "Osmium"
| Ir  -> "Iridium"
| Pt  -> "Platinum"
| Au  -> "Gold"
| Hg  -> "Mercury"
| Tl  -> "Thallium"
| Pb  -> "Lead"
| Bi  -> "Bismuth"
| Po  -> "Polonium"
| At  -> "Astatine"
| Rn  -> "Radon"
| Fr  -> "Francium"
| Ra  -> "Radium"
| Ac  -> "Actinium"
| Th  -> "Thorium"
| Pa  -> "Protactinium"
| U   -> "Uranium"
| Np  -> "Neptunium"
| Pu  -> "Plutonium"
| Am  -> "Americium"
| Cm  -> "Curium"
| Bk  -> "Berkelium"
| Cf  -> "Californium"
| Es  -> "Einsteinium"
| Fm  -> "Fermium"
| Md  -> "Mendelevium"
| No  -> "Nobelium"
| Lr  -> "Lawrencium"
| Rf  -> "Rutherfordium"
| Db  -> "Dubnium"
| Sg  -> "Seaborgium"
| Bh  -> "Bohrium"
| Hs  -> "Hassium"
| Mt  -> "Meitnerium"
| Ds  -> "Darmstadtium"
| Rg  -> "Roentgenium"
| Cn  -> "Copernicium"
| Nh  -> "Nihonium"
| Fl  -> "Flerovium"
| Mc  -> "Moscovium"
| Lv  -> "Livermorium"
| Ts  -> "Tennessine"
| Og  -> "Oganesson"

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
  | Cs  -> 55
  | Ba  -> 56
  | La  -> 57
  | Ce  -> 58
  | Pr  -> 59
  | Nd  -> 60
  | Pm  -> 61
  | Sm  -> 62
  | Eu  -> 63
  | Gd  -> 64
  | Tb  -> 65
  | Dy  -> 66
  | Ho  -> 67
  | Er  -> 68
  | Tm  -> 69
  | Yb  -> 70
  | Lu  -> 71
  | Hf  -> 72
  | Ta  -> 73
  | W   -> 74
  | Re  -> 75
  | Os  -> 76
  | Ir  -> 77
  | Pt  -> 78
  | Au  -> 79
  | Hg  -> 80
  | Tl  -> 81
  | Pb  -> 82
  | Bi  -> 83
  | Po  -> 84
  | At  -> 85
  | Rn  -> 86
  | Fr  -> 87
  | Ra  -> 88
  | Ac  -> 89
  | Th  -> 90
  | Pa  -> 91
  | U   -> 92
  | Np  -> 93
  | Pu  -> 94
  | Am  -> 95
  | Cm  -> 96
  | Bk  -> 97
  | Cf  -> 98
  | Es  -> 99
  | Fm  -> 100
  | Md  -> 101
  | No  -> 102
  | Lr  -> 103
  | Rf  -> 104
  | Db  -> 105
  | Sg  -> 106
  | Bh  -> 107
  | Hs  -> 108
  | Mt  -> 109
  | Ds  -> 110
  | Rg  -> 111
  | Cn  -> 112
  | Nh  -> 113
  | Fl  -> 114
  | Mc  -> 115
  | Lv  -> 116
  | Ts  -> 117
  | Og  -> 118
  in Charge.of_int result


let of_charge c = match (Charge.to_int c) with
|  0    -> X
|  1    -> H
|  2    -> He
|  3    -> Li
|  4    -> Be
|  5    -> B
|  6    -> C
|  7    -> N
|  8    -> O
|  9    -> F
|  10   -> Ne
|  11   -> Na
|  12   -> Mg
|  13   -> Al
|  14   -> Si
|  15   -> P
|  16   -> S
|  17   -> Cl
|  18   -> Ar
|  19   -> K
|  20   -> Ca
|  21   -> Sc
|  22   -> Ti
|  23   -> V
|  24   -> Cr
|  25   -> Mn
|  26   -> Fe
|  27   -> Co
|  28   -> Ni
|  29   -> Cu
|  30   -> Zn
|  31   -> Ga
|  32   -> Ge
|  33   -> As
|  34   -> Se
|  35   -> Br
|  36   -> Kr
|  37   -> Rb
|  38   -> Sr
|  39   -> Y
|  40   -> Zr
|  41   -> Nb
|  42   -> Mo
|  43   -> Tc
|  44   -> Ru
|  45   -> Rh
|  46   -> Pd
|  47   -> Ag
|  48   -> Cd
|  49   -> In
|  50   -> Sn
|  51   -> Sb
|  52   -> Te
|  53   -> I
|  54   -> Xe
|  55   -> Cs
|  56   -> Ba
|  57   -> La
|  58   -> Ce
|  59   -> Pr
|  60   -> Nd
|  61   -> Pm
|  62   -> Sm
|  63   -> Eu
|  64   -> Gd
|  65   -> Tb
|  66   -> Dy
|  67   -> Ho
|  68   -> Er
|  69   -> Tm
|  70   -> Yb
|  71   -> Lu
|  72   -> Hf
|  73   -> Ta
|  74   -> W
|  75   -> Re
|  76   -> Os
|  77   -> Ir
|  78   -> Pt
|  79   -> Au
|  80   -> Hg
|  81   -> Tl
|  82   -> Pb
|  83   -> Bi
|  84   -> Po
|  85   -> At
|  86   -> Rn
|  87   -> Fr
|  88   -> Ra
|  89   -> Ac
|  90   -> Th
|  91   -> Pa
|  92   -> U
|  93   -> Np
|  94   -> Pu
|  95   -> Am
|  96   -> Cm
|  97   -> Bk
|  98   -> Cf
|  99   -> Es
|  100  -> Fm
|  101  -> Md
|  102  -> No
|  103  -> Lr
|  104  -> Rf
|  105  -> Db
|  106  -> Sg
|  107  -> Bh
|  108  -> Hs
|  109  -> Mt
|  110  -> Ds
|  111  -> Rg
|  112  -> Cn
|  113  -> Nh
|  114  -> Fl
|  115  -> Mc
|  116  -> Lv
|  117  -> Ts
|  118  -> Og
| x -> raise (ElementError ("Element of charge "^(string_of_int x)^" unknown"))


let covalent_radius x =
  let result = function
  | X   -> 0.
  | H   -> 0.31
  | He  -> 0.28
  | Li  -> 1.28
  | Be  -> 0.96
  | B   -> 0.85
  | C   -> 0.76
  | N   -> 0.71
  | O   -> 0.66
  | F   -> 0.57
  | Ne  -> 0.58
  | Na  -> 1.66
  | Mg  -> 1.41
  | Al  -> 1.21
  | Si  -> 1.11
  | P   -> 1.07
  | S   -> 1.05
  | Cl  -> 1.02
  | Ar  -> 1.06
  | K   -> 2.03
  | Ca  -> 1.76
  | Sc  -> 1.70
  | Ti  -> 1.60
  | V   -> 1.53
  | Cr  -> 1.39
  | Mn  -> 1.39
  | Fe  -> 1.32
  | Co  -> 1.26
  | Ni  -> 1.24
  | Cu  -> 1.32
  | Zn  -> 1.22
  | Ga  -> 1.22
  | Ge  -> 1.20
  | As  -> 1.19
  | Se  -> 1.20
  | Br  -> 1.20
  | Kr  -> 1.16
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
  | Cs  -> 2.44
  | Ba  -> 2.15
  | La  -> 2.07
  | Ce  -> 2.04
  | Pr  -> 2.03
  | Nd  -> 2.01
  | Pm  -> 1.99
  | Sm  -> 1.98
  | Eu  -> 1.98
  | Gd  -> 1.96
  | Tb  -> 1.94
  | Dy  -> 1.92
  | Ho  -> 1.92
  | Er  -> 1.89
  | Tm  -> 1.90
  | Yb  -> 1.87
  | Lu  -> 1.87
  | Hf  -> 1.75
  | Ta  -> 1.70
  | W   -> 1.62
  | Re  -> 1.51
  | Os  -> 1.44
  | Ir  -> 1.41
  | Pt  -> 1.36
  | Au  -> 1.36
  | Hg  -> 1.32
  | Tl  -> 1.45
  | Pb  -> 1.46
  | Bi  -> 1.48
  | Po  -> 1.40
  | At  -> 1.50
  | Rn  -> 1.50
  | Fr  -> 2.60
  | Ra  -> 2.21
  | Ac  -> 2.15
  | Th  -> 2.06
  | Pa  -> 2.00
  | U   -> 1.96
  | Np  -> 1.90
  | Pu  -> 1.87
  | Am  -> 1.80
  | Cm  -> 1.69
  | Bk  -> raise (ElementError "Covalent radius not defined for Bk")
  | Cf  -> raise (ElementError "Covalent radius not defined for Cf")
  | Es  -> raise (ElementError "Covalent radius not defined for Es")
  | Fm  -> raise (ElementError "Covalent radius not defined for Fm")
  | Md  -> raise (ElementError "Covalent radius not defined for Md")
  | No  -> raise (ElementError "Covalent radius not defined for No")
  | Lr  -> raise (ElementError "Covalent radius not defined for Lr")
  | Rf  -> raise (ElementError "Covalent radius not defined for Rf")
  | Db  -> raise (ElementError "Covalent radius not defined for Db")
  | Sg  -> raise (ElementError "Covalent radius not defined for Sg")
  | Bh  -> raise (ElementError "Covalent radius not defined for Bh")
  | Hs  -> raise (ElementError "Covalent radius not defined for Hs")
  | Mt  -> raise (ElementError "Covalent radius not defined for Mt")
  | Ds  -> raise (ElementError "Covalent radius not defined for Ds")
  | Rg  -> raise (ElementError "Covalent radius not defined for Rg")
  | Cn  -> raise (ElementError "Covalent radius not defined for Cn")
  | Nh  -> raise (ElementError "Covalent radius not defined for Nh")
  | Fl  -> raise (ElementError "Covalent radius not defined for Fl")
  | Mc  -> raise (ElementError "Covalent radius not defined for Mc")
  | Lv  -> raise (ElementError "Covalent radius not defined for Lv")
  | Ts  -> raise (ElementError "Covalent radius not defined for Ts")
  | Og  -> raise (ElementError "Covalent radius not defined for Og")
  in
  Units.angstrom_to_bohr *. (result x)
  |> Positive_float.of_float


let vdw_radius x =
  let result = function
  | X   -> Some 0.
  | H   -> Some 1.20
  | He  -> Some 1.40
  | Li  -> Some 1.82
  | Be  -> None
  | B   -> None
  | C   -> Some 1.70
  | N   -> Some 1.55
  | O   -> Some 1.52
  | F   -> Some 1.47
  | Ne  -> Some 1.54
  | Na  -> Some 2.27
  | Mg  -> Some 1.73
  | Al  -> Some 1.94
  | Si  -> Some 2.10
  | P   -> Some 1.80
  | S   -> Some 1.80
  | Cl  -> Some 1.75
  | Ar  -> Some 1.88
  | K   -> Some 2.75
  | Ca  -> None
  | Sc  -> None
  | Ti  -> None
  | V   -> Some 1.98
  | Cr  -> Some 1.94
  | Mn  -> Some 1.93
  | Fe  -> Some 1.93
  | Co  -> Some 1.92
  | Ni  -> Some 1.63
  | Cu  -> Some 1.40
  | Zn  -> Some 1.39
  | Ga  -> Some 1.87
  | Ge  -> None
  | As  -> Some 1.85
  | Se  -> Some 1.90
  | Br  -> Some 1.85
  | Kr  -> Some 2.02
  | Rb  -> Some 3.03
  | Sr  -> Some 2.49
  | Y   -> None
  | Zr  -> None
  | Nb  -> None
  | Mo  -> None
  | Tc  -> None
  | Ru  -> None
  | Rh  -> None
  | Pd  -> Some 1.63
  | Ag  -> Some 1.72
  | Cd  -> Some 1.58
  | In  -> Some 1.93
  | Sn  -> Some 2.17
  | Sb  -> Some 2.06
  | Te  -> Some 2.06
  | I   -> Some 1.98
  | Xe  -> Some 2.16
  | Cs  -> None
  | Ba  -> None
  | La  -> None
  | Ce  -> None
  | Pr  -> None
  | Nd  -> None
  | Pm  -> None
  | Sm  -> None
  | Eu  -> None
  | Gd  -> None
  | Tb  -> None
  | Dy  -> None
  | Ho  -> None
  | Er  -> None
  | Tm  -> None
  | Yb  -> None
  | Lu  -> None
  | Hf  -> None
  | Ta  -> None
  | W   -> None
  | Re  -> None
  | Os  -> None
  | Ir  -> None
  | Pt  -> Some 1.75
  | Au  -> Some 1.66
  | Hg  -> Some 1.55
  | Tl  -> Some 1.96
  | Pb  -> Some 2.02
  | Bi  -> None
  | Po  -> None
  | At  -> None
  | Rn  -> None
  | Fr  -> None
  | Ra  -> None
  | Ac  -> None
  | Th  -> None
  | Pa  -> None
  | U   -> Some 1.86
  | Np  -> None
  | Pu  -> None
  | Am  -> None
  | Cm  -> None
  | Bk  -> None
  | Cf  -> None
  | Es  -> None
  | Fm  -> None
  | Md  -> None
  | No  -> None
  | Lr  -> None
  | Rf  -> None
  | Db  -> None
  | Sg  -> None
  | Bh  -> None
  | Hs  -> None
  | Mt  -> None
  | Ds  -> None
  | Rg  -> None
  | Cn  -> None
  | Nh  -> None
  | Fl  -> None
  | Mc  -> None
  | Lv  -> None
  | Ts  -> None
  | Og  -> None
  in
  match result x with
  | Some y -> Some (Positive_float.of_float @@ Units.angstrom_to_bohr *. y )
  | None -> None


let mass x =
  let result = function
  | X   ->  0.
  | H   ->  1.0079
  | He  ->  4.002602
  | Li  ->  6.941
  | Be  ->  9.0121831
  | B   -> 10.81
  | C   -> 12.011
  | N   -> 14.0067
  | O   -> 15.9994
  | F   -> 18.998403163
  | Ne  -> 20.1797
  | Na  -> 22.98976928
  | Mg  -> 24.305
  | Al  -> 26.9815385
  | Si  -> 28.0855
  | P   -> 30.973761998
  | S   -> 32.06
  | Cl  -> 35.453
  | Ar  -> 39.948
  | K   -> 39.0983
  | Ca  -> 40.078
  | Sc  -> 44.955908
  | Ti  -> 47.867
  | V   -> 50.9415
  | Cr  -> 51.9961
  | Mn  -> 54.938044
  | Fe  -> 55.845
  | Co  -> 58.933194
  | Ni  -> 58.6934
  | Cu  -> 63.546
  | Zn  -> 65.38
  | Ga  -> 69.723
  | Ge  -> 72.630
  | As  -> 74.921595
  | Se  -> 78.971
  | Br  -> 79.904
  | Kr  -> 83.798
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
  | Cs  -> 132.90545196
  | Ba  -> 137.327
  | La  -> 138.90547
  | Ce  -> 140.116
  | Pr  -> 140.90766
  | Nd  -> 144.242
  | Pm  -> 145.
  | Sm  -> 150.36
  | Eu  -> 151.964
  | Gd  -> 157.25
  | Tb  -> 158.92535
  | Dy  -> 162.500
  | Ho  -> 164.93033
  | Er  -> 167.259
  | Tm  -> 168.93422
  | Yb  -> 173.045
  | Lu  -> 174.9668
  | Hf  -> 178.49
  | Ta  -> 180.94788
  | W   -> 183.84
  | Re  -> 186.207
  | Os  -> 190.23
  | Ir  -> 192.217
  | Pt  -> 195.084
  | Au  -> 196.966569
  | Hg  -> 200.592
  | Tl  -> 204.38
  | Pb  -> 207.2
  | Bi  -> 208.98040
  | Po  -> 209.
  | At  -> 210.
  | Rn  -> 222.
  | Fr  -> 223.
  | Ra  -> 226.
  | Ac  -> 227.
  | Th  -> 232.0377
  | Pa  -> 231.03588
  | U   -> 238.02891
  | Np  -> 237.
  | Pu  -> 244.
  | Am  -> 243.
  | Cm  -> 247.
  | Bk  -> 247.
  | Cf  -> 251.
  | Es  -> 252.
  | Fm  -> 257.
  | Md  -> 258.
  | No  -> 259.
  | Lr  -> 262.
  | Rf  -> 267.
  | Db  -> 270.
  | Sg  -> 269.
  | Bh  -> 270.
  | Hs  -> 270.
  | Mt  -> 278.
  | Ds  -> 281.
  | Rg  -> 281.
  | Cn  -> 285.
  | Nh  -> 286.
  | Fl  -> 289.
  | Mc  -> 289.
  | Lv  -> 293.
  | Ts  -> 293.
  | Og  -> 294.
  in
  result x
  |> Positive_float.of_float


