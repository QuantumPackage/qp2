file=$1
grep "Ndet,E,E+PT2,pt2_minus,pt2_plus,pt2_abs=" $file | cut -d "=" -f 2 > ${file}.conv_fci_tc
