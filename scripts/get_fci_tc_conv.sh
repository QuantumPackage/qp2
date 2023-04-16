file=$1
grep "Ndet,E,E+PT2,E+RPT2,|PT2|=" $file | cut -d "=" -f 2 > ${file}.conv_fci_tc
