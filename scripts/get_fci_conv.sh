file=$1
grep "# Ndet,E,E+PT2,E+RPT2,absPT2=" $file | tail -n 1 | cut -d "=" -f 2 | awk '{print $1, $3}' | column -s ' ' -t > $file.conv_fci
