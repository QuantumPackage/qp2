file=$1
grep "N_det             =" $1 | cut -d "=" -f 2 > N_det_tmp
grep "E               =" $file | cut -d "=" -f 2 > E_tmp
grep "E+PT2           =" $file | cut -d "=" -f 2 | cut -d "+" -f 1 > E+PT2_tmp
grep "E+rPT2          =" $file | cut -d "=" -f 2 | cut -d "+" -f 1 > E+rPT2_tmp
paste N_det_tmp E_tmp E+PT2_tmp E+rPT2_tmp | column -s ' ' -t > $file.conv_fci
rm N_det_tmp E_tmp E+PT2_tmp E+rPT2_tmp
