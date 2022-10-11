
input=h2o
basis=dz
EZFIO=${input}_${basis}_bi_ortho
file=${EZFIO}.tc_fci.out 
grep "Ndet,E,E+PT2,E+RPT2,|PT2|=" ${file} | cut -d "=" -f 2 > data_${EZFIO}
file=${EZFIO}.tc_fci_normal_order.out
grep "Ndet,E,E+PT2,E+RPT2=" ${file} | cut -d "=" -f 2 > data_${EZFIO}_normal

#EZFIO=${input}_${basis}_ortho
#file=${EZFIO}.tc_fci.out 
#grep "Ndet, E_tc, E+PT2 =" ${file} | cut -d "=" -f 2 > data_${EZFIO}
#file=${EZFIO}.tc_fci_normal_order.out
#grep "Ndet, E_tc, E+PT2 =" ${file} | cut -d "=" -f 2 > data_${EZFIO}_normal

#zip data_${input}_${basis}.zip data*
