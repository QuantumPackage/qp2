#!/bin/bash
#
# This script runs a CIPSI calculation as a sequence of single CIPSI iterations.
# After each iteration, the EZFIO directory is saved.
#
# Usage: cipsi_save [EZFIO_FILE] [NDET]
#
# Example: cipsi_save file.ezfio 10000

EZ=$1
NDETMAX=$2

qp set_file ${EZ}
qp reset -d
qp set determinants read_wf true
declare -i NDET
NDET=1
while [[ ${NDET} -lt ${NDETMAX} ]]
do
        NDET=$(($NDET + $NDET))
        qp set determinants n_det_max $NDET
        qp run fci > ${EZ}.out
        NDET=$(qp get determinants n_det)
	mv ${EZ}.out ${EZ}.${NDET}.out
        cp -r ${EZ} ${EZ}.${NDET}
done

