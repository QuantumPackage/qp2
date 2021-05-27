#!/bin/bash
# This script is used by the determinants module, and should not be used by users.
# It copies the EZFIO/determinants directory to EZFIO/save/determinants/determinants.n_det


if [[ -z ${QP_ROOT} ]] ; then
  print "The QP_ROOT environment variable is not set."
  print "Please reload the quantum_package.rc file."
  exit -1
fi

EZFIO="$1"

if [[ -z "${EZFIO}" ]] ; then
  echo "Error in $0"
  exit 1
fi

NDET=$(head -1 "${EZFIO}/determinants/n_det" | xargs printf "%09d") #xargs trims the result
NSTATE=$(head -1 "${EZFIO}/determinants/n_states" | xargs) #xargs trims the result
DESTINATION="save/determinants/determinants.${NSTATE}.${NDET}.tar"

cd "${EZFIO}"

BACKUP="${DESTINATION}.old"
if [[ -f "${BACKUP}" ]] ; then
        rm -f "${BACKUP}"
fi

if [[ -f "${DESTINATION}" ]] ; then
        mv "${DESTINATION}" "${BACKUP}"
fi
if [[ ! -d save/determinants ]] ; then
        mkdir -p save/determinants
fi

tar cf ${DESTINATION} determinants
#tar cf dets_tmp.tar determinants
#mv --backup=t dets_tmp.tar ${DESTINATION}


