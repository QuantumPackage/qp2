#!/bin/bash
# This script is used by the MOs module, and should not be used by users.
# It copies the EZFIO/mo_basis directory in the save/EZFIO/mo_basis/xxx
# directory, where xxx is the corresponding mo_label.
# Wed Apr  2 14:35:15 CEST 2014


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

if [[ ! -f "${EZFIO}/mo_basis/mo_label" ]] ; then
  LABEL='no_label'
else
  LABEL=$(head -1 "${EZFIO}/mo_basis/mo_label" | xargs) #xargs trims the result
fi

DESTINATION="save/mo_basis/${LABEL}"

cd "${EZFIO}"

if [[ ! -d save/mo_basis ]] ; then
        mkdir -p save/mo_basis
fi

BACKUP="${DESTINATION}.old"
if [[ -d "${BACKUP}" ]] ; then
        rm -rf "${BACKUP}"
fi

if [[ -d "${DESTINATION}" ]] ; then
        mv "${DESTINATION}" "${BACKUP}"
fi

cp -r mo_basis "${DESTINATION}"

