#!/bin/bash
# This is a sample PBS script
# temps CPU a ajuster au calcul
   #PBS -l cput=2000:00:00
   #PBS -l nodes=1:ppn=16
# memoire a ajuster au calcul
   #PBS -l vmem=100gb
# a changer
# Pour savoir sur quel noeud on est
#echo $HOSTNAME
# Startdir = ou sont les fichiers d'input, par defaut HOMEDIR
#
StartDir=$PBS_O_WORKDIR
echo $StartDir
#
# SCRATCHDIR = espace temporaire (local au noeud et a vider apres le calcul)
# NE PAS MODIFIER
ulimit -s unlimited
export SCRATCHDIR=/scratch/$USER/$PBS_JOBID
#
cd $StartDir


############################################################################
#### EXAMPLE OF SCRIPT TO RUN A CIPSI CALCULATION ON 5 STATES ON THE Ne^+ CATION
#### USING NATURAL ORBITALS OF A SMALL CIPSI AS MOS 
#### ALL STATES WILL HAVE THE SAME SPIN SIMETRY : A DOUBLET 

####### YOU PUT THE PATH TO YOUR 
QP_ROOT=/home_lct/eginer/programs/qp2
source ${QP_ROOT}/quantum_package.rc 
####### YOU LOAD SOME LIBRARIES 
alias python3='/programmes/installation/Python/3.7.1/bin/python3'
type -a python3

export OMP_NUM_THREADS=16

module load intel2016_OMPI-V2

source ~/programs/qp2/quantum_package.rc
./script.sh h2o dz O 1
