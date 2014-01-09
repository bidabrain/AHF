#!/bin/sh
#$ -N ahf
#$ -M ds381@sussex.ac.uk
#$ -m bea
#$ -cwd
##$ -pe openmpi 16
#$ -q mps_amd.q
#$ -S /bin/bash
# source modules environment:
module add sge
module add openmpi/gcc/64/1.6.5

DIR="/home/d/ds/ds381/Source/AHF-Github/bin/"
##DIR="/home/d/ds/ds381/apps/AHF/bin/"

##mpirun -np 16 ${DIR}/AHF-v1.0-075 $1
$DIR/AHF-v1.0-075 $1