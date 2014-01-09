#!/bin/sh
#$ -N ramses2gadget
#$ -M ds381@sussex.ac.uk
#$ -m bea
#$ -cwd
#$ -pe openmpi 10
#$ -q mps_gpu.q
#$ -S /bin/bash
# source modules environment:
module add sge
module add gcc/4.8.1
module add openmpi/gcc/64/1.6.5

echo "Running ramses2gadget from $1"

AHF_CONVERT_DIR="/home/d/ds/ds381/Source/AHF-Github/ahf-v1.0-075/convert"

AHF_BIN_DIR="/home/d/ds/ds381/Source/AHF-Github/ahf-v1.0-075/bin"

##$DIR/ramses2gadget -g /mnt/lustre/scratch/phys/ds381/test10_rhd/output_00181

## Convert ramses to gadget format
mpirun -np 10 $AHF_CONVERT_DIR/ramses2gadget -g $1

## Run AHF
mpirun -np 10 $AHF_BIN_DIR/AHF-v1.0-075 $2/ahf.input

## Delete Gadget data
## cd $1
## rm ramses2gadget_*