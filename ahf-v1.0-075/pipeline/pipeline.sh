#!/bin/bash

# source modules environment:
module add sge

cd $1

echo "Running from $PWD"

AHF_OUTPUT_DIR="/mnt/lustre/scratch/phys/ds381/AHF_Output"

AHF_PIPELINE_DIR="/home/d/ds/ds381/Source/AHF-Github/ahf-v1.0-075/pipeline"

#DIRS="output_00*"

DIR="output_00001"

FOO="$PWD/ramses2gadget_$DIR."

#trap "exit" INT
#for DIR in $DIRS
#do

echo "Starting job $DIR"
mkdir $AHF_OUTPUT_DIR/job_$DIR
cp $AHF_PIPELINE_DIR/ahf.input $AHF_OUTPUT_DIR/job_$DIR/

## Sub foo for bar in AHF input file

sed -i 's,foo,'"$FOO"',g' $AHF_OUTPUT_DIR/job_$DIR/ahf.input
sed -i 's,bar,'"$AHF_OUTPUT_DIR/job_$DIR/ahf_"',g' $AHF_OUTPUT_DIR/job_$DIR/ahf.input

echo "Running ramses2gadget on /mnt/lustre/scratch/phys/ds381/test10_rhd/output_00001"

## Submit job
qsub $AHF_PIPELINE_DIR/ramses2gadget.q "$PWD/$DIR"


echo "Done"