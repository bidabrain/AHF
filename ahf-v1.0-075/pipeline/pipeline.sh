#!/bin/sh
#$ -N ramses2gadget
#$ -M ds381@sussex.ac.uk
#$ -m bea
#$ -cwd
#$ -pe openmpi 32
#$ -q mps_gpu.q
#$ -S /bin/bash
# source modules environment:
module add sge
module add gcc/4.8.1
module add openmpi/gcc/64/1.6.5

## Serial pipeline for running AHF on ramses data

################################FUNCTIONS#########################

usage(){
	echo "Usage: $0 filename"
	exit 1
}
 
# define is_file_exits function 
# $f -> store argument passed to the script
is_file_exits(){
	local f="$1"
	[[ -f "$f" ]] && return 0 || return 1
}

run_ramses2gadget()
{

	local RUN_DIR=$1

	if $verbose  ;  then
	echo "Running ramses2gadget on $RUN_DIR with $NSLOTS mpi nodes"
	fi

	## Convert ramses to gadget format
	if ! $debug  ;  then
		mpirun -np $NSLOTS $AHF_CONVERT_DIR/ramses2gadget -g "$RUN_DIR"
	fi

	echo "File conversion done"
}

run_ahf()
{

	local OUT_DIR=$1

	if $verbose  ;  then
	echo "Running AHF with $NSLOTS mpi nodes, output at $OUT_DIR"
	fi

	## Run AHF
	if ! $debug  ;  then
		mpirun -np $NSLOTS $AHF_BIN_DIR/AHF-v1.0-075 $OUT_DIR/ahf.input
	fi

	echo "AHF Done"
}

################################END#########################

# invoke  usage
# call usage() function if filename not supplied
[[ $# -eq 0 ]] && usage

verbose=true
debug=false

if $debug  ;  then
	echo "Debug flag on, won't run MPI tasks"
fi

cd $1

echo "Running from $PWD"

# Declare variables

AHF_CONVERT_DIR="/home/d/ds/ds381/Source/AHF-Github/ahf-v1.0-075/convert"

AHF_BIN_DIR="/home/d/ds/ds381/Source/AHF-Github/ahf-v1.0-075/bin"

AHF_OUTPUT_DIR="/mnt/lustre/scratch/phys/ds381/AHF_Output"

AHF_PIPELINE_DIR="/home/d/ds/ds381/Source/AHF-Github/ahf-v1.0-075/pipeline"

GADGET_DIR="/mnt/lustre/scratch/phys/ds381/scratch"

DIRS="output_00*"

#DIR="output_00001"

# End declare variables

# Loop over all snapshots

trap "exit" INT
for DIR in $DIRS
do

	echo "Starting job $DIR"

	#Convert ramses -> gadget
	run_ramses2gadget $PWD/$DIR

	#Run AHF
	# ic_filename to be set in ahf.input
	FOO="$GADGET_DIR/ramses2gadget_$DIR."

	mkdir $AHF_OUTPUT_DIR/job_$DIR
	cp $AHF_PIPELINE_DIR/ahf.input $AHF_OUTPUT_DIR/job_$DIR/

	## Set ic_filename and outfile_prefix and write changes to ahf.input

	sed -i 's,foo,'"$FOO"',g' $AHF_OUTPUT_DIR/job_$DIR/ahf.input
	sed -i 's,bar,'"$AHF_OUTPUT_DIR/job_$DIR/ahf_"',g' $AHF_OUTPUT_DIR/job_$DIR/ahf.input

	run_ahf $AHF_OUTPUT_DIR/job_$DIR

	## Delete Gadget data
	echo "Deleting ramses2gadget_$DIR.* from $GADGET_DIR/"
	rm $GADGET_DIR/ramses2gadget_output_00*

done

echo "Done"
