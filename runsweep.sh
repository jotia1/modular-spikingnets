#! /bin/bash

#SBATCH --job-name=paramsweep
#SBATCH --partition=batch
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=joshua.arnold1@uqconnect.edu.au

JOBNAME="b7"

ARGS="--error=$JOBNAME.err --output=$JOBNAME.out"
srun -n1 $ARGS paramsweep.sh > $JOBNAME.log


wait
