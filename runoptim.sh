#! /bin/bash

#SBATCH --job-name=GAtrace
#SBATCH --partition=batch
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=joshua.arnold1@uqconnect.edu.au
#SBATCH --cpus-per-task=12

# Alternatives; 
#       1:SSDVL BAYESOPT 
#       2:SDVL BAYESOPT 
#       3:SSDVL GA
SCRIPT="3"

NOTES="'Do a reasonable GA run to see if logging population matters'"

ALPHA_MIN=1
ALPHA_MAX=25
BETA_MIN=1
BETA_MAX=25
FGI_MIN=0
FGI_MAX=10
ETA_MIN="0.001"
ETA_MAX="0.1"
SIM_TIME=150

CPUS=$SLURM_CPUS_PER_TASK
NAME="'$SLURM_JOB_NAME'"
SLURMTIME=$((1 * 9000)) # days x hours x minutes x seconds

ARGS="$SCRIPT $ALPHA_MIN $ALPHA_MAX $BETA_MIN $BETA_MAX $ETA_MIN $ETA_MAX $FGI_MIN $FGI_MAX $SIM_TIME $CPUS $NAME $SLURM_JOB_ID $NOTES $SLURMTIME"
CARGS="${ARGS//[[:blank:]]/,}"

echo $CARGS

JOBLOG="joblog.txt"
echo "-------------------------------------------------------- $SLURM_JOB_ID" >> joblog.txt
echo $NOTES >> $JOBLOG
echo $CARGS >> $JOBLOG

srun matlab -nosplash -nodisplay -nodesktop -r "try, cd('experiments'); executeoptim($CARGS), catch me, fprintf('%s / %s\n',me.identifier,me.message), end, exit" >> slurm-$SLURM_JOB_ID.out



wait
