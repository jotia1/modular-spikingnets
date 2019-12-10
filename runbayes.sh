#! /bin/bash

#SBATCH --job-name=bayes
#SBATCH --partition=kaleen
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=joshua.arnold1@uqconnect.edu.au

#SBATCH --cpus-per-task=24

# Alternatives; 
#       1:SSDVL BAYESOPT 
#       2:SDVL BAYESOPT 
#       3:SSDVL GA
SCRIPT="1"

ALPHA_MIN=1
ALPHA_MAX=25
BETA_MIN=1
BETA_MAX=25
FGI_MIN=0
FGI_MAX=10
ETA_MIN="0.001"
ETA_MAX="0.1"
SIM_TIME=150

#echo $SLURM_JOB_ID

echo "matlab $SCRIPT $ALPHA_MIN $ALPHA_MAX $BETA_MIN $BETA_MAX $ETA_MIN $ETA_MAX $FGI_MIN $FGI_MAX $SIM_TIME"

ARGS="$SCRIPT $ALPHA_MIN $ALPHA_MAX $BETA_MIN $BETA_MAX $ETA_MIN $ETA_MAX $FGI_MIN $FGI_MAX $SIM_TIME"
CARGS="${ARGS//[[:blank:]]/,}"

srun matlab -nosplash -nodisplay -nodesktop -r "try, cd('experiments'); executebayes($CARGS), catch me, fprintf('%s / %s\n',me.identifier,me.message), end, exit" >> slurm-$SLURM_JOB_ID.out



wait
