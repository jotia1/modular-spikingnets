#! /bin/bash

#SBATCH --job-name=grid
#SBATCH --partition=coursework
##SBATCH --nodelist=r730-1
#SBATCH --mail-type=end
#SBATCH --mail-user=joshua.arnold1@uqconnect.edu.au
#SBATCH --cpus-per-task=16

# Learning Rule;
#       SDVL
#       SSDVL
LEARNINGRULE="'SSDVL'"

# Task;
#       JITTER
#       FREQUENCY
#       NUMAFFERENTS
#       DROPOUT
#       FGI
#       GRID
TASK="'GRID'"
TASKSTART="1"
TASKSTEP="1"
TASKEND="625"

NOTES="'Running actual experiments for paper, corrected variable names for drp and frq'"

ALPHA1=6
ALPHA2=6
BETA1=16
BETA2=16
FGI="0.0228"
ETAMEAN="0.04"
ETAVAR="0.04"
SIMTIME=150
REPEATS=3

CPUS=$SLURM_CPUS_PER_TASK
NAME="'$SLURM_JOB_NAME'"

ARGS="$NAME $LEARNINGRULE $CPUS $SLURM_JOB_ID $TASK $TASKSTART $TASKSTEP $TASKEND $ALPHA1 $ALPHA2 $BETA1 $BETA2 $FGI $ETAMEAN $ETAVAR $SIMTIME $REPEATS $NOTES"
CARGS="${ARGS//[[:blank:]]/,}"

echo $CARGS

JOBLOG="joblog.txt"
echo "-------------------------------------------------------- $SLURM_JOB_ID" >> joblog.txt
echo $NOTES >> $JOBLOG
echo $CARGS >> $JOBLOG

srun matlab -nosplash -nodisplay -nodesktop -r "try, cd('experiments'); runtask($CARGS), catch me, fprintf('%s / %s\n',me.identifier,me.message), end, exit" >> slurm-$SLURM_JOB_ID.out

wait
