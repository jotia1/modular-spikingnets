#! /bin/bash

#SBATCH --job-name=timesd
#SBATCH --partition=batch
##SBATCH --nodelist=r730-1
#SBATCH --mail-type=end
#SBATCH --mail-user=joshua.arnold1@uqconnect.edu.au
#SBATCH --cpus-per-task=24

# Learning Rule;
#       SDVL
#       SSDVL
LEARNINGRULE="'SDVL'"

# Task;
#       JITTER          FREQUENCY
#       NUMAFFERENTS    DROPOUT
#       FGI             GRID
#       TIME            BAYES
#       GENALGO
TASK="'TIME'"
TASKSTART="150"
TASKSTEP="75"
TASKEND="450"

NOTES="'Test what impact simulation time has on network performance.'"

FGI="0.0228"
ETAMEAN="0.04"
SIMTIME=150
REPEATS=24

if [ $LEARNINGRULE = "'SDVL'" ] ; then
    ALPHA1=3
    ALPHA2=5
    BETA1=5
    BETA2=5
    ETAVAR="0.05"
else
    ALPHA1=6
    ALPHA2=6
    BETA1=16
    BETA2=16
    ETAVAR="0.04"
fi

CPUS=$SLURM_CPUS_PER_TASK
NAME="'$SLURM_JOB_NAME'"

if [ $TASK = "'BAYES'" ] || [ $TASK = "'GENALGO'" ] ; then
    ALPHA_MIN=1
    ALPHA_MAX=25
    BETA_MIN=1
    BETA_MAX=25
    FGI_MIN=0
    FGI_MAX=10
    ETA_MIN="0.001"
    ETA_MAX="0.1"
    SLURMTIME=$((1 * 9000)) # days x hours x minutes x seconds
    ARGS="$NAME $LEARNINGRULE $CPUS $SLURM_JOB_ID $TASK $ALPHA_MIN $ALPHA_MAX $BETA_MIN $BETA_MAX $ETA_MIN $ETA_MAX $FGI_MIN $FGI_MAX $SIMTIME $NOTES $SLURMTIME"
    MATLABSCRIPT="executeoptim"

else
    ARGS="$NAME $LEARNINGRULE $CPUS $SLURM_JOB_ID $TASK $TASKSTART $TASKSTEP $TASKEND $ALPHA1 $ALPHA2 $BETA1 $BETA2 $FGI $ETAMEAN $ETAVAR $SIMTIME $REPEATS $NOTES"
    MATLABSCRIPT="runtask"
fi

CARGS="${ARGS//[[:blank:]]/,}"
echo $CARGS

JOBLOG="joblog.txt"
echo "-------------------------------------------------------- $SLURM_JOB_ID" >> joblog.txt
echo $NOTES >> $JOBLOG
echo $CARGS >> $JOBLOG

srun matlab -nosplash -nodisplay -nodesktop -r "try, cd('experiments'); $MATLABSCRIPT($CARGS), catch me, fprintf('%s / %s\n',me.identifier,me.message), end, exit" >> slurm-$SLURM_JOB_ID.out

wait
