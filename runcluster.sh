#! /bin/bash

#SBATCH --job-name=IPdrp
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
#       JITTER          FREQUENCY
#       NUMAFFERENTS    DROPOUT
#       FGI             GRID
#       TIME            BAYES
#       GENALGO         CUSTOM
TASK="'DROPOUT'"
TASKSTART="0.021"
TASKSTEP="0.01"
TASKEND="0.022"

NOTES="'Repeat experiments with new IP'"

FGI="0.0228"
ETAMEAN="0.04"
SIMTIME=150
REPEATS=16

if [ $LEARNINGRULE = "'SDVL'" ] ; then
    ALPHA1=3
    ALPHA2=5
    BETA1=5
    BETA2=5
    ETAVAR="0.05"
else  # If sSDVL 
    ALPHA1=3
    ALPHA2=20
    BETA1=20
    BETA2=20
    ETAVAR="0.04"
fi

CPUS=$SLURM_CPUS_PER_TASK
NAME="'$SLURM_JOB_NAME'"

if [ $TASK = "'JITTER'" ] ; then
    TASKSTART="0"
    TASKSTEP="2"
    TASKEND="20"
elif [ $TASK = "'NUMAFFERENTS'" ] ; then
    TASKSTART="0"
    TASKSTEP="100"
    TASKEND="1000"
elif [ $TASK = "'FREQUENCY'" ] ; then
    TASKSTART="1"
    TASKSTEP="1"
    TASKEND="10"
elif [ $TASK = "'DROPOUT'" ] ; then
    TASKSTART="0.0"
    TASKSTEP="0.1"
    TASKEND="1.0"
elif [ $TASK = "'FGI'" ] ; then
    TASKSTART="0.0222"
    TASKSTEP="0.0001"
    TASKEND="0.0234"
fi

if [ $TASK = "'BAYES'" ] || [ $TASK = "'GENALGO'" ] ; then
    ALPHA_MIN=1
    ALPHA_MAX=25
    BETA_MIN=1
    BETA_MAX=25
    FGI_MIN=0
    FGI_MAX=10
    ETA_MIN="0.001"
    ETA_MAX="0.1"
    SLURMTIME=$((4 * 24 * 60 * 60)) # days x hours x minutes x seconds
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
