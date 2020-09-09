#! /bin/bash

#SBATCH --job-name=stdp
#SBATCH --partition=batch
##SBATCH --nodelist=r730-1
#SBATCH --mail-type=end
#SBATCH --mail-user=joshua.arnold1@uqconnect.edu.au
#SBATCH --cpus-per-task=1

# Task;
# Task is the name of the variable to be set before each experiment.
#       jit             Pf
#       naf             dropout
#       fgi             GRID
#       sim_time_sec    BAYES
#       GENALGO         CUSTOM
#       THRESHOLDLR     STDP
#       Any other field in the net structure...
TASK="'w_init'"
TASKSTART="1"
TASKSTEP="1"
TASKEND="6"

NOTES="'See how the value of a1 influences the performance, full range between 1-20 for 300 sec'"

REPEATS=16    #$SLURM_CPUS_PER_TASK

#VARSTOSET="{'v_thres_to_save', [2001], 'dynamic_threshold', true, 'thres_freq', 5}" 
#VARSTOSET="{'thres_lr', 0.021, 'ip_decay', 0.1, 'avg_period', 1000, 'v_thres_to_save', [2001], 'sim_time_sec', 150}"
VARSTOSET="{'thres_lr', 0.0, 'sim_time_sec', 300}"

VARSTOOPTIMISE="{}"
LBS="[]"
UBS="[]"

##    STDP vars to optimise
#VARSTOOPTIMISE="{'Apre', 'Apost', 'taupre', 'taupost', 'fgi'}"
#LBS="[0.1, 0.1, 1, 1, 0.0220]"
#UBS="[2.0, 2.0, 40, 40, 0.0236]"

CPUS=$SLURM_CPUS_PER_TASK
NAME="'$SLURM_JOB_NAME'"

if [ $TASK = "'xjit'" ] ; then
    TASKSTART="0"
    TASKSTEP="2"
    TASKEND="20"
elif [ $TASK = "'xNp'" ] ; then
    TASKSTART="0"
    TASKSTEP="100"
    TASKEND="1000"
elif [ $TASK = "'xPf'" ] ; then
    TASKSTART="1"
    TASKSTEP="1"
    TASKEND="10"
elif [ $TASK = "'xdropout'" ] ; then
    TASKSTART="0.0"
    TASKSTEP="0.1"
    TASKEND="1.0"
elif [ $TASK = "'xfgi'" ] ; then
    TASKSTART="0.0222"
    TASKSTEP="0.0002"
    TASKEND="0.0234"
#elif [ $TASK = "'thres_lr'" ] ; then
#    TASKSTART="0.001"
#    TASKSTEP="0.005"
#    TASKEND="0.041"
fi

if [ $TASK = "'BAYES'" ] || [ $TASK = "'GENALGO'" ] ; then
    SLURMTIME=$((4 * 24 * 60 * 60)) # days x hours x minutes x seconds
    ARGS="$NAME, $CPUS, $SLURM_JOB_ID, $TASK, $VARSTOSET, $VARSTOOPTIMISE, $LBS, $UBS, $NOTES, $SLURMTIME" 
    MATLABSCRIPT="executeoptim"
else
    ARGS="$NAME, $CPUS, $SLURM_JOB_ID, $TASK, $TASKSTART, $TASKSTEP, $TASKEND, $REPEATS, $VARSTOSET, $NOTES"
    MATLABSCRIPT="runtask"
fi

echo "$MATLABSCRIPT $ARGS"

printf -v date '%(%Y-%m-%d %H:%M:%S)T\n' -1
JOBLOG="joblog.txt"
echo "------------------------- $date--------------------------------------------- $SLURM_JOB_ID" >> joblog.txt
echo "$MATLABSCRIPT $ARGS" >> $JOBLOG

srun matlab -nosplash -nodisplay -nodesktop -r "try, cd('experiments'); $MATLABSCRIPT($ARGS), catch me, fprintf('%s / %s\n',me.identifier,me.message), end, exit" >> slurm-$SLURM_JOB_ID.out

wait
