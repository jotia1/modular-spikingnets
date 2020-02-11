#! /bin/bash

#SBATCH --job-name=optimwork
#SBATCH --partition=batch
##SBATCH --nodelist=r730-0
#SBATCH --mail-type=end
#SBATCH --mail-user=joshua.arnold1@uqconnect.edu.au
#SBATCH --cpus-per-task=4

# Task;
# Task is the name of the variable to be set before each experiment.
#       jit             Pf
#       naf             dropout
#       fgi             GRID
#       sim_time_sec    BAYES
#       GENALGO         CUSTOM
#       THRESHOLDLR     STDP
#       Any other field in the net structure...
TASK="'GENALGO'"
TASKSTART="0.021"
TASKSTEP="0.005"
TASKEND="0.03"

NOTES="'Test flexible pipeline to GA'"

REPEATS=4

#VARSTOSET="{'v_thres_to_save', [2001], 'dynamic_threshold', true, 'thres_freq', 5}" 
VARSTOSET="{'sim_time_sec', 10, 'nu', 0, 'nv', 0, 'fgi', 0.0228}"

VARSTOOPTIMISE="{'Apre', 'Apost', 'taupre', 'taupost'}"
LBS="[0.1, 0.1, 1, 1]"
UBS="[2.0, 2.0, 40, 40]"

CPUS=$SLURM_CPUS_PER_TASK
NAME="'$SLURM_JOB_NAME'"

if [ $TASK = "'jit'" ] ; then
    TASKSTART="0"
    TASKSTEP="2"
    TASKEND="20"
elif [ $TASK = "'naf'" ] ; then
    TASKSTART="0"
    TASKSTEP="100"
    TASKEND="1000"
elif [ $TASK = "'Pf'" ] ; then
    TASKSTART="1"
    TASKSTEP="1"
    TASKEND="10"
elif [ $TASK = "'dropout'" ] ; then
    TASKSTART="0.0"
    TASKSTEP="0.1"
    TASKEND="1.0"
elif [ $TASK = "'fgi'" ] ; then
    TASKSTART="0.0222"
    TASKSTEP="0.0002"
    TASKEND="0.0234"
elif [ $TASK = "'thres_lr'" ] ; then
    TASKSTART="0.001"
    TASKSTEP="0.005"
    TASKEND="0.041"
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
