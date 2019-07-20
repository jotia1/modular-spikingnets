#! /bin/bash

#SBATCH --job-name=patsize
#SBATCH --partition=batch
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=joshua.arnold1@uqconnect.edu.au

srun -n1 matlab -nosplash -nodisplay -nodesktop -r "try, run('experiments/paper/runcluster.m'), catch me, fprintf(    '%s / %s\n',me.identifier,me.message), end, exit"


wait
