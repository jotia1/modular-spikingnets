#! /bin/bash

#SBATCH --job-name=tpxtn
#SBATCH --partition=batch
#SBATCH --nodelist=r730-0
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=joshua.arnold1@uqconnect.edu.au

#SBATCH --cpus-per-task=3

srun matlab -nosplash -nodisplay -nodesktop -r "try, run('analysis/calculatebatchmetric.m'), catch me, fprintf(    '%s / %s\n',me.identifier,me.message), end, exit"


wait
