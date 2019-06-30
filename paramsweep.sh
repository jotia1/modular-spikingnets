#! /bin/bash
## Usage: ./paramsweep.sh
# A link in the chain to run a parameter sweep 

#MATLAB_LINE='experiments//cpsSDVLMNIST'

matlab -nosplash -nodisplay -nodesktop -r "try, run('experiments/cpsSDVLMNIST.m'), catch me, fprintf('%s / %s\n',me.identifier,me.message), end, exit" > matlab.out


