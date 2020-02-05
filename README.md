# modular-spikingnets

## basic structure
The core file is simulator/spikingnet.m this file forms the simulator for which all other files are support structures. 
A network is defined as a net struct (see networks folder for examples) and then passed to spikingnet.m to run a single simulation.
Multiple simulations forming an experiment can be defined using experiments/runtask.m which is essentially is just a for loop around a call to spikingnet.m responsible for keeping a set of experiments together. 
When a simulation is done spikingnet.m returns a struct called out which can be analysed further using the functions in 
analysis/ or visualised with the functions in visualisation/.
Currently only randomly generated patterns embedded in poisson noise is supported but all associated code and future datasources will be stored in the datasources/ folder. 
Miscellaneous other support tools can be found in tools.
To run this experiment on a slurm enable compute cluster see runcluster.sh. 
Some additional optimisation algorithms and bits and pieces can be found in experiments. 

## Running a simulation
To run a simulation several pieces are needed, a network definition, and a datasource at the bare minimum.
These can be found in networks/ and datasources/ respectively, it is suggested you review the file networks/networks.md for additional information on the network structs. 
Additionally reading the code in experiments/runtask.m will show example usage for a full pipeline. 
