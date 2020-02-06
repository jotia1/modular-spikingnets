# How to define networks
Networks are defined in a struct, typically called net, which is passed to the spiking net simulator which produces a new strcut, typically called out, containing output information from running that network.
Typically for an experiment a base network will be stored in a getter function and a variable of interest will be manipulated.
Currently the simulator only supports single layer networks with arbitrary inputs.
In this documentation only units in the output layer are refered to as neurons, the input layer are referred to as inputs.
Much of the early testing has been done with only a single neuron and whilst effort has been taken to ensure more neurons can be added the functionality of multiple neurons is not guaranteed at this time.

For parameters listed below, a value in brackets before the description is a suggested value that may (or may not) help set things in the right range.


## Setting weights, delays, and variances
Weights, variances and delays are set using the `net.w`, `net.variance`, and `net.delays` matrices respectively.
Each matrix is `NxN` where `N` is the number of inputs and neurons in a simulation.
The weight from input 3 and neuron 4 is specified at location `net.w(3, 4)`, that is, the row specifies which inputs a connection comes from and the column specifies at which neuron the connection terminates.
This is consistent for `net.w`, `net.delays`, and `net.variances`, a value of 0 indicates no connection.
The user is required to ensure all three matrices connections agree when setting values.

## Groups and dependencies
Several variables exist together and influence similar parts of the simulation.
Some variables explicitly gate whether another will be used and should be considered.


### SDVL fields
To turn SDVL off set `net.nu` and `net.nv` to 0, otherwise SDVL will occur.
The variables of SDVL learning are set through the expected:
- `net.a1` - (3) Alpha 1.
- `net.a2` - (20) Alpha 2.
- `net.b1` - (20) Beta 1.
- `net.b2` - (20) Beta 2.
- `net.nu` - (0.04) Mean learning rate.
- `net.nv` - (0.04) Variance learning rate.
- `net.fgi` - (0.0228) The Fixed Global Integral to use.
- `net.fixed_integrals` - TODO ; what is this?

For more details of the parameters see [1]. Additionally SDVL requires some less obious parameters, these include:
- `net.variance_min` - (0.1) The minimum allowed variance.
- `net.variance_max` - (10) The maximum allowed varince.
- `net.neuron_tau` - (20) TODO : what is this?
- `net.delay_max` - (20) The maximum delay allowed.


### STDP fields
To turn STDP off set `net.Apre` and `net.Apost` to 0, otherwise STDP will occur. 
The variables of STDP are set through:
- `net.w_max` - (10) The max weight allowed.
- `net.taupre` - (20) TODO - Need to check which this corresponds to in STDP paper.
- `net.taupost` - (20) TODO - As per above.
- `net.Apre` - (0.1) TODO - As per above.
- `net.Apost` - (-0.1) TODO - As per above.

For additional information on STDP see [2].


### IP fields
Intrinsic plasticity is still **under active development** and the **API is volatile**.
Intrinsic plasticity can be enables by setting `net.dynamic_threshold = true;`
Additional fields include:
- `net.thres_frq` - (5) The frequency IP should try and maintain.
- `net.thres_lr` - (0.021) The learning rate IP should used (if applicable).
- `net.thres_rise` - (nan) Used for direct control of the threshold, part of an old rule, **currently deprecated**.


### Simulated Annealing fields
An implementation to mirror the original work [1], set `net.use_simulated_annealing = true;` to use or `false` to turn off. 
Parameters are:
- `net.If` - (0.0222) The final FGI to use.
- `net.Tf` - (30) The number of seconds over which to anneal (from starting FGI to `net.If`)

Simualted annealing seems to just be linear decaying from the initial FGI to the final FGI `net.If`, not FGI is a current hence why 'I' is used.


### Neuron types
Most experiments have been conducted using LIF neurons although very early work used Izhikevich RS neurons, there is no reason they should not still work but functionality is not gauranteed.


#### LIF fields
LIF neurons are the default are are simulted using the equation ... TODO.
To use LIF neurons set `net.use_izhikevich = false;`. 
Parameters are:
- `net.v_rest` - (-65) the resting voltage of the neuron.
- `net.v_reset` - (-70) the voltage to reset a neuron to after firing (allows for a soft refactory period).
- `v_thres` - (-55) The voltage at which a spike (action potential) occurs. 


#### Izhikevich fields 
To use Izhikevich neurons set `net.use_izhikevich = true;` else will use LIF.
Parameters for Izhikevich are set inside the spikingnet.m file. At the moment only RS units are implemnted. 

## Pattern / data generation fields
TODO section. 

## Traces / data logging
TODO section.

## Other important fields
Some other important fields that exist seperate to the groups below include:
- `net.sim_time_secs` - (150) The number of sections a simulation should run for. 
- `net.print_progress` - (true) Whether the simulator should print progress information (timing information).
- `net.group_sizes` - ([2000, 1]) The size of each layer in the network, currently only a single layer of inputs and a single layer of neruons are supported. 
- `net.N` - (2001) The total number of units (input and neurons) in the simulation. Usually calculated as `sum(net.group_sizes)`.
- `net.rand_seed` - (-1) the random seed to use when running the network. A value of -1 means the random seed has already been set external to the simulation and should not be set again. Any other value will be used as the random seed.
- `net.supervising` - (false) Used in the replication of [1]. Causes the (assumes 1) output neuron to fire after a pattern presentation if the neuron has not fired recently. 


Planned changes:
- [ ] Change switching between LIF and Izhikevich so LIF seems like the default
- [ ] Clean up / simplify section on IP once rule is established
- [ ] Write down correspondance between STDP fields and the parameters metioned in [2].
- [ ] Find out what `net.neuron_tau` is... the membrane voltage of an LIF??.
- [ ] Find out what `net.fixed_integrals` does.
- [ ] Do the pattern and data generation section. 
- [ ] Do the traces / data logging section.

# References
[1] - Wright, P. W., and J. Wiles. “Learning Transmission Delays in Spiking Neural Networks: A Novel Approach to Sequence Learning Based on Spike Delay Variance.” In The 2012 International Joint Conference on Neural Networks (IJCNN), 1–8, 2012. https://doi.org/10.1109/IJCNN.2012.6252371.
[2] - Masquelier, Timothée, Rudy Guyonneau, and Simon J. Thorpe. “Spike Timing Dependent Plasticity Finds the Start of Repeating Patterns in Continuous Spike Trains.” PLOS ONE 3, no. 1 (January 2, 2008): e1377. https://doi.org/10.1371/journal.pone.0001377.

