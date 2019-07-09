
rng(1);
SIM_TIME = 60;
[orig_inp, orig_ts, labels] = mnist2input(SIM_TIME);

fitnessfunction = @(x) optimisenetwork(x, orig_inp, orig_ts, labels, SIM_TIME );
numvariables = 9; 

opts = optimoptions(@ga,'PlotFcn',{@gaplotbestf,@gaplotstopping});
opts.PopulationSize = 15;

lb = [1; 0.1; 0.1; 1; 1; 1; 1; 1; 1];
ub = [10; 0.5; 0.5; 9; 9; 9; 9; 10; Inf];
[x,Fval,exitFlag,Output] = ga(fitnessfunction,numvariables,[],[],[],[], lb, ub,[],opts);



% FitnessFunction = @rastriginsfcn;
% numberOfVariables = 2;
% 
% opts = optimoptions(@ga,'PlotFcn',{@gaplotbestf,@gaplotstopping}, 'SelectionFcn',@selectiontournament, ...
%                         'FitnessScalingFcn',@fitscalingprop);
% opts.PopulationSize = 10;
% opts.InitialPopulationRange = [-1 0; 1 2];
% 
% opts = optimoptions(opts,'MaxGenerations',150,'MaxStallGenerations', 100);
% 
% [x,Fval,exitFlag,Output] = ga(FitnessFunction,numberOfVariables,[],[],[],[],[],[],[],opts);
% 

fprintf('The number of generations was : %d\n', Output.generations);
fprintf('The number of function evaluations was : %d\n', Output.funccount);
fprintf('The best function value found was : %g\n', Fval);


function [ err ] = optimisenetwork( x, orig_inp, orig_ts, labels, SIM_TIME )

    %rng(1);  %Set earlier
    net = getcpsSDVLmnist(SIM_TIME);
    net.rand_seed = -1;
    net.inp = orig_inp;
    net.ts = orig_ts;
    
    net.fgi = x(1);
    net.a1 = x(4); net.a2 = x(5);
    net.b1 = x(6); net.b2 = x(7);
    net.nu = x(2); net.nv = x(3);
    net.If = x(8);
    net.Tf = x(9);
    
    out = spikingnet(net);

    err = 1 - percentagecorrect(net, out, labels);
    %acc = -calcAccuracy(net, out);
    %acc = -totalspikes(net, out);

end


function [ net ] = getcpsSDVLmnist(sim_time_sec)

net = struct();

net.print_progress = false;
net.group_sizes = [28, 2];
net.N = sum(net.group_sizes);
net.rand_seed = 1;

net.delays = zeros(net.N);
net.delays(1:28, 29:30) = 5 + randi(4, net.group_sizes) - 2;
net.variance = zeros(net.N);
net.variance(1:28, 29:30) = 2  + (rand(net.group_sizes) - 0.5) * 2;
net.w = zeros(net.N);
net.w(1:28, 29:30) = 1;

net.fgi = 7.75;
net.nu = 0.0315;
net.nv = 0.0218;
net.a1 = 2;%-33.2;
net.a2 = 1; %35.6;
net.b1 = 5; %24.9;
net.b2 = 5; %82.3;
net.variance_min = 0.1;
net.variance_max = 10;
net.neuron_tau = 20;
net.delay_max = 20;

net.v_rest = -65;
net.v_reset = -70;
net.v_thres = -55;
net.w_max = 10;
net.taupre = 20;
net.taupost = 20;
net.Apost = 0;%0.1;
net.Apre = 0;%-0.1;

net.use_izhikevich = true;
net.fixed_integrals = true;
net.dynamic_threshold = false;
net.thres_rise = 10; %[mV]
net.thres_freq = 1;  %[Hz]
net.lateral_inhibition = [29, 30];

net.use_sim_annealing = true;
net.If = 3.226; %4.026;
net.Tf = 514.6;

net.sim_time_sec = sim_time_sec;
net.seq = [];
net.seq_freq_ms = 500;
net.supervised_seconds = 0;
net.supervising = false;

    
net.voltages_to_save = [1 : net.N];
net.variance_to_save = [1 : net.N];
net.delays_to_save = [1 : net.N];
net.v_thres_to_save = [];

end
