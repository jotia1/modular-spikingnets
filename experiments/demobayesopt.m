%% play with optimisation


%
SIM_TIME=60;
fgivar = optimizableVariable('fgi', [7.0, 15.0], 'Type', 'real');
nuvar = optimizableVariable('nu', [0.02, 0.05], 'Type', 'real');
nvvar = optimizableVariable('nv', [0.02, 0.05], 'Type', 'real');
a1var = optimizableVariable('a1', [1, 9], 'Type', 'integer');
a2var = optimizableVariable('a2', [1, 9], 'Type', 'integer');
b1var = optimizableVariable('b1', [1, 9], 'Type', 'integer');
b2var = optimizableVariable('b2', [1, 9], 'Type', 'integer');
[orig_inp, orig_ts, labels] = mnist2input(SIM_TIME);

fun = @(x) optimisenetwork(x, orig_inp, orig_ts, labels, SIM_TIME);
results = bayesopt(fun, [fgivar, a1var, a2var, b1var, b2var, nuvar, nvvar],'Verbose',1,...
    'AcquisitionFunctionName','expected-improvement-plus', ...
    'NumSeedPoints', 40, 'MaxObjectiveEvaluations',30, 'UseParallel',true )
%'UseParallel',true 


% load ionosphere
% rng default
% num = optimizableVariable('n',[1,30],'Type','integer');
% dst = optimizableVariable('dst',{'chebychev','euclidean','minkowski'},'Type','categorical');
% c = cvpartition(351,'Kfold',5);
% fun = @(x)kfoldLoss(fitcknn(X,Y,'CVPartition',c,'NumNeighbors',x.n,...
%     'Distance',char(x.dst),'NSMethod','exhaustive'));
% results = bayesopt(fun,[num,dst],'Verbose',0,...
%     'AcquisitionFunctionName','expected-improvement-plus')



function [ acc ] = optimisenetwork( x, orig_inp, orig_ts, labels, SIM_TIME )

    rng(1);
    net = getcpsSDVLmnist(SIM_TIME);
    net.rand_seed = -1;
    net.inp = orig_inp;
    net.ts = orig_ts;
    
    net.fgi = x.fgi;
    net.a1 = x.a1; net.a2 = x.a2;
    net.b1 = x.b1; net.b2 = x.b2;
    net.nu = x.nu; net.nv = x.nv;
    
    out = spikingnet(net);

    acc = -percentagecorrect(net, out, labels);
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

net.fixed_integrals = false;
net.dynamic_threshold = false;
net.thres_rise = 10; %[mV]
net.thres_freq = 1;  %[Hz]
net.lateral_inhibition = [29, 30];

net.sim_time_sec = sim_time_sec;
net.seq = [];
net.seq_freq_ms = 500;
net.supervised_seconds = 0;
net.supervising = false;

    
net.voltages_to_save = [];
net.variance_to_save = [];
net.delays_to_save = [];

end