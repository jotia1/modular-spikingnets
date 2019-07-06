function [ net, out ] = condition1()
%% CONDITION1 - Replication of wright, wiles 2012 Condition 1
%   See paper for additional details.

rng(1); 
net = getCond1Net();
net.rand_seed = -1;

out = spikingnet(net);
plot3inp(out);
calcAccuracy(net, out)

figure;
rasterspiketimes(out.spike_time_trace, 3, 2);

disp('break');

end

function [ net ] = getCond1Net()

net = struct();

net.print_progress = true;
net.group_sizes = [3, 2];
net.N = sum(net.group_sizes);
net.rand_seed = 1;

net.delays = zeros(net.N);
net.delays(1:3, 4:5) = 5;
net.variance = zeros(net.N);
net.variance(1:3, 4:5) = 2;
net.w = zeros(net.N);
net.w(1:3, 4:5) = 1;

net.fgi = 13;
net.nu = 0.03;
net.nv = 0.01;
net.a1 = 3;%-33.2;
net.a2 = 2; %35.6;
net.b1 = 5; %24.9;
net.b2 = 5; %82.3;
net.variance_min = 0.1;
net.variance_max = 10;
net.neuron_tau = 20;
net.delay_max = 15;

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

net.sim_time_sec = 60;
net.seq = [0, 3, 7];
net.seq_freq_ms = 500;
net.supervised_seconds = 15;
net.supervising = true;
net.lateral_inhibition = [];
 
%[net.inp, net.ts] = createTwoPatternInput(net.seq, fliplr(net.seq), net.seq_freq_ms, net.sim_time_sec, net.supervised_seconds, 13);
% [p1inp, p1ts] = createTwoPatternInput(net.seq, fliplr(net.seq), net.seq_freq_ms, 120, 60, 13);
% [p2inp, p2ts] = createTwoPatternInput(fliplr(net.seq), net.seq, net.seq_freq_ms, 30, 0, 13);
% [p3inp, p3ts] = createTwoPatternInput(fliplr(net.seq), net.seq, net.seq_freq_ms, 120, 60, 13);
% [p4inp, p4ts] = createTwoPatternInput(net.seq, fliplr(net.seq), net.seq_freq_ms, 230, 0, 13);
% 
% net.inp = [p1inp, p2inp, p3inp, p4inp];
% net.ts = [p1ts, p2ts + (120 * 1000), p3ts + (150 * 1000), p4ts + (270 * 1000)];

[ inp1, ts1 ] = createInput(net.seq, net.seq_freq_ms, 25, net.supervised_seconds, 13);
[ inp2, ts2 ] = createInput(fliplr(net.seq), net.seq_freq_ms, 5, 0, 13);
[ inp3, ts3 ] = createInput(fliplr(net.seq), net.seq_freq_ms, 25, net.supervised_seconds, 13);

net.inp = [inp1, inp2, inp3];
net.ts = [ts1, ts2 + (25 * 1000), ts3 + (30 * 1000)];

net.voltages_to_save = [1 : net.N];
net.variance_to_save = [1 : net.N];
net.delays_to_save = [1 : net.N];
net.v_thres_to_save = [];

end