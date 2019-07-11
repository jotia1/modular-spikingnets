function [ net, out ] = condition1()
%% CONDITION1 - Replication of wright, wiles 2012 Condition 1
%   See paper for additional details.
%
%   NOTE: If you just look at the raster plot of this experiment it looks
%   good. It is important to also consider the output response time (which
%   paints a grim picture (response time increasing, learning is diverging
%   for second pattern?). Wright, Wiles (2012) describes condition1 as only
%   being 60 seconds (as replicated here) but fig2 which is described as if
%   it is condition1 last for 300 seconds and the supervision times etc.
%   are much longer (which will likely allow learning to be properly
%   cemented). See the matlab script fig2.m for an attempt at replication
%   of figure2 from the paper. 

rng(1); 
net = getCond1Net();
net.rand_seed = -1;

out = spikingnet(net);
%plot3inp(out);
figure;
subplot(2, 2, 1);
plottrace(out.delayst, [4]);
title('delays trace');

subplot(2, 2, 2);
plottrace(out.vart, [4]);
title('variance trace');

subplot(2, 2, 3);
plottrace(out.vt, [4]);
title('voltage');

subplot(2, 2, 4);
rasterspiketimes(out.spike_time_trace, 3, 1);
title('Neuron raster plot (compare with reponse time plot)');

figure;
plotresponsetime(net, out, [4]);
title('Neuron 4s response time');

end

function [ net ] = getCond1Net()

net = struct();

net.print_progress = true;
net.group_sizes = [3, 2];
net.N = sum(net.group_sizes);
net.rand_seed = 1;

net.delays = zeros(net.N);
net.delays(1:3, 4) = 5;
net.variance = zeros(net.N);
net.variance(1:3, 4) = 2;
net.w = zeros(net.N);
net.w(1:3, 4) = 1;

net.fgi = 13;
net.nu = 0.03;
net.nv = 0.01;
net.a1 = 3;
net.a2 = 2;
net.b1 = 5;
net.b2 = 5; 
net.variance_min = 0.1;
net.variance_max = 10;
net.neuron_tau = 20;
net.delay_max = 15;

net.v_rest = -70;
net.v_reset = -65;
net.v_thres = 35;
net.w_max = 10;
net.taupre = 20;
net.taupost = 20;
net.Apost = 0;
net.Apre = 0;

net.use_izhikevich = true;
net.fixed_integrals = true;
net.dynamic_threshold = false;
net.thres_rise = 10; %[mV]
net.thres_freq = 1;  %[Hz]

net.sim_time_sec = 60;
net.seq = [0, 3, 7];
net.seq_freq_ms = 500;
net.supervised_seconds = 15;
net.supervising = true;
net.lateral_inhibition = [];
net.use_simulated_annealing = false;
net.If = 1;
net.Tf = 1;

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