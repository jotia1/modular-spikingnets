function [ net, out ] = fig2()
%% FIG2 - Replication of wright, wiles 2012 Figure 2
%   See paper for additional details.
%
%   Produces a plot similar to figure 2 from the paper, some minor
%   differences that might be worth further investigation. Two major points
%   the bumps in the delay and variance plots, what causes these? and why
%   does the output keep firing at times 121+ (when the pattern reverses).
%   I suspect these are both coming from the same issue, that the fgi is
%   too high and causes the output to fire despite a full reversal of the
%   pattern (then as it shifts delays it looses this). With a lower fgi the
%   plot is reproduced very faithfully. fgi = 12 is left as the default
%   here (13 used in paper).


rng(1); 
net = getfig2Net();
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

figure;
subplot(2, 1, 1);
plot(out.debug(:, 1:3));
title('peaks');

subplot(2, 1, 2);
plot(out.debug(:, 4));
title('Iapp');

end

function [ net ] = getfig2Net()

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

net.fgi = 12;
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

net.use_simulated_annealing = false;
net.If = 5;
net.Tf = 10;

net.sim_time_sec = 300;
net.seq = [0, 3, 7];
net.seq_freq_ms = 500;
net.supervised_seconds = 60;
net.supervising = true;
net.lateral_inhibition = [];

[ inp1, ts1 ] = createInput(net.seq, net.seq_freq_ms, 60, net.supervised_seconds, 13);
[ inp2, ts2 ] = createInput(net.seq, net.seq_freq_ms, 60, 0, 13);
[ inp3, ts3 ] = createInput(fliplr(net.seq), net.seq_freq_ms, 30, 0, 13);
[ inp4, ts4 ] = createInput(fliplr(net.seq), net.seq_freq_ms, 60, net.supervised_seconds, 13);
[ inp5, ts5 ] = createInput(fliplr(net.seq), net.seq_freq_ms, 60, 0, 13);
[ inp6, ts6 ] = createInput(net.seq, net.seq_freq_ms, 30, 0, 13);

net.inp = [inp1, inp2, inp3, inp4, inp5, inp6];
net.ts = [ts1, ts2 + (60 * 1000), ...
                ts3 + (120 * 1000), ...
                ts4 + (150 * 1000), ...
                ts5 + (210 * 1000), ...
                ts6 + (270 * 1000)];

net.voltages_to_save = [1 : net.N];
net.variance_to_save = [1 : net.N];
net.delays_to_save = [1 : net.N];
net.v_thres_to_save = [];

end