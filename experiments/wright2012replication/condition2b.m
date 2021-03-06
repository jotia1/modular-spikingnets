function [ net, out ] = condition2b()
%% CONDITION2B - Replication of wright, wiles 2012 Condition 2B
%   See paper for additional details.

% NOTES : Successful instantiations can be created. A technically correct
% set of parameters (for fgi and If) have been set. Original values are in
% the comments (but do not work for this simulator). The technically
% correct solution is diverging and so whilst correct (zoom in on raster
% plot for evidence), it seems this solution is not actually stable.

% TODO:
%   - Generate two patterns with a given distance
%   DONE - Reimplement simulated annealing
%   - Each distance is tried 30 times and averaged
%   


rng(1); 
net = getcondition2bNet();
net.rand_seed = -1;

out = spikingnet(net);
%plot3inp(out);

neurons_to_plot = [11, 12];

figure;
subplot(2, 2, 1);
plottrace(out.delayst, neurons_to_plot);
title('delays trace');

subplot(2, 2, 2);
plottrace(out.vart, neurons_to_plot);
title('variance trace');

subplot(2, 2, 3);
plottrace(out.vt, neurons_to_plot);
title('voltage');

subplot(2, 2, 4);
rasterspiketimes(out.spike_time_trace, net.group_sizes(1), net.group_sizes(2));
title('Neuron raster plot (compare with reponse time plot)');

figure;
plotresponsetime(net, out, neurons_to_plot);
title('Neuron response time');

% figure;
% subplot(2, 1, 1);
% plot(out.debug(:, 1:3));
% title('peaks');
% 
% subplot(2, 1, 2);
% plot(out.debug(:, 4));
% title('Iapp');

end

function [ net ] = getcondition2bNet()

net = struct();

net.print_progress = true;
net.group_sizes = [10, 2];
net.N = sum(net.group_sizes);
net.rand_seed = 1;

n_inp = net.group_sizes(1);
net.delays = zeros(net.N);
net.delays(1:n_inp, n_inp+1:end) = rand(n_inp, net.group_sizes(2)) * 8 + 2;  % U(2, 10)
net.variance = zeros(net.N);
net.variance(1:n_inp, n_inp+1:end) = rand(n_inp, net.group_sizes(2)) * 2 + 1;  % U(2, 10)
net.w = zeros(net.N);
net.w(1:n_inp, n_inp+1:end) = 1;

net.fgi = 4.552; %6.752;
net.nu = 0.0338;
net.nv = 0.0357;
net.a1 = 1;
net.a2 = 3;
net.b1 = 8;
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

net.sim_time_sec = ceil(100 * 2 / 5);
net.seq = [];
net.seq_freq_ms = 200;
net.supervised_seconds = 0;
net.supervising = false;
net.lateral_inhibition = [];

net.use_simulated_annealing = true;
net.If = 3.226; %4.026;
net.Tf = 31.4;

% +1 to account for [0, 12], "No zeros please" - Matlab (2019)
seqs = generatesequences(10, 2, 12) + 1; 
s1 = seqs(1, :);
s2 = seqs(2, :);

[net.inp, net.ts] = createTwoPatternInput(s1, s2, net.seq_freq_ms, net.sim_time_sec, net.supervised_seconds, 14);


net.voltages_to_save = [1 : net.N];
net.variance_to_save = [1 : net.N];
net.delays_to_save = [1 : net.N];
net.v_thres_to_save = [];

end