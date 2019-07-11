function [ net, out ] = condition3()
%% CONDITION3 - Replication of wright, wiles 2012 Condition 3
%   See paper for additional details.

% NOTES : Not complete -> Long run times, other things to do.

% TODO:
%   DONE - Generate 40 patterns and order correctly
%   - Each distance is tried 30 times and averaged
%   


rng(1); 
net = getcondition3Net();
net.rand_seed = -1;

out = spikingnet(net);
%plot3inp(out);

neurons_to_plot = [11 : 50];

figure;
% subplot(2, 2, 1);
% plottrace(out.delayst, neurons_to_plot);
% title('delays trace');
% 
% subplot(2, 2, 2);
% plottrace(out.vart, neurons_to_plot);
% title('variance trace');

subplot(2, 1, 1);
plottrace(out.vt, neurons_to_plot);
title('voltage');

subplot(2, 1, 2);
rasterspiketimes(out.spike_time_trace, net.group_sizes(1), net.group_sizes(2));
title('Neuron raster plot (compare with reponse time plot)');
% 
% figure;
% plotresponsetime(net, out, neurons_to_plot);
% title('Neuron response time');

% figure;
% subplot(2, 1, 1);
% plot(out.debug(:, 1:3));
% title('peaks');
% 
% subplot(2, 1, 2);
% plot(out.debug(:, 4));
% title('Iapp');

end

function [ net ] = getcondition3Net()

net = struct();

net.print_progress = true;
net.group_sizes = [10, 40];
net.N = sum(net.group_sizes);
net.rand_seed = 1;

n_inp = net.group_sizes(1);
net.delays = zeros(net.N);
net.delays(1:n_inp, n_inp+1:end) = rand(n_inp, net.group_sizes(2)) * 8 + 2;  % U(2, 10)
net.variance = zeros(net.N);
net.variance(1:n_inp, n_inp+1:end) = rand(n_inp, net.group_sizes(2)) * 2 + 1;  % U(2, 10)
net.w = zeros(net.N);
net.w(1:n_inp, n_inp+1:end) = 1;

net.fgi = 6.754; %6.754;
net.nu = 0.0315;
net.nv = 0.0218;
net.a1 = 4;
net.a2 = 2;
net.b1 = 1;
net.b2 = 1; 
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

net.sim_time_sec = ceil(80 * 40 / 5);
net.seq = [];
net.seq_freq_ms = 200;
net.supervised_seconds = 0;
net.supervising = false;
net.lateral_inhibition = [11:50];

net.use_simulated_annealing = true;
net.If = 3.226; %4.026;
net.Tf = 514.6;

% +1 to account for [0, 12], "No zeros please" - Matlab (2019)
num_patts = net.group_sizes(2);
seqs = generatesequences(10, num_patts, 12) + 1; 
gap = net.seq_freq_ms * num_patts;
inp = [];
ts = [];
for i = 1 : num_patts
    seq = seqs(i, :);
    [seq_inp, seq_ts ] = createInput(seq, gap, net.sim_time_sec, 0, 14);
    inp = [inp, seq_inp];
    ts = [ts, seq_ts + (i - 1) * net.seq_freq_ms];
end
[net.inp, net.ts] = sortspiketimes(inp, ts);

net.voltages_to_save = [1 : net.N];
net.variance_to_save = [1:3];
net.delays_to_save = [1:3];
net.v_thres_to_save = [];

end