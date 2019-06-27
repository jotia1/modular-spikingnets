function [ net, out ] = SDVLmnist()
%% SDVLmnist - This file runs Matsubaras MNIST task with SDVL 
%
%   See Matsubaras  paper to understand the task 


rng(1); 
net = getSDVLmnist();
net.rand_seed = -1;

% Need to break symmetry in network
%net.delays(1:3, 4:5) = 5 + (rand(3, 2) - 0.5);
%net.variance(1:3, 4:5) = 2 + (rand(3, 2) - 0.5);


out = spikingnet(net);

figure
rasterspiketimes(out.spike_time_trace, 28, 2);
%subplot(2, 2, 1);
%plot(out.vt');

%subplot(2, 2, 2);
%plot(squeeze(out.vart(1, :, :))');
%plot(squeeze(out.vart(2, :, :))');
%title('Variance');


%subplot(2, 2, 3);
%plot(squeeze(out.delayst(1, :, :))');
%plot(squeeze(out.delayst(2, :, :))');
%title('Delays')


disp('done');

%plot3inp(out);
%calcAccuracy(net, out)

end

function [ net ] = getSDVLmnist()

net = struct();

net.print_progress = true;
net.group_sizes = [28, 2];
net.N = sum(net.group_sizes);
net.rand_seed = 1;

net.delays = zeros(net.N);
net.delays(1:28, 29:30) = 5;
net.variance = zeros(net.N);
net.variance(1:28, 29:30) = 2;
net.w = zeros(net.N);
net.w(1:28, 29:30) = 1;

net.fgi = 7.75;
net.nu = 0.075;
net.nv = 0.1;
net.a1 = 2;%-33.2;
net.a2 = 1; %35.6;
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

net.sim_time_sec = 300;
net.seq = [0, 3, 7];
net.seq_freq_ms = 500;
net.supervised_seconds = 0;
net.supervising = false;

 
%[net.inp, net.ts] = createTwoPatternInput(net.seq, fliplr(net.seq), net.seq_freq_ms, net.sim_time_sec, net.supervised_seconds, 13);
%[p1inp, p1ts] = createTwoPatternInput(net.seq, fliplr(net.seq), net.seq_freq_ms, 120, 60, 13);
%[p2inp, p2ts] = createTwoPatternInput(fliplr(net.seq), net.seq, net.seq_freq_ms, 30, 0, 13);
%[p3inp, p3ts] = createTwoPatternInput(fliplr(net.seq), net.seq, net.seq_freq_ms, 120, 60, 13);
%[p4inp, p4ts] = createTwoPatternInput(net.seq, fliplr(net.seq), net.seq_freq_ms, 230, 0, 13);

%net.inp = [p1inp, p2inp, p3inp, p4inp];
%net.ts = [p1ts, p2ts + (120 * 1000), p3ts + (150 * 1000), p4ts + (270 * 1000)];

[net.inp, net.ts] = mnist2input(net.sim_time_sec);

    
net.voltages_to_save = [];
net.variance_to_save = [];
net.delays_to_save = [];

end