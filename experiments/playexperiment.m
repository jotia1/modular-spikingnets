function [ net, out ] = playexperiment()
%% PLAYEXPERIMENT is a script to trial running the full network
%   No guarentees of this script exist, it is purely for debugging and
%   testing purposes.

%% Currently set up for: 
%       FIXING INTEGRAL PROBLEM
% Fixed global integral (fgi) is not in fact fixed as it depends on when
% the gaussian is sampled as to how much current is applied. For a sharp
% peek perhaps only one sample (at say 60) would get sampled then 1 ms
% later the sample would be 0.001 (for 60.001 current applied). Conversely 
% a spread gaussian might have many peaks of lower value but averaging to 
% ~100 current applied. 

rng(1); 
net = getplaynet();
net.rand_seed = -1;

out = spikingnet(net);
%figure;
subplot(2, 1, 1);
plot(out.vt(4:end, 1:40)');
title('Example of fixed gloabl integral failing to account for short delays. Top is voltage response.');
legend({'N4', 'N5', 'N6'});

subplot(2, 1, 2);
% NOTE: This depends on out.debug being set in the simulator (which is
% unusual). debug should be the Iapp used. 
plot(out.debug(1:40, 4:end))
title('Iapp with an analytically calculated FGI.')
legend({'N4 - delay 1ms', 'N5- delay 2ms', 'N6 - delay 10ms'});



end

function [ net ] = getplaynet()

net = struct();

net.print_progress = true;
net.group_sizes = [3, 3];
net.N = sum(net.group_sizes);
net.rand_seed = 1;

net.delays = zeros(net.N);
net.delays(1, 4) = 1;
net.delays(2, 5) = 2;
net.delays(3, 6) = 19;
%net.delays(1:3, 4:5) = 5;
net.variance = zeros(net.N);
net.variance(1, 4) = 5;
net.variance(2, 5) = 5;
net.variance(3, 6) = 10;
%net.variance(1:3, 4:5) = 2;
net.w = zeros(net.N);
%net.w(1:3, 4:5) = 1;
net.w(1, 4) = 1;
net.w(2, 5) = 1;
net.w(3, 6) = 1;

net.fgi = 13;
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

net.dynamic_threshold = false;
net.thres_rise = 10; %[mV]
net.thres_freq = 1;  %[Hz]

net.sim_time_sec = 1;
net.seq = [0, 3, 7];
net.seq_freq_ms = 500;
net.supervised_seconds = 0;
net.supervising = false;
 
%[net.inp, net.ts] = createTwoPatternInput(net.seq, fliplr(net.seq), net.seq_freq_ms, net.sim_time_sec, net.supervised_seconds, 13);
% Set up specific data to see issue
net.inp = [1, 2, 3];
net.ts = [3, 3, 3];

net.voltages_to_save = [1 : net.N];
net.variance_to_save = [1 : net.N];
net.delays_to_save = [1 : net.N];

end