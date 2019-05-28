function [ value ] = gaplay( params )

net = getgaplaynet();
net.rand_seed = -1;

net.a1 = params(1);
net.a2 = params(2);
net.b1 = params(3);
net.b2 = params(4);

% Need to break symmetry in network
%net.delays(1:3, 4:5) = 5 + (rand(3, 2) - 0.5);
%net.variance(1:3, 4:5) = 2 + (rand(3, 2) - 0.5);

out = spikingnet(net);

value = - calcAccuracy(net, out);

end



function [ net ] = getgaplaynet()

net = struct();

net.print_progress = false;
net.group_sizes = [3, 2];
net.N = sum(net.group_sizes);
net.rand_seed = 1;

net.delays = zeros(net.N);
net.delays(1:3, 4:5) = 5;
net.variance = zeros(net.N);
net.variance(1:3, 4:5) = 2;
net.w = zeros(net.N);
net.w(1:3, 4:5) = 1;

net.fgi = 65;
net.nu = 0.04;
net.nv = 0.05;
net.a1 = 3;
net.a2 = 1;
net.b1 = 5;
net.b2 = 5;
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

net.sim_time_sec = 150;
net.seq = [0, 3, 7];
net.seq_freq_ms = 1000;
net.supervised_seconds = 0;
 
[net.inp, net.ts] = createTwoPatternInput(net.seq, fliplr(net.seq), net.seq_freq_ms, net.sim_time_sec, net.supervised_seconds, 13);
    
net.voltages_to_save = [1 : net.N];
net.variance_to_save = [1 : net.N];
net.delays_to_save = [1 : net.N];

end