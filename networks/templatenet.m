function [ net ] = templatenet(  )
%GETDEFAULTNET Return a default network specification as an example
%   Returns a network definition that can be passed to spikingnet.m to run
%   a spiking net simulation. The network structure and parameters are
%   specified here.

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

net.fgi = 65;
net.nu = 0.03;
net.nv = 0.01;
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

net.dynamic_threshold = true;
net.thres_rise = 10; %[mV]
net.thres_freq = 1;  %[Hz]

net.sim_time_sec = 100;
net.seq = [0, 3, 7];
net.seq_freq_ms = 500;
net.supervised_seconds = 100;
 
[net.inp, net.ts] = createInput(net.seq, net.seq_freq_ms, net.sim_time_sec, net.supervised_seconds, 13);
    
net.voltages_to_save = [1 : net.N];
net.variance_to_save = [1 : net.N];
net.delays_to_save = [1 : net.N];



end

