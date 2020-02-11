function [ net ] = defaultnetwork( )

net = struct();

net.print_progress = true;
net.group_sizes = [2000, 1];
N_inp = net.group_sizes(1);
net.N = sum(net.group_sizes);
net.rand_seed = -1;

net.delays = zeros(net.N);
net.delays(1:N_inp, N_inp+1:end) = 5;
net.variance = zeros(net.N);
net.variance(1:N_inp, N_inp+1:end) = 2;
net.w = zeros(net.N);
net.w(1:N_inp, N_inp+1:end) = 1;

net.fgi = 0.0228;
net.nu = 0.04;
net.nv = 0.04;
net.a1 = 3;
net.a2 = 20;
net.b1 = 20;
net.b2 = 20;
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

net.use_izhikevich = false;
net.fixed_integrals = true;
net.dynamic_threshold = false;
net.thres_rise = 10; %[mV]
net.thres_freq = 1;  %[Hz]
net.lateral_inhibition = [];
net.use_simulated_annealing = false;
net.If = 3.226; %4.026;
net.Tf = 514.6;
net.dv_mem_max = 1;
net.thres_lr = 0.021;

net.sim_time_sec = 150;
net.test_seconds = 50;
net.supervising = false;
 
net.Tp = 50;
net.Np = 500;
net.Df = 10;
net.Pf = 5;
net.dropout = 0.0;
net.pattfun = [];
[net.pinp, net.pts] =  generateuniformpattern( net.Tp, net.Np );
net.data_generator = @() balancedpoisson(net.Tp, net.Df, net.group_sizes(1), net.Np, net.Pf, net.pinp, net.pts, net.pattfun, net.dropout);

net.voltages_to_save = [2001];
net.variance_to_save = [];
net.delays_to_save = [];
net.v_thres_to_save = [];
net.iapp_to_save = [2001];

end
