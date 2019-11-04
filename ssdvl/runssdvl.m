


net = defaultpapernetwork();

net.sim_time_sec = 150;
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

%net.fgi = 0.035;
net.fgi = 0.0236;
net.nu = 0.04;
net.nv = net.nu;
net.a1 = 3;
%net.a2 = 1;
net.b1 = 5;
%net.b2 = 5;

net.Tp = 50;
net.Np = 25;
net.Df = 10;
net.Pf = 5;
net.pattfun = [];
[net.pinp, net.pts] =  generateuniformpattern( net.Tp, net.Np );
net.data_generator = @() balancedpoisson(net.Tp, net.Df, N_inp, net.Np, net.Pf, net.pinp, net.pts, net.pattfun);

net.voltages_to_save = [net.N];
net.variance_to_save = [];
net.delays_to_save = [];
net.v_thres_to_save = [];
net.iapp_to_save = [net.N];


out = ssdvl(net);

value = offsetaccuracy( net, out, 50, net.sim_time_sec)

plotall(net, out);