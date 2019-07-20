%% cpsSDVLMNIST - Cluster Parameter Sweep with SDVL and MNIST data
addpath(genpath('../'));

%% Constants and outputs
SIM_TIME = 300;

output_folder = newoutputfolder();

%% Parameters to sweep
r_fgi = [ 8.2, 8.3, 8.4, 8.5 ];  
r_a1 =  [ 1, 3, 5, 7 ];
r_a2 = [ 1, 3, 5, 7 ];
r_b1 = [ 1 ]    %, 1, 3, 5, 7 ];
r_b2 = [ 1, 3, 5, 7 ];

[orig_inp, orig_ts] = mnist2input(SIM_TIME);


%% Run exp
sweep_tic = tic;
for i_b2 = 1 : numel(r_b2)
    b2 = r_b2(i_b2);
    
    for i_b1 = 1 : numel(r_b1)
        b1 = r_b1(i_b1);
        
        for i_a2 = 1 : numel(r_a2)
            a2 = r_a2(i_a2);
            
            for i_a1 = 1 : numel(r_a1)
                a1 = r_b2(i_a1);
                
                for i_fgi = 1 : numel(r_fgi)
                    fgi = r_fgi(i_fgi);
                      
                    %% Save results
                    cvt_fgi = sprintf('%.1f', fgi);
                    cvt_fgi(2) = '-';
                    filename = sprintf('%s/%d_%d_%d_%d_%s', output_folder, a1, a2, b1, b2, cvt_fgi);
                    filename_old = sprintf('%d_%d_%d_%d_%s', a1, a2, b1, b2, cvt_fgi);
                    if exist(filename_old) == 2
                        continue;
                    end

                    %% Build network
                    rng(1);
                    net = getcpsSDVLmnist(SIM_TIME);
                    net.rand_seed = -1;
                    net.inp = orig_inp;
                    net.ts = orig_ts;
                    
                    net.fgi = fgi;
                    net.a1 = a1; net.a2 = a2;
                    net.b1 = b1; net.b2 = b2;
                    
                    %% Run network 
                    exp_tic = tic;
                    out = spikingnet(net);
                    
                    save(filename, 'net', 'out');
                    
                    %% log
                    exp_num = i_b2 * (5 ^ 4) + i_b1 * (5 ^ 3) + i_a2 * (5 ^ 2)  + i_a1 * 5 + i_fgi;
                    exp_num = exp_num - (5 ^ 4) - (5 ^ 3) - (5 ^ 2) - (5 ^ 1);
                    fprintf('This exp: %.3f, exp avg.: %.3f, exp num: %d, progress: %.3f \n', ...
                        toc(exp_tic), toc(sweep_tic) / exp_num, exp_num, exp_num / (5^5));
                    
                end
            end
        end
    end
end



function [ net ] = getcpsSDVLmnist(sim_time_sec)

net = struct();

net.print_progress = false;
net.group_sizes = [28, 2];
net.N = sum(net.group_sizes);
net.rand_seed = 1;

net.delays = zeros(net.N);
net.delays(1:28, 29:30) = 5 + randi(4, net.group_sizes) - 2;
net.variance = zeros(net.N);
net.variance(1:28, 29:30) = 2  + (rand(net.group_sizes) - 0.5) * 2;
net.w = zeros(net.N);
net.w(1:28, 29:30) = 1;

net.fgi = 7.75;
net.nu = 0.0315;
net.nv = 0.0218;
net.a1 = 2;%-33.2;
net.a2 = 1; %35.6;
net.b1 = 5; %24.9;
net.b2 = 5; %82.3;
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

net.fixed_integrals = false;
net.dynamic_threshold = false;
net.thres_rise = 10; %[mV]
net.thres_freq = 1;  %[Hz]
net.lateral_inhibition = [29, 30];

net.sim_time_sec = sim_time_sec;
net.seq = [];
net.seq_freq_ms = 500;
net.supervised_seconds = 0;
net.supervising = false;

    
net.voltages_to_save = [];
net.variance_to_save = [];
net.delays_to_save = [];

end
