function [ net, out ] = sdvlpattsize( params )
%% SDVLPATTSIZE - function will plot the effect of pattern size on SDVL
%   The Performance of SDVL with a varying pattern size is explored.
%
%   Example usage:
%       [ net, out ] = sdvlpattsize( [3, 1, 5, 5, 0.0236] )

net = defaultpapernetwork();

net.a1 = params(1);
net.a2 = params(2);
net.b1 = params(3);
net.b2 = params(4);
net.fgi = params(5);

%% Parameters
net.run_date = datestr(datetime);
net.sim_time_sec = 150;
net.Tp = 50;
net.Df = 10;
net.num_repeats = 50;
net.Np = 500;
net.Pf = 5;
net.fgi = 0.0234;
net.dropout = 0.0;
N_inp = net.group_sizes(1);
net.var_range = 0.0220 : 0.0001 : 0.0228;
var_range = net.var_range;
%% Simulated annealing params
net.use_simulated_annealing = false;
net.If = 0.0222;
net.Tf = 30;
net.test_seconds = 50;
duplicate_num = 3;

exp_name = sprintf('150jits8', duplicate_num);
output_folder = newoutputfolder(exp_name);
net.output_folder = output_folder;
values = zeros(numel(net.var_range), net.num_repeats);
count = 1;


for repeat = 1 : net.num_repeats
    count = 1;
    for var = net.var_range

        net.preset_seed = duplicate_num * 17 + repeat * 19 + count * 23; 
        rng(net.preset_seed);
        net.rand_seed = -1;  % don't set inside simulator.
        
        % Set experiment variable
        %net.Pf = var;
        %net.dropout = var;
        %net.Np = var;  % Num aff
        net.jit = var;
        %net.fgi = var;

        %pvariances = rand(1, net.Np) * net.pvar_max;
        %pattfun = @(pinp, pts) mixedvariancefunc(pinp, pts, pvariances);
        %pattfun = @(pinp, pts) patterndropoutfunc(pinp, pts, net.dropout);
        net.pattfun = @(pinp, pts) gaussianjitter(pinp, pts, net.jit);
        %pattfun = [];
        
        [net.pinp, net.pts] = generateuniformpattern( net.Tp, net.Np );
        net.data_generator = @() balancedpoisson(net.Tp, net.Df, N_inp, net.Np, net.Pf, net.pinp, net.pts, net.pattfun, net.dropout);
        net.repeat = repeat;
        net.count = count;
        
        out = spikingnet(net);

        %value = detectionrate(net, out)
        value = offsetaccuracy(net, out, net.Tp, net.test_seconds)
        
        values(count, repeat) = value;
        fprintf('%s: count: %d, repeat: %d \n', output_folder, count, repeat);
        try 
            filename = sprintf('%s/%s_%d_%d', output_folder, exp_name, count, repeat);
            % decrease file sizes...
            out.spike_time_trace = out.spike_time_trace(out.spike_time_trace(:, 2) == 2001, :);
            save(filename, 'net', 'out', 'count', 'repeat', '-v7.3');
        catch exception
            err = lasterror;
            fprintf('Failed to write: %s\n%s\n\n%s\n\n', filename, getReport(exception), err.message);
        end

        try
            filename = sprintf('%s/res_%s_%d_%d', output_folder, exp_name, count, repeat);
            save(filename, 'values', 'var_range', '-v7.3');
        catch exception
            fprintf('Failed to write results: %s\n%s\n\n', filename, getReport(exception));
        end
        count = count + 1;
    end
end

filename = sprintf('%s/res_%s_final', output_folder, exp_name);
save(filename, 'values', 'var_range', '-v7.3');

end
