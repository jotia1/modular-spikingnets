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
net.sim_time_sec = 50;
net.Tp = 50;
net.Df = 10;
net.num_repeats = 50;
net.Np = 1000;
net.Pf = 2;
net.fgi = 0.0236;
N_inp = net.group_sizes(1);
net.var_range = 0.0215 : 0.0001 : 0.0226;
%% Simulated annealing params
net.use_simulated_annealing = true;
net.If = 0.0222;
net.Tf = 30;
net.test_seconds = 20;

exp_name = 'salIf';
output_folder = newoutputfolder(exp_name);
net.output_folder = output_folder;
values = zeros(numel(net.var_range), net.num_repeats);
count = 1;


for repeat = 1 : net.num_repeats
    count = 1;
    for If = net.var_range
        
        net.pattfun = [];
        [net.pinp, net.pts] = generatenoise(net.Np, net.Df, net.Tp);
        net.data_generator = @() poisspattfreq(net.Tp, net.Df, N_inp, net.Np, net.Pf, net.pinp, net.pts, net.pattfun);
        net.repeat = repeat;
        net.count = count;
        
        net.rand_seed = repeat * 19 + count * 23; 
        net.If = If;
        out = spikingnet(net);

        %value = detectionrate(net, out)
        value = offsetaccuracy(net, out, net.Tp, net.test_seconds)
        
        values(count, repeat) = value;
        fprintf('%s: count: %d, repeat: %d \n', output_folder, count, repeat);
        try 
            filename = sprintf('%s/%s_%d_%d', output_folder, exp_name, count, repeat);
            save(filename, 'net', 'out', 'count', 'repeat');
        catch exception
            fprintf('Failed to write: %s\n%s\n\n', filename, getReport(exception));
        end

        try
            filename = sprintf('%s/res_%s_%d_%d', output_folder, exp_name, count, repeat);
            save(filename, 'values');
        catch exception
            fprintf('Failed to write results: %s\n%s\n\n', filename, getReport(exception));
        end
        count = count + 1;
    end
end

filename = sprintf('%s/res_%s_final', output_folder, exp_name);
save(filename, 'values');

end
