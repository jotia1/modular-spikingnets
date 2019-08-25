%% RUNPAREXP - Run a parallel experiment 
%   Run an experiment with repeats running in parallel
%
%   Note : Saving values inside parfor not yet supported.

addpath(genpath('../../'));
parpool(24); 

duplicate_num = 3;
num_repeats = 50;

exp_name = sprintf('low450', duplicate_num);
output_folder = newoutputfolder(exp_name);
var_range = 0.0218 : 0.0002 : 0.0228;
num_range = numel(var_range);
values = zeros(num_range, num_repeats);


parfor repeat = 1 : num_repeats
    
    net = defaultpapernetwork();

    %% Parameters
    net.run_date = datestr(datetime);
    net.sim_time_sec = 450;
    net.test_seconds = 50;
    net.var_range = var_range;
    
    net.Tp = 50;
    net.Df = 10;
    net.num_repeats = num_repeats;
    net.Np = 500;
    net.Pf = 5;
    net.fgi = 0.0226;
    net.dropout = 0.0;
    N_inp = net.group_sizes(1);
    net.output_folder = output_folder;

    
    %% Simulated annealing params
    net.use_simulated_annealing = false;
    net.If = 0.0222;
    net.Tf = 30;
    
    for count = 1 : num_range
        var = net.var_range(count);

        net.preset_seed = duplicate_num * 17 + repeat * 19 + count * 23; 
        rng(net.preset_seed);
        net.rand_seed = -1;  % don't set inside simulator.
        
        % Set experiment variable
        %net.Pf = var;
        %net.dropout = var;
        %net.Np = var;  % Num aff
        %net.jit = var;
        net.fgi = var;

        net.pattfun = [];
        %pvariances = rand(1, net.Np) * net.pvar_max;
        %net.pattfun = @(pinp, pts) mixedvariancefunc(pinp, pts, pvariances);
        %net.pattfun = @(pinp, pts) patterndropoutfunc(pinp, pts, net.dropout);
        %net.pattfun = @(pinp, pts) gaussianjitter(pinp, pts, net.jit);
        
        [net.pinp, net.pts] = generateuniformpattern( net.Tp, net.Np );
        net.data_generator = @() balancedpoisson(net.Tp, net.Df, N_inp, net.Np, net.Pf, net.pinp, net.pts, net.pattfun, net.dropout);
        net.repeat = repeat;
        net.count = count;
        
        out = spikingnet(net);

        %value = detectionrate(net, out)
        value = offsetaccuracy(net, out, net.Tp, net.test_seconds)
        out.accuracy = value;
        
        values(count, repeat) = value;
        fprintf('progress %s: count: %d, repeat: %d \n', output_folder, count, repeat);
        try 
            filename = sprintf('%s/%s_%d_%d', output_folder, exp_name, count, repeat);
            % decrease file sizes...
            out.spike_time_trace = out.spike_time_trace(out.spike_time_trace(:, 2) == 2001, :);
            prog_save(filename, net, out, count, repeat);
        catch exception
            err = lasterror;
            fprintf('Failed to write: %s\n%s\n\n%s\n\n', filename, getReport(exception), err.message);
        end

        try
            filename = sprintf('%s/res_%s_%d_%d', output_folder, exp_name, count, repeat);
            %parfor_save(filename, values, var_range);
        catch exception
            fprintf('Failed to write results: %s\n%s\n\n', filename, getReport(exception));
        end
    end
end

filename = sprintf('%s/res_%s_final', output_folder, exp_name);
save(filename, 'values', 'var_range', '-v7.3');

function [] = prog_save(filename, net, out, count, repeat)
    save(filename, 'net', 'out', 'count', 'repeat', '-v7.3');
end



function [] = parfor_save(varargin)
% Taken from:
% https://au.mathworks.com/matlabcentral/answers/135285-how-do-i-use-save-with-a-parfor-loop-using-parallel-computing-toolbox
    fname=varargin{1};    
    for i=2:nargin
       eval([inputname(i),'=varargin{i};']);  
       if i==2
            save('-mat', '-v7.3', fname,inputname(i));
       else
           save('-mat','-v7.3', fname,inputname(i),'-append');
       end        
    end
end
