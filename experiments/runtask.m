function [] = runtask(exp_name, cpus, slurm_id, task, taskstart, taskstep, taskend, repeats, vars_to_set, notes)

    assert(mod(numel(vars_to_set),2)==0, 'vars_to_set must be key value pairs, got an odd number.');

    addpath(genpath('../'));
    parpool(cpus); 

    duplicate_num = mod(slurm_id, 100);
    num_repeats = repeats;

    output_folder = sprintf('%s_%d', exp_name, slurm_id);
    mkdir(output_folder); 

    var_range = [taskstart : taskstep : taskend];
    num_range = numel(var_range);
    values = zeros(num_range, num_repeats);
    pocs = zeros(num_range, num_repeats);
    css = zeros(num_range, num_repeats);
    nmos = zeros(num_range, num_repeats);
    iss = zeros(num_range, num_repeats);
    tpxtn = zeros(num_range, num_repeats);


    parfor repeat = 1 : num_repeats

        for count = 1 : num_range
        
            net = defaultnetwork();
            net.task = task;

            %% log info
            net.run_date = datestr(datetime);
            net.var_range = var_range;
            
            net.num_repeats = num_repeats;
            net.output_folder = output_folder;
            
            try
                var = net.var_range(count);
                
                for i = 1 : 2 : numel(vars_to_set)
                    net.(vars_to_set{i}) = vars_to_set{i + 1};
                end

                net.preset_seed = duplicate_num * 17 + repeat * 19 + count * 23; 
                rng(net.preset_seed);
                net.rand_seed = -1;  % don't set inside simulator.
                
                % Set experiment variable
                net.(task) = var;
                
                net.pattfun = [];
                if strcmp(task, 'jit')
                    net.pattfun = @(pinp, pts) gaussianjitter(pinp, pts, net.jit);
                elseif strcmp(task, 'dropout')
                    net.pattfun = @(pinp, pts) patterndropoutfunc(pinp, pts, net.dropout);                
                elseif strcmp(task, 'GRID')
                    grid_size = 25;
                    [alpha, beta] = ind2sub([grid_size, grid_size], var);
                    net.a1 = alpha; net.a2 = alpha;
                    net.b1 = beta; net.b2 = beta;
                    net.nv = net.nu;
                end     
                
                [net.pinp, net.pts] = generateuniformpattern( net.Tp, net.Np );
                net.data_generator = @() balancedpoisson(net.Tp, net.Df, net.group_sizes(1), net.Np, net.Pf, net.pinp, net.pts, net.pattfun, net.dropout);
                net.repeat = repeat;
                net.count = count;


                out = spikingnet(net);


                value = trueposxtrueneg(net, out);
                out.accuracy = value;

                %  LOG multiple metrics.
                pocs(count, repeat) = percentoffsetscorrect(net, out);
                css(count, repeat) = correctspikes(net, out);
                nmos(count, repeat) = missedoffsets(net, out);
                iss(count, repeat) = incorrectspikes(net, out);
                tpxtn(count, repeat) = value;
                
                values(count, repeat) = value;
                fprintf('progress %s: count: %d, repeat: %d \n', output_folder, count, repeat);

                % Log info
                out.slurm_id = slurm_id;
                out.notes = notes;
                out.name = exp_name;

            catch exception
                err = lasterror;
                fprintf('-----------------       ===============       ERRORRR: progress %s: count: %d, repeat: %d \n\n\n%s\n\n%s\n\n', output_folder, count, repeat, getReport(exception), err.message);
            end

            try 
                filename = sprintf('%s/%s_%d_%d', output_folder, exp_name, count, repeat);

                % decrease file sizes...
                out.spike_time_trace = out.spike_time_trace(out.spike_time_trace(:, 2) == 2001, :);
                %out.w = ;
                out.timing_info = [];
                out.variance = sparse(out.variance);
                out.delays = sparse(out.delays);

                
                % Sparsify net
                net.w = sparse(net.w);
                net.delays = sparse(net.delays);
                net.variance = sparse(net.variance);

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
    save(filename, 'values', 'var_range', 'pocs', 'css', 'iss', 'nmos', 'tpxtn', '-v7.3');

end %% function runtask

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
