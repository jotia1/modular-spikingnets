function [] = runtask(exp_name, learning_rule, cpus, slurm_id, task, taskstart, taskstep, taskend, alpha1, alpha2, beta1, beta2, fgi, etamean, etavar, sim_time_secs, repeats, notes)



    addpath(genpath('../'));
    parpool(cpus); 

    duplicate_num = mod(slurm_id, 10);
    num_repeats = repeats;
    grid_size = 25;

    %   duplicate number,
    %   parpool, var_range
    %   patt_fun, ??? = var
    %   exp_name

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
        
            net = defaultpapernetwork();

            %% log info
            net.run_date = datestr(datetime);
            net.sim_time_sec = sim_time_secs;
            net.test_seconds = 50;
            net.var_range = var_range;

            
            net.Tp = 50;
            net.Df = 10;
            net.num_repeats = num_repeats;
            net.Np = 500;
            net.Pf = 5;
            net.fgi = fgi;
            net.dropout = 0.0;
            N_inp = net.group_sizes(1);
            net.output_folder = output_folder;

            
            %% Simulated annealing params
            net.use_simulated_annealing = false;
            net.If = 0.0222;
            net.Tf = 30;

            % Set alpha beta eta
            net.a1 = alpha1; net.a2 = alpha2;
            net.b1 = beta1; net.b2 = beta2;
            net.nu = etamean; net.nv = etavar;

        
            try
                var = net.var_range(count);

                net.preset_seed = duplicate_num * 17 + repeat * 19 + count * 23; 
                rng(net.preset_seed);
                net.rand_seed = -1;  % don't set inside simulator.
                
                % Set experiment variable
                net.pattfun = [];
                if strcmp(task, 'FGI')
                    net.fgi = var;
                elseif strcmp(task, 'JITTER')
                    net.jit = var;
                    net.pattfun = @(pinp, pts) gaussianjitter(pinp, pts, net.jit);
                elseif strcmp(task, 'FREQUENCY')
                    net.Pf = var;
                elseif strcmp(task, 'DROPOUT')
                    net.dropout = var;
                    net.pattfun = @(pinp, pts) patterndropoutfunc(pinp, pts, net.dropout);
                elseif strcmp(task, 'NUMAFFERENTS')
                    net.Np = var;
                elseif strcmp(task, 'GRID')
                    tvar = net.var_range(count);
                    [alpha, beta] = ind2sub([grid_size, grid_size], tvar);
                    net.a1 = alpha; net.a2 = alpha;
                    net.b1 = beta; net.b2 = beta;
                    net.nv = net.nu;
                elseif strcmp(task, 'TIME')
                    net.sim_time_sec = var;
                else
                    fprintf('INVALID TASK: %s\n', task); 
                    assert false
                end

                [net.pinp, net.pts] = generateuniformpattern( net.Tp, net.Np );
                net.data_generator = @() balancedpoisson(net.Tp, net.Df, N_inp, net.Np, net.Pf, net.pinp, net.pts, net.pattfun, net.dropout);
                net.repeat = repeat;
                net.count = count;

                if strcmp(learning_rule, 'SDVL')
                    out = spikingnet(net);
                elseif strcmp(learning_rule, 'SSDVL');
                    net.a2 = net.a1;
                    net.b2 = net.b1;
                net.nu = net.nv
                out = ssdvl(net);
                else
                    fprintf('INVALID LEARNING RULE: %s\n', learning_rule);
                end

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
                out.learning_rule = learning_rule;

            catch exception
                err = lasterror;
                fprintf('-----------------       ===============       ERRORRR: progress %s: count: %d, repeat: %d \n\n\n%s\n\n%s\n\n', output_folder, count, repeat, getReport    (exception), err.message);
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
