%% RUNPARSSDVL - Run a parallel SSDVL experiment
%   Same as runparexp but for ssdvl

addpath(genpath('../'));
parpool(16); 

duplicate_num = 3;
num_repeats = 16;


%%    GENERATE PETE DIAGONAL
%final_list = [];
%for betas = 17 : 40
%    for alphas = 1 : 13
%        alpha = betas + alphas - 7;
%        ind = sub2ind([46, 46], alpha, beta);
%        final_list(end + 1) = ind;
%    end
%end


%%     GENERATE TOP RIGHT
final_list = [];
for beta = 25 : -1 : 1;
    for alpha = 1 : beta
        %[alpha, beta]
        final_list = [final_list, sub2ind([25, 25], alpha, beta)];
    end
end

%%%    GENERATE BOT LEFT
%final_list = [];                                                                
%for beta = 20 : -1 : 17;                                                        
%    for alpha = 1 : beta                                                        
%        %[alpha, beta]                                                          
%        final_list = [final_list, sub2ind([25, 25], alpha, beta)];              
%    end                                                                         
%end 


exp_name = sprintf('ssgridhalf', duplicate_num);
notes = 'Rerun grid with fgi set to 0.0228, net not being reinitialised error fixed, tpxtn as default and last 50 sec frozen for evaluation.';
output_folder = newoutputfolder(exp_name);
grid_size = 25;
var_range = final_list; %1 : 1 : grid_size * grid_size;
num_range = numel(var_range);
values = zeros(num_range, num_repeats);
pocs = zeros(num_range, num_repeats);
css = zeros(num_range, num_repeats);
nmos = zeros(num_range, num_repeats);
iss = zeros(num_range, num_repeats);
tpxtn = zeros(num_range, num_repeats);


parfor repeat = 1 : num_repeats

    for count = 1 : num_range
	    try
            net = defaultpapernetwork();

            %% Parameters
            net.run_date = datestr(datetime);
            net.sim_time_sec = 150;
            net.test_seconds = 50;
            net.var_range = var_range;

            net.Tp = 50;
            net.Df = 10;
            net.num_repeats = num_repeats;
            net.Np = 500;
            net.Pf = 5;
            net.fgi = 0.0228;
            net.dropout = 0.0;
            N_inp = net.group_sizes(1);
            net.output_folder = output_folder;


            %% Simulated annealing params
            net.use_simulated_annealing = false;
            net.If = 0.0222;
            net.Tf = 30;
            tvar = net.var_range(count);
            [alpha, beta] = ind2sub([grid_size, grid_size], tvar);

            net.preset_seed = duplicate_num * 17 + repeat * 19 + count * 23; 
            rng(net.preset_seed);
            net.rand_seed = -1;  % don't set inside simulator.

            % Set experiment variable
            net.a1 = alpha; 
            net.a2 = net.a1;
            net.b1 = beta; 
            net.b2 = net.b1;
            net.nu = net.nv;

            net.pattfun = [];
            [net.pinp, net.pts] = generateuniformpattern( net.Tp, net.Np );
            net.data_generator = @() balancedpoisson(net.Tp, net.Df, N_inp, net.Np, net.Pf, net.pinp, net.pts, net.pattfun, net.dropout);
            net.repeat = repeat;
            net.count = count;
            net.beta = beta;

            out = ssdvl(net);

            %value = detectionrate(net, out)
            value = trueposxtrueneg(net, out);
            
            %  LOG multiple metrics.
            tpxtn(count, repeat) = value;
            pocs(count, repeat) = percentoffsetscorrect(net, out);
            css(count, repeat) = correctspikes(net, out);
            nmos(count, repeat) = missedoffsets(net, out);
            iss(count, repeat) = incorrectspikes(net, out);
            out.accuracy = value;

            values(count, repeat) = value;
            fprintf('progress %s: alpha: %d, beta: %d, repeat: %d\n', output_folder, net.a1, net.b1, repeat);
        catch exception
            err = lasterror;
            fprintf('------------- ERRROR with simulator---------------------: \n%s\n\n%s\n\n', getReport(exception), err.message);
        end


        try 
            filename = sprintf('%s/%s_%d_%d_%d', output_folder, exp_name, net.a1, net.b1, repeat);
            % decrease file sizes...
            out.spike_time_trace = out.spike_time_trace(out.spike_time_trace(:, 2) == net.N, :);
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
            filename = sprintf('%s/res_%s_%d_%d_%d', output_folder, exp_name, net.a1, net.b1, repeat);
            %parfor_save(filename, values, var_range);
        catch exception
            fprintf('Failed to write results: %s\n%s\n\n', filename, getReport(exception));
        end
    end
end

filename = sprintf('%s/res_%s_final', output_folder, exp_name);
save(filename, 'values', 'var_range', 'pocs', 'css', 'iss', 'nmos', 'tpxtn', '-v7.3');

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
