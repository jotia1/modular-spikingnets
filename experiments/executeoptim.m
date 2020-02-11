function [ ] = executeoptim(exp_name, cpus, slurm_id, task, vars_to_set, vars_to_optimise, lbs, ubs, exp_notes, slurm_time)

    addpath(genpath('../'));
    
    % Try stopping 30 mins early so running trials can end before getting slurmed
    slurm_time = slurm_time - (30 * 60);  % 30 minutes x 60 seconds
    rng(1);

    parpool(cpus);
    output_folder = sprintf('%s_%d', exp_name, slurm_id); 
    mkdir(output_folder);
    
    if strcmp(task, 'BAYES')
        result = runbayes(vars_to_set, vars_to_optimise, lbs, ubs, output_folder, slurm_time);
    elseif strcmp(task, 'GENALGO')
        result = runga(vars_to_set, vars_to_optimise, lbs, ubs, output_folder, slurm_time);
    else
        fprintf('ERROR: Unrecognised task: %s\n', task);
        return;
    end

    filename = sprintf('%s/opt_%s_final', output_folder, exp_name);
    save(filename, 'result', 'slurm_id', 'exp_notes', '-v7.3');

end

function [results] = runga(vars_to_set, vars_to_optimise, lbs, ubs, output_folder, slurm_time)
    fprintf('Start optimising %s with GA\n', output_folder);

%     if strcmp(learning_rule, 'SSDVL')
%         lb = [alphamin, betamin, etamin, fgimin];
%         ub = [alphamax, betamax, etamax, fgimax];
%         fun = @(x) ssdvlnetwork(x, sim_time_sec);
%     elseif strcmp(learning_rule, 'SDVL')
%         lb = [alphamin, alphamin, betamin, betamin, etamin, etamin, fgimin];
%         ub = [alphamax, alphamax, betamax, betamax, etamax, etamax, fgimax];
%         fun = @(x) sdvlnetwork(x, sim_time_sec);
%     else
%         fprintf('ERROR: runga - Unknown learning rule: %s\n', learning_rule);
%         return;
%     end
    
    fun = @(x) sdvlnetwork(x, vars_to_set, vars_to_optimise);
        
    numvariables = numel(lbs); 
    popsize = 12;

    opts = optimoptions(@ga, ...
        'OutputFcn', @(options, state, flag) logpop(options, state, flag, output_folder), ...
        'PopulationSize', popsize, ...
        'MaxGenerations', 10, ...
        'MaxStallGenerations', 7,...
        'MaxTime', slurm_time, ...
        'FunctionTolerance', 1e-2, ...
        'FitnessScalingFcn', 'fitscalingrank', ...
        'SelectionFcn', 'selectionstochunif', ...
        'EliteCount', ceil(0.05 * popsize), ...
        'MutationFcn', 'mutationadaptfeasible', ...
        'CrossoverFraction', 0.8, ...
        'CrossoverFcn', 'crossoverscattered', ...
        'UseParallel',true);

    [x,Fval,exitFlag,Output,population,scores] = ga(fun, numvariables, [],[],[],[], lbs, ubs,[],opts);

    results = struct();
    results.x = x; results.exitFlag = exitFlag; 
    results.Fval = Fval; results.Output = Output;
    results.population = population; results.scores = scores;

    fprintf('The number of generations was : %d\n', Output.generations);            
    fprintf('The number of function evaluations was : %d\n', Output.funccount);     
    fprintf('The best function value found was : %g\n', Fval);

end


function [results] = runbayes(learning_rule, alphamin, alphamax, betamin, betamax, etamin, etamax, fgimin, fgimax, sim_time_sec, foldername, slurm_time)
    if strcmp(learning_rule, 'SSDVL')

        fgivar = optimizableVariable('fgi', [fgimin, fgimax], 'Type', 'integer');
        nuvar = optimizableVariable('nu', [etamin, etamax], 'Type', 'real');
        a1var = optimizableVariable('a1', [alphamin, alphamax], 'Type', 'integer');
        b1var = optimizableVariable('b1', [betamin, betamax], 'Type', 'integer');
        opt_vars = [fgivar, a1var, b1var, nuvar];
        fun = @(x) ssdvlnetwork(x, sim_time_sec);

    elseif strcmp(learning_rule, 'SDVL')

        fgivar = optimizableVariable('fgi', [fgimin, fgimax], 'Type', 'integer');
        nuvar = optimizableVariable('nu', [etamin, etamax], 'Type', 'real');
        nvvar = optimizableVariable('nv', [etamin, etamax], 'Type', 'real');
        a1var = optimizableVariable('a1', [alphamin, alphamax], 'Type', 'integer');
        a2var = optimizableVariable('a2', [alphamin, alphamax], 'Type', 'integer');
        b1var = optimizableVariable('b1', [betamin, betamax], 'Type', 'integer');
        b2var = optimizableVariable('b2', [betamin, betamax], 'Type', 'integer');
        opt_vars = [fgivar, a1var, a2var, b1var, b2var, nuvar, nvvar];
        fun = @(x) sdvlnetwork(x, sim_time_sec);

    else
        fprintf('ERROR: runbayes - Unknown learning rule: %s\n', learning_rule);
        return;
    end

    fprintf('Start optimising %s with BAYES\n', learning_rule);

    savefilename = sprintf('%s/progress.mat', foldername);
    results = bayesopt(fun, opt_vars, 'Verbose', 2,...
        'AcquisitionFunctionName','expected-improvement-plus', ...
        'NumSeedPoints', 10, 'MaxObjectiveEvaluations',100, ...
        'SaveFileName', savefilename, 'OutputFcn', @saveToFile, ...
        'MaxTime', slurm_time, ...
        'UseParallel',true);
end



% function [ acc ] = ssdvlnetwork( x, SIM_TIME )
% 
%     %rng(1);
%     net = defaultpapernetwork();
%     net.rand_seed = -1;
%     net.run_date = datestr(datetime);
%     net.sim_time_sec = SIM_TIME;
% 
%     net.Tp = 50;
%     net.Df = 10;
%     net.Np = 500;
%     net.Pf = 5;
%     net.dropout = 0.0;
%     net.test_seconds = 50;
% 
%     try
%     if istable(x)  %% BAYES OPT
%         net.fgi = (220 + x.fgi) / 1e4; % Hack to avoid need for many ptables
%         net.a1 = x.a1; net.a2 = x.a1;
%         net.b1 = x.b1; net.b2 = x.b1;
%         net.nu = x.nu; net.nv = x.nu;
%     else    %% GA OPT
%         net.fgi = (220 + x(4)) / 1e4;
%         net.a1 = x(1); net.a2 = x(1);
%         net.b1 = x(2); net.b2 = x(2);
%         net.nu = x(3); net.nv = x(3);
%     end
%     catch err
%         fid = fopen('errorFile','a+');
%         fprintf(fid, '%s', err.getReport('extended', 'hyperlinks','off'))
%         fclose(fid)
%         acc = 1.1;
%         return
%     end
% 
%     net.use_simulated_annealing = false;
%     net.If = 0.0222;
%     net.Tf = 30;
% 
%     % Set up pattern
%     net.pattfun = [];
%     [net.pinp, net.pts] = generateuniformpattern( net.Tp, net.Np );         
%     net.data_generator = @() balancedpoisson(net.Tp, net.Df, net.group_sizes(1), net.Np, net.Pf, net.pinp, net.pts, net.pattfun, net.dropout);
% 
% 
%     try    
%         out = ssdvl(net);
%     catch err
%         fid = fopen('errorFile','a+');
%         fprintf(fid, '%s', err.getReport('extended', 'hyperlinks','off'))
%         fclose(fid)
%         acc = 1.1;
%         return
%     end
% 
%     acc = 1 - trueposxtrueneg(net, out);
%     %%%   DO NOT USE BELOW METRICS - MIGHT HAVE ERRORS
%     %acc = -percentagecorrect(net, out, labels);
%     %acc = -calcAccuracy(net, out);
%     %acc = -totalspikes(net, out);
% 
% end

function [ acc ] = sdvlnetwork( x, vars_to_set, vars_to_optimise )
    
    %rng(1);
    net = defaultnetwork();
    net.rand_seed = -1;
    net.run_date = datestr(datetime);

    net.Tp = 50;
    net.Df = 10;
    net.Np = 500;
    net.Pf = 5;
    net.dropout = 0.0;
    net.test_seconds = 50;
    
    if istable(x)  %% BAYES OPT
        net.fgi = (220 + x.fgi) / 1e4; % Hack to avoid need for many ptables
        net.a1 = x.a1; net.a2 = x.a1;
        net.b1 = x.b1; net.b2 = x.b1;
        net.nu = x.nu; net.nv = x.nu;
    else    %% GA OPT
        %net.fgi = (220 + x(7)) / 1e4;
        %net.a1 = x(1); net.a2 = x(2);
        %net.b1 = x(3); net.b2 = x(4);
        %net.nu = x(5); net.nv = x(6);
        % Set up any variables that need optimising.
        for i = 1 : numel(vars_to_optimise)
            if strcmp(vars_to_optimise{i}, 'fgi')
                net.('fgi') = (220 + x(i)) / 1e4;
                continue;
            end
            net.(vars_to_optimise{i}) = x(i);
        end
        
    end
    
    % TODO : remove after standardising default network
    net.use_simulated_annealing = false;
    net.If = 0.0222;
    net.Tf = 30;
    
    % Set any variables that need setting
    assert(mod(numel(vars_to_set),2)==0, 'vars_to_set must be key value pairs, got an odd number.');
    for i = 1 : 2 : numel(vars_to_set)
        net.(vars_to_set{i}) = vars_to_set{i + 1};
    end

    % Set up pattern
    net.pattfun = [];
    [net.pinp, net.pts] = generateuniformpattern( net.Tp, net.Np );         
    net.data_generator = @() balancedpoisson(net.Tp, net.Df, net.group_sizes(1), net.Np, net.Pf, net.pinp, net.pts, net.pattfun, net.dropout);
    
    try    
        out = spikingnet(net);
    catch err
        fid = fopen('errorFile','a+');
        fprintf(fid, '%s', err.getReport('extended', 'hyperlinks','off'))
        fprintf('%s\n\nSee experiments/errorFile for more.', err.getReport('extended', 'hyperlinks','off'));
        fclose(fid)
        acc = 1.1;
        return
    end

    acc = 1 - trueposxtrueneg(net, out);
    %%%   DO NOT USE BELOW METRICS - MIGHT HAVE ERRORS
    %acc = -percentagecorrect(net, out, labels);
    %acc = -calcAccuracy(net, out);
    %acc = -totalspikes(net, out);

end

function [state, options, optchanged] = logpop(options, state, flag, outfolder)

    filename = sprintf('%s/poplog.mat', outfolder);
    fprintf('Try to save pop for Generation: %d\nFilename: %s, folder: %s\n', state.Generation, filename, outfolder);

    if exist(filename) == 0 % No file yet
        poplog = [state.Population, state.Score];
    else
        load(filename, 'poplog');
        poplog(:, :, end+1) = [state.Population, state.Score];
    end

    save(filename, 'poplog', '-v7.3');

    % Also force re-evaluation of elites and stop duplicates in population
    state.EvalElites = true;
    state.HaveDuplicates = false;

    optchanged = false;
end
