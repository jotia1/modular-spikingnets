function [] = executeoptim(exp_name, learning_rule, cpus, slurm_id, task, alphamin, alphamax, betamin, betamax, etamin, etamax, fgimin, fgimax, sim_time_sec, exp_notes, slurm_time)
    %% play with optimisation
    addpath(genpath('../'));

    % Try stopping 30 mins early so running trials can end before getting slurmed
    slurm_time = slurm_time - (30 * 60);  % 30 minutes x 60 seconds
    rng(1);
    SIM_TIME=sim_time_sec;
    parpool(cpus);
    output_folder = sprintf('%s_%d', exp_name, slurm_id); %newoutputfolder(exp_name);
    
    if strcmp(task, 'BAYES')
        if strcmp(learning_rule, 'SDVL')
            result = optimisesdvl(alphamin, alphamax, betamin, betamax, etamin, etamax, fgimin, fgimax, sim_time_sec, output_folder, slurm_time);
        elseif stcmp(learning_rule, 'SSDVL')
            result = optimisessdvl(alphamin, alphamax, betamin, betamax, etamin, etamax, fgimin, fgimax, sim_time_sec, output_folder, slurm_time);
        else
            fprintf('ERROR: Unrecognised learning_rule: %s\n', learning_rule);
            return;
        end
    elseif strcmp(task, 'GENALGO')
        if strcmp(learning_rule, 'SDVL')
            fprintf('ERROR: Not yet implemented, %s for %s\n', task, learning_rule);
        elseif stcmp(learning_rule, 'SSDVL')
            result = gassdvl(alphamin, alphamax, betamin, betamax, etamin, etamax, fgimin, fgimax, sim_time_sec, output_folder, slurm_time);

        else
            fprintf('ERROR: Unrecognised learning_rule: %s\n', learning_rule);
            return;
        end
    else
        fprintf('ERROR: Unrecognised task: %s\n', task);
        return;
    end


    filename = sprintf('%s/opt_%s_final', output_folder, exp_name);
    save(filename, 'result', 'slurm_id', 'exp_notes', '-v7.3');

end

function [results] = gassdvl(alphamin, alphamax, betamin, betamax, etamin, etamax, fgimin, fgimax, sim_time_sec, foldername, slurm_time)
    fprintf('start optimising with GA');

    fun = @(x) ssdvlnetwork(x, sim_time_sec);
    numvariables = 4; 

    opts = optimoptions(@ga, ...
        'OutputFcn', @(options, state, flag) logpop(options, state, flag, foldername), ...
        'PopulationSize', 12, ...
        'MaxGenerations', 10, ...
        'MaxStallGenerations', 7,...
        'MaxTime', slurm_time, ...
        'FunctionTolerance', 1e-2, ...
        'SelectionFcn', 'selectionstochunif', ...
        'MutationFcn', 'mutationadaptfeasible', ...
        'CrossoverFcn', 'crossoverscattered', ...
        'UseParallel',true);

    % Bounds:
    lb = [alphamin, betamin, etamin, fgimin];
    ub = [alphamax, betamax, etamax, fgimax];
    %lb = [1; 0.1; 0.1; 1; 1; 1; 1; 1; 1];
    %ub = [10; 0.5; 0.5; 9; 9; 9; 9; 10; Inf];

    [x,Fval,exitFlag,Output,population,scores] = ga(fun, numvariables, [],[],[],[], lb, ub,[],opts);

    results = struct();
    results.x = x; results.exitFlag = exitFlag; 
    results.Fval = Fval; results.Output = Output;
    results.population = population; results.scores = scores;

    fprintf('The number of generations was : %d\n', Output.generations);            
    fprintf('The number of function evaluations was : %d\n', Output.funccount);     
    fprintf('The best function value found was : %g\n', Fval);

end


function [results] = optimisessdvl(alphamin, alphamax, betamin, betamax, etamin, etamax, fgimin, fgimax, sim_time_sec, foldername, slurm_time)
    fprintf('Start optimising ssdvl');
    fgivar = optimizableVariable('fgi', [fgimin, fgimax], 'Type', 'integer');
    nuvar = optimizableVariable('nu', [etamin, etamax], 'Type', 'real');
    a1var = optimizableVariable('a1', [alphamin, alphamax], 'Type', 'integer');
    b1var = optimizableVariable('b1', [betamin, betamax], 'Type', 'integer');
    savefilename = sprintf('%s/prog.mat', foldername);

    fun = @(x) ssdvlnetwork(x, sim_time_sec);
    results = bayesopt(fun, [fgivar, a1var, b1var, nuvar],'Verbose',1,...
        'AcquisitionFunctionName','expected-improvement-plus', ...
        'NumSeedPoints', 100, 'MaxObjectiveEvaluations',1000, ...
        'SaveFileName', savefilename, 'OutputFcn', @saveToFile, ...
        'MaxTime', slurm_time, ... 
        'UseParallel',true);
end


function [results] = optimisesdvl(alphamin, alphamax, betamin, betamax, etamin, etamax, fgimin, fgimax, sim_time_sec, foldername, slurm_time)
    fprintf('Start optimising SDVL');
    fgivar = optimizableVariable('fgi', [fgimin, fgimax], 'Type', 'integer');
    nuvar = optimizableVariable('nu', [etamin, etamax], 'Type', 'real');
    nvvar = optimizableVariable('nv', [etamin, etamax], 'Type', 'real');
    a1var = optimizableVariable('a1', [alphamin, alphamax], 'Type', 'integer');
    a2var = optimizableVariable('a2', [alphamin, alphamax], 'Type', 'integer');
    b1var = optimizableVariable('b1', [betamin, betamax], 'Type', 'integer');
    b2var = optimizableVariable('b2', [betamin, betamax], 'Type', 'integer');
    savefilename = sprintf('%s/prog.mat', foldername);

    fun = @(x) sdvlnetwork(x, sim_time_sec);
    results = bayesopt(fun, [fgivar, a1var, a2var, b1var, b2var, nuvar, nvvar],'Verbose',1,...
        'AcquisitionFunctionName','expected-improvement-plus', ...
        'NumSeedPoints', 100, 'MaxObjectiveEvaluations',1000, ...
        'SaveFileName', savefilename, 'OutputFcn', @saveToFile, ...
        'MaxTime', slurm_time, ...
        'UseParallel',true);
end



function [ acc ] = ssdvlnetwork( x, SIM_TIME )

    %rng(1);
    net = defaultpapernetwork();
    net.rand_seed = -1;
    net.run_date = datestr(datetime);
    net.sim_time_sec = SIM_TIME;

    net.Tp = 50;
    net.Df = 10;
    net.Np = 500;
    net.Pf = 5;
    net.dropout = 0.0;
    net.test_seconds = 50;

    try
    if istable(x)  %% BAYES OPT
        net.fgi = (220 + x.fgi) / 1e4; % Hack to avoid need for many ptables
        net.a1 = x.a1; net.a2 = x.a1;
        net.b1 = x.b1; net.b2 = x.b1;
        net.nu = x.nu; net.nv = x.nu;
    else    %% GA OPT
        net.fgi = (220 + x(4)) / 1e4;
        net.a1 = x(1); net.a2 = x(1);
        net.b1 = x(2); net.b2 = x(2);
        net.nu = x(3); net.nv = x(3);
    end
    catch err
        fid = fopen('errorFile','a+');
        fprintf(fid, '%s', err.getReport('extended', 'hyperlinks','off'))
        fclose(fid)
        acc = 1.1;
        return
    end

    net.use_simulated_annealing = false;
    net.If = 0.0222;
    net.Tf = 30;

    % Set up pattern
    net.pattfun = [];
    [net.pinp, net.pts] = generateuniformpattern( net.Tp, net.Np );         
    net.data_generator = @() balancedpoisson(net.Tp, net.Df, net.group_sizes(1), net.Np, net.Pf, net.pinp, net.pts, net.pattfun, net.dropout);


    try    
        out = ssdvl(net);
    catch err
        fid = fopen('errorFile','a+');
        fprintf(fid, '%s', err.getReport('extended', 'hyperlinks','off'))
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

function [ acc ] = sdvlnetwork( x, SIM_TIME )

    %rng(1);
    net = defaultpapernetwork();
    net.rand_seed = -1;
    net.run_date = datestr(datetime);
    net.sim_time_sec = SIM_TIME;

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
        net.fgi = (220 + x(7)) / 1e4;
        net.a1 = x(1); net.a2 = x(2);
        net.b1 = x(3); net.b2 = x(4);
        net.nu = x(5); net.nv = x(6);
    end


    net.use_simulated_annealing = false;
    net.If = 0.0222;
    net.Tf = 30;

    % Set up pattern
    net.pattfun = [];
    [net.pinp, net.pts] = generateuniformpattern( net.Tp, net.Np );         
    net.data_generator = @() balancedpoisson(net.Tp, net.Df, net.group_sizes(1), net.Np, net.Pf, net.pinp, net.pts, net.pattfun, net.dropout);
    
    try    
        out = spikingnet(net);
    catch err
        fid = fopen('errorFile','a+');
        fprintf(fid, '%s', err.getReport('extended', 'hyperlinks','off'))
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

    optchanged = false;
end
