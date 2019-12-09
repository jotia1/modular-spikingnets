function [] = executebayes(script, alphamin, alphamax, betamin, betamax, etamin, etamax, fgimin, fgimax, sim_time_sec)
%% play with optimisation
addpath(genpath('../'));

%
rng(1);
SIM_TIME=sim_time_sec;
fgivar = optimizableVariable('fgi', [fgimin, fgimax], 'Type', 'real');
nuvar = optimizableVariable('nu', [etamin, etamax], 'Type', 'real');
a1var = optimizableVariable('a1', [alphamin, alphamax], 'Type', 'integer');
b1var = optimizableVariable('b1', [betamin, betamax], 'Type', 'integer');
%[orig_inp, orig_ts, labels] = mnist2input(SIM_TIME);

fun = @(x) optimisenetwork(x, SIM_TIME);
results = bayesopt(fun, [fgivar, a1var, b1var, nuvar],'Verbose',1,...
    'AcquisitionFunctionName','expected-improvement-plus', ...
    'NumSeedPoints', 50, 'MaxObjectiveEvaluations',1000, 'UseParallel',true)


% load ionosphere
% rng default
% num = optimizableVariable('n',[1,30],'Type','integer');
% dst = optimizableVariable('dst',{'chebychev','euclidean','minkowski'},'Type','categorical');
% c = cvpartition(351,'Kfold',5);
% fun = @(x)kfoldLoss(fitcknn(X,Y,'CVPartition',c,'NumNeighbors',x.n,...
%     'Distance',char(x.dst),'NSMethod','exhaustive'));
% results = bayesopt(fun,[num,dst],'Verbose',0,...
%     'AcquisitionFunctionName','expected-improvement-plus')



function [ acc ] = optimisenetwork( x, SIM_TIME )

    %rng(1);
    net = defaultpapernetwork();
    net.rand_seed = -1;
    net.run_date = datestr(datetime);
    net.sim_time_sec = SIM_TIME

    net.Tp = 50;
    net.Df = 10;
    net.Np = 500;
    net.Pf = 5;
    net.dropout = 0.0;
    net.test_seconds = 50;
    
    net.fgi = x.fgi;
    net.a1 = x.a1; net.a2 = x.a1;
    net.b1 = x.b1; net.b2 = x.b1;
    net.nu = x.eta; net.nv = x.eta;

    net.use_simulated_annealing = false;
    net.If = 0.0222;
    net.Tf = 30;

    % Set up pattern
    net.pattfun = [];
    [net.pinp, net.pts] = generateuniformpattern( net.Tp, net.Np );         
    net.data_generator = @() balancedpoisson(net.Tp, net.Df, N_inp, net.Np, net.Pf, net.pinp, net.pts, net.pattfun, net.dropout);


    
    out = ssdvl(net);

    acc = trueposxtrueneg(net, out);
    %%%   DO NOT USE BELOW METRICS - MIGHT HAVE ERRORS
    %acc = -percentagecorrect(net, out, labels);
    %acc = -calcAccuracy(net, out);
    %acc = -totalspikes(net, out);

end

end % function 'executebayes'
