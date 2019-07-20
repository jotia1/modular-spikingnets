function [ net, out ] = sdvlpattsize( params )
%% SDVLPATTSIZE - function will plot the effect of pattern size on SDVL
%   The Performance of SDVL with a varying pattern size is explored.
%
%   Example usage:
%       [ net, out ] = sdvlpattsize( [3, 1, 5, 5, 0.0236] )

% Needed for running as main on cluster
addpath(genpath('../../'));

net = defaultpapernetwork();

net.a1 = params(1);
net.a2 = params(2);
net.b1 = params(3);
net.b2 = params(4);
net.fgi = params(5);

net.sim_time_sec = 15;

Tp = 50;
Df = 10;

output_folder = newoutputfolder();

for Np = 50: 50 : 1000
    [pinp, pts] = generateuniformpattern(Tp, Np);
    net.data_generator = @() repeatingrefinedpattern(Tp, Df, net.group_sizes(1), pinp, pts);
    
    out = spikingnet(net);
    
    filename = sprintf('%s/exp1_%d', output_folder, Np);
    save(filename, 'net', 'out');
    
end


end
