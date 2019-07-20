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

net.sim_time_sec = 15;

Tp = 50;
Df = 10;

output_folder = newoutputfolder();
values = zeros(20, 5);
count = 1;

for Np = 50: 50 : 1000
    for repeat = 1 : 5
        [pinp, pts] = generateuniformpattern(Tp, Np);
        net.data_generator = @() repeatingrefinedpattern(Tp, Df, net.group_sizes(1), pinp, pts);
    
        out = spikingnet(net);

        value = detectionrate(net, out);
        
        values(count, repeat) = value;
    
        filename = sprintf('%s/exp1_%d', output_folder, Tp);
        save(filename, 'net', 'out');
    end
    count = count + 1;
end

filename = sprintf('%s/results', output_folder);
save(filename, 'values');


end
