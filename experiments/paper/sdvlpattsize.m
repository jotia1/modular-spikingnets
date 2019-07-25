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
values = zeros(20, 10);
count = 1;

rng(56);

for repeat = 1 : 10
    count = 1;
    for Np = 50: 50 : 1000
        [pinp, pts] = generateuniformpattern(Tp, Np);
        net.data_generator = @() repeatingrefinedpattern(Tp, Df, net.group_sizes(1), pinp, pts);
    
        out = spikingnet(net);

        value = detectionrate(net, out)
        
        values(count, repeat) = value;
        try 
            filename = sprintf('%s/exp1_%d_%d', output_folder, Np, repeat);
            save(filename, 'net', 'out');
        catch exception
            fprintf('Failed to write: %s', filename);
        end

        try
            filename = sprintf('%s/results_%d_%d', output_folder, Np, repeat);
            save(filename, 'values');
        catch exception
            fprintf('Failed to write results: %s', filename);
        end
        count = count + 1;
    end
end

filename = sprintf('%s/results_%d', output_folder, Np);
save(filename, 'values');

end
