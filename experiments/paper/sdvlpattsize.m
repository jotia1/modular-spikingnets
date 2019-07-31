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

net.sim_time_sec = 50;

Tp = 50;
Df = 10;
num_repeats = 50;
Np = 1000;
Pf = 2;
N_inp = net.group_sizes(1);

exp_name = 'relfrq';
output_folder = newoutputfolder(exp_name);
values = zeros(10, num_repeats);
count = 1;

rng(60);

for repeat = 1 : num_repeats
    count = 1;
%    for Np = 50: 50 : 1000
%    for fgi = 0.0230 : 0.0001 : 0.0238
    for Pf =  0.5 : 0.5 : 5.0
        %[pinp, pts] = generateuniformpattern(Tp, Np);
        %net.data_generator = @() repeatingrefinedpattern(Tp, Df, net.group_sizes(1), Np, pinp, pts);
        [~, ~, net.pinp, net.pts, ~] = poisspattfreq(Tp, Df, N_inp, Np, Pf);
        net.data_generator = @() poisspattfreq(Tp, Df, N_inp, Np, Pf, net.pinp, net.pts);
        
%        net.fgi = fgi;
        out = spikingnet(net);

        %value = detectionrate(net, out)
        value = offsetaccuracy(net, out, Tp)
        
        values(count, repeat) = value;
        fprintf('%s: count: %d, repeat: %d \n', output_folder, count, repeat);
        try 
            filename = sprintf('%s/%s_%d_%d', output_folder, exp_name, count, repeat);
            idxs = out.spike_time_trace(:, 2) == 2001;
            trace = out.spike_time_trace(idxs, :);
            offsets = out.offsets;
            save(filename, 'trace', 'offsets', 'count', 'repeat');
        catch exception
            fprintf('Failed to write: %s\n', filename);
        end
%
        try
            filename = sprintf('%s/res_%s_%d_%d', output_folder, exp_name, count, repeat);
            save(filename, 'values');
        catch exception
            fprintf('Failed to write results: %s\n', filename);
        end
        count = count + 1;
    end
end

filename = sprintf('%s/res_%s_final', output_folder, exp_name);
save(filename, 'values');

end
