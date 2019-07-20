function [ ] = plotvardels(net, out, neurons, asraster)
%% PLOTVARDELS - Plot variance and delays as either a trace or final values
%   Given a network and output, plot the delays and variances as subplots
%   in a new figure. By default plot the traces of the 2nd neuron that was
%   traced. Optionally add asraster = true to plot the distribution of the 
%   final delay and variance values. 
%
%   Example usage:
%       plotvardels(net, out, [2], true)

if ~exist('neurons', 'var')
    neurons = [2];
end

if ~exist('asraster', 'var')
    neurons = false;
end

if asraster
    
figure;
subplot(2, 1, 1);
if asraster
    plot(out.delays(:, 2001), '.');
else
    plottrace(out.delayst, neurons);
end
title('Delays');

subplot(2, 1, 2);
if asraster
    plot(out.variance(:, 2001), '.');
else
    plottrace(out.vart, neurons);
end
title('Variance');


end