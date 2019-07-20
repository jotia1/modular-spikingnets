function [ percentage ] = detectionrate(net, out, Pf, Tp)
%% DETECTIONRATE - Calculate the percentage of spikes in correct windows
%   Calculate the percentage of spikes that occur after the start of a
%   pattern but before the end of the pattern plus delay max. 
%
%   Parameters:
%       net - the network struct
%       out - the network output struct
%       Pf - the Pattern frequency. Default 2Hz
%       Tp - Time of the pattern (length). Default 50ms
%
%   Assumes there is only 1 output neuron. 

if ~exist('Pf', 'var')
    Pf = 2; 
end

if ~exist('Tp', 'var')
    Tp = 50;
end

pattern_freq_ms = 1000 / Pf;

assert(net.group_sizes(2) == 1, 'More than 1 output neuron, invalidates assumptions');

spike_times = out.spike_time_trace(out.spike_time_trace(:, 2) == sum(net.group_sizes), 1);

spike_offsets = mod(spike_times, pattern_freq_ms);

correctly_timed = sum(spike_offsets < (Tp + net.delay_max));

percentage = correctly_timed / numel(spike_offsets);

end