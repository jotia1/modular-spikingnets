function [ accuracy ] = offsetaccuracy(net, out, Tp, testing_seconds)
%% OFFSETACCURACY - Calculate the number of spikes within offsets
%   Only calculates the accuracy in the last half of the simulation which
%   is considered after training by default, else the last number of
%   seconds can be passed in as an argument. 
%
%   Parameters:
%       net - the network struct
%       out - the network output struct
%       Tp - Time of the pattern (length). Default 50ms
%       testing_seconds - Number of seconds from end to include in
%           calculation.
%
%   Assumes there is only 1 output neuron. 

if ~exist('testing_seconds', 'var')
    testing_seconds = net.sim_time_sec / 2;
end
testing_ms = testing_seconds * 1000;

N = sum(net.group_sizes);
filter = (out.spike_time_trace(:, 2) == N) & (out.spike_time_trace(:, 1) > testing_ms);
output_spike_times = out.spike_time_trace(filter, 1);

test_offsets = out.offsets(out.offsets >= testing_ms);

correct_spikes = 0;
for i = 1 : numel(test_offsets)
    offset = test_offsets(i);
    
    offset_spikes = sum(output_spike_times > offset & output_spike_times < (offset + Tp + net.delay_max));
    correct_spikes = correct_spikes + offset_spikes;
    
end

accuracy = correct_spikes / numel(output_spike_times);
if numel(output_spike_times) == 0
    accuracy = 0;
end
end