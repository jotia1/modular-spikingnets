function [ accuracy ] = offsetaccuracy(net, out, Tp)
%% OFFSETACCURACY - Calculate the number of spikes within offsets
%
%   Parameters:
%       net - the network struct
%       out - the network output struct
%       Pf - the Pattern frequency. Default 2Hz
%       Tp - Time of the pattern (length). Default 50ms
%
%   Assumes there is only 1 output neuron. 

N = sum(net.group_sizes);
output_spike_times = out.spike_time_trace(out.spike_time_trace(:, 2) == N, 1);

correct_spikes = 0;
for i = 1 : numel(out.offsets)
    offset = out.offsets(i);
    
    
    offset_spikes = sum(output_spike_times > offset & output_spike_times < (offset + Tp + net.delay_max));
    correct_spikes = correct_spikes + offset_spikes;
    
end

accuracy = correct_spikes / numel(output_spike_times);
end