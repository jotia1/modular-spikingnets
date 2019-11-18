function [ result ] = correctspikes(net, out)
%% CORRECTSPIKES - Calculate the number of spikes within offsets
%   Count the number of spikes occuring during pattern periods. NOTE this
%   includes duplicates in the same period.
%
%   Parameters:
%       net - the network struct
%       out - the network output struct
%
%   Assumes there is only 1 output neuron. 

testing_seconds = net.test_seconds;
training_ms = (net.sim_time_sec - testing_seconds) * 1000;

N = sum(net.group_sizes);
filter = (out.spike_time_trace(:, 2) == N) & (out.spike_time_trace(:, 1) > training_ms);
output_spike_times = out.spike_time_trace(filter, 1);

test_offsets = out.offsets(out.offsets >= training_ms);

correct_spikes = 0;
for i = 1 : numel(test_offsets)
    offset = test_offsets(i);
    
    offset_spikes = sum(output_spike_times > offset & output_spike_times < (offset + net.Tp + net.delay_max));
    correct_spikes = correct_spikes + offset_spikes;
    
end

result = correct_spikes;

end