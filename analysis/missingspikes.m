function [ result ] = missedoffsets(net, out)
%% MISSEDOFFSETS - Calculate the number of offsets without a spike
%   TODO
%
%   Parameters:
%       net - the network struct
%       out - the network output struct
%
%   Assumptions:
%       - there is only 1 output neuron. 

testing_seconds = net.test_seconds;
training_ms = (net.sim_time_sec - testing_seconds) * 1000;

N = sum(net.group_sizes);
filter = (out.spike_time_trace(:, 2) == N) & (out.spike_time_trace(:, 1) > training_ms);
output_spike_times = out.spike_time_trace(filter, 1);

test_offsets = out.offsets(out.offsets >= training_ms);

missing_spikes = 0;
for i = 1 : numel(test_offsets)
    offset = test_offsets(i);
    
    offset_spikes = sum(output_spike_times > offset & output_spike_times < (offset + net.Tp + net.delay_max));
    if offset_spikes == 0
        missing_spikes = missing_spikes + 1;
    end
    
end

result = missing_spikes;

end