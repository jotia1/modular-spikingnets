function [ result ] = incorrectspikes(net, out)
%% INCORRECTSPIKES - Calculate the number of spikes not within offsets
%   Count the number of spikes not occuring during pattern periods.
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

result = numel(output_spike_times) - correctspikes(net, out);

end