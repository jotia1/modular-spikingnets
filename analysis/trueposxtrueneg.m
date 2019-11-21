function [ result ] = trueposxtrueneg(net, out)
%% TRUEPOSXTRUENEG - True positives multiplied by the True negatives
%   As the name suggests the nunber of true positives (defined as the
%   percentage of patterns with at least one spike) multiplied by the True
%   negatives (defined as the percentage of gaps between patterns with at
%   least one spike).
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

%
last_offset = -net.Tp - net.delay_max;

true_positives = 0;
true_negatives = 0;
num_true_neg_slots = 0;
for i = 1 : numel(test_offsets)
    offset = test_offsets(i);
    
    % Have we had any spikes since the end of the last pattern?
    pre_offset_spikes = sum(output_spike_times >= last_offset + net.Tp + net.delay_max & output_spike_times <= offset);
    
    % Count the true neg slots as any slots with at least Tp + Dm space
    if offset > last_offset + net.Tp + net.delay_max  % if there has been time between slots
        num_true_neg_slots = num_true_neg_slots + 1;
        if pre_offset_spikes == 0 % no spikes occured since last offset
            true_negatives = true_negatives + 1;
        end
    end
        
    last_offset = offset;
    
    % Deal with pattern offsets
    offset_spikes = sum(output_spike_times > offset & output_spike_times < (offset + net.Tp + net.delay_max));
    if offset_spikes > 0   %% at least one spike occurs during pattern 
        true_positives = true_positives + 1;
    end
    
end

prop_TP = true_positives / numel(test_offsets);
prop_TN = true_negatives / num_true_neg_slots;

result = prop_TP * prop_TN;

end