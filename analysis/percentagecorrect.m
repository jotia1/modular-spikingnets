function [ acc ] = percentagecorrect(net, out, labels)
%% PERCENTAGECORRECT - Given a networks output, calculate the percentage correct
%   Given a network and its output calculate the percentage of correct
%   classifications. Assumes presentations every 500ms, assumes a correct
%   classification in one in which only the allocated neuron fires (one or
%   many times) in response to a pattern. Patterns are 28ms long, the
%   longest delay usually <= 20ms, so only the window of 48ms after
%   presentation are considered in the classification.

% TODO : Currently considered any spike in the full 500ms (should only
% consider the 48ms). 

num_inputs = net.group_sizes(1);

output_spikes = out.spike_time_trace(out.spike_time_trace(:, 2) > num_inputs, :);
output_spikes_presentation_nums = ceil(output_spikes(:, 1) / 500);

total_presentations = ceil(net.sim_time_sec * 1000 / 500);

%% Get binary vectors for which outputs fired at which presentations
output29 = zeros(total_presentations, 1);
idxs = output_spikes(:, 2) == 29;
output_spikes_presentation_nums(idxs)
output29(output_spikes_presentation_nums(idxs)) = 1;

output30 = zeros(total_presentations, 1);
idxs = output_spikes(:, 2) == 30;
output_spikes_presentation_nums(idxs)
output30(output_spikes_presentation_nums(idxs)) = 1;

%% Now filter by only correct outputs for each presentation
uniquespikes = xor(output29, output30);
filtered29as0 = uniquespikes & output29 & labels == 1; 
filtered30as8 = uniquespikes & output30 & labels == 9;
filteredop1 = or(filtered29as0, filtered30as8);

% In case the neurons are responding in opposite
filtered29as8 = uniquespikes & output29 & labels == 9; 
filtered30as0 = uniquespikes & output30 & labels == 1;
filteredop2 = or(filtered29as8, filtered30as0);

%% Calculate final accuracy as a ratio of correct unique spikes
acc1 = sum(filteredop1(:)) / total_presentations;
acc2 = sum(filteredop2(:)) / total_presentations;

acc = max(acc1, acc2);


end