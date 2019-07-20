function [] = rasterspiketimes( spike_times_trace, N_inp, N_hid )
%% RASTERSPIKETIMES - Draw a rater plot of the spike times
%   Given a (s, 2) sized matrix with s spikes, draw a raster plot with the
%   input spikes in black and the output spikes in red. The supplied data
%   can be an (s, 2) sized matrix where the first column is the neuron
%   index and the second column is the time (in ms) of firing or the struct
%   output with a field called 'spike_time_trace'. 
%   
%   NOTE: will offset the raster plot such that the first spike happens at
%   1 ms. 


if isfield(spike_times_trace, 'spike_time_trace')  % if this is a struct of 'out'
    spike_times_trace = spike_times_trace.spike_time_trace;
end

if ~exist('N_inp', 'var') || ~exist('N_hid', 'var')
    N_inp = max(spike_times_trace(:, 1));
    N_hid = 0;
end

thousands = round(spike_times_trace(1, 1) / 1000);
offset = thousands * 1000; 
filter = find(spike_times_trace(:,2) <= N_inp + N_hid & spike_times_trace(:,2) > N_inp & spike_times_trace(:,1) > offset);
l2_spike_idxs = spike_times_trace(filter, 2);
l2_spike_times = spike_times_trace(filter, 1);
filter = find(spike_times_trace(:,2) < N_inp + N_hid & spike_times_trace(:, 1) > offset);
l1_idxs = spike_times_trace(filter, 2);
l1_times = spike_times_trace(filter, 1);

% Note this is the SPIKE TIME (not arrival time)
plot(l1_times - offset, l1_idxs, '.k', 'MarkerSize', 8);
hold on 
plot(l2_spike_times - offset, l2_spike_idxs, '.r', 'MarkerSize', 8)
ax = gca;
i = 0;

%axis([0, 30, -30 N_v + 50]);
while i < numel(l2_spike_times)
    i = i + 1;
    pos = l2_spike_times(i);
    c_idx = mod(l2_spike_idxs(i) - N_inp - 1, size(ax.ColorOrder, 1)) + 1;
    colour = ax.ColorOrder(c_idx, :);
    plot( [pos pos] - offset, get( gca, 'Ylim' ), '--', 'Color', colour, 'LineWidth',2);
end

% highlight pattern
do_highlights = true;
if do_highlights
    D_MAX = 20;
    patt_locs = 0:500:max(spike_times_trace(:, 1));
    p_length = 20;
    hold on
    for i = 1:numel(patt_locs)
        start_loc = patt_locs(i);
        xpos = start_loc + p_length + D_MAX;
        patch([start_loc xpos xpos start_loc ], [ 0 0 2002 2002  ], [0.6 0.4 0.9], 'FaceAlpha',0.5, 'EdgeColor','none')
    end
    drawnow;
    hold off
end

%axis([0, max(spike_times_trace(:, 1)), 0, 2102]);
title(sprintf('Time offset by %d thousands', thousands));
xlabel('Time (ms)');
ylabel('Neuron number');


end