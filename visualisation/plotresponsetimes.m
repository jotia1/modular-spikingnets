function [ ] = plotresponsetimes(net, out)

offsets = out.offsets;
N = sum(net.group_sizes);
filter = (out.spike_time_trace(:, 2) == N);
output_spike_times = out.spike_time_trace(filter, 1);

times = [];

for i = 1 : numel(offsets)
    offset = offsets(i);
    next_spikes = output_spike_times(output_spike_times > offset & output_spike_times <= (offset + net.Tp + net.delay_max));
    if numel(next_spikes) > 0
        times = [times; offset, next_spikes(1) - offset];
    end
end

plot(times(:, 1), times(:, 2), '.');
title('Time of first spike within pattern (if any)');

end