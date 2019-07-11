function [] = plottrace( trace, neurons, addlegend )

if ~exist('addlegend', 'var')
    addlegend = false;
end

leg = cell(1, numel(neurons));
for n = 1 : numel(neurons)

    neuron = neurons(n);
    plot(squeeze(trace(neuron, :, :))');
    hold on
    leg{n} = ['N', num2str(neuron)]; 
    
end
%plot(out.spike_time_trace(out.spike_time_trace(:, 2) == 5, 1), mod(out.spike_time_trace(out.spike_time_trace(:, 2) == 5, 1), 500), '.')
hold off
%axis([-Inf Inf 5 14]);
title('Output trace');
if addlegend
    legend(leg);
end

end