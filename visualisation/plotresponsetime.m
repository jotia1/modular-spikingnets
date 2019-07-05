function [] = plotresponsetime(net, out, neurons)

leg = cell(1, numel(neurons));
for n = 1 : numel(neurons)

    neuron = neurons(n);
    
    plot(out.spike_time_trace(out.spike_time_trace(:, 2) == neuron, 1), mod(out.spike_time_trace(out.spike_time_trace(:, 2) == neuron, 1), 500), '.')
    hold on
    plot(out.spike_time_trace(out.spike_time_trace(:, 2) == neuron, 1), movmean(mod(out.spike_time_trace(out.spike_time_trace(:, 2) == neuron, 1), 500), 20));
    
    leg{n} = ['N', num2str(neuron)]; 
    leg{n} = ['N', num2str(neuron), ' movmean'];
    
end
%plot(out.spike_time_trace(out.spike_time_trace(:, 2) == 5, 1), mod(out.spike_time_trace(out.spike_time_trace(:, 2) == 5, 1), 500), '.')
hold off
%axis([-Inf Inf 5 14]);
title('Output spike time');
legend(leg);


end