function [ num_spikes ] = totalspikes(net, out)

spikes = out.spike_time_trace(:, 2) == 29 | out.spike_time_trace(:, 2) == 30;
                    
num_spikes = sum(spikes(:));



end
