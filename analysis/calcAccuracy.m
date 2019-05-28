function [ accuracy ] = calcAccuracy( net, out )

    num_inputs = net.group_sizes(1);
    
    output_spikes = out.spike_time_trace(out.spike_time_trace(:, 2) > num_inputs, :);
    output_spikes_presentation_nums = ceil(output_spikes(:, 1) / 500);
    
    total_presentations = ceil(net.sim_time_sec * 1000 / 500);
    res = zeros(total_presentations, 1);
    even_idx = -1; odd_idx = -1;  % which neuron respresents even presentations and which odd
    
    for p = 1 : total_presentations
        presentation_spikes = output_spikes(output_spikes_presentation_nums == p, :);
         
        single_spike_occured = numel(unique(presentation_spikes(:, 2))) == 1; % xor 4 and 5 fired
        if ~single_spike_occured
            continue;
        end
        output_idx = presentation_spikes(1, 2); % Must be at least one spike
        
        if mod(p, 2) == 0
            if even_idx == -1 && output_idx ~= odd_idx
                even_idx = output_idx;
                %fprintf('Set Even_idx: %d\n', output_idx);
            end
            res(p) = output_idx == even_idx;
        else
            if odd_idx == -1 && output_idx ~= even_idx
                odd_idx = output_idx;
                %fprintf('Set Odd_idx: %d\n', odd_idx);
            end
            res(p) = output_idx == odd_idx;
        end
        
        
        
    end

% 
%     output_within_20ms_idxs = mod(out.spike_time_trace(:, 1), 500) < 20 & out.spike_time_trace(:, 2) > num_inputs;
%     presentation_number = ceil(out.spike_time_trace(:, 1) / 500);
%     
%     spikes_per_presentation = accumarray( presentation_number, output_within_20ms_idxs);
%     total_presentations = net.sim_time_sec * 1000 / 500;
%     assert(total_presentations >= max(presentation_number) & min(presentation_number) > 0, 'Sometime failed dividing up presentations');
%     
%     odd_presentations = 1 : 2 : total_presentations;
%     presentation_result = zeros(total_presentations, 1);
%     even_output = -1; odd_output = -1;
%      
%     for p = 1 : total_presentations
%         outputs_fired = output_within_20ms_idxs(presentation_number == p);
%         
%         presentation_result(p) = numel(unique(outputs_fired)) == 1;  % If only 4 xor 5 fired.
%         if mod(p, 2) == 0 && presentation_result(p) % even
%             if even_output == -1  % First even 
%                 even_output = unique(outputs_fired);
%             end
%             presentation_result(p) = presentation_result(p) & unique(outputs_fired) == even_output;
%         else if presentation_result(p) % odd
%             if odd_output == -1  % First even 
%                 odd_output = unique(outputs_fired);
%             end
%             presentation_result(p) = presentation_result(p) & unique(outputs_fired) == odd_output;
%         end
%     end
    
    %plot(presentation_result);
    num_correct = sum(res(:));
    
    accuracy = num_correct / total_presentations;

end