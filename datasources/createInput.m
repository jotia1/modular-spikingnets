function [ inp, ts ] = createInput(seq, freq, time_sec, sprvsn_time_sec, sprvsn_spike_time_ms)
% CREATEINPUT - Given a sequence create network input
%   Given a sequence of firing times (assume each represents a single
%   neuron firing at that time)  create a valid network input that is this
%   firing sequence repeated every 'freq' ms's. The input will be
%   'time_sec' long and will have supervision for the first
%   'sprvsn_time_sec' seconds. Supervision spike will occur after
%   'spvsn_spike_time_ms' ms's.
%
%   Assumptions:
%       - seq is a zero indexed pattern
%       -        
%
%   Example
%       [inp, ts] = createInput([0 1 7], 500, 100, 50, 13)

N = numel(seq);
ms_per_sec = 1000;
num_presentations = ceil(time_sec * (ms_per_sec / freq));
inp = repmat(1:N, 1, num_presentations); 
ts = reshape(repmat(((0:num_presentations -1 )' * freq)', N, 1), 1, num_presentations * N) + ...
        repmat(seq, 1, num_presentations) + 1; 

% Set up supervision
num_supervisions = ceil(sprvsn_time_sec * (ms_per_sec / freq));
inp = [inp, ones(1, num_supervisions) * length(seq) + 1];
ts = [ts, (0:freq:((sprvsn_time_sec * ms_per_sec) - 1)) + sprvsn_spike_time_ms + 1];

[ts, idxs] = sort(ts);
inp = inp(idxs);

end