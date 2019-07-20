function [ inp, ts ] = createTwoPatternInput( p1, p2, freq, time_sec, sprvsn_time_sec, sprvsn_spike_time_ms )
%CREATETWOPATTERNINPUT Helper wrapper function around create input
%   Assumes 2 Hz rate
%   Assumes neurons 4 and 5 are outputs for supervision
%   Assumes 500ms between patterns
%
%

assert(length(p1) == length(p2), 'Patterns must be the same length');

[p1inp, p1ts] = createInput(p1, freq * 2, time_sec, sprvsn_time_sec, sprvsn_spike_time_ms);
[p2inp, p2ts] = createInput(p2, freq * 2, time_sec, sprvsn_time_sec, sprvsn_spike_time_ms);
p2inp(p2inp == length(p1) + 1) = length(p1) + 2;

inp = [p1inp, p2inp];
ts = [p1ts, p2ts + freq];

[ts, idxs] = sort(ts);
inp = inp(idxs);

end

