function [ inp, ts, patt_inp, patt_ts ] = repeatingcrudepattern( gaps_ms, N_total, N_patt, patt_inp, patt_ts, patt_length_ms, total_length_ms )
%REPEATINGCRUDEPATTERN - Return a repeating pattern in noise every 500ms
%   Generates a repeating pattern (pattern itself is noiseless) embedded in
%   a noisy background.

if ~exist('gaps_ms', 'var') || isempty(gaps_ms)
    gaps_ms = 5;
end

offset = 0;
inp = [];
ts = [];
while offset < total_length_ms
    [ninp, nts, patt_inp, patt_ts] = crudepattinnoise(N_total, N_patt, patt_inp, patt_ts, patt_length_ms, 500);
    [inp, ts] = joineventstreams(inp, ts, ninp, nts + offset);
    offset = offset + gaps_ms
end

end

