function [ inp, ts ] = joineventstreams( inp1, ts1, inp2, ts2 )
%JOINEVENTSTREAMS - Join torgether two event streams and sort by time
%   Sort two event streams given the firing neurons index (inp1/inp2) and
%   the time at which those neurons fired (ts1/ts2). Any duplicates are
%   removed.
%
%   Example usage:
%       [ inp, ts ] = joineventstreams( [1; 2], [3; 5], [2; 3], [8; 4] )

assert((isempty(inp1) || size(inp1, 2) == 1) && (isempty(inp2) || size(inp2, 2) == 1), 'inps must have a single column');
assert((isempty(ts1) || size(ts1, 2) == 1) && (isempty(ts2) || size(ts2, 2) == 1), 'tss must have a single column');

inp = [inp1; inp2];
ts = [ts1; ts2];

% Sort and remove duplicates using unique
data = unique([ts, inp], 'rows');
inp = data(:, 2);
ts = data(:, 1);

end

