function [ inp, ts ] = patterndropoutfunc(inp, ts, percent)
%% PATTERNDROPOUTFUNC - Drop out a random percentage of the pattern
%
%
%   Parameters:
%       percent - the percentage of the pattern to randomly delete
%

keepers = rand(size(inp)) > percent;
inp = inp(keepers);
ts = ts(keepers);

end