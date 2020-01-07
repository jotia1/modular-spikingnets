function [inp , ts ] = gaussianjitter(inp, ts, jitter)
%% GAUSSIANJITTER - Designed to implement jitter from Masquelier paper
%   Add Gaussian noise to a signal with mean 0 and standard deviation
%   jitter
%
%   Parameters
%       inp, ts - spike train 
%       jitter - integer number of ms noise to add
%
%   Usage: 
%       [inp , ts ] = gaussianjitter(inp, ts, 1)

% TODO : Hardcoded wrap around 50ms pattern
% Need to subtract 1 then do wrapping, then add one to account for non-zero
% indexing of matlab and now event streams - Otherwise 50ms wraps round to
% 0 ms. (which is bad, should wrap from 51 to 1).
ts = mod(ts + round(randn(size(ts)) .* jitter) - 1, 50) + 1;

end