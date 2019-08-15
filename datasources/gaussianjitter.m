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

ts = ts + round(rand(size(ts)) * jitter);

end