function [ inp, ts ] = mixedvariancefunc(inp, ts, variances)
%% MIXEDVARIANCEFUNC - Modify a pattern to have mixed variance
%
%
%   Parameters:
%       variances - max variance for each input where variances(i) is for
%           inp(i) and numel(variances) == Np
%
%   NOTE: Within a given presentation a given neurons spikes are all
%   shifted in the same way, That is if the neurons presentations variance
%   is -2 ms and that neuron has two spikes in the pattern, BOTH will be
%   moved by -2. Variance is calculated per neuron per presentation not per
%   spike per presentation. 

inst_var = round(randn(size(variances)) .* variances);
tweaks = inst_var(inp);
%tweaks = round(min(max(randn(size(ts)) * max_variance, -max_variance), max_variance));
ts = ts + tweaks;
    


end