function [ patt_inp, patt_ts ] = generateuniformpattern( patt_length_ms, num_neurons )
%GENERATEUNIFORMPATTERN - Generate spike times for 'num_nurons' neurons
%   Generate a firing time for each of a number of neurons with a uniform
%   distribution.
%
%   Example usage:
%       [ patt_inp, patt_ts ] = generateuniformpattern( 100, 10 )

    patt_inp = (1:num_neurons);
    patt_ts = sort(randi([1, patt_length_ms], 1, num_neurons));
    
end

