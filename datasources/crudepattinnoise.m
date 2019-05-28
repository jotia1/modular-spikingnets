function [ inp, ts, patt_inp, patt_ts ] = crudepattinnoise( N_total, N_patt, patt_inp, patt_ts, patt_length_ms, total_length_ms )
%CRUDEPATTINNOISE - Insert a pattern on top of noise ignoring issues
%   Given a pattern embed the pattern on top of background noise, ignore
%   any issues to do with neuron firing rates and generating data with good
%   statistics. Defaults to generating a uniform pattern in 1 second of
%   noise. 
%   
%   [inp, ts] = crudepattinnoise( 10, 3 )
%
%   Example usage:
%   

if ~exist('N_total', 'var') || isempty(N_total)
    N_total = 3;
end

if ~exist('N_patt', 'var') || isempty(N_patt)
    N_patt = 3;
end

if ~exist('patt_length_ms', 'var') || isempty(patt_length_ms)
    patt_length_ms = 15;
end

if ~exist('total_length_ms', 'var') || isempty(total_length_ms)
    total_length_ms = 1000;
end

if ~exist('patt_inp', 'var') || ~exist('patt_ts', 'var') || ...
    isempty(patt_inp) || isempty(patt_ts)
    [patt_inp, patt_ts] = generateuniformpattern(patt_length_ms, N_patt);
end

noise_freq = 5;     % [Hz]
[noise_inp, noise_ts] = generatenoise(N_total, noise_freq, total_length_ms);

[inp, ts] = joineventstreams(patt_inp, patt_ts, noise_inp, noise_ts);

end

